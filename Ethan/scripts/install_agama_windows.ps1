# Build the supplied AGAMA source in the same virtual environment as the N-Baddies GUI.
# Run from PowerShell: powershell -ExecutionPolicy Bypass -File .\scripts\install_agama_windows.ps1
$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$python = Join-Path $root '.venv\Scripts\python.exe'
$bat = Join-Path $PSScriptRoot 'install_agama_windows.bat'

if (-not (Test-Path $python)) {
    throw "Project virtual environment missing. First run start_gui.bat in $root"
}
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path $vswhere)) {
    throw 'Visual Studio Build Tools not found. Install the Desktop development with C++ workload (VS 2022 or newer).'
}
$vs = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vs) {
    throw 'MSVC C++ tools not found. Open Visual Studio Installer and add Desktop development with C++.'
}
$vcvars = Join-Path $vs 'VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) {
    throw "Cannot find MSVC environment setup: $vcvars"
}

Write-Host "AGAMA build with compiler from $vs"
Write-Host "Python environment: $python"
& $python -c "import platform,sys; print('Python:',sys.version); print('Architecture:',platform.architecture()[0])"

# AGAMA's Windows setup.py searches sys.exec_prefix\libs. In a virtual environment
# sys.exec_prefix is the venv root, which normally lacks python312.lib. Link the
# compiled extension against the matching base interpreter's import library.
$baseLib = & $python -c "import os,sys; print(os.path.join(sys.base_prefix, 'libs', 'python%d%d.lib' % sys.version_info[:2]))"
if ($LASTEXITCODE -ne 0) { throw 'Could not determine base Python library path.' }
if (-not (Test-Path -LiteralPath $baseLib)) {
    throw "Python import library not found at $baseLib. The Microsoft Store Python installation may lack development libraries. Install the x64 Python 3.12 distribution from python.org, recreate the .venv with it, and retry."
}
$venvLibs = Join-Path $root '.venv\libs'
New-Item -ItemType Directory -Path $venvLibs -Force | Out-Null
Copy-Item -LiteralPath $baseLib -Destination (Join-Path $venvLibs (Split-Path $baseLib -Leaf)) -Force
Write-Host "Provided Python import library: $baseLib -> $venvLibs"

# AGAMA upstream setup.py has an inconsistent Windows venv check:
# it reports sys.base_prefix from the compiled extension, then compares it
# against sys.prefix. Patch ONLY that check in an extracted copy of the source.
$sourceZip = Join-Path $root 'vendor\Agama-master.zip'
if (-not (Test-Path -LiteralPath $sourceZip)) {
    throw "Cannot find AGAMA source archive at $sourceZip"
}
$sourceTemp = Join-Path $env:TEMP "nbody-agama-venv-fix-$PID"
New-Item -ItemType Directory -Path $sourceTemp -Force | Out-Null
try {
    Expand-Archive -LiteralPath $sourceZip -DestinationPath $sourceTemp -Force
    $agamaSource = Join-Path $sourceTemp 'Agama-master'
    $setupPath = Join-Path $agamaSource 'setup.py'
    if (-not (Test-Path -LiteralPath $setupPath)) {
        throw "AGAMA setup.py missing from extracted source: $setupPath"
    }
    $setup = [System.IO.File]::ReadAllText($setupPath)
    $old = "sysprefix = 'sys.prefix=' + os.path.realpath(sys.prefix)"
    $new = "sysprefix = 'sys.prefix=' + os.path.realpath(sys.base_prefix)"
    if (-not $setup.Contains($old)) {
        throw 'The uploaded AGAMA version has changed; cannot safely apply Windows venv fix.'
    }
    $setup = $setup.Replace($old, $new)
    [System.IO.File]::WriteAllText($setupPath, $setup, ([System.Text.UTF8Encoding]::new($false)))
    Write-Host 'Patched AGAMA Windows virtual-environment prefix check in temporary build source.'

    # The supplied source uses the standard C++ alternative token "or" in one
    # condition, but the MSVC 2019 default options do not enable that spelling.
    # Replace it with the identical || operator only in our temporary copy.
    $interfacePath = Join-Path $agamaSource 'src\interface_python.cpp'
    $interface = [System.IO.File]::ReadAllText($interfacePath)
    $oldOperator = "if(type_str[14]=='t' or type_str[14]=='T')"
    $newOperator = "if(type_str[14]=='t' || type_str[14]=='T')"
    if (-not $interface.Contains($oldOperator)) {
        throw 'AGAMA source has changed; cannot safely apply MSVC alternative-operator fix.'
    }
    $interface = $interface.Replace($oldOperator, $newOperator)
    [System.IO.File]::WriteAllText($interfacePath, $interface, ([System.Text.UTF8Encoding]::new($false)))
    Write-Host 'Patched AGAMA MSVC alternative-operator compatibility in temporary build source.'


    # AGAMA's math_cubature.h is missing the Windows-specific header guard
    # already used in its other C++ files. MSVC has <malloc.h> rather than
    # POSIX <alloca.h>. Patch only the extracted copy, never the vendor ZIP.
    $cubaturePath = Join-Path $agamaSource 'src\math_cubature.h'
    if (-not (Test-Path -LiteralPath $cubaturePath)) {
        throw "AGAMA source missing expected header: $cubaturePath"
    }
    $cubature = [System.IO.File]::ReadAllText($cubaturePath)
    $oldAllocaInclude = '#include <alloca.h>'
    $newAllocaInclude = @'
#ifndef _MSC_VER
#include <alloca.h>
#else
#include <malloc.h>
#endif
'@
    if (-not $cubature.Contains($oldAllocaInclude)) {
        throw 'AGAMA source has changed; cannot safely apply MSVC alloca header fix.'
    }
    $cubature = $cubature.Replace($oldAllocaInclude, $newAllocaInclude)
    [System.IO.File]::WriteAllText($cubaturePath, $cubature, ([System.Text.UTF8Encoding]::new($false)))
    Write-Host 'Patched AGAMA MSVC alloca.h compatibility in temporary build source.'

    # MSVC environment variables must persist inside the same cmd.exe process as the build.
    $tmp = Join-Path $env:TEMP "nbody-agama-$PID.cmd"
    @(
        '@echo off',
        "call `"$vcvars`"",
        'if errorlevel 1 exit /b 1',
        "call `"$bat`" `"$agamaSource`"",
        'exit /b %ERRORLEVEL%'
    ) | Set-Content $tmp -Encoding ASCII
    try {
        & $env:ComSpec /d /c $tmp
        if ($LASTEXITCODE -ne 0) {
            throw "AGAMA build failed (exit code $LASTEXITCODE). Read compiler/GSL errors above."
        }
    } finally {
        Remove-Item $tmp -ErrorAction SilentlyContinue
    }
} finally {
    Remove-Item -LiteralPath $sourceTemp -Recurse -Force -ErrorAction SilentlyContinue
}
Write-Host 'Verifying the AGAMA extension in the GUI environment:'
& $python -c 'import agama,sys; print("AGAMA:",agama.__file__); print("G:",agama.G); print("Interpreter:",sys.executable)'
if ($LASTEXITCODE -ne 0) { throw 'AGAMA import failed in the project environment.' }
Write-Host 'Next: .\.venv\Scripts\python.exe -m experiments.agama_selftest'
