AGAMA Windows MSVC alloca.h fix
=================================

This patch retains the existing Python import-library, virtual-environment
prefix, and alternative C++ operator fixes. It also corrects the missing
MSVC-specific header guard in AGAMA src/math_cubature.h.

Use:
1. Extract ZIP INTO your existing PHYS332-FP\Ethan project directory.
   The file scripts\install_agama_windows.ps1 must REPLACE the existing one.
2. Open PowerShell inside the Ethan project directory.
3. Run:
   powershell.exe -ExecutionPolicy Bypass -File .\scripts\install_agama_windows.ps1
4. If AGAMA asks to download GSL or Eigen, answer Y.
   CVXOPT is optional for your galaxy-generation use case.
5. After a successful install, verify:
   .\.venv\Scripts\python.exe -c "import agama; print(agama.__file__)"
   .\.venv\Scripts\python.exe -m experiments.agama_selftest

Only the Windows installer script is replaced. The vendor AGAMA ZIP,
application physics code and results are not modified. A new temporary
copy of AGAMA is created and patched for each installation attempt.

This patch is source-inspected and verified for its transforms, but has not
been compiled with MSVC in this Linux environment.
