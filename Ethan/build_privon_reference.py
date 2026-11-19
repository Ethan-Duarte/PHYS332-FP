import math
import os
import subprocess
from datetime import datetime
from pathlib import Path

root = Path.home() / "nbody-tools/zeno-reference"
out = root / (
    "privon-build-" + datetime.now().strftime("%Y%m%d-%H%M%S")
)
out.mkdir()

env = dict(os.environ, ZENOPATH=str(root))
log = out / "generation.log"

# Correct the radial derivative term in a separate source copy.
# Appendix B13 requires a squared denominator.
source = root / "src/nbody/gspmodels/gspdisk.c"
text = source.read_text()
old = "r_mu/(r+r_mu)*(r+r_mu)"

if text.count(old) != 2:
    raise RuntimeError(
        "Unexpected disk source layout. Stop and send this output."
    )

patched = out / "gspdisk_b13.c"
patched.write_text(
    text.replace(old, "r_mu/((r+r_mu)*(r+r_mu))")
)

subprocess.run(
    [
        "gcc", "-std=gnu99", "-DLINUX",
        "-include", "unistd.h",
        "-Wno-error=implicit-function-declaration",
        "-I" + str(root / "inc"),
        "-I" + str(root / "src/nbody/gspmodels"),
        "-L" + str(root / "lib"),
        "-o", str(out / "gspdisk_b13"),
        str(patched),
        str(root / "lib/libGSP.a"),
        "-lNBody", "-lClib", "-lgsl", "-lgslcblas", "-lm",
    ],
    check=True,
)


def run(name, *arguments):
    executable = (
        out / "gspdisk_b13"
        if name == "gspdisk"
        else root / "bin" / name
    )
    command = [str(executable), *map(str, arguments)]
    print("Running:", name, flush=True)

    with log.open("a") as stream:
        stream.write("\n" + " ".join(command) + "\n")
        stream.flush()
        subprocess.run(
            command,
            env=env,
            stdout=stream,
            stderr=stream,
            check=True,
        )


# Tapered Hernquist bulge: normalize including its outer taper.
a, b, target_mass = 0.02, 4.0, 0.0625
rho_b = a / (2 * math.pi * b * (a + b) ** 3)
inner_mass = b * b / (a + b) ** 2
normalization = target_mass / (
    inner_mass + 2 * math.pi * b**3 * rho_b
)

table = out / "bulge.txt"
with table.open("w") as stream:
    stream.write("# tapered Hernquist bulge\n")
    stream.write("# radius density mass density_copy\n")

    for i in range(2049):
        r = 1e-5 * math.exp(math.log(1e6) * i / 2048)

        if r <= b:
            rho = normalization * a / (
                2 * math.pi * r * (a + r) ** 3
            )
            mass = normalization * r*r / (a + r) ** 2
        else:
            e = math.exp(-2 * (r - b) / b)
            rho = normalization * rho_b * (b/r)**2 * e
            mass = normalization * (
                inner_mass
                + 2 * math.pi * b**3 * rho_b * (1-e)
            )

        stream.write(
            f"{r:.16e} {rho:.16e} {mass:.16e} {rho:.16e}\n"
        )

run(
    "tablegsp",
    f"in={table}", f"out={out}/bulge.gsp",
    "smooth=false", "mtot=.0625",
)

run(
    "halogsp",
    f"out={out}/halo.gsp",
    "m_a=.16", "a=.25", "b=.98015",
    "taper=sw", "npoint=2049", "rrange=0.00001:10",
)

run(
    "expdgsp",
    f"out={out}/disk_spherical.gsp",
    "mtot=.1875", "alpha=12", "zdisk=0",
    "npoint=2049", "rrange=0.00001:10",
    "smartrange=true",
)

run(
    "gspadd",
    f"in={out}/bulge.gsp,{out}/halo.gsp",
    f"out={out}/spheroid.gsp",
)

run(
    "gspadd",
    f"in={out}/spheroid.gsp,{out}/disk_spherical.gsp",
    f"out={out}/total.gsp",
)

run(
    "gspsmooth",
    f"in={out}/total.gsp",
    f"out={out}/total_smooth.gsp",
    "eps=.0115", "kappa=1.975",
)

run(
    "gspsmooth",
    f"in={out}/spheroid.gsp",
    f"out={out}/spheroid_smooth.gsp",
    "eps=.0115", "kappa=2.025",
)

# First generator check: quarter of the paper's follow-up count.
for component, count, seed in [
    ("bulge", 2048, 42),
    ("halo", 8192, 43),
]:
    run(
        "gsprealize",
        f"gsp={out}/{component}.gsp",
        f"grav={out}/total_smooth.gsp",
        f"out={out}/{component}.snap",
        f"nbody={count}", f"seed={seed}",
    )

run(
    "gspdisk",
    f"grav={out}/spheroid_smooth.gsp",
    f"out={out}/disk.snap",
    "mdisk=.1875", "alpha=12", "zdisk=.0075",
    "epsilon=.010606601717798213",
    "mu=2", "r_mu=.075", "eta=.25",
    "rcut=1", "ndisk=6144", "seed=44",
)

print("\nGenerated components:", out)
print("Generation log:", log)

for line in log.read_text().splitlines():
    if line.startswith("0.2500"):
        print("At R=3Rd; columns begin radius, circular speed:")
        print(line)
    if "WARNING" in line:
        print(line)

print("Generation finished. No evolution launched.")