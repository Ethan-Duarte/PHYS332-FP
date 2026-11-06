# Optional AGAMA backend: status and setup

- Reference: https://github.com/GalacticDynamics-Oxford/Agama
- Reference workflow: https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_self_consistent_model3.py
- Paper: https://doi.org/10.1093/mnras/sty2672

`nbody/galaxies/sources/agama.py` implements a three-component model based on the documented API: a vertically sech² exponential disk, Hernquist (`gamma=1`) spheroid and cored Dehnen (`gamma=0`) halo; construct initial density potential, derive a quasi-isothermal disk DF and quasi-spherical bulge/halo DFs, iterate a self-consistent potential, then sample particles from each DF.

AGAMA includes a compiled extension and is **not included in `requirements*.txt`**. Consult its official install instructions and a compatible Python version / compiler. On Windows, a supported Linux/WSL setup may be easier than compiling extensions directly. Once available in the same Python environment, verify with:

```python
import agama
print(agama.__file__)
```

Then choose **AGAMA (equilibrium candidate)** in the GUI; it will not silently switch to native preview on failure. Save the returned snapshot, inspect the component masses and measured density profiles, check the gravitational constant / AGAMA units against `G=1`, and run isolated stability for several dynamical times. Investigate any mass mismatch between the generated DF and the requested input configuration before treating it as equilibrium.

AGAMA has its own random sampler; the adapter calls `agama.setRandomSeed(seed)` per the documented API. Even with a seed, save the exact particle snapshot for platform-independent reproducibility.

**Environment note:** AGAMA was unavailable in the build environment and its installation could not be completed there; the adapter has not been run against the extension. Do not cite the adapter as a validated galaxy generator yet.
