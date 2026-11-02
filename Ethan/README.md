# N-Baddies: N-body Galaxy Simulator — Version 1

A small, fresh Python project. The direct force calculation, leapfrog, and octree descend from our earlier student implementations, but the file layout is new.

## Status (important)

- **Working/tested:** softened direct gravity, Barnes–Hut tree gravity, leapfrog, particle structure, saved snapshots, exact isotropic Plummer sphere, analytic density profiles, illustrative disk–bulge–halo sampling, simple parabolic encounter, basic diagnostics.
- **Illustrative only:** `mice_preview()` generates the published-style *spatial* components, but uses heuristic velocities. Its galaxy is **not** a validated collisionless equilibrium; resulting tidal shapes cannot yet be presented as a Mice reconstruction.
- **Unverified optional backend:** `galaxies/sources/agama.py` follows AGAMA's documented self-consistent DF process, but AGAMA is not installed in the test environment. Its API behavior, equilibrium quality, sampling reproducibility, and unit consistency require testing on a local AGAMA installation.
- **Still needed for real validation:** isolated galaxy stability for disk–bulge–halo, adequate resolution, velocity profiles, convergence, physically matched orbital/viewing parameters, observational comparisons.

## Quick start

Run from this folder:

```bash
pip install -r requirements.txt
python -m pytest -q
python -m experiments.two_body
python -m experiments.plummer_validation
python -m experiments.mice_preview
python -m experiments.force_benchmark
```

Mice preview writes `results/mice_preview.png` and `results/mice_preview.npz`.

Optional AGAMA is independently installed (`pip install agama`, which compiles a C++ extension and may require platform-specific dependencies); use:

```python
import json
from nbody.galaxies.generate import galaxy
cfg = json.load(open('presets/mice.json'))
model = galaxy(cfg, backend='agama')
```

**AGAMA seed warning:** the `seed` parameter governs native NumPy sampling only. AGAMA's internal sampler isn't established here to use that seed; save the resulting `.npz` for exact replay.

## Directory map

- `nbody/particles.py`: common particle arrays `(x,v,m,component,galaxy)`, combine and recenter.
- `nbody/gravity/`: direct and octree-based Barnes–Hut acceleration.
- `nbody/integrate.py`: kick–drift–kick timestep and periodic snapshots.
- `nbody/galaxies/profiles.py`: Hernquist/Dehnen and exponential disk analytic profiles.
- `nbody/galaxies/sources/native.py`: Plummer exact equilibrium and approximate Mice preview.
- `nbody/galaxies/sources/agama.py`: separate, optional multicomponent equilibrium backend.
- `nbody/galaxies/generate.py`: generation entry point.
- `nbody/encounters.py`: orientations, offsets, Keplerian parabolic initial orbit.
- `nbody/analysis.py`: energy, momentum, radius and core separation.
- `nbody/io.py`: portable snapshot writing/loading.
- `presets/mice.json`: Barnes-inspired mass and spatial scales, developer-level particle counts.
- `experiments/`: independently executable scripts.
- `tests/`: automated mathematical and plumbing checks.

## Units and model details

Internal gravitational constant is **G=1**. Barnes (2004) provides a physical mapping for his reference Mice fit: 1 length unit ≈35.7 kpc, 1 speed unit ≈180 km/s, 1 time unit ≈194 Myr, 1 mass unit ≈2.69×10¹¹ solar masses. These are *model-matching scales*, not fundamental conversion constants.

The `mice.json` reference uses Barnes's idealized collisionless density shapes: stellar disk ∝ exp(-R/Rd) sech²(z/zd), Hernquist bulge (Dehnen gamma=1), cored Dehnen halo (gamma=0). It does **not** reproduce Barnes's SPH gas component. Truncated spatial sampling is conditionally normalized to the specified component mass, so this is not an exact untruncated Barnes realization.

`encounters.parabolic_pair` places the galaxy COMs on an initial *point-mass* parabolic Kepler orbit. The live, extended halos alter the subsequent orbit; measure the actual pericenter from the evolution. Disk `rotation` adopts the explicit convention `Rz(node_angle) @ Rx(inclination)`; it is not a published sky-viewing transformation.

**Softening**: both gravity solvers calculate acceleration proportional to `(r² + eps²)^(-3/2)` with the same `eps`. The diagnostic energy uses the matching potential `-Gmimj/sqrt(r² + eps²)`.

**Barnes–Hut performance**: this clean Python implementation is useful for learning and force comparison. It may be slower than NumPy-vectorized direct summation for modest N; it needs a separate performance/optimization study before large simulations.

See `docs/WALKTHROUGH.md` for a conceptual lesson, `docs/VALIDATION.md` for the tests required before making scientific claims, and `docs/BIBLIOGRAPHY.md` for the running references.

## References

- Barnes, J. E. (2004), *Shock-induced star formation in a model of the Mice*: https://doi.org/10.1111/j.1365-2966.2004.07725.x
- Hernquist, L. (1990), analytic spherical galaxy: https://doi.org/10.1086/168845
- Dehnen, W. (1993), spherical gamma-family profiles: https://doi.org/10.1093/mnras/265.1.250
- Barnes & Hut (1986), tree algorithm: https://doi.org/10.1038/324446a0
- Vasiliev (2019), AGAMA: https://doi.org/10.1093/mnras/sty2672
- AGAMA three-component worked example: https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_self_consistent_model_bdh.py
