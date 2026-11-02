# Walkthrough: What the code is doing, and why

## The data flow

`mice.json` -> `galaxies/generate.py` -> `Particles` -> `encounters.parabolic_pair()` -> `integrate.evolve()` using `gravity/direct.py` or `gravity/barnes_hut.py` -> `analysis.py` -> `io.save()`.

Think of this as **three distinct physics operations**: create a self-consistent galaxy, place it on an external orbit, then evolve all live particles under gravity. We already implement the data plumbing, but the preview doesn't yet provide equilibrium disk/halo velocities.

## 1. particles.py (data container)

`x[i]` is the three-dimensional position of particle i; `v[i]` is its three-dimensional velocity; `m[i]` is its positive mass. Two integer labels track which physical component (disk/bulge/halo) and which original galaxy a particle belongs to, even after the galaxies interpenetrate.

Example: `x.shape == (500,3)` means 500 particles, each with x,y,z. `m.shape == (500,)` means one mass per particle. The force calculation only needs x and m; the leapfrog uses v too; the tail analysis needs the origin labels.

## 2. gravity/direct.py (the reference force)

For target particle i, the acceleration is

`a_i = G sum_(j != i) m_j (x_j-x_i) / (|x_j-x_i|^2 + eps^2)^(3/2)`.

The `d[i,j,:]` array is the vector from i to j; we add squared softening and make the diagonal infinite to remove self-interaction. This NumPy implementation is concise and fast for modest N but stores an N×N separation matrix.

## 3. gravity/barnes_hut.py (approximation)

The octree recursively divides a cube into up to eight children. Each node records total mass, center of mass and physical size. If a group is far enough from the target (`size / distance < theta`) and does **not** contain the target, the group is approximated by its total mass at its center of mass. Otherwise the tree opens the node and sums its children. Leaves calculate pairwise forces explicitly.

A smaller theta is normally more accurate and more expensive. `theta=0` opens all nodes and should approach the same softened force as direct summation. Our original Python-recursion implementation may be slower than the NumPy direct version for small N.

## 4. integrate.py (kick–drift–kick)

At each step of length dt:

`v_half = v + (dt/2) a(x)`;
`x_new = x + dt v_half`;
`v_new = v_half + (dt/2) a(x_new)`.

The acceleration is calculated through a callable `force(x,m)`, so the integrator doesn't need to know whether the force came from direct summation or Barnes–Hut. It yields independent snapshots instead of allocating a huge history automatically. Leapfrog is time-reversible and symplectic for a fixed-step separable Hamiltonian under an exact consistent force; a Barnes–Hut approximation can weaken some of the ideal conservation properties.

## 5. galaxies/profiles.py (mathematical definitions)

These functions encode target *spatial distributions*. For Dehnen spheres, `M(<r) = M [r/(r+a)]^(3-gamma)`; gamma=1 is Hernquist and gamma=0 is the cored halo. For an exponential disk, `M(<R) = M[1-(1+R/Rd)exp(-R/Rd)]`.

These formulae say **where** particles should be, not **how fast they should move**. Density cannot by itself specify a general equilibrium velocity distribution.

## 6. galaxies/sources/native.py (sampling)

The *Plummer* generator is a physically useful test. Its enclosed mass gives an invertible radial CDF. For its isotropic distribution function `f(E) proportional E^(7/2)`, the squared speed as a fraction of escape speed follows a `Beta(3/2,9/2)` distribution. We use these laws to sample positions and speeds, and isotropic random directions for each. This is exact in the continuous unsoftened Plummer model, not a claim that the finite-N softened realization will remain perfectly motionless.

The **Mice preview**, by contrast, samples an exponential/sech² disk and Hernquist/Dehnen spheres correctly in space (with a finite cutoff) but uses a spherical circular-speed approximation plus heuristic velocity dispersions. This is useful to test the data flow, NOT a physically equilibrium-validated galaxy.

## 7. galaxies/sources/agama.py (optional equilibrium generator)

AGAMA builds a combined gravitational potential, uses a disk quasi-isothermal DF and quasi-spherical DFs for bulge/halo, recomputes the component density fields, and iterates. It then samples phase-space coordinates. In our adapter, a negative disk scaleHeight is **half** the Barnes sech² denominator because AGAMA writes the vertical law as sech²(z/(2h)).

The adapter is based on the upstream examples but **has not been executed in our environment**, where AGAMA isn't installed. Do not equate code written with tested equilibrium generation. AGAMA's own internal sampling RNG is also not guaranteed here to follow our NumPy seed.

## 8. encounters.py (moving whole galaxies)

Once a model galaxy is generated, rigidly rotate positions and velocity vectors, translate the galaxy center, and add a bulk velocity. No internal galaxy physics should be changed by that placement. For a parabolic two-point-mass reference orbit, specific orbital energy is zero:

`v^2 = 2G(M1+M2)/r` and `h = sqrt(2G(M1+M2) r_peri)`.

Therefore at separation r, `v_t = h/r` and incoming `v_r = -sqrt(v^2-v_t^2)`. These formulas specify an *initial reference orbit*; live extended halo interaction changes it. Actual pericenter must be measured.

## 9. analysis.py (diagnostics)

The energy is kinetic `0.5 sum m v²` plus the Plummer-softened pair potential `-sum_(i<j) G m_i m_j/sqrt(r_ij²+eps²)`. The potential calculation is O(N²), so use it on modest snapshots or sampled particles. Momentum and angular momentum check for unexpected numerical drift.

`core_separation` uses the center of mass of each bulge's sampled particles, not the exact parameterized galaxy COM. At low N these differ by statistical noise; the observed pair core separation also depends on how an observer defines the nuclear centroid.

## 10. io.py (provenance)

Saved `.npz` files include position/velocity/mass arrays, component and galaxy IDs, and JSON metadata. The config and implementation versions should be recorded for every result you later report, because a seed alone does not guarantee the same output after generator implementation changes.

## The scientific gate

Passing `pytest` demonstrates specific implementation checks. It does NOT validate the disk–bulge–halo equilibrium. Next: run AGAMA on a supported installation, check its output and initial requested mass profiles, evolve the resulting galaxy in isolation at multiple resolutions, and only then conduct the Mice observational comparison.
