# Reproduced results — preliminary investigation, not a stable-galaxy release

## What changed our diagnosis

Fixing the softened equilibrium was necessary for mean-force consistency,
but was NOT sufficient to stop the live disk thickening. The corrected model
still grows in height by 59.85% over t=0.5. The earlier emphasis on softening
as the main explanation was too strong. The old and new runs have different
DFs/realizations, so their 55.04% and 59.85% growth cannot rank the methods.

The corrected mean odd vertical force at R=Rd,z=zd is 1.0009 times its own
smooth reference. The previous audit compared a softened particle field to
an unsoftened reference and obtained about 0.66. The field-consistency error
is corrected at these probes; remaining growth must not be labelled solved.

## Matched controls

| Experiment | Duration | Height change | Vertical dispersion change |
|---|---:|---:|---:|
| Corrected fixed smooth field, disk tracers | 0.5 | -0.061% | -0.655% |
| Corrected fully live 12,000 particles | 0.5 | +59.847% | +26.308% |
| Same live disk, fixed smooth spheroids | 0.5 | +6.474% | +3.390% |

The last control removes both spheroid discreteness AND its dynamical response.
It is a diagnostic, not a substitute merger model: it omits halo wakes and
self-consistent dynamical friction. It cannot identify numerical halo heating
alone.

| Same corrected DF, t=0.125 | Disk N | Bulge N | Halo N | Height change |
|---|---:|---:|---:|---:|
| Baseline | 4800 | 1200 | 6000 | +12.041% |
| Disk resolution x4 | 19200 | 1200 | 6000 | +10.502% |
| Halo resolution x4 | 4800 | 1200 | 24000 | +6.490% |

The disk-resolution run copies the baseline spheroid particles exactly.
The halo-resolution run copies the baseline disk and bulge exactly. These
copies were checked with array_equal for position, velocity and mass.
Only the changed component is sampled anew from the reconstructed fixed DF.
The static smooth potential is copied verbatim. Seed/noise uncertainty remains;
one realization is not a resolution-convergence demonstration.

Together the controls implicate live spheroids, with halo resolution showing
a larger benefit than stellar resolution in this short test. This supports
spurious halo heating/noise-driven response as a leading explanation, not a
unique proof excluding collective modes or all bulge effects. A force-scatter
budget is not itself a heating budget.

The 512-step fully live energy change is +1.5874e-6 of the initial absolute
energy. This demonstrates good global energy behavior, NOT equilibrium or
adequate spatial resolution. Integration refinement on the corrected model
remains to be done.

## Checks actually completed

- Kernel versus independent adaptive spherical integral; gradient/symmetry;
  small-softening limit: PASS.
- Parallel direct versus original serial direct: bitwise equality in random
  particle tests at two softenings, 1/2/4 threads, and on the supplied 12k snapshot.
- Local 12k benchmark: serial 0.655 s, parallel four-thread 0.213 s per force.
  Timing is machine-specific and includes contention from other local jobs.
- New SCM: seven iterations; last sampled force/density change 0.2265%.
- Quadrature and expansion force checks at disk probes: <=0.636% vertical,
  <=0.258% radial, compared with the 36/18/36 refined quadrature used for
  the saved candidate. Current source uses an even finer 48/24/48 check.
- Fresh non-interpolated Staeckel action finder density checks: max 2.617%,
  2.193%, 1.245% for disk/bulge/halo, using the stated density floor.
  Non-interpolated does not mean mathematically exact actions.
- DF versus interpolated density integrated masses: max discrepancy 1.212%.
- Matched smooth/live runs, both separate component-resolution controls,
  and the live-disk/static-spheroid control: completed; JSON histories included.

Provisional numerical construction gates (2% force,3% density,2% mass) are
engineering screens, not literature-derived stability criteria.

An exploratory epsilon=0.0125 construction with the original coarse quadrature
was REJECTED by its numerical gate (vertical discrepancy about 3.19%). No live
result from it is claimed. Current code chooses quadrature by an upfront
refinement check; the larger-epsilon live experiment has not been reproduced.

## Next work, in order

1. Refine dt on the corrected halo-24k snapshot at the same t=0.125.
2. Repeat independent seeds to measure variation, then extend the halo-24k
   isolation control over the intended encounter duration and measured orbital
   periods. Do not infer its t=0.5 behavior by multiplying the short result.
3. Test larger softening only with regenerated, numerically checked DFs;
   compare force resolution and the measured physical profiles as well as drift.
4. If adequate halo mass resolution makes runs expensive, compile the tree
   construction as well as traversal and verify it against direct gravity.
   The current Python-built tree was slower at the measured N values.
5. Choose a defensible production resolution/softening from these comparisons,
   then enable that generator in the GUI/library. Existing cached galaxies
   remain unvalidated and were not overwritten by this experimental patch.

The corrected construction should be retained, but no included live snapshot
is certified as a stable Mice progenitor or a validated Mice encounter.
