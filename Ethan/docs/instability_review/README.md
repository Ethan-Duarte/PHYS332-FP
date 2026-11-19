# Softened equilibrium investigation

Baseline inspected: PHYS332-FP commit 793b7286108942bee375f0ebf9d2915938c46ff2.
Input: user-uploaded instability_diagnostic_bundle.zip, November 18, 2026.

## Evidence and limits

The saved new smooth-field control has r50 -1.14%, fitted-plane RMS height
-1.08%, sigma_z +0.16%, over t=0.5. The corresponding live 12,000-particle
control has r50 +1.11%, height +55.04%, sigma_z +8.67%. These are one-seed
screens, not stability certification. Timestep refinement on earlier cached
models makes little difference to their thickening; it does not automatically
establish timestep convergence for a new softening/model.

Direct convolution of the exported smooth potential with a Plummer kernel at
epsilon=0.00625 gives vertical restoring force 0.6697 of the unsoftened
reference at R=Rd,z=zd (36 radial,18 polar,36 azimuth quadrature).
The 24,000-halo particle audit gives 0.658 there. Thus a systematic difference
exists even before finite-particle fluctuations. This is not a unique causal
explanation of all disk thickening.

Component noise was recomputed using 96 azimuths and the odd vertical force
at +/-zd. Values are spatial standard deviations divided by the original
smooth restoring force, NOT standard errors or measured heating rates.
At R=2Rd, epsilon=0.0015625, the shared 4,800-particle disk contributes 30.8%,
the 1,200-particle bulge 0.6%, and the 24,000-particle halo 14.2%.
At epsilon=0.00625 these become 11.5%,0.5%,4.1%. Covariances mean the
component standard deviations must not be summed as scalar variances without
checking correlations. Only one disk/bulge realization and two independent
halo realizations are available. The odd-force diagnostic suppresses even
vertical offsets and does not characterize every bending perturbation.

Requested component masses are present in the snapshots. Untruncated sampled
bulge/halo tails extend thousands of code-length units from the center. Global
center-of-mass recentering is translation invariant for a fully isolated live
run, but using it to align particles with an external potential or to place
visible galaxy centers on encounter orbits can be misleading. The supplied
control avoided that recentering.

Compiled direct is measured faster than the current tree at 12k and 30k.
Use direct for these controls, with benchmarked forecasts and checkpoints.
This says nothing universal about the asymptotic efficiency of Barnes-Hut.

## Candidate correction

`agama_softened_control.py` is a separate experimental generator. The original
GUI adapter and all existing snapshots remain unchanged.

1. Construct the initial smooth field and convolve it with the Plummer kernel.
2. Define action-based disk/bulge/halo DFs using that softened initial field.
3. Recompute each DF's unsmoothed particle density in the current softened field.
4. Construct its Newtonian field, then convolve that field with the same kernel.
5. Repeat until sampled forces/densities converge. Refresh the ActionFinder
   whenever the potential changes. Keep each DF fixed as a function of actions
   during this loop; do not reinitialize it on every iteration.
6. Sample the unsmoothed DF density, not the smoothed gravitational source.
   Normalize sampled component masses to the DF integrals, not unrelated
   original analytic target masses.
7. Check quadrature refinement and potential-expansion errors separately.

The kernel is W_e(r)=3e^2/[4*pi*(r^2+e^2)^(5/2)]. Convolving the Newtonian
potential with W_e produces the smooth field corresponding to the live
Plummer pair interaction. This does not remove finite-N noise or suppress
physical collective modes. The original density parameters are starting
targets; the actual final density, thickness, and rotation must be measured.
This is a collisionless Mice-style approximation, not a gas simulation or a
reproduction of every Barnes (2004) modelling choice.

## Validation protocol

- Kernel: independent adaptive spherical integral, potential/force gradient,
  symmetry, and small-epsilon limit.
- Field construction: refined convolution and expansion checks. The provisional
  2% numerical gate is an engineering screen, not a literature stability bound.
- Sampled field: component-resolved exact-force audit; interpret spatial noise
  separately from mean mismatch.
- Short screen: matched smooth and live runs, default 32 steps. A short pass
  cannot establish long-term stability.
- Isolation: t=0.5 and timestep refinement, then several measured orbital
  periods covering the intended encounter duration. Inspect time histories,
  annular plane-removed thickness, sigma_z, radial size, tilt, A2 and energy.
- Repeat seeds and increase stellar/halo resolution independently. Genuine
  collective instability and numerical heating require different remedies.
- After a candidate is acceptable, test tree versus direct from identical
  snapshots, using a separate force/trajectory convergence study.

Do not evolve old snapshots under a new epsilon and call them corrected ICs.
Do not globally rescale vertical velocities, erase vertical forces, or damp
an isolated run until it appears quiet and call that a validated equilibrium.

## Source inventory

- Barnes (2004), Mice model: https://arxiv.org/abs/astro-ph/0402248
- Vasiliev (2019), AGAMA: https://arxiv.org/abs/1802.08239
- AGAMA reference/API: https://github.com/GalacticDynamics-Oxford/Agama/blob/master/doc/reference.tex
- Barnes (2012), softening as smoothing: https://arxiv.org/abs/1205.2729
- Ludlow et al. (2021), halo-driven spurious heating: https://arxiv.org/abs/2105.03561
- Yurin & Springel (2014), GALIC: https://arxiv.org/abs/1402.1623
- Hernquist (1990), bulge profile: https://adsabs.harvard.edu/pdf/1990ApJ...356..359H
- Toomre (1964), local disk stability: https://articles.adsabs.harvard.edu/pdf/1964ApJ...139.1217T

References support methods and possible mechanisms; their results are not
substitutes for measurements on these particular candidates.
