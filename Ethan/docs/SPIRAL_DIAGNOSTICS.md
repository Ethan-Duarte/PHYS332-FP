# V1.4.1 — Interpreting transient spiral structure

**Purpose:** distinguish an interesting evolving density pattern from an equilibrium/instability claim. These diagnostics do **not** establish a physical spiral mode or validate the native galaxy generator.

## Fourier measurements (the mathematical core)

At a saved time, recenter the stellar disk (`component == 0`) and put its stars in radial annuli. For annulus `k`, define

\[
C_m(R_k,t)=\frac{\sum_{j\in k}m_j e^{im\phi_j}}{\sum_{j\in k}m_j},\qquad
A_m=|C_m|,\quad m=1,2,3,4.
\]

The weights are actual particle masses; the underlying halo and bulge particles influence the *evolution*, but not these stellar Fourier sums. `spiral_modes.csv` reports the mass-weighted mean of the **absolute annular** strengths, rather than the global complex sum. That matters because winding spirals have different phases at different radii; a global sum can cancel a real pattern.

- `m=1`: lopsidedness; `m=2`: two-fold pattern (bar **or** two-armed spiral); `m=3,4`: higher multiplicity.
- `phase2=arg(C2)/2`, defined modulo π; a nearly constant phase with radius is *suggestive* of a bar, while radial phase variation is compatible with a spiral. Both require more than a single plot to interpret.
- For independently sampled, azimuthally uniform stars, the **RMS** noise in `A_m` is approximately `sqrt(sum(m_i²))/sum(m_i)=1/sqrt(N_eff)` per annulus. This is a rough reference, not a p-value. It is especially important at low N.
- Fixed annular edges are `[0.5,1,1.5,2,2.5,3,4,5] R_d`. Annuli with fewer than 20 stellar particles are omitted. At low resolution, the measured amplitudes are noisy and some bins may be missing.

## Disk heating and the approximate Toomre Q

We measure stellar radial and azimuthal velocity dispersions after subtracting the disk's bulk velocity, and vertical RMS thickness after recentering it. The structural plot reports changes relative to their initial values.

The locally axisymmetric stellar Toomre parameter is

\[
Q_\star(R)=\frac{\sigma_R(R)\,\kappa(R)}{3.36\,G\,\Sigma_\star(R)}.
\]

We estimate `Sigma_star` from the mass of disk particles per annular area. We estimate the circular speed from the azimuthally averaged radial gravitational acceleration of **all live particles**, using the **same softened Newtonian interaction** as the simulation. We then obtain

\[
\Omega^2=\frac{v_c^2}{R^2},\qquad
\kappa^2=\frac{1}{R^3}\frac{d}{dR}\left(R^4\Omega^2\right).
\]

This is an **approximate diagnostic**. `Q` is meaningful for a thin, nearly axisymmetric, roughly equilibrated stellar disk; those assumptions are poor for the rapidly relaxing native preview. The finite radial sampling, softening, asymmetric mass distribution, and particle noise also affect `kappa` and `Q`. Nonpositive estimated `kappa²` is reported as NaN, not forced to a physical number. Do not label `Q < 1` in a preview run proof of an unstable astrophysical disk.

## Outputs

Running an *isolated disk–bulge–halo* experiment automatically adds:

- `spiral_modes.csv`: time, mean annular `A1`–`A4`, annular noise RMS, radial dispersion, disk thickness.
- `spiral_annuli.csv`: per-annulus counts, masses, mode strengths, complex `m=2` phase, dispersions, approximate `kappa`, `Q_star`, surface density.
- `spiral_diagnostics.png`: four plots (modes, heating, radial A2, approximate Q).
- `spiral_phase.png`: radial mode-amplitude and phase profiles at initial, peak and final times.
- `spiral_report.json`: peak A2 time/strength and changes in velocity dispersion and disk thickness, with scientific caveats.

The GUI's **Analyze disk spirals** tab can regenerate these from a saved *isolated disk* trajectory without rerunning the gravitational simulation. The 3D viewer remains unchanged. For isolated Plummer spheres, the disk diagnostic is intentionally omitted.

## Suggested experimental comparison

Run a native disk in isolation with a single fixed configuration, then change (a) the random seed, (b) stellar and halo particle counts independently, and (c) timestep and softening. Plot `A2(t)` beside `sigma_R(t)` and thickness. Note whether the peak's timing and amplitude are larger than the expected finite-N noise and stable against resolution. A transient `m=2` enhancement alone is not enough to distinguish spiral arms from bar formation or disequilibrium.

For V1.5 the key goal remains a genuinely self-consistent disk DF, preferably generated using AGAMA and then tested in isolation. All measurements remain exploratory until that benchmark is passed.

## References

- Toomre (1964), stellar disk Q and axisymmetric stability: https://articles.adsabs.harvard.edu/pdf/1964ApJ...139.1217T
- Fujii et al. (2011), time-dependent spiral structure, finite-N sensitivity, disk heating: https://doi.org/10.1088/0004-637X/730/2/109
- Sellwood & Masters (2022), review of spiral mechanisms: https://doi.org/10.1146/annurev-astro-052920-104505
- Barnes (2004), our Mice-style reference: https://doi.org/10.1111/j.1365-2966.2004.07725.x
