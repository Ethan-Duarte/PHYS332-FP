# inspect: MEASURED

No automatic certification of physical or observational validity. Passing applies only to named tests and declared tolerances.

| Check | Status | Value | Limit |
|---|---|---|---|
| Array shapes, finite values, positive mass and increasing times | PASS |  |  |
| max_com_residual_over_initial_rms_radius | MEASURED | 2.529790717500557e-16 | None |
| max_momentum_change_over_initial_scalar_internal_momentum | MEASURED | 1.4379520465109815e-16 | None |
| max_L_vector_change_over_initial_L_norm | MEASURED | 5.531946339487211e-15 | None |
| max_L_vector_change_over_initial_scalar_L | MEASURED | 3.769190433497327e-16 | None |
| max_sampled_energy_relative_change | MEASURED | 1.8904619634839422e-06 | None |
| g0 all_disk r50 evolution | MEASURED | 0.0026371224429697173 | None |
| g0 all_disk r90 evolution | MEASURED | 0.0007392043781750962 | None |
| g0 all_disk height evolution | MEASURED | 0.003784734318499562 | None |
| g0 all_disk sigma_z evolution | MEASURED | 0.0011163748486874925 | None |
| g0 all_disk sigma_R evolution | MEASURED | 0.0024064199289542463 | None |
| g0 initial_R_le_2R50 r50 evolution | MEASURED | 0.000183026367037144 | None |
| g0 initial_R_le_2R50 r90 evolution | MEASURED | 0.003326475878229323 | None |
| g0 initial_R_le_2R50 height evolution | MEASURED | 0.004380224386150845 | None |
| g0 initial_R_le_2R50 sigma_z evolution | MEASURED | 0.000783055356514474 | None |
| g0 initial_R_le_2R50 sigma_R evolution | MEASURED | 0.0011509619637377771 | None |

## Notes

- Energy is calculated with direct pair sums at saved frames only; it does not bound all-step drift.
- Energy, if requested, assumes fixed-mass self-gravity only; no external, cosmological or thermodynamic energy terms.
- Fitted planes remove rigid tilt, not all warping. Radial cell residuals separate some coherent geometry but are not exact heating estimators.
- Cell residuals cover only cells with >=10 particles; inspect cell_mass_coverage.
- Initial_R_le_2R50 is a fixed set of initial particle IDs; annuli use instantaneous radii and fixed initial radial edges.
- A1/A2 and rotation sign are coordinate conventions; mean_rotation_abs avoids arbitrary eigenvector signs.
- Structural changes in encounters may be physical. No stability gate is applied to an encounter.
- Component labels: 0=disk, 1=bulge, 2=halo. Float32 replay precision limits very small residuals.
