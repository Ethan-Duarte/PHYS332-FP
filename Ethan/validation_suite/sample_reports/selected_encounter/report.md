# inspect: INCOMPLETE

No automatic certification of physical or observational validity. Passing applies only to named tests and declared tolerances.

| Check | Status | Value | Limit |
|---|---|---|---|
| Array shapes, finite values, positive mass and increasing times | PASS |  |  |
| max_com_residual_over_initial_rms_radius | MEASURED | 1.8712037416664022e-05 | None |
| max_momentum_change_over_initial_scalar_internal_momentum | MEASURED | 2.9596985971283536e-05 | None |
| max_L_vector_change_over_initial_L_norm | MEASURED | 3.2905642238078e-05 | None |
| max_L_vector_change_over_initial_scalar_L | MEASURED | 7.00897582567849e-06 | None |
| max_sampled_energy_relative_change | MEASURED | 0.00021088869923709935 | None |
| Isolated progenitor stability | MISSING |  |  |

## Notes

- Energy is calculated with direct pair sums at saved frames only; it does not bound all-step drift.
- Fitted planes remove rigid tilt, not all warping. Radial cell residuals separate some coherent geometry but are not exact heating estimators.
- Cell residuals cover only cells with >=10 particles; inspect cell_mass_coverage.
- Initial_R_le_2R50 is a fixed set of initial particle IDs; annuli use instantaneous radii and fixed initial radial edges.
- A1/A2 and rotation sign are coordinate conventions; mean_rotation_abs avoids arbitrary eigenvector signs.
- Structural changes in encounters may be physical. No stability gate is applied to an encounter.
- Component labels: 0=disk, 1=bulge, 2=halo. Float32 replay precision limits very small residuals.
