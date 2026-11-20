# inspect: INCOMPLETE

No automatic certification of physical or observational validity. Passing applies only to named tests and declared tolerances.

| Check | Status | Value | Limit |
|---|---|---|---|
| Array shapes, finite values, positive mass and increasing times | PASS |  |  |
| max_com_residual_over_initial_rms_radius | MEASURED | 2.084778730233072e-05 | 0.0001 |
| max_momentum_change_over_initial_scalar_internal_momentum | MEASURED | 4.375677433085716e-05 | 0.0001 |
| max_L_vector_change_over_initial_L_norm | MEASURED | 0.0001351189891426141 | 0.001 |
| max_L_vector_change_over_initial_scalar_L | MEASURED | 9.206329314938604e-06 | None |
| Energy conservation | MISSING |  |  |
| g0 all_disk r50 evolution | MEASURED | 0.07800376658538884 | 0.05 |
| g0 all_disk r90 evolution | MEASURED | 0.04968574798288086 | None |
| g0 all_disk height evolution | MEASURED | 0.43034043911995873 | 0.05 |
| g0 all_disk sigma_z evolution | MEASURED | 0.18629504466794433 | 0.05 |
| g0 all_disk sigma_R evolution | MEASURED | 0.5085704261142543 | None |
| g0 initial_R_le_2R50 r50 evolution | MEASURED | 0.06253344304743355 | 0.05 |
| g0 initial_R_le_2R50 r90 evolution | MEASURED | 0.1485814977342775 | None |
| g0 initial_R_le_2R50 height evolution | MEASURED | 0.3818947655460325 | 0.05 |
| g0 initial_R_le_2R50 sigma_z evolution | MEASURED | 0.18195070423805015 | 0.05 |
| g0 initial_R_le_2R50 sigma_R evolution | MEASURED | 0.48366351249961514 | None |

## Notes

- Energy is calculated with direct pair sums at saved frames only; it does not bound all-step drift.
- Fitted planes remove rigid tilt, not all warping. Radial cell residuals separate some coherent geometry but are not exact heating estimators.
- Cell residuals cover only cells with >=10 particles; inspect cell_mass_coverage.
- Initial_R_le_2R50 is a fixed set of initial particle IDs; annuli use instantaneous radii and fixed initial radial edges.
- A1/A2 and rotation sign are coordinate conventions; mean_rotation_abs avoids arbitrary eigenvector signs.
- Structural changes in encounters may be physical. No stability gate is applied to an encounter.
- Component labels: 0=disk, 1=bulge, 2=halo. Float32 replay precision limits very small residuals.
