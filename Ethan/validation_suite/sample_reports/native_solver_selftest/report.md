# selftest: PASS_FOR_SPECIFIED_CHECKS

No automatic certification of physical or observational validity. Passing applies only to named tests and declared tolerances.

| Check | Status | Value | Limit |
|---|---|---|---|
| Direct force vs independent vectorized sum | PASS | 4.600591911373043e-16 | 1e-10 |
| Softened circular orbit: endpoint RMS / separation | PASS | 0.0007476390367758474 | 0.01 |
| Leapfrog second-order convergence on fixture | PASS | 3.9973420871096206 |  |

## Notes

- Acceptance bounds are declared for this numerical fixture, not galaxy or observational tolerances.
- This tests the native project direct solver and leapfrog. It does not test Zeno or GADGET binaries.
- Two orbital periods; independent DOP853 rtol=1e-12. Passing this fixture is not full galaxy validation.
