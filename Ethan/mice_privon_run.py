"""Privon-inspired orientation experiment using the existing encounter runner."""
import runpy
import sys
from pathlib import Path

import nbody.encounters as encounters

original_pair = encounters.parabolic_pair


def privon_pair(galaxy_a, galaxy_b, separation, pericenter, orientations):
    # The runner already converts from the pericenter reference frame
    # to its simulation coordinates. Preserve that conversion.
    reference_angle = orientations[0][1] + 30.0

    # Modify the same list so the runner records the actual angles used.
    orientations[:] = [
        [15.0, -35.0 + reference_angle],  # NGC 4676A
        [25.0, 200.0 + reference_angle],  # NGC 4676B
    ]

    print("Privon-inspired disk orientations relative to pericenter:")
    print("  A: inclination=15 deg, node=-35 deg")
    print("  B: inclination=25 deg, node=200 deg")
    print("Simulation-coordinate orientations:", orientations)

    return original_pair(
        galaxy_a, galaxy_b, separation, pericenter, orientations
    )


encounters.parabolic_pair = privon_pair

runner = Path(__file__).resolve().with_name("mice_candidate_run.py")
if not runner.is_file():
    raise FileNotFoundError(f"Encounter runner missing: {runner}")

sys.argv[0] = str(runner)
runpy.run_path(str(runner), run_name="__main__")