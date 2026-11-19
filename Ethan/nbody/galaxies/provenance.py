"""Separate construction-force compatibility from evidence of equilibrium."""
import numpy as np


def check_evolution_softening(metadata, epsilon):
    """Reject explicit mismatches; return limitations for legacy/unknown ICs.

    A legacy top-level epsilon only records a declared value. It is not proof
    the original DF was constructed in a softened field or passed isolation.
    """
    epsilon = float(epsilon)
    if not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError('Evolution epsilon must be finite and positive')
    force = metadata.get('construction_force_model') or {}
    kernel = force.get('kernel')
    if kernel and kernel not in ('plummer','newtonian'):
        raise ValueError(f'Unsupported construction kernel: {kernel}')
    expected = force.get('epsilon') if kernel=='plummer' else metadata.get('epsilon')
    if expected is not None:
        expected = float(expected)
        if not np.isfinite(expected) or expected <= 0:
            raise ValueError('Invalid source softening metadata')
        if not np.isclose(epsilon,expected,rtol=0,atol=1e-15):
            raise ValueError(f'Snapshot declares epsilon={expected:g}; requested {epsilon:g}. Regenerate or explicitly construct a new candidate for a different force law.')
    if kernel=='newtonian':
        return 'Newtonian construction evolved with Plummer softening; equilibrium mismatch must be measured.'
    if kernel=='plummer':
        return 'Construction softening matches; live equilibrium and resolution convergence remain separate tests.'
    return 'Construction force law is unrecorded; an epsilon annotation does not establish equilibrium compatibility.'
