"""Select generator explicitly. Do not mistake preview kinematics for equilibrium."""
from .sources import native


def galaxy(config, backend='preview'):
    if backend == 'agama':
        from .sources.agama import generate
        return generate(config)
    if backend == 'preview':
        return native.mice_preview(config)
    raise ValueError(f'Unknown galaxy backend: {backend}')
