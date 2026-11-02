"""Fast pipeline demo, not a validated model of NGC 4676.

Run: python -m experiments.mice_preview
"""
import json
from pathlib import Path
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from nbody.galaxies.generate import galaxy
from nbody.encounters import parabolic_pair
from nbody.integrate import evolve
from nbody.gravity.direct import acceleration
from nbody.analysis import core_separation
from nbody.io import save

root = Path(__file__).resolve().parents[1]
config = json.loads((root/'presets/mice.json').read_text())
g1 = galaxy(config, 'preview')
config2 = dict(config, seed=config['seed']+1)
g2 = galaxy(config2, 'preview')
e = config['encounter']
p = parabolic_pair(g1, g2, e['separation'], e['pericenter'], e['orientations'])
force = partial(acceleration, eps=0.01)
out = root/'results'
out.mkdir(exist_ok=True)
snaps = list(evolve(p, force, dt=0.01, steps=40, save_every=10))
for t, state in snaps:
    print(f't={t:.2f}: bulge core separation = {core_separation(state):.4f}')
last = snaps[-1][1]
save(out/'mice_preview.npz', last, {'config': config, 'backend': 'native-preview', 'force': 'direct', 'eps': 0.01, 'dt': 0.01, 'steps': 40, 'time': snaps[-1][0]})
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, state, title in [(axes[0], p, 'Initial'), (axes[1], last, 'Final')]:
    for gid in (0, 1):
        mask = (state.galaxy == gid) & (state.component == 0)
        ax.scatter(state.x[mask, 0], state.x[mask, 1], s=2, label=f'Galaxy {gid+1}')
    ax.set(xlim=(-0.9, 0.9), ylim=(-0.9, 0.9), xlabel='x', ylabel='y', title=title)
    ax.set_aspect('equal')
axes[0].legend(loc='upper right', fontsize=7)
fig.tight_layout()
fig.savefig(out/'mice_preview.png', dpi=160)
print('Saved results/mice_preview.npz and results/mice_preview.png')
