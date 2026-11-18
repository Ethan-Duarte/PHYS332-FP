"""Isolated dynamical screening from a saved particle realization.

Usage:
python -m experiments.validate_galaxy --snapshot results/galaxies/agama_smoke_900.npz
python -m experiments.validate_galaxy --plummer --steps 100 --dt 0.005
"""
import argparse
import json
from datetime import datetime
from pathlib import Path
from nbody.io import load
from nbody.galaxies.sources.native import plummer
from nbody.validation import validate_initial_state

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--snapshot', type=Path)
parser.add_argument('--plummer', action='store_true')
parser.add_argument('--dt', type=float, default=0.003)
parser.add_argument('--steps', type=int, default=80)
parser.add_argument('--eps', type=float, default=0.0125)
parser.add_argument('--save-every', type=int, default=10)
parser.add_argument('--solver', choices=['Fast direct','Barnes–Hut'], default='Fast direct')
args = parser.parse_args()
if args.plummer:
    p = plummer(300, mass=1.0, a=1.0, seed=42)
    backend = 'analytic Plummer DF'
else:
    p, meta = load(args.snapshot)
    backend = meta.get('backend', 'saved snapshot')
out = root/'results'/'validation'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
result = validate_initial_state(p, out, args.dt, args.steps, args.eps, args.save_every, backend, args.solver)
print(json.dumps({k:v for k,v in result.items() if k not in ('base','half_dt')}, indent=2))
print('Reports and trajectories:', out)
