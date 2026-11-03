"""Run from project root: python -m experiments.compare_solvers"""
import json
from pathlib import Path
from nbody.galaxies.generate import galaxy
from nbody.encounters import parabolic_pair
from nbody.benchmarks import compare

root = Path(__file__).resolve().parents[1]
config = json.loads((root/'presets/mice.json').read_text())
config['particles'] = {'disk':300, 'bulge':60, 'halo':300}
a = galaxy(config, 'preview')
b = galaxy({**config,'seed':config['seed']+1}, 'preview')
p = parabolic_pair(a,b,1.,.25,config['encounter']['orientations'])
print('Force benchmark for',len(p.m),'live particles (single snapshot)')
print('Solver | Seconds | Median relative error % | 95% error % | Max absolute acceleration error')
for row in compare(p.x,p.m):
    print(' | '.join(map(str,row)))
