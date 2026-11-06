"""Run same initial Mice-style particles with each force algorithm."""
import argparse
import json
from pathlib import Path
from nbody.galaxies.generate import galaxy
from nbody.encounters import parabolic_pair
from nbody.solver_study import compare_trajectories

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--particles-per-galaxy',type=int,default=200)
    ap.add_argument('--steps',type=int,default=30)
    ap.add_argument('--theta',type=float,default=.5)
    ap.add_argument('--out',default='results/solver_study')
    a=ap.parse_args()
    c=json.loads((Path(__file__).resolve().parents[1]/'presets/mice.json').read_text())
    c['particles']={k:max(10,round(a.particles_per_galaxy*w)) for k,w in [('disk',.58),('bulge',.08),('halo',.34)]}
    one=galaxy(c,'preview');two=galaxy({**c,'seed':c['seed']+1},'preview')
    p=parabolic_pair(one,two,1.,.25,c['encounter']['orientations'])
    results,plot=compare_trajectories(p,a.out,dt=.006,steps=a.steps,theta=a.theta)
    for name,row in results.items():
        print(name,'seconds=',round(row['seconds'],3),'p95 force=',round(100*row['force_p95_error'],3),'%')
    print('Plot:',plot)

if __name__=='__main__':main()
