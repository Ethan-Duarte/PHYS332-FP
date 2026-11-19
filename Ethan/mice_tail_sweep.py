"""Five controlled 16384-particle encounter trials using installed WSL Zeno."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from nbody.io import load

JOBS=[('control_s1p25_rp0p25','1p25',.25),
      ('size_s1p50_rp0p25','1p50',.25),
      ('size_s1p75_rp0p25','1p75',.25),
      ('orbit_s1p25_rp0p20','1p25',.20),
      ('orbit_s1p50_rp0p20','1p50',.20)]


def summarize(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    after=[.125,.25,.375,.5,.625]
    fig,axes=plt.subplots(len(JOBS),len(after),figsize=(16,13),facecolor='#08090e')
    az=np.radians(120);B=np.array([[-np.sin(az),np.cos(az),0],[0,0,1.]])
    summaries=[]
    for row,(name,size,rp) in enumerate(JOBS):
        with np.load(root/name/'run_trajectory.npz') as data:
            t,x,c,g=data['t'],data['x'],data['component'],data['galaxy']
        centers=np.stack([np.median(x[:,(c==1)&(g==k)],axis=1) for k in (0,1)],axis=1)
        separation=np.linalg.norm(centers[:,1]-centers[:,0],axis=1)
        peri=t[separation.argmin()]
        summaries.append(dict(name=name,scale=float(size.replace('p','.')),nominal_pericenter=rp,
                              closest_saved_time=float(peri),closest_saved_separation=float(separation.min())))
        for col,a in enumerate(after):
            requested=peri+a
            if requested>t[-1]+1e-10:raise ValueError(f'{name}: requested post-passage time exceeds saved run')
            fi=int(np.argmin(abs(t-requested)))
            core=centers[fi]@B.T;vec=core[0]-core[1];dist=np.linalg.norm(vec)
            if dist<1e-12:raise ValueError('Degenerate projected core separation')
            # One rotation, translation, and uniform scale for the whole system.
            ang=-np.arctan2(vec[1],vec[0]);R=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])
            q=((x[fi]@B.T)-core.mean(0))@R.T/dist
            ax=axes[row,col];ax.set_facecolor('#08090e')
            for k,color in [(0,'#62d5e8'),(1,'#f465b5')]:
                mask=(c==0)&(g==k);ax.scatter(q[mask,0],q[mask,1],s=.35,alpha=.55,color=color,linewidths=0)
            ax.set(xlim=(-4.,4.),ylim=(-3.,3.));ax.set_aspect('equal')
            ax.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)
            for spine in ax.spines.values():spine.set_visible(False)
            ax.set_title(f't={t[fi]:.3f}; after={t[fi]-peri:.3f}',color='#bbbccc',fontsize=9)
            if col==0:ax.set_ylabel(f'Size ×{size.replace("p",".")}\nnominal rp={rp:g}',color='white',fontsize=11)
    fig.suptitle('Tail search • azimuth 120°, elevation 0° • 16,384 particles per pair',color='white',fontsize=16)
    fig.text(.5,.02,'Each entire projection is aligned and uniformly scaled to core separation = 1. Shared plot limits; extended debris may lie outside the window.',ha='center',color='#bbbccc',fontsize=9)
    fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig(root/'tail_comparison.png',dpi=150,facecolor=fig.get_facecolor());plt.close(fig)
    (root/'encounter_summary.json').write_text(json.dumps(summaries,indent=2))
    print(f'Saved comparison: {root / "tail_comparison.png"}',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--zeno-root',default='/home/whisky/nbody-tools/zeno-reference')
    ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args()
    source=args.source.resolve();p,meta=load(source)
    if len(p.m)!=65536 or meta.get('model_id')!='zeno_reference_bdh_v1':ap.error('Use the full original jeans_candidate.npz')
    root=args.out.resolve();root.mkdir(parents=True,exist_ok=False)
    scripts=Path(__file__).resolve().parent
    subprocess.run([sys.executable,str(scripts/'prepare_mice_size_screen.py'),str(source),
                    '--out',str(root/'candidates'),'--disk','3072','--bulge','1024','--halo','4096',
                    '--scales','1.25','1.5','1.75'],check=True)
    plan=dict(particles_per_pair=16384,disk_per_galaxy=3072,bulge_per_galaxy=1024,halo_per_galaxy=4096,
              source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),duration=1.25,dt=1/512,theta=.7,
              note='Whole-galaxy homology; no disk-only stretching. Same subsample across trials. No equilibrium or resolution-convergence certification.',
              jobs=[dict(name=n,candidate=f'candidate_s{s}.npz',nominal_pericenter=r) for n,s,r in JOBS])
    (root/'sweep_plan.json').write_text(json.dumps(plan,indent=2))
    if args.prepare_only:
        print('Prepared candidates and five-run plan; no evolution launched.',flush=True);return
    for i,(name,size,rp) in enumerate(JOBS,1):
        print(f'\nTRIAL {i}/{len(JOBS)}: {name}',flush=True)
        subprocess.run([sys.executable,str(scripts/'zeno_mice_run.py'),str(root/'candidates'/f'candidate_s{size}.npz'),
                        '--out',str(root/name),'--zeno-root',args.zeno_root,'--duration','1.25',
                        '--pericenter',str(rp)],check=True)
    summarize(root)
    print('All five trials complete. Upload tail_comparison.png and encounter_summary.json first.',flush=True)

if __name__=='__main__':main()
