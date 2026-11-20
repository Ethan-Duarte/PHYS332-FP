"""Reproduce prepared final state from the exact prior uploaded parent states.
This is the preparation audit; users can run the supplied prepared NPZ directly.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rotation(i,node):
 i,node=np.deg2rad([i,node]);c,s=np.cos(node),np.sin(node)
 return np.array([[c,-s,0],[s,c,0],[0,0,1]])@np.array([[1,0,0],[0,np.cos(i),-np.sin(i)],[0,np.sin(i),np.cos(i)]])

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
 args=ap.parse_args()
 if args.out.exists():raise ValueError('Output exists')
 with np.load(args.baseline,allow_pickle=False) as d:
  x=d['x'][0].astype(float);v=d['v'][0].astype(float);m=d['m'].copy();c=d['component'].copy();meta=json.loads(d['metadata'].item())
 with np.load(args.parent,allow_pickle=False) as d:
  X=d['x'][0].astype(float);V=d['v'][0].astype(float);M=d['m'].copy();C=d['component'].copy();G=d['galaxy'].copy();pm=json.loads(d['metadata'].item())
 expected=meta['source']['size_screen']['source_sha256']
 if pm['source_sha256']!=expected:raise ValueError('Parent candidate hash mismatch')
 if [int(np.sum(c==j)) for j in range(3)]!=[3072,1024,32768]:raise ValueError('Wrong baseline component sampling')
 s=meta['source']['size_screen']['scale'];R=rotation(*meta['orientations'][0]);P=rotation(*pm['orientations'][0])
 disk=c==0;pdisk=(G==0)&(C==0)
 dc=np.average(x[disk],axis=0,weights=m[disk]);dv=np.average(v[disk],axis=0,weights=m[disk])
 pc=np.average(X[pdisk],axis=0,weights=M[pdisk]);pv=np.average(V[pdisk],axis=0,weights=M[pdisk])
 mapped_x=(((X-pc)@P)*s-np.array(meta['removed_disk_position']))@R.T+dc
 mapped_v=(((V-pv)@P)/np.sqrt(s)-np.array(meta['removed_disk_velocity']))@R.T+dv
 rng=np.random.default_rng(20261120);xx=[];vv=[];mm=[];cc=[];reports=[]
 for component,count in [(0,12288),(1,4096)]:
  old=c==component;pool=np.flatnonzero((G==0)&(C==component));tree=cKDTree(mapped_x[pool]);dist,indices=tree.query(x[old])
  if len(np.unique(indices))!=old.sum() or dist.max()>2e-6:raise ValueError('Original sample cannot be uniquely mapped into parent')
  velerr=float(np.max(np.linalg.norm(mapped_v[pool[indices]]-v[old],axis=1)))
  if velerr>2e-5:raise ValueError('Parent velocity mapping inconsistent')
  others=np.setdiff1d(np.arange(len(pool)),indices);added=rng.choice(others,size=count-int(old.sum()),replace=False)
  # Preserve the original stellar coordinates exactly; add distinct parent particles.
  cx=np.vstack([x[old],mapped_x[pool[added]]]);cv=np.vstack([v[old],mapped_v[pool[added]]])
  xx.append(cx);vv.append(cv);mm.append(np.full(count,m[old].sum()/count));cc.append(np.full(count,component))
  reports.append(dict(component=component,old_count=int(old.sum()),new_count=count,retained_original_particles=int(old.sum()),new_distinct_parent_particles=len(added),max_parent_position_match_error=float(dist.max()),max_parent_velocity_match_error=velerr))
 halo=c==2;xx.append(x[halo]);vv.append(v[halo]);mm.append(m[halo]);cc.append(c[halo])
 xx=np.concatenate(xx);vv=np.concatenate(vv);mm=np.concatenate(mm);cc=np.concatenate(cc);gg=np.zeros(len(mm),int)
 assert np.array_equal(xx[cc==2],x[halo]) and np.array_equal(vv[cc==2],v[halo]) and np.array_equal(mm[cc==2],m[halo])
 for k in range(3):assert np.isclose(mm[cc==k].sum(),m[c==k].sum(),rtol=1e-13)
 report=dict(baseline_sha256=sha(args.baseline),parent_trajectory_sha256=sha(args.parent),parent_candidate_sha256=expected,seed=20261120,components=reports,halo_state_unchanged_before_common_recenter=True,particle_counts=[int(np.sum(cc==k)) for k in range(3)],eps=float(meta['eps']),dt=1/512,theta=float(meta['theta']),duration=1.25,save_every=8,
  interpretation='Single matched-parent stellar sampling comparison. Not a seed ensemble, exact DF, or automatic equilibrium certification.')
 # Preserve earlier metadata as history, explicitly update current sampling.
 meta=dict(meta,final_stellar_resolution_control=report,dt=1/512,duration=1.25,save_every=8,
  sampling='Nested distinct parent stellar sampling; original halo positions, velocities and masses retained before common recentering.',
  current_component_counts={'disk':12288,'bulge':4096,'halo':32768},validation_mode='isolate',
  note='Final planned isolated validation run; no orbital parameter change. Historical source counts describe ancestors, not this current snapshot.')
 np.savez_compressed(args.out,t=np.array([0.]),x=xx[None],v=vv[None],m=mm,component=cc,galaxy=gg,metadata=json.dumps(meta))
 args.out.with_suffix('.preparation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

if __name__=='__main__':main()
