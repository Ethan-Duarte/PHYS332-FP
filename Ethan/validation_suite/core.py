"""Independent measurements; no observational claim is inferred from conservation."""
import csv,hashlib,json
from dataclasses import dataclass
from pathlib import Path
import numpy as np

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csv_write(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def fraction(a,b):return float(a/b) if abs(b)>1e-30 else None
def quantile(x,m,q):
    ids=np.argsort(x);return float(x[ids[min(np.searchsorted(np.cumsum(m[ids]),q*m.sum()),len(ids)-1)]])

@dataclass
class Run:
    path:Path
    t:np.ndarray
    x:np.ndarray
    v:np.ndarray
    m:np.ndarray
    c:np.ndarray
    g:np.ndarray
    meta:dict
    labels:bool

def load(path,snapshot=False):
    path=Path(path)
    with np.load(path,allow_pickle=False) as d:
        x=np.asarray(d['x'],dtype=float);v=np.asarray(d['v'],dtype=float);m=np.asarray(d['m'],dtype=float)
        if x.ndim==2 and snapshot:x=x[None];v=v[None]
        t=np.asarray(d['t'],dtype=float) if 't' in d else (np.array([0.]) if snapshot else None)
        if t is None:raise ValueError('Trajectory requires t; use initial for a snapshot')
        labels='component' in d and 'galaxy' in d
        c=np.asarray(d['component']) if 'component' in d else np.full(len(m),-1)
        g=np.asarray(d['galaxy']) if 'galaxy' in d else np.zeros(len(m),dtype=int)
        meta=json.loads(str(d['metadata'].item())) if 'metadata' in d else {}
    if len(m)==0 or m.ndim!=1 or x.shape!=(len(t),len(m),3) or v.shape!=x.shape:raise ValueError('Expected x,v [frames,N,3]; m [N]')
    if t.ndim!=1 or len(t)==0 or np.any(np.diff(t)<=0):raise ValueError('Times must be strictly increasing')
    if c.shape!=m.shape or g.shape!=m.shape:raise ValueError('Labels must have shape [N]')
    if not all(np.isfinite(a).all() for a in [t,x,v,m,c,g]):raise ValueError('Nonfinite data')
    if np.any(m<=0):raise ValueError('Masses must be positive')
    if not np.equal(c,np.floor(c)).all() or not np.equal(g,np.floor(g)).all():raise ValueError('Labels must be integers')
    return Run(path,t,x,v,m,c.astype(int),g.astype(int),meta,labels)

def law(run,epsilon=None,G=None,kernel=None):
    model=run.meta.get('construction_force_model',{})
    src=run.meta.get('source',{});model=src.get('construction_force_model',model)
    recorded_eps=run.meta.get('eps',run.meta.get('epsilon',model.get('epsilon')))
    recorded_G=run.meta.get('G',model.get('G'))
    # Metadata is recorded evidence, not inferred from which package generated it.
    if epsilon is not None and recorded_eps is not None and not np.isclose(epsilon,recorded_eps,rtol=1e-10,atol=1e-14):raise ValueError('Explicit epsilon disagrees with metadata')
    if G is not None and recorded_G is not None and not np.isclose(G,recorded_G):raise ValueError('Explicit G disagrees with metadata')
    eps=recorded_eps if epsilon is None else epsilon;gg=recorded_G if G is None else G
    kk=kernel or run.meta.get('evolution_kernel') or model.get('kernel')
    if eps is None or gg is None or kk is None:raise ValueError('Need recorded force law or explicit --epsilon, --G and --kernel')
    if kk.lower()!='plummer':raise ValueError('Only Plummer kernel implemented; no silent substitution')
    if eps<0 or gg<=0:raise ValueError('Require epsilon >=0 and G>0')
    return float(eps),float(gg),'plummer'

def reference_force(x,m,ids,epsilon,G):
    result=[]
    for i in ids:
        delta=x-x[i];r2=np.sum(delta*delta,axis=1)+epsilon**2
        r2[i]=np.inf
        result.append(np.sum((G*m/r2**1.5)[:,None]*delta,axis=0))
    return np.asarray(result)

try:
    from numba import njit,prange
    @njit(parallel=True)
    def _potential(x,m,eps,G):
        rows=np.zeros(len(m))
        for i in prange(len(m)):
            v=0.
            for j in range(i):
                r2=eps*eps
                for k in range(3):r2+=(x[i,k]-x[j,k])**2
                v-=G*m[i]*m[j]/np.sqrt(r2)
            rows[i]=v
        return rows.sum()
except ImportError:
    def _potential(x,m,eps,G):
        if len(m)>3000:raise ValueError('Install numba for energy on N>3000; inspect without energy still works')
        return sum(float(np.sum(-G*m[i]*m[:i]/np.sqrt(np.sum((x[:i]-x[i])**2,axis=1)+eps**2))) for i in range(len(m)))

def plane(x,v,m):
    x=x-np.average(x,axis=0,weights=m);v=v-np.average(v,axis=0,weights=m)
    values,axes=np.linalg.eigh((x*m[:,None]).T@x/m.sum())
    normal=axes[:,0];z=x@normal;vz=v@normal;q=x@axes[:,1:];R=np.linalg.norm(q,axis=1)
    phi=np.arctan2(q[:,1],q[:,0]);vr=np.sum(q*(v@axes[:,1:]),axis=1)/np.maximum(R,1e-30)
    tang=(q[:,0]*(v@axes[:,1:])[:,1]-q[:,1]*(v@axes[:,1:])[:,0])/np.maximum(R,1e-30)
    return R,z,vz,phi,vr,tang,normal

def structure(x,v,m):
    R,z,vz,phi,vr,vt,n=plane(x,v,m)
    meanvz=float(np.average(vz,weights=m));meanvr=float(np.average(vr,weights=m))
    return dict(r50=quantile(R,m,.5),r90=quantile(R,m,.9),
        height=float(np.sqrt(np.average(z*z,weights=m))),
        sigma_z=float(np.sqrt(np.average((vz-meanvz)**2,weights=m))),
        sigma_R=float(np.sqrt(np.average((vr-meanvr)**2,weights=m))),
        mean_rotation_abs=float(abs(np.average(vt,weights=m))),
        A1=float(abs(np.average(np.exp(1j*phi),weights=m))),
        A2=float(abs(np.average(np.exp(2j*phi),weights=m))),
        tilt_deg=float(np.degrees(np.arccos(np.clip(abs(n[2]),0,1)))))

def structure_history(run):
    rows=[];radial=[]
    if not run.labels:return rows,radial
    for gal in np.unique(run.g):
        ids=np.flatnonzero((run.g==gal)&(run.c==0))
        if len(ids)<10:continue
        R0,*_=plane(run.x[0,ids],run.v[0,ids],run.m[ids]);r50=quantile(R0,run.m[ids],.5)
        inner=ids[R0<=2*r50];edges=r50*np.array([0,.5,1,2,4,np.inf])
        for i,t in enumerate(run.t):
            rows.append(dict(frame=i,time=float(t),galaxy=int(gal),aperture='all_disk',**structure(run.x[i,ids],run.v[i,ids],run.m[ids])))
            if len(inner)>=10:rows.append(dict(frame=i,time=float(t),galaxy=int(gal),aperture='initial_R_le_2R50',**structure(run.x[i,inner],run.v[i,inner],run.m[inner])))
            R,z,vz,phi,vr,vt,n=plane(run.x[i,ids],run.v[i,ids],run.m[ids]);mm=run.m[ids]
            for b in range(len(edges)-1):
                mask=(R>=edges[b])&(R<edges[b+1]);count=int(mask.sum())
                if count<10:continue
                w=mm[mask];zz=z[mask];vv=vz[mask];mz=float(np.average(zz,weights=w));mv=float(np.average(vv,weights=w))
                # Within-cell and between-cell vertical second moments; not a full mode decomposition.
                cellmeans=[];cellvars=[];cellmasses=[]
                for sector in range(8):
                    ss=mask&(phi>=-np.pi+sector*np.pi/4)&(phi<-np.pi+(sector+1)*np.pi/4)
                    if ss.sum()<10:continue
                    mass=mm[ss];mean=float(np.average(z[ss],weights=mass));cellmasses.append(float(mass.sum()));cellmeans.append(mean);cellvars.append(float(np.average((z[ss]-mean)**2,weights=mass)))
                within=float(np.sqrt(np.average(cellvars,weights=cellmasses))) if cellmasses else None
                coherent=float(np.sqrt(np.average(np.asarray(cellmeans)**2,weights=cellmasses))) if cellmasses else None
                radial.append(dict(time=float(t),galaxy=int(gal),bin=b,R_lo=float(edges[b]),R_hi=None if np.isinf(edges[b+1]) else float(edges[b+1]),count=count,mass=float(w.sum()),mean_z=mz,mean_vz=mv,
                    height_about_annulus_mean=float(np.sqrt(np.average((zz-mz)**2,weights=w))),sigma_z_about_annulus_mean=float(np.sqrt(np.average((vv-mv)**2,weights=w))),
                    cell_residual_height=within,cell_mean_height_rms=coherent,cell_mass_coverage=float(sum(cellmasses)/w.sum()),
                    mean_rotation_abs=float(abs(np.average(vt[mask],weights=w)))))
    return rows,radial

def conserved(run,energy_samples=0,epsilon=None,G=None,kernel=None):
    m=run.m;x=run.x;v=run.v;M=m.sum()
    com=np.einsum('tnj,n->tj',x,m)/M;P=np.einsum('tnj,n->tj',v,m)
    r=x-com[:,None];u=v-P[:,None]/M;L=np.sum(np.cross(r,u)*m[None,:,None],axis=1)
    R0=np.sqrt(np.sum(m*np.sum(r[0]**2,axis=1))/M)
    pscale=np.sum(m*np.linalg.norm(u[0],axis=1));lscale=np.sum(m*np.linalg.norm(np.cross(r[0],u[0]),axis=1))
    drift=com-com[0]-(run.t-run.t[0])[:,None]*P[0]/M
    summary=dict(max_com_residual_over_initial_rms_radius=fraction(np.linalg.norm(drift,axis=1).max(),R0),
        max_momentum_change_over_initial_scalar_internal_momentum=fraction(np.linalg.norm(P-P[0],axis=1).max(),pscale),
        max_L_vector_change_over_initial_L_norm=fraction(np.linalg.norm(L-L[0],axis=1).max(),np.linalg.norm(L[0])),
        max_L_vector_change_over_initial_scalar_L=fraction(np.linalg.norm(L-L[0],axis=1).max(),lscale))
    rows=[]
    for i,t in enumerate(run.t):rows.append(dict(time=float(t),com_x=float(com[i,0]),com_y=float(com[i,1]),com_z=float(com[i,2]),momentum_change=float(np.linalg.norm(P[i]-P[0])),L_vector_change=float(np.linalg.norm(L[i]-L[0])),com_ballistic_residual=float(np.linalg.norm(drift[i]))))
    energies=[]
    if energy_samples:
        eps,gg,kk=law(run,epsilon,G,kernel)
        for i in np.unique(np.linspace(0,len(run.t)-1,min(energy_samples,len(run.t))).astype(int)):
            K=float(.5*np.sum(m[:,None]*v[i]**2));U=float(_potential(x[i],m,eps,gg))
            energies.append(dict(frame=int(i),time=float(run.t[i]),kinetic=K,potential=U,energy=K+U))
        e0=energies[0]['energy']
        for e in energies:e['relative_change']=fraction(e['energy']-e0,abs(e0))
        good=[abs(e['relative_change']) for e in energies if e['relative_change'] is not None]
        summary['max_sampled_energy_relative_change']=max(good) if good else None
        summary['energy_samples']=len(energies);summary['energy_force_law']=dict(epsilon=eps,G=gg,kernel=kk)
    return summary,rows,energies

def criterion(name,value,policy,key,notes=''):
    limit=policy.get('limits',{}).get(key)
    active=policy.get('reviewed') is True
    status='MISSING' if value is None else ('MEASURED' if limit is None or not active else ('PASS' if abs(value)<=limit else 'FAIL'))
    return dict(name=name,status=status,value=value,limit=limit,policy_reviewed=active,notes=notes)

def camera(az,el):
    a,b=np.radians([az,el]);return np.array([[-np.sin(a),np.cos(a),0],[np.sin(b)*np.cos(a),np.sin(b)*np.sin(a),-np.cos(b)],[-np.cos(b)*np.cos(a),-np.cos(b)*np.sin(a),-np.sin(b)]])

def comparison(a,b,kind):
    if a.labels!=b.labels:raise ValueError('Label availability differs')
    epsa,Ga,_=law(a);epsb,Gb,_=law(b)
    if not np.isclose(epsa,epsb) or not np.isclose(Ga,Gb):raise ValueError('Force laws differ; not a controlled comparison')
    provenance=[]
    if kind in ['timestep','solver']:
        if not np.array_equal(a.m,b.m) or not np.array_equal(a.c,b.c) or not np.array_equal(a.g,b.g):raise ValueError('Particle masses, order or labels differ')
        # Permit the uniform bulk subtraction used by isolation controls, not arbitrary rotations or scaling.
        ax=a.x[0]-np.average(a.x[0],axis=0,weights=a.m);bx=b.x[0]-np.average(b.x[0],axis=0,weights=b.m)
        av=a.v[0]-np.average(a.v[0],axis=0,weights=a.m);bv=b.v[0]-np.average(b.v[0],axis=0,weights=b.m)
        if not np.allclose(ax,bx,rtol=2e-6,atol=2e-6) or not np.allclose(av,bv,rtol=2e-6,atol=2e-6):raise ValueError('Initial particle states differ beyond replay precision')
        if kind=='timestep':
            if a.meta.get('theta')!=b.meta.get('theta'):raise ValueError('Theta changed in timestep comparison')
            if a.meta.get('dt') is None or b.meta.get('dt') is None or a.meta['dt']==b.meta['dt']:raise ValueError('Distinct recorded timesteps required')
            if a.meta.get('solver_sha256') is None or b.meta.get('solver_sha256') is None:provenance.append('Solver binary identity missing')
            elif a.meta['solver_sha256']!=b.meta['solver_sha256']:raise ValueError('Solver binaries differ in timestep comparison')
        if kind=='solver' and a.meta.get('dt')!=b.meta.get('dt'):raise ValueError('Timestep changed in solver comparison')
    else:
        for gal in np.union1d(a.g,b.g):
            for comp in np.union1d(a.c,b.c):
                ma=a.m[(a.g==gal)&(a.c==comp)].sum();mb=b.m[(b.g==gal)&(b.c==comp)].sum()
                if not np.isclose(ma,mb,rtol=1e-8):raise ValueError('Component masses differ')
        if kind=='resolution' and len(a.m)==len(b.m):raise ValueError('Resolution comparison needs different particle counts')
        if a.meta.get('dt')!=b.meta.get('dt') or a.meta.get('theta')!=b.meta.get('theta'):provenance.append('dt/theta differ or are absent; resolution effect not isolated')
        pa=a.meta.get('validation_parent_model_id');pb=b.meta.get('validation_parent_model_id')
        if pa is None or pb is None or pa!=pb:provenance.append('Explicit matching validation_parent_model_id missing; common DF not established')
    pairs=[]
    for i,t in enumerate(a.t):
        j=int(np.argmin(abs(b.t-t)))
        if abs(b.t[j]-t)<=1e-8*max(1,abs(t)):pairs.append((i,j))
    if len(pairs)<2:raise ValueError('Need at least two common saved times; no time interpolation is performed')
    sa,_=structure_history(a);sb,_=structure_history(b)
    aa={(r['frame'],r['galaxy'],r['aperture']):r for r in sa};bb={(r['frame'],r['galaxy'],r['aperture']):r for r in sb}
    rows=[]
    for i,j in pairs:
        for gal in np.unique(a.g):
            for aperture in ['all_disk','initial_R_le_2R50']:
                r=aa.get((i,int(gal),aperture));s=bb.get((j,int(gal),aperture))
                if not r or not s:continue
                for key in ['r50','r90','height','sigma_z','sigma_R','A2']:
                    rows.append(dict(time=float(a.t[i]),galaxy=int(gal),aperture=aperture,metric=key,a=r[key],b=s[key],relative_difference=fraction(s[key]-r[key],r[key])))
    # Compare sampled centre orbits as an additional independent observable.
    if len(np.unique(a.g))==2 and len(np.unique(b.g))==2 and a.labels and b.labels:
        for i,j in pairs:
            vals=[]
            for run,index in [(a,i),(b,j)]:
                centers=[np.median(run.x[index,(run.g==g)&(run.c==1)],axis=0) for g in sorted(np.unique(run.g))]
                vals.append(float(np.linalg.norm(centers[1]-centers[0])))
            rows.append(dict(time=float(a.t[i]),galaxy=-1,aperture='pair',metric='bulge_separation',a=vals[0],b=vals[1],relative_difference=fraction(vals[1]-vals[0],vals[0])))
    return rows,provenance,pairs
