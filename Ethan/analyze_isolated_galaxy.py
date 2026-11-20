"""Analyze saved isolated disk/bulge/halo NPZ histories; no dynamics rerun.

python analyze_isolated_galaxy.py --find-results results --out analysis_output
python analyze_isolated_galaxy.py --trajectory path/to/run_trajectory.npz --out analysis_output
Requires NumPy; Matplotlib is optional for plots.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def quantile(a,w,q):
    order=np.argsort(a)
    return float(np.interp(q,np.cumsum(w[order])/w.sum(),a[order]))


def plane(x,v,w,reference=None):
    center=np.average(x,axis=0,weights=w);bulk=np.average(v,axis=0,weights=w)
    x=x-center;v=v-bulk
    values,basis=np.linalg.eigh((x*w[:,None]).T@x/w.sum())
    n=basis[:,0]
    spin=np.sum(w[:,None]*np.cross(x,v),axis=0)
    if reference is None:
        if n@spin<0:n=-n
        reference=np.array([1.,0.,0.])
        if abs(reference@n)>.9:reference=np.array([0.,1.,0.])
    else:
        if n@reference[1]<0:n=-n
        reference=reference[0]
    e1=reference-n*(reference@n);e1/=np.linalg.norm(e1);e2=np.cross(n,e1)
    xx=x@e1;yy=x@e2;z=x@n;vz=v@n
    r=np.hypot(xx,yy);phi=np.arctan2(yy,xx)
    vr=np.divide(xx*(v@e1)+yy*(v@e2),r,out=np.zeros_like(r),where=r>0)
    return r,phi,z,vz,vr,n,e1,center,bulk


def fit_modes(phi,z,vz,w):
    a=np.column_stack([np.ones(len(phi)),np.cos(phi),np.sin(phi),np.cos(2*phi),np.sin(2*phi)])
    root=np.sqrt(w)
    beta=np.linalg.lstsq(a*root[:,None],z*root,rcond=None)[0]
    residual=z-a@beta
    # Remove low-order azimuthal bulk velocity, then a linear breathing term.
    b=np.column_stack([a,z])
    vb=np.linalg.lstsq(b*root[:,None],vz*root,rcond=None)[0]
    vres=vz-b@vb
    return dict(height=float(np.sqrt(np.average(z*z,weights=w))),
        residual_height=float(np.sqrt(np.average(residual**2,weights=w))),
        m1_bending_amplitude=float(np.hypot(beta[1],beta[2])),
        m2_bending_amplitude=float(np.hypot(beta[3],beta[4])),
        sigma_z=float(np.sqrt(np.average((vz-np.average(vz,weights=w))**2,weights=w))),
        residual_sigma_z=float(np.sqrt(np.average(vres**2,weights=w))),
        coherent_height_variance_share=float(np.average((a@beta)**2,weights=w)/max(np.average(z*z,weights=w),1e-30)),
        breathing_gradient=float(vb[-1]))


def find_run(root):
    matches=[];errors=[]
    for f in sorted(root.rglob('run_trajectory.npz')):
        try:
            # Inspect short time/label arrays before decompressing phase space.
            with np.load(f,allow_pickle=False) as d:
                if not {'t','component','galaxy','metadata'}.issubset(d.files):continue
                t=d['t'];c=d['component'];g=d['galaxy'];meta=json.loads(d['metadata'].item())
                if len(t)<2 or not np.isclose(t[0],0) or not np.isclose(t[-1],1.25):continue
                if len(np.unique(g))!=1:continue
                if [int(np.sum(c==k)) for k in range(3)]!=[3072,1024,32768]:continue
                if not np.isclose(meta.get('eps',-1),.009375):continue
                matches.append(f)
        except (OSError,ValueError,KeyError) as e:errors.append(f'{f}: {e}')
    if errors:print('Unreadable candidates:\n'+'\n'.join(errors))
    if len(matches)!=1:
        listing='\n'.join(str(f) for f in matches) or '(none)'
        raise ValueError(f'Expected one completed matching halo-refined run; found {len(matches)}:\n{listing}\nUse --trajectory with the exact path if there are multiple copies or different sampling.')
    return matches[0]


def analyze(path,out,rd_override=None):
    with np.load(path,allow_pickle=False) as d:
        t=d['t'].astype(float);x=d['x'];v=d['v'];m=d['m'].astype(float);c=d['component'];g=d['galaxy']
        meta=json.loads(d['metadata'].item()) if 'metadata' in d else {}
    if x.ndim!=3 or v.shape!=x.shape or x.shape!=(len(t),len(m),3):raise ValueError('Expected x,v shaped (frames,particles,3)')
    if len(np.unique(g))!=1:raise ValueError('This analysis expects one isolated galaxy, not an encounter')
    if len(t)<2 or np.any(np.diff(t)<=0) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(v)) or not np.all(np.isfinite(m)) or np.any(m<=0):raise ValueError('Invalid time, positions, velocities or masses')
    if not all(np.any(c==k) for k in [0,1,2]):raise ValueError('Expected component 0=disk, 1=bulge, 2=halo')
    source=meta.get('source',{})
    if rd_override is not None:rd=rd_override;rd_basis='explicit --rd'
    elif source.get('model_id')=='zeno_reference_bdh_size_screen_v1':
        rd=float(source['size_screen']['scale'])/12;rd_basis='recorded whole-galaxy size scaling times reference Rd=1/12'
    elif source.get('model_id')=='zeno_reference_bdh_v1':rd=1/12;rd_basis='reference model Rd=1/12'
    else:raise ValueError('Unknown disk scale length: supply --rd explicitly')
    if rd<=0:raise ValueError('Rd must be positive')
    out.mkdir(parents=True,exist_ok=False)
    disk=c==0;wd=m[disk];reference=None;fixed_inner=None
    global_rows=[];annuli=[];components=[];profiles=[]
    edges=rd*np.geomspace(.05,100,61)
    first_profile={};bins=[(0,.5),(.5,1),(1,1.5),(1.5,2),(2,3),(3,4),(4,5)]
    for fi,time in enumerate(t):
        r,phi,z,vz,vr,n,e1,dc,dv=plane(x[fi,disk],v[fi,disk],wd,reference)
        if reference is None:reference=(e1,n);fixed_inner=r<3*rd
        inner=r<3*rd
        a2=abs(np.average(np.exp(2j*phi[inner]),weights=wd[inner]))
        row=dict(time=float(time),height=float(np.sqrt(np.average(z*z,weights=wd))),
            median_abs_z=quantile(abs(z),wd,.5),r50=quantile(r,wd,.5),
            sigma_z=float(np.sqrt(np.average(vz*vz,weights=wd))),
            sigma_R=float(np.sqrt(np.average((vr-np.average(vr,weights=wd))**2,weights=wd))),
            inner_A2=float(a2),inner_count=int(inner.sum()),inner_effective_N=float(wd[inner].sum()**2/np.sum(wd[inner]**2)),
            tilt_from_initial_deg=float(np.degrees(np.arccos(np.clip(n@reference[1],-1,1)))))
        row['A2_sampling_scale']=float(1/np.sqrt(row['inner_effective_N']))
        ii=plane(x[fi,disk][fixed_inner],v[fi,disk][fixed_inner],wd[fixed_inner])
        row['fixed_initial_inner_height']=float(np.sqrt(np.average(ii[2]**2,weights=wd[fixed_inner])))
        global_rows.append(row)
        for lo,hi in bins:
            mask=(r>=lo*rd)&(r<hi*rd);count=int(mask.sum())
            if count<50:continue
            w=wd[mask];a2b=abs(np.average(np.exp(2j*phi[mask]),weights=w))
            annuli.append(dict(time=float(time),R_lo_Rd=lo,R_hi_Rd=hi,count=count,
                A2=float(a2b),**fit_modes(phi[mask],z[mask],vz[mask],w)))
        for k,name in [(0,'disk'),(1,'bulge'),(2,'halo')]:
            use=c==k;w=m[use];xx=x[fi,use].astype(float);vv=v[fi,use].astype(float)
            cc=np.average(xx,axis=0,weights=w);cv=np.average(vv,axis=0,weights=w)
            rad=np.linalg.norm(xx-cc,axis=1);diskrad=np.linalg.norm(xx-dc,axis=1)
            offsets=cc-dc
            components.append(dict(time=float(time),component=name,count=int(use.sum()),mass=float(w.sum()),
                spherical_r10=quantile(rad,w,.1),spherical_r50=quantile(rad,w,.5),spherical_r90=quantile(rad,w,.9),
                offset_from_disk=float(np.linalg.norm(offsets)),vertical_offset_from_disk=float(offsets@n),
                relative_bulk_speed=float(np.linalg.norm(cv-dv))))
            shell=np.histogram(rad,bins=edges,weights=w)[0]
            disk_cumulative=np.array([w[diskrad<edge].sum()/w.sum() for edge in edges[1:]])
            self_cumulative=np.array([w[rad<edge].sum()/w.sum() for edge in edges[1:]])
            if fi==0:first_profile[k]=(disk_cumulative.copy(),self_cumulative.copy())
            for j in range(len(shell)):
                profiles.append(dict(time=float(time),component=name,r_inner_Rd=float(edges[j]/rd),r_outer_Rd=float(edges[j+1]/rd),
                    spherical_density_self_center=float(shell[j]/(4*np.pi/3*(edges[j+1]**3-edges[j]**3))),
                    enclosed_fraction_self_center=float(self_cumulative[j]),enclosed_fraction_disk_center=float(disk_cumulative[j]),
                    enclosed_fraction_change_self_center=float(self_cumulative[j]-first_profile[k][1][j]),
                    enclosed_fraction_change_disk_center=float(disk_cumulative[j]-first_profile[k][0][j])))
        if fi%16==0 or fi==len(t)-1:print(f'Analyzed frame {fi+1}/{len(t)}, t={time:g}',flush=True)
    write_csv(out/'disk_modes.csv',global_rows);write_csv(out/'annular_bending.csv',annuli)
    write_csv(out/'component_structure.csv',components);write_csv(out/'component_profiles.csv',profiles)
    first,last=global_rows[0],global_rows[-1]
    summary=dict(source_path=str(path.resolve()),source_sha256=sha256(path),frames=len(t),particles=len(m),
        time_span=[float(t[0]),float(t[-1])],disk_scale_length=rd,scale_length_basis=rd_basis,
        component_counts={name:int(np.sum(c==k)) for k,name in enumerate(['disk','bulge','halo'])},
        initial=first,final=last,relative_changes={key:last[key]/first[key]-1 for key in ['height','median_abs_z','r50','sigma_z','sigma_R','fixed_initial_inner_height']},
        peak_inner_A2=max(global_rows,key=lambda q:q['inner_A2']),
        component_radius_changes={name:{key:next(q for q in reversed(components) if q['component']==name)[key]/next(q for q in components if q['component']==name)[key]-1 for key in ['spherical_r10','spherical_r50','spherical_r90']} for name in ['disk','bulge','halo']},
        caveats=['Measurements only; no automatic physical or observational validation.',
        'Annuli use instantaneous cylindrical radius; radial migration changes their membership.',
        'Low-order bending/breathing subtraction is descriptive, not a unique separation of random heating.',
        'Inner A2 includes bars and two-armed spirals; A2 alone does not distinguish them.',
        'A2 sampling scale is an order-of-magnitude Poisson reference, not a confidence interval.',
        'Spherical component profiles and offsets do not measure the full local restoring force.'],metadata=meta)
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    try:make_plots(out,global_rows,annuli,components,profiles,rd)
    except ImportError:print('Matplotlib unavailable; CSV/JSON analysis completed without plots.')
    print('Finished:',out.resolve());print('Send summary.json and the four CSV files, or zip this output folder.')


def make_plots(out,rows,annuli,components,profiles,rd):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'figure.facecolor':'#08090e','axes.facecolor':'#08090e','text.color':'white','axes.labelcolor':'#ccccda','xtick.color':'#ccccda','ytick.color':'#ccccda','axes.edgecolor':'#555565','legend.facecolor':'#08090e','legend.labelcolor':'white'})
    fig,axs=plt.subplots(2,2,figsize=(11,8));tt=[q['time'] for q in rows]
    for key,label,color in [('height','Whole disk RMS height','#62d5e8'),('fixed_initial_inner_height','Fixed initial inner particles','#f465b5')]:
        yy=np.array([q[key] for q in rows]);axs[0,0].plot(tt,yy/yy[0],label=label,color=color)
    axs[0,0].set_title('Thickness / initial');axs[0,0].legend(fontsize=8)
    for key,label,color in [('sigma_R','Radial','#f465b5'),('sigma_z','Vertical','#62d5e8')]:
        yy=np.array([q[key] for q in rows]);axs[0,1].plot(tt,yy/yy[0],label=label,color=color)
    axs[0,1].set_title('Velocity dispersion / initial');axs[0,1].legend()
    axs[1,0].plot(tt,[q['inner_A2'] for q in rows],color='#ffbd73',label='Inner A2')
    axs[1,0].plot(tt,[q['A2_sampling_scale'] for q in rows],color='#aaaabc',ls=':',label='1/sqrt(effective N), indicative only')
    axs[1,0].set_title('Inner two-fold structure');axs[1,0].legend(fontsize=8)
    for name,color in [('bulge','#ffbd73'),('halo','#f465b5')]:
        qq=[q for q in components if q['component']==name]
        axs[1,1].plot([q['time'] for q in qq],[q['offset_from_disk']/rd for q in qq],label=name,color=color)
    axs[1,1].set_title('Component center offset from disk / Rd');axs[1,1].legend()
    for ax in axs.flat:ax.set_xlabel('Code time')
    fig.tight_layout();fig.savefig(out/'evolution.png',dpi=160);plt.close(fig)
    ranges=sorted(set((q['R_lo_Rd'],q['R_hi_Rd']) for q in annuli));fig,axs=plt.subplots(2,4,figsize=(15,7))
    for ax,(lo,hi) in zip(axs.flat,ranges):
        qq=[q for q in annuli if q['R_lo_Rd']==lo and q['R_hi_Rd']==hi]
        for key,color,label in [('height','#f465b5','Raw height'),('residual_height','#62d5e8','After fitted bending')]:
            initial=qq[0][key];ax.plot([q['time'] for q in qq],[q[key]/initial for q in qq],color=color,label=label)
        ax.set_title(f'{lo:g}–{hi:g} Rd');ax.set_xlabel('Code time')
    axs.flat[0].legend(fontsize=8)
    for ax in list(axs.flat)[len(ranges):]:ax.set_visible(False)
    fig.suptitle('Annular thickness / its own initial value; evolving radial membership')
    fig.tight_layout();fig.savefig(out/'bending.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    times=sorted(set(q['time'] for q in profiles));chosen=sorted(set([times[0],times[len(times)//2],times[-1]]))
    for ax,name in zip(axs,['disk','bulge','halo']):
        for time,color in zip(chosen,['#62d5e8','#ffbd73','#f465b5']):
            qq=[q for q in profiles if q['component']==name and q['time']==time]
            ax.plot([q['r_outer_Rd'] for q in qq],[q['enclosed_fraction_disk_center'] for q in qq],label=f't={time:g}',color=color)
        ax.set_xscale('log');ax.set_title(name);ax.set_xlabel('Spherical radius about disk center / Rd');ax.set_ylabel('Enclosed mass fraction');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'component_profiles.png',dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    src=p.add_mutually_exclusive_group(required=True);src.add_argument('--trajectory',type=Path);src.add_argument('--find-results',type=Path)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--rd',type=float)
    a=p.parse_args()
    try:
        source=a.trajectory if a.trajectory is not None else find_run(a.find_results)
        print('Using trajectory:',source.resolve(),flush=True)
        analyze(source,a.out,a.rd)
    except (OSError,ValueError,KeyError) as e:p.exit(1,f'Analysis stopped: {e}\n')


if __name__=='__main__':main()
