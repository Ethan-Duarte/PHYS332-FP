"""CLI for reproducible, explicitly scoped N-body validation measurements."""
import argparse,csv,json,subprocess,sys,traceback
import platform,importlib.metadata
from pathlib import Path
import numpy as np
from core import *

def output_dir(path):
    path=Path(path);path.mkdir(parents=True,exist_ok=False);return path
def policy_for(a):
    policy=json.loads(a.policy.read_text()) if getattr(a,'policy',None) else {}
    if policy and policy.get('reviewed') is True and not policy.get('rationale'):raise ValueError('Reviewed policy requires a rationale')
    for key,value in policy.get('limits',{}).items():
        if not isinstance(value,(int,float)) or not np.isfinite(value) or value<0:raise ValueError('Policy limits must be finite nonnegative numbers: '+key)
    return policy
def finish(out,command,inputs,checks,measurements,notes,policy=None):
    statuses=[x['status'] for x in checks]
    versions={}
    for name in ['numpy','scipy','numba','matplotlib']:
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]=None
    suite_hashes={p.name:digest(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
    result=dict(schema_version=1,command=command,suite_source_sha256=suite_hashes,environment=dict(python=sys.version,platform=platform.platform(),packages=versions),inputs=[dict(path=str(Path(p).resolve()),sha256=digest(p)) for p in inputs],
        status='FAIL' if 'FAIL' in statuses else ('INCOMPLETE' if any(s in ['MISSING','INCONCLUSIVE'] for s in statuses) else ('PASS_FOR_SPECIFIED_CHECKS' if statuses and all(s=='PASS' for s in statuses) else 'MEASURED')),
        checks=checks,measurements=measurements,notes=notes,policy=policy or {},
        claim='No automatic certification of physical or observational validity. Passing applies only to named tests and declared tolerances.')
    write_json(out/'report.json',result)
    lines=[f'# {command}: {result["status"]}','',result['claim'],'','| Check | Status | Value | Limit |','|---|---|---|---|']
    for c in checks:lines.append(f"| {c['name']} | {c['status']} | {c.get('value','')} | {c.get('limit','')} |")
    lines+=['','## Notes','']+['- '+x for x in notes]
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],report=str(out/'report.json'),checks=checks),indent=2))
    return 2 if result['status']=='FAIL' else 0

def inspect(a):
    out=output_dir(a.out);r=load(a.input);policy=policy_for(a)
    if a.context=='isolated' and len(np.unique(r.g))!=1:raise ValueError('Isolated context needs a single galaxy; pair cannot certify isolated stability')
    if a.energy_samples and not a.self_gravity_only and r.meta.get('dynamics_model')!='self_gravitating_collisionless':raise ValueError('Energy requires --self-gravity-only or dynamics_model=self_gravitating_collisionless metadata; external fields/gas need different energy terms')
    summary,conservation,energy=conserved(r,a.energy_samples,a.epsilon,a.G,a.kernel)
    sh,radial=structure_history(r)
    csv_write(out/'conservation.csv',conservation);csv_write(out/'sampled_energy.csv',energy)
    csv_write(out/'structure.csv',sh);csv_write(out/'radial_structure.csv',radial)
    checks=[dict(name='Array shapes, finite values, positive mass and increasing times',status='PASS')]
    for key,value in summary.items():
        if isinstance(value,(int,float)) or value is None:
            if key!='energy_samples':checks.append(criterion(key,value,policy,key))
    if not energy:checks.append(dict(name='Energy conservation',status='MISSING',notes='Use --energy-samples >=2 with known force law'))
    changes=[]
    for g in np.unique(r.g):
        for aperture in ['all_disk','initial_R_le_2R50']:
            group=[x for x in sh if x['galaxy']==g and x['aperture']==aperture]
            if not group:continue
            for key in ['r50','r90','height','sigma_z','sigma_R','A2']:
                vals=[fraction(x[key]-group[0][key],group[0][key]) for x in group]
                measured=max(abs(x) for x in vals if x is not None) if any(x is not None for x in vals) else None
                changes.append(dict(galaxy=int(g),aperture=aperture,metric=key,endpoint_relative_change=vals[-1],max_abs_relative_change=measured))
                if a.context=='isolated' and key not in ['A2']:
                    checks.append(criterion(f'g{g} {aperture} {key} evolution',measured,policy,'isolated_'+key))
    if not r.labels:checks.append(dict(name='Disk/component diagnostics',status='MISSING',notes='Need component and galaxy labels'))
    if a.context!='isolated':checks.append(dict(name='Isolated progenitor stability',status='MISSING',notes='Requires a separate isolated trajectory'))
    if len(r.t)<2:checks.append(dict(name='Evolution evidence',status='MISSING'))
    notes=['Energy is calculated with direct pair sums at saved frames only; it does not bound all-step drift.',
        'Energy, if requested, assumes fixed-mass self-gravity only; no external, cosmological or thermodynamic energy terms.',
        'Fitted planes remove rigid tilt, not all warping. Radial cell residuals separate some coherent geometry but are not exact heating estimators.',
        'Cell residuals cover only cells with >=10 particles; inspect cell_mass_coverage.',
        'Initial_R_le_2R50 is a fixed set of initial particle IDs; annuli use instantaneous radii and fixed initial radial edges.',
        'A1/A2 and rotation sign are coordinate conventions; mean_rotation_abs avoids arbitrary eigenvector signs.',
        'Structural changes in encounters may be physical. No stability gate is applied to an encounter.',
        'Component labels: 0=disk, 1=bulge, 2=halo. Float32 replay precision limits very small residuals.']
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,3,figsize=(13,3.5))
        for gal in np.unique(r.g):
            for aperture,ls in [('all_disk','-'),('initial_R_le_2R50','--')]:
                group=[z for z in sh if z['galaxy']==gal and z['aperture']==aperture]
                if group:
                    for ax,key in zip(axes[:2],['height','sigma_z']):
                        if group[0][key]>0:ax.plot([z['time'] for z in group],[z[key]/group[0][key] for z in group],ls,label=f'g{gal} {aperture}')
        if sh:axes[0].legend(fontsize=6)
        for ax,title in zip(axes,['Relative fitted-plane height','Relative vertical dispersion','Sampled relative energy change']):ax.set(xlabel='Code time',title=title)
        if energy:axes[2].plot([e['time'] for e in energy],[e['relative_change'] for e in energy])
        fig.tight_layout();fig.savefig(out/'diagnostics.png',dpi=140);plt.close(fig)
    except ImportError:notes.append('matplotlib unavailable; numerical outputs were still written')
    return finish(out,'inspect',[a.input]+([a.policy] if a.policy else []),checks,dict(conservation=summary,structural_changes=changes,metadata=r.meta),notes,policy)

def compare(a):
    out=output_dir(a.out);left=load(a.a);right=load(a.b);policy=policy_for(a)
    rows,issues,pairs=comparison(left,right,a.kind);csv_write(out/'comparison.csv',rows)
    checks=[dict(name='Compatible controlled inputs',status='PASS' if not issues else 'INCONCLUSIVE',notes='; '.join(issues)),dict(name='Common saved times',status='PASS',value=len(pairs))]
    metrics=[]
    for key in sorted(set(x['metric'] for x in rows)):
        vals=[abs(x['relative_difference']) for x in rows if x['metric']==key and x['relative_difference'] is not None]
        val=max(vals) if vals else None;metrics.append(dict(metric=key,max_relative_difference=val))
        checks.append(criterion(key+' agreement',val,policy,'compare_'+key))
    if not rows:checks.append(dict(name='Structural/orbit comparison',status='MISSING',notes='Requires labels and disk or bulge data'))
    return finish(out,'compare',[a.a,a.b]+([a.policy] if a.policy else []),checks,dict(kind=a.kind,metrics=metrics,common_time_start=float(left.t[pairs[0][0]]),common_time_end=float(left.t[pairs[-1][0]])),
        ['Only common saved times are compared; no interpolation or time stretching.',
         'Agreement concerns ensemble structures/orbit, not indefinite particle-by-particle equality.',
         'A single resolution or seed comparison does not establish a convergence rate.',*issues],policy)

def forces(a):
    out=output_dir(a.out);r=load(a.input,snapshot=True);policy=policy_for(a);eps,G,k=law(r,a.epsilon,a.G,a.kernel)
    if a.frame<0 or a.frame>=len(r.t):raise ValueError('Frame index outside file')
    with np.load(a.acceleration,allow_pickle=False) as d:
        ids=d['indices'];acc=np.asarray(d['a'],dtype=float)
        provenance=json.loads(d['metadata'].item()) if 'metadata' in d else {}
    if ids.ndim!=1 or not np.equal(ids,np.floor(ids)).all():raise ValueError('indices must be integer vector')
    ids=ids.astype(int)
    if len(ids)==0 or len(np.unique(ids))!=len(ids) or np.any(ids<0) or np.any(ids>=len(r.m)) or acc.shape!=(len(ids),3) or not np.isfinite(acc).all():raise ValueError('Invalid acceleration sample')
    issues=[]
    for key,want in [('trajectory_sha256',digest(a.input)),('frame',a.frame),('epsilon',eps),('G',G),('kernel','plummer')]:
        if key not in provenance:issues.append('Acceleration provenance missing '+key)
        elif provenance[key]!=want:raise ValueError('Acceleration sample provenance mismatch: '+key)
    ref=reference_force(r.x[a.frame],r.m,ids,eps,G);err=np.linalg.norm(acc-ref,axis=1);norm=np.linalg.norm(ref,axis=1)
    if not np.isfinite(ref).all():raise ValueError('Reference force is singular/nonfinite; check coincident unsoftened particles')
    floor=max(float(np.median(norm))*1e-12,1e-30);relative=err/np.maximum(norm,floor)
    csv_write(out/'force_errors.csv',[dict(index=int(i),reference_norm=float(n),absolute_error=float(e),relative_error=float(re),near_zero_reference=bool(n<floor)) for i,n,e,re in zip(ids,norm,err,relative)])
    values=dict(median_relative_force_error=float(np.median(relative)),p95_relative_force_error=float(np.quantile(relative,.95)),max_relative_force_error=float(relative.max()),sample_size=len(ids),relative_denominator_floor=floor)
    checks=[dict(name='Acceleration input provenance',status='PASS' if not issues else 'INCONCLUSIVE',notes='; '.join(issues))]
    checks += [criterion(key,value,policy,key) for key,value in values.items() if key.endswith('error')]
    return finish(out,'forces',[a.input,a.acceleration]+([a.policy] if a.policy else []),checks,values,
        ['Independent vectorized direct Plummer reference; self-force explicitly excluded.',
         'Errors measured at supplied target indices only. Include inner, outer and low-density targets; a random sample may miss difficult locations.',*issues],policy)

def initial(a):
    out=output_dir(a.out);r=load(a.input,snapshot=True);policy=policy_for(a);spec=json.loads(a.targets.read_text());checks=[];results=[]
    input_files=[a.input,a.targets]+([a.policy] if a.policy else [])
    if not r.labels:raise ValueError('Initial-condition targets need component and galaxy labels')
    for target in spec['components']:
        gal=int(target['galaxy']);comp=int(target['component']);mask=(r.g==gal)&(r.c==comp)
        if not mask.any():checks.append(dict(name=f'g{gal} c{comp} exists',status='FAIL'));continue
        x=r.x[0,mask];m=r.m[mask];x-=np.average(x,axis=0,weights=m)
        if target.get('coordinate','spherical_radius')=='disk_radius':radius=plane(x,r.v[0,mask],m)[0]
        elif target.get('coordinate','spherical_radius')=='spherical_radius':radius=np.linalg.norm(x,axis=1)
        else:raise ValueError('Unknown target coordinate')
        entry=dict(galaxy=gal,component=comp,mass=float(m.sum()),count=int(mask.sum()))
        if 'mass' in target:
            if target['mass']<=0:raise ValueError('Target mass must be positive')
            error=fraction(m.sum()-target['mass'],target['mass']);entry['mass_relative_error']=error
            checks.append(criterion(f'g{gal} c{comp} target mass',error,policy,'initial_mass_relative_error'))
        if 'cdf_csv' in target:
            path=a.targets.parent/target['cdf_csv'];data=np.genfromtxt(path,delimiter=',',names=True)
            input_files.append(path)
            rr=np.atleast_1d(data['radius']);ff=np.atleast_1d(data['cumulative_fraction'])
            if len(rr)<2 or not np.isfinite(rr).all() or not np.isfinite(ff).all() or np.any(np.diff(rr)<=0) or np.any(np.diff(ff)<0) or ff.min()<0 or ff.max()>1:raise ValueError('Invalid cumulative profile table')
            ids=np.argsort(radius);sortedr=radius[ids];cm=np.cumsum(m[ids])/m.sum();indices=np.searchsorted(sortedr,rr,side='right')
            empirical=np.array([cm[i-1] if i else 0 for i in indices]);error=float(np.max(abs(empirical-ff)))
            entry['max_cdf_difference_on_table_grid']=error;checks.append(criterion(f'g{gal} c{comp} radial CDF',error,policy,'initial_cdf_difference'))
            csv_write(out/f'cdf_g{gal}_c{comp}.csv',[dict(radius=float(z),target=float(t),sample=float(e)) for z,t,e in zip(rr,ff,empirical)])
        if 'moments_csv' in target:
            path=a.targets.parent/target['moments_csv'];input_files.append(path)
            with path.open(newline='') as f:bins=list(csv.DictReader(f))
            if not bins:raise ValueError('Velocity-profile table is empty')
            RR,zz,vz,phi,vr,vt,normal=plane(r.x[0,mask],r.v[0,mask],m)
            moments=[]
            for j,b in enumerate(bins):
                lo=float(b['R_lo']);hi=float(b['R_hi'])
                if not 0<=lo<hi:raise ValueError('Invalid velocity-profile bin')
                use=(RR>=lo)&(RR<hi)
                if use.sum()<10:
                    checks.append(dict(name=f'g{gal} c{comp} moment bin {j}',status='MISSING',notes='Fewer than 10 sampled particles'));continue
                weights=m[use]
                measured=dict(rotation_abs=float(abs(np.average(vt[use],weights=weights))),sigma_R=float(np.sqrt(np.average((vr[use]-np.average(vr[use],weights=weights))**2,weights=weights))),sigma_z=float(np.sqrt(np.average((vz[use]-np.average(vz[use],weights=weights))**2,weights=weights))))
                for key in measured:
                    if b.get(key,'').strip():
                        desired=float(b[key]);err=fraction(measured[key]-desired,desired)
                        checks.append(criterion(f'g{gal} c{comp} bin{j} {key} target',err,policy,'initial_'+key+'_relative_error'))
                        moments.append(dict(bin=j,R_lo=lo,R_hi=hi,metric=key,measured=measured[key],target=desired,relative_error=err))
            csv_write(out/f'moments_g{gal}_c{comp}.csv',moments)
        results.append(entry)
    checks.append(dict(name='Full phase-space equilibrium',status='MISSING',notes='Spatial targets do not establish velocity support or stationarity; use velocity profiles and matched isolation controls'))
    return finish(out,'initial',input_files,checks,dict(components=results),
        ['CDF error is evaluated at supplied table radii, not an exact KS statistic; no p-value is inferred.',
         'Profiles use each component COM. Use target tables that include the actual truncation/taper.',
         'A profile inherited from a different galaxy construction must not be treated as the target for this one.'],policy)

def observations(a):
    out=output_dir(a.out);r=load(a.input);policy=policy_for(a);s=json.loads(a.setup.read_text())
    if not r.labels:raise ValueError('Observation comparison requires labels')
    if not r.t[0]<=s['time']<=r.t[-1]:raise ValueError('Observation time outside trajectory')
    if s.get('length_unit')!='kpc' or s.get('velocity_unit')!='km/s':raise ValueError('Setup must explicitly use kpc and km/s')
    L=float(s['length_scale']);V=float(s['velocity_scale'])
    if L<=0 or V<=0:raise ValueError('Positive physical scales required')
    fi=int(np.argmin(abs(r.t-s['time'])));B=camera(s['azimuth'],s['elevation']);roll=np.radians(s.get('roll',0));R=np.array([[np.cos(roll),-np.sin(roll)],[np.sin(roll),np.cos(roll)]])
    # Origin is the entire system's instantaneous mass COM, with bulk velocity removed.
    xx=r.x[fi]-np.average(r.x[fi],axis=0,weights=r.m);vv=r.v[fi]-np.average(r.v[fi],axis=0,weights=r.m)
    xy=(xx@B[:2].T)@R.T*L+np.array(s['sky_offset']);rv=vv@B[2]*V+s['systemic_velocity']
    constraints=a.setup.parent/s['constraints_csv']
    with constraints.open(newline='') as f:regions=list(csv.DictReader(f))
    if not regions:raise ValueError('No observational regions')
    output=[]
    for row in regions:
        if row['role'] not in ['fit','holdout']:raise ValueError('Every region needs role=fit or holdout')
        if row['kind'] not in ['tail','nucleus']:raise ValueError('kind must be tail or nucleus')
        gal=int(row['galaxy']);component=0 if row['kind']=='tail' else 1
        pool=(r.c==component)&((r.g==gal) if gal>=0 else True)
        sx,sy,rad=[float(row[k]) for k in ['x','y','radius']]
        if rad<=0 or not np.isfinite([sx,sy,rad]).all():raise ValueError('Invalid observation aperture')
        use=pool&(np.sum((xy-[sx,sy])**2,axis=1)<=rad*rad);count=int(use.sum())
        result=dict(id=row['id'],role=row['role'],kind=row['kind'],galaxy=gal,count=count,
            region_mass_fraction=fraction(r.m[use].sum(),r.m[pool].sum()),x=sx,y=sy,radius=rad,
            model_velocity_median=None,model_v16=None,model_v84=None,observed_vlos=None,observed_sigma_v=None,velocity_residual=None,residual_over_measurement_sigma=None)
        if row.get('vlos','').strip():
            obs=float(row['vlos']);sigma=float(row['sigma_v'])
            if sigma<=0 or not np.isfinite([obs,sigma]).all():raise ValueError('Velocity observations need finite positive sigma_v')
            result.update(observed_vlos=obs,observed_sigma_v=sigma)
            if count:
                med=quantile(rv[use],r.m[use],.5);result.update(model_velocity_median=med,model_v16=quantile(rv[use],r.m[use],.16),model_v84=quantile(rv[use],r.m[use],.84),velocity_residual=med-obs,residual_over_measurement_sigma=(med-obs)/sigma)
        output.append(result)
    csv_write(out/'observational_regions.csv',output)
    np.savez_compressed(out/'sky_pv.npz',xy_kpc=xy,rv_kms=rv,m=r.m,component=r.c,galaxy=r.g)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        fig,axes=plt.subplots(1,3,figsize=(14,4))
        for j,gal in enumerate(np.unique(r.g)):
            mask=(r.c==0)&(r.g==gal);color=['#179bb2','#cc428a'][j%2]
            axes[0].scatter(*xy[mask].T,s=.5,alpha=.35,c=color)
            axes[1].scatter(xy[mask,0],rv[mask],s=.5,alpha=.35,c=color)
            axes[2].scatter(rv[mask],xy[mask,1],s=.5,alpha=.35,c=color)
        for row in output:
            color='#df8f00' if row['role']=='fit' else '#16813b'
            axes[0].add_patch(Circle((row['x'],row['y']),row['radius'],fill=False,color=color,lw=1))
            axes[0].annotate(row['id'],(row['x'],row['y']),fontsize=6,color=color)
            if row['observed_vlos'] is not None:
                axes[1].errorbar(row['x'],row['observed_vlos'],yerr=row['observed_sigma_v'],fmt='o',color=color,ms=3)
                axes[2].errorbar(row['observed_vlos'],row['y'],xerr=row['observed_sigma_v'],fmt='o',color=color,ms=3)
        axes[0].set(xlabel='X (kpc)',ylabel='Y down (kpc)',title='Model and observation apertures',aspect='equal');axes[0].invert_yaxis()
        axes[1].set(xlabel='X (kpc)',ylabel='Receding velocity (km/s)',title='X-velocity')
        axes[2].set(xlabel='Receding velocity (km/s)',ylabel='Y down (kpc)',title='Velocity-Y');axes[2].invert_yaxis()
        fig.suptitle('Supplied constraints: orange=fit, green=holdout. Aperture radii are not positional error bars.',fontsize=10)
        fig.tight_layout();fig.savefig(out/'observational_comparison.png',dpi=150);plt.close(fig)
    except ImportError:pass
    checks=[]
    for row in output:
        checks.append(dict(name=row['id']+' populated',status='PASS' if row['count'] else 'FAIL',value=row['count']))
        if row['observed_vlos'] is not None:checks.append(criterion(row['id']+' velocity residual / observational sigma',row['residual_over_measurement_sigma'],policy,'observation_velocity_sigma'))
    holdout=[z for z in output if z['role']=='holdout' and z['observed_vlos'] is not None]
    checks.append(dict(name='Independent velocity holdout supplied',status='PASS' if holdout else 'MISSING',value=len(holdout)))
    return finish(out,'observations',[a.input,a.setup,constraints]+([a.policy] if a.policy else []),checks,
        dict(frame=fi,time=float(r.t[fi]),setup=s,regions=output,time_unit_Myr=977.79222168*L/V,mass_unit_Msun=V*V*L/4.30091727003628e-6),
        ['Positive LOS velocity is away from observer; X is camera right and Y is camera down after the specified roll.',
         'No camera, time or scale is fitted by this command. Supply one registered setup for all constraints.',
         'Region occupancy is not image-intensity agreement. Model velocity width is reported separately; residual/sigma is not a calibrated chi-square.',
         'Fit rows are fitted constraints, not independent validation. Holdout rows must genuinely have been withheld.',
         'No automatic gas, dust, star-formation or optical-emissivity prediction.'],policy)

def batch(a):
    out=output_dir(a.out);spec=json.loads(a.manifest.read_text());results=[]
    if not spec['jobs']:raise ValueError('Batch must contain at least one job')
    allowed={'inspect','compare','forces','initial','observations','selftest'}
    for i,job in enumerate(spec['jobs']):
        args=job['args']
        if not isinstance(args,list) or not args or args[0] not in allowed or '--out' in args:raise ValueError('Batch args must start with allowed command and omit --out')
        dest=(out/f'{i:02d}_{args[0]}').resolve()
        done=subprocess.run([sys.executable,str(Path(__file__).resolve()),*map(str,args),'--out',str(dest)],cwd=a.manifest.resolve().parent,text=True,capture_output=True)
        (out/f'{i:02d}_stdout.txt').write_text(done.stdout);(out/f'{i:02d}_stderr.txt').write_text(done.stderr)
        report=json.loads((dest/'report.json').read_text()) if (dest/'report.json').exists() else {}
        results.append(dict(job=i,command=args[0],exit_code=done.returncode,status=report.get('status','ERROR'),report=str(dest/'report.json')))
    write_json(out/'batch.json',dict(jobs=results,claim='Results remain separate; no automatic whole-project validation verdict.'))
    print(json.dumps(results,indent=2));return 2 if any(x['exit_code'] for x in results) else 0

def parser():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    def common(name):
        s=sub.add_parser(name);s.add_argument('--out',type=Path,required=True);s.add_argument('--policy',type=Path);return s
    def gravity(s):
        s.add_argument('--epsilon',type=float);s.add_argument('--G',type=float);s.add_argument('--kernel',choices=['plummer'])
    s=common('inspect');s.add_argument('input',type=Path);s.add_argument('--context',choices=['isolated','encounter','unspecified'],default='unspecified');s.add_argument('--energy-samples',type=int,default=0);s.add_argument('--self-gravity-only',action='store_true');gravity(s)
    s=common('compare');s.add_argument('a',type=Path);s.add_argument('b',type=Path);s.add_argument('--kind',choices=['timestep','solver','resolution','seed'],required=True)
    s=common('forces');s.add_argument('input',type=Path);s.add_argument('acceleration',type=Path);s.add_argument('--frame',type=int,default=0);gravity(s)
    s=common('initial');s.add_argument('input',type=Path);s.add_argument('targets',type=Path)
    s=common('observations');s.add_argument('input',type=Path);s.add_argument('setup',type=Path)
    s=common('selftest');s.add_argument('--repo-root',type=Path,default=Path('.'))
    s=sub.add_parser('batch');s.add_argument('manifest',type=Path);s.add_argument('--out',type=Path,required=True)
    return p

def main():
    a=parser().parse_args()
    existed=Path(a.out).exists()
    try:
        if a.command=='inspect' and (a.energy_samples<0 or a.energy_samples==1):raise ValueError('Energy samples must be 0 or >=2')
        if a.command=='selftest':
            from solver_selftest import run
            return run(a,finish,output_dir)
        return globals()[a.command](a)
    except Exception as e:
        # Preserve an explicit machine-readable invalid-input result; never overwrite another run.
        print(type(e).__name__+': '+str(e),file=sys.stderr)
        if not existed and Path(a.out).exists() and not (Path(a.out)/'report.json').exists():
            write_json(Path(a.out)/'error.json',dict(status='ERROR',error=type(e).__name__,message=str(e)))
        return 3

if __name__=='__main__':sys.exit(main())
