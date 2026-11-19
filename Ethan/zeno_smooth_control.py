"""Fixed smooth-field diagnostic for the existing Zeno reference disk.

Requires NumPy and SciPy. Reads the original seven-column initial.txt;
does not regenerate, recenter, damp, or change any initial velocities.
The first 24,576 rows must be disk particles from the documented stack.
G=1; exact Plummer kernel epsilon=.0075. This is a diagnostic control,
not a new initial-condition generator or a validated Mice reconstruction.
"""
import argparse
import csv
import json
import hashlib
import time
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre, ellipk, ellipe
from scipy.interpolate import RectBivariateSpline, CubicSpline
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid


def vertical_jeans_variance(z, restoring_force, height):
    """Solve rho*sigma_z^2 = integral_z^infinity rho*g_z dz.

    Assumes an axisymmetric stationary target with zero R-z stress. Matching
    this moment is necessary under that assumption, not a full DF solution.
    The supplied grid must reach negligible target density at its upper edge.
    """
    z = np.asarray(z, float)
    force = np.asarray(restoring_force, float)
    if (height <= 0 or z.ndim != 1 or len(z) < 3 or z[0] != 0
            or np.any(np.diff(z) <= 0) or force.shape[-1] != len(z)
            or not np.isfinite(force).all() or np.any(force < -1e-10)):
        raise ValueError('Invalid vertical Jeans grid or restoring force')
    rho = 1 / np.cosh(z / height)**2
    if rho[-1] > 1e-10:
        raise ValueError('Vertical grid must extend into negligible target density')
    integrand = rho * np.maximum(force, 0)
    pressure = -cumulative_trapezoid(integrand[..., ::-1], z[::-1], axis=-1,
                                     initial=0)[..., ::-1]
    return pressure / rho

EPS = .0075
ZD = .0075
RD = 1/12


def quad_disk(order=8, vertical=64):
    t, w = roots_legendre(order)
    edges = np.linspace(0, 1, 51)
    radius = ((edges[:-1,None]+edges[1:,None])/2+.01*t).ravel()
    rw = np.tile(.01*w, 50)
    u, uw = roots_legendre(vertical)
    heights = ZD*np.arctanh(u)
    mass = (.1875*144*radius*np.exp(-12*radius)*rw)[:,None]*(uw/2)[None,:]
    mass /= 1-13*np.exp(-12)
    return radius[:,None], heights[None,:], mass


def disk_field(points, quadrature, acceleration=False):
    rp, zp, mass = quadrature
    values = []
    for start in range(0, len(points), 12):
        q = points[start:start+12]
        R = q[:,0,None,None]
        dz = q[:,1,None,None]-zp
        ap = (R+rp)**2+dz*dz+EPS*EPS
        am = (R-rp)**2+dz*dz+EPS*EPS
        m = 4*R*rp/ap
        K = ellipk(m)
        if acceleration:
            E = ellipe(m)
            az = -2/np.pi*np.sum(mass*dz*E/(am*np.sqrt(ap)),axis=(1,2))
            ar = -np.sum(mass*(K-(rp*rp-R*R+dz*dz+EPS*EPS)*E/am)/
                         (np.pi*R*np.sqrt(ap)),axis=(1,2))
            values.append(np.column_stack([ar,az]))
        else:
            values.append(-2/np.pi*np.sum(mass*K/np.sqrt(ap),axis=(1,2)))
    return np.concatenate(values)


def sphere_quadrature():
    # Composite logarithmic quadrature with explicit taper breakpoints.
    edges = np.unique(np.r_[0, np.geomspace(1e-7, 80, 301), .98015, 4.])
    t,w = roots_legendre(8)
    s = ((edges[:-1,None]+edges[1:,None])/2+
         (edges[1:,None]-edges[:-1,None])*t/2).ravel()
    sw = ((edges[1:,None]-edges[:-1,None])*w/2).ravel()
    a,b = .02,4.
    rb = a/(2*np.pi*b*(a+b)**3)
    norm = .0625/(b*b/(b+a)**2+2*np.pi*b**3*rb)
    bulge = np.where(s<=b, norm*a/(2*np.pi*s*(s+a)**3),
                     norm*rb*(b/s)**2*np.exp(-2*(s-b)/b))
    a,b = .25,.98015
    mu = .16/(np.log(2)-.5)
    rb = mu/(4*np.pi*b*(a+b)**2)
    gamma = b/a-(a+3*b)/(a+b)
    halo = np.where(s<=b, mu/(4*np.pi*s*(a+s)**2),
                    rb*(s/b)**gamma*np.exp(-(s-b)/a))
    return s, 4*np.pi*s*s*(bulge+halo)*sw


def sphere_field(radii, quadrature, acceleration=False):
    s, mass = quadrature
    result=[]
    for start in range(0,len(radii),128):
        r=np.asarray(radii[start:start+128])[:,None]
        hp=np.sqrt((r+s)**2+EPS*EPS)
        hm=np.sqrt((r-s)**2+EPS*EPS)
        # Rationalized shell potential avoids cancellation near r=0.
        potential=-2/(hp+hm)
        if acceleration:
            derivative=-2*((r+s)/hp+(r-s)/hm)/(hp+hm)**2
            result.append(np.sum(mass*derivative,axis=1))
        else:
            result.append(np.sum(mass*potential,axis=1))
    return np.concatenate(result)


def measure(x,v):
    p=x-x.mean(0);vel=v-v.mean(0)
    _,axes=np.linalg.eigh(p.T@p)
    n=axes[:,0]
    if n[2]<0:n=-n
    z=p@n;vz=vel@n
    r=np.sqrt(np.maximum((p*p).sum(1)-z*z,0))
    row=dict(height=float(np.sqrt(np.mean(z*z))),
             height90=float(np.quantile(abs(z),.9)),
             sigma_z=float(np.std(vz)),R50=float(np.median(r)),
             R90=float(np.quantile(r,.9)))
    for lo,hi in [(0,1),(1,2),(2,3),(3,4)]:
        k=(r>=lo*RD)&(r<hi*RD)
        row[f'ring_{lo}_{hi}_height']=float(np.std(z[k]))
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--initial',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('results/zeno_smooth_control'))
    p.add_argument('--duration',type=float,default=1.5)
    p.add_argument('--disk',type=int,default=24576)
    p.add_argument('--nr',type=int,default=200)
    p.add_argument('--nz',type=int,default=145)
    p.add_argument('--vertical-jeans', action='store_true',
                   help='Create and screen an experimental local-moment-adjusted candidate; not a full DF equilibrium')
    a=p.parse_args()
    if a.duration<=0 or a.disk<100 or a.nr<80 or a.nz<80:
        p.error('Positive duration, >=100 disk particles, and grids >=80 required.')
    data=np.loadtxt(a.initial)
    if data.shape!=(65536,7) or not np.isfinite(data).all() or a.disk!=24576:
        raise SystemExit('Expected original 65,536-row, seven-column reference file and 24,576 disk particles.')
    if not np.allclose(data[:a.disk,0],.1875/a.disk,rtol=1e-8,atol=0):
        raise SystemExit('Disk masses do not match the reference stack.')
    if not np.isclose(data[a.disk:32768,0].sum(),.0625,rtol=1e-7):
        raise SystemExit('Bulge slice does not match the reference stack.')
    a.out.mkdir(parents=True,exist_ok=False)
    print('Building the fixed smooth field. No particles are being evolved yet.',flush=True)
    start=time.perf_counter()
    dq=quad_disk();sq=sphere_quadrature()
    R=np.r_[0,np.geomspace(1e-5,3,a.nr-1)]
    z=ZD*np.sinh(np.linspace(0,np.arcsinh(1/ZD),a.nz))
    rr,zz=np.meshgrid(R,z,indexing='ij')
    points=np.column_stack([rr.ravel(),zz.ravel()])
    phi=disk_field(points,dq).reshape(rr.shape)
    # Tabulate spherical potential separately so its steep central curvature
    # does not degrade the flattened disk interpolation.
    radial=np.r_[0,np.geomspace(1e-6,4,1800)]
    sph_phi=sphere_field(radial,sq)
    sph=CubicSpline(radial,sph_phi,bc_type=((1,0),'not-a-knot'))
    np.savez_compressed(a.out/'smooth_spheroid_field.npz',radius=radial,
                        potential=sph_phi,epsilon=EPS,G=1.,
                        model_id='zeno_reference_bdh_v1',
                        purpose='fixed reference bulge+halo diagnostic; no disk field')
    # Even extensions enforce reflection symmetry of the interpolated potential.
    RR=np.r_[-R[:0:-1],R];ZZ=np.r_[-z[:0:-1],z]
    pp=np.concatenate([phi[:0:-1],phi],axis=0)
    pp=np.concatenate([pp[:,:0:-1],pp],axis=1)
    disk=RectBivariateSpline(RR,ZZ,pp,kx=3,ky=3,s=0)

    def force(x):
        r=np.hypot(x[:,0],x[:,1]);s=np.linalg.norm(x,axis=1)
        if np.any(r>R[-1]) or np.any(abs(x[:,2])>z[-1]) or np.any(s>radial[-1]):
            raise RuntimeError('A tracer left the tabulated field. No extrapolation is allowed.')
        ar=-disk.ev(r,x[:,2],dx=1)
        az=-disk.ev(r,x[:,2],dy=1)
        inv=np.divide(1,r,out=np.zeros_like(r),where=r>0)
        result=np.column_stack([ar*x[:,0]*inv,ar*x[:,1]*inv,az])
        gs=-sph(s,1)
        result+=x*np.divide(gs,s,out=np.zeros_like(s),where=s>0)[:,None]
        return result

    def energy(x,v):
        r=np.hypot(x[:,0],x[:,1]);s=np.linalg.norm(x,axis=1)
        return .5*(v*v).sum(1)+disk.ev(r,x[:,2])+sph(s)

    rng=np.random.default_rng(723)
    probes=np.column_stack([rng.uniform(.005,.65,96),
                            rng.choice([-1,1],96)*rng.uniform(.001,.035,96)])
    actual=disk_field(probes,dq,True)
    ref=disk_field(probes,quad_disk(16,128),True)
    radii=np.linalg.norm(probes,axis=1)
    ss=sphere_field(radii,sq,True)
    ref+=probes*(ss/radii)[:,None]
    actual+=probes*(ss/radii)[:,None]
    interp=force(np.column_stack([probes[:,0],np.zeros(len(probes)),probes[:,1]]))[:,[0,2]]
    relative=np.linalg.norm(interp-ref,axis=1)/np.linalg.norm(ref,axis=1)
    vertical=abs(interp[:,1]-ref[:,1])/np.maximum(abs(ref[:,1]),1e-10)
    quaderror=np.linalg.norm(actual-ref,axis=1)/np.linalg.norm(ref,axis=1)
    checks=dict(force_error_p95=float(np.quantile(relative,.95)),
                force_error_max=float(max(relative)),
                vertical_force_error_max=float(max(vertical)),
                disk_quadrature_error_max=float(max(quaderror)),
                field_seconds=time.perf_counter()-start,
                note='Provisional numerical gates for this control, not published stability thresholds.')
    (a.out/'field_checks.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2),flush=True)
    if max(relative)>.005 or max(vertical)>.01 or max(quaderror)>.002:
        raise SystemExit('Field check failed. Send field_checks.json; no evolution was launched.')

    adjustment = None
    if a.vertical_jeans:
        # Use a fine integration grid independently of the potential grid.
        jz = ZD*np.sinh(np.linspace(0,np.arcsinh(.3/ZD),1000))
        jr,jh = np.meshgrid(R,jz,indexing='ij')
        points3 = np.column_stack([jr.ravel(),np.zeros(jr.size),jh.ravel()])
        gz = -force(points3)[:,2].reshape(jr.shape)
        variance = vertical_jeans_variance(jz,gz,ZD)
        table = RegularGridInterpolator((R,jz),variance,bounds_error=True)
        r = np.hypot(data[:a.disk,1],data[:a.disk,2])
        target = table(np.column_stack([r,abs(data[:a.disk,3])]))
        previous = np.pi*ZD*.1875*144/(2*np.pi)*np.exp(-12*r)
        factors = np.sqrt(target/previous)
        if not np.isfinite(factors).all() or np.any(factors<=0):
            raise RuntimeError('Non-positive Jeans correction; no candidate written')
        # Only the vertical random velocities change. Positions, in-plane
        # velocities, masses and component ordering remain exactly unchanged.
        data[:a.disk,6] *= factors
        from nbody.io import save
        from nbody.particles import Particles
        candidate = Particles(data[:,1:4],data[:,4:7],data[:,0],
                              np.repeat([0,1,2],[24576,8192,32768]),
                              np.zeros(len(data),dtype=int))
        adjustment = dict(method='vertical Jeans second moment in softened reference field',
                          assumption='zero R-z stress; separable initial vertical sampling',
                          factor_min=float(factors.min()),factor_max=float(factors.max()),
                          status='experimental candidate; full DF and live stability unverified')
        save(a.out/'jeans_candidate.npz',candidate,dict(
            backend='Zeno reference, experimental moment adjustment',epsilon=EPS,
            construction_force_model=dict(kernel='plummer',epsilon=EPS,G=1.),
            model_id='zeno_reference_bdh_v1',
            equilibrium_status='unvalidated',adjustment=adjustment,
            source_sha256=hashlib.sha256(a.initial.read_bytes()).hexdigest(),
            frame='unchanged model origin; no recentering',
            note='Spheroid particles unchanged; not a reproduction of the original paper'))
        print(json.dumps(adjustment,indent=2),flush=True)

    histories={}
    ends={}
    for label,dt in [('base',1/512),('half_dt',1/1024)]:
        steps=round(a.duration/dt)
        if not np.isclose(steps*dt,a.duration):
            raise SystemExit('Duration must be a multiple of both timesteps.')
        x=data[:a.disk,1:4].copy();v=data[:a.disk,4:7].copy()
        initial_energy=energy(x,v)
        acc=force(x);rows=[]
        every=round((1/32)/dt)
        for step in range(steps+1):
            if step%every==0 or step==steps:
                en=energy(x,v);row={'time':step*dt,**measure(x,v)}
                drift=abs(en-initial_energy)/np.maximum(abs(initial_energy),1e-12)
                row['orbit_energy_error_p95']=float(np.quantile(drift,.95))
                row['orbit_energy_error_max']=float(max(drift))
                rows.append(row)
                print(f'{label}: t={step*dt:.5f}, height={row["height"]:.7f}',flush=True)
            if step==steps:break
            v+=.5*dt*acc;x+=dt*v
            acc=force(x);v+=.5*dt*acc
        with (a.out/f'{label}_history.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
        np.savez_compressed(a.out/f'{label}_final.npz',x=x,v=v)
        histories[label]=rows
        ends[label]=(x.copy(),v.copy())
    report={'model':'Fixed continuum density; exact Plummer kernel, interpolated potential',
            'initial_file':str(a.initial),'duration':a.duration,'field_checks':checks,
            'initial_velocity_adjustment':adjustment,
            'changes_percent':{},'limits':'No live self-gravity, particle scattering, or growing collective field. Not a stable-galaxy certification.'}
    for label,rows in histories.items():
        report['changes_percent'][label]={k:100*(rows[-1][k]/rows[0][k]-1)
                                         for k in ['height','height90','sigma_z','R50','R90']}
    report['timestep_endpoint_position_rms']=float(np.sqrt(np.mean(np.sum((ends['base'][0]-ends['half_dt'][0])**2,axis=1))))
    report['thickness_timestep_difference_percent']=100*(histories['base'][-1]['height']/histories['half_dt'][-1]['height']-1)
    (a.out/'summary.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    print(f'Results: {a.out.resolve()}',flush=True)


if __name__=='__main__':
    main()
