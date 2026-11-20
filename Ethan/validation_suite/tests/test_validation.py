import unittest,tempfile,json,sys,subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core import load,law,reference_force,_potential,structure,comparison,criterion,camera
CLI=Path(__file__).resolve().parents[1]/'validate.py'

class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        rng=np.random.default_rng(52);self.x=rng.normal(size=(32,3));self.x[:,2]*=.03;self.v=rng.normal(size=(32,3))*.1
        self.m=np.ones(32)/32;self.meta=dict(eps=.01,G=1,evolution_kernel='plummer',dt=.01,theta=.5,solver_sha256='fixture')
    def tearDown(self):self.tmp.cleanup()
    def fixture(self,name='a.npz',**updates):
        d=dict(t=np.array([0.,.1]),x=np.stack([self.x,self.x]),v=np.stack([self.v,self.v]),m=self.m,component=np.zeros(32,dtype=int),galaxy=np.zeros(32,dtype=int),metadata=json.dumps(self.meta));d.update(updates)
        p=self.root/name;np.savez_compressed(p,**d);return p
    def cli(self,*args):return subprocess.run([sys.executable,str(CLI),*map(str,args)],capture_output=True,text=True)
    def test_force_has_no_self_term_and_correct_mass_sign(self):
        x=np.array([[0.,0,0],[2.,0,0]]);m=np.array([2.,3.]);a=reference_force(x,m,[0,1],.5,1)
        self.assertAlmostEqual(a[0,0],6/(4.25**1.5));self.assertTrue(np.allclose(m[0]*a[0]+m[1]*a[1],0))
        self.assertAlmostEqual(_potential(x,m,.5,1),-6/np.sqrt(4.25))
    def test_force_law_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):law(load(self.fixture()),epsilon=.02)
    def test_unknown_force_law_not_silently_assumed(self):
        with self.assertRaises(ValueError):law(load(self.fixture(metadata='{}')))
    def test_rigid_tilt_does_not_become_heating(self):
        B=camera(30,45);a=structure(self.x,self.v,self.m);b=structure(self.x@B.T+17,self.v@B.T-2,self.m)
        for key in ['height','r50','sigma_z','sigma_R']:self.assertAlmostEqual(a[key],b[key],places=10)
    def test_bad_initial_match_is_rejected(self):
        a=load(self.fixture());meta=dict(self.meta,dt=.005)
        b=load(self.fixture('b.npz',x=np.stack([self.x*2,self.x]),metadata=json.dumps(meta)))
        with self.assertRaises(ValueError):comparison(a,b,'timestep')
    def test_changed_epsilon_not_called_convergence(self):
        a=load(self.fixture());b=load(self.fixture('b.npz',metadata=json.dumps(dict(self.meta,eps=.02,dt=.005))))
        with self.assertRaises(ValueError):comparison(a,b,'timestep')
    def test_shared_times_and_halfstep_input_accepted(self):
        a=load(self.fixture());b=load(self.fixture('b.npz',metadata=json.dumps(dict(self.meta,dt=.005))))
        rows,issues,pairs=comparison(a,b,'timestep');self.assertEqual(len(pairs),2);self.assertFalse(issues);self.assertTrue(rows)
    def test_unreviewed_policy_is_not_a_pass(self):
        r=criterion('test',.01,dict(reviewed=False,limits={'x':.1}),'x');self.assertEqual(r['status'],'MEASURED')
        r=criterion('test',.2,dict(reviewed=True,limits={'x':.1}),'x');self.assertEqual(r['status'],'FAIL')
    def test_nonfinite_and_bad_time_rejected(self):
        xx=np.stack([self.x,self.x]);xx[1,0,0]=np.nan
        with self.assertRaises(ValueError):load(self.fixture(x=xx))
        with self.assertRaises(ValueError):load(self.fixture('b.npz',t=np.array([0.,0.])))
    def test_inspect_does_not_certify_missing_evidence(self):
        p=self.fixture();out=self.root/'result';r=self.cli('inspect',p,'--context','isolated','--out',out)
        self.assertEqual(r.returncode,0,r.stderr);report=json.loads((out/'report.json').read_text());self.assertEqual(report['status'],'INCOMPLETE')
        self.assertTrue(any(z['name']=='Energy conservation' and z['status']=='MISSING' for z in report['checks']))
    def test_existing_output_never_modified(self):
        p=self.fixture();out=self.root/'existing';out.mkdir();(out/'marker').write_text('keep')
        r=self.cli('inspect',p,'--out',out);self.assertEqual(r.returncode,3);self.assertEqual(sorted(x.name for x in out.iterdir()),['marker'])
    def test_unknown_external_dynamics_cannot_get_pair_energy_verdict(self):
        p=self.fixture();r=self.cli('inspect',p,'--energy-samples','2','--out',self.root/'energyout');self.assertEqual(r.returncode,3)
    def test_observation_table_no_data_rejected(self):
        p=self.fixture();(self.root/'obs.csv').write_text('id,role,kind,galaxy,x,y,radius,vlos,sigma_v\n')
        setup=dict(time=0.,azimuth=120,elevation=0,length_unit='kpc',velocity_unit='km/s',length_scale=1,velocity_scale=1,sky_offset=[0,0],systemic_velocity=0,constraints_csv='obs.csv')
        sp=self.root/'setup.json';sp.write_text(json.dumps(setup));r=self.cli('observations',p,sp,'--out',self.root/'obsout');self.assertEqual(r.returncode,3)
    def test_synthetic_velocity_holdout_and_fit_kept_separate(self):
        p=self.fixture();(self.root/'obs.csv').write_text('id,role,kind,galaxy,x,y,radius,vlos,sigma_v\nsynthetic,holdout,tail,0,0,0,100,0,1\n')
        sp=self.root/'setup.json';sp.write_text(json.dumps(dict(time=0,azimuth=120,elevation=0,length_unit='kpc',velocity_unit='km/s',length_scale=1,velocity_scale=1,sky_offset=[0,0],systemic_velocity=0,constraints_csv='obs.csv')))
        r=self.cli('observations',p,sp,'--out',self.root/'obsout');self.assertEqual(r.returncode,0,r.stderr)
        report=json.loads((self.root/'obsout/report.json').read_text());self.assertEqual(report['status'],'MEASURED')
    def test_forged_acceleration_provenance_rejected(self):
        p=self.fixture();ap=self.root/'acc.npz';np.savez(ap,indices=np.array([0]),a=np.zeros((1,3)),metadata=json.dumps(dict(frame=0,trajectory_sha256='wrong')))
        r=self.cli('forces',p,ap,'--out',self.root/'forceout');self.assertEqual(r.returncode,3)
    def test_target_profiles_are_hashed_and_do_not_certify_equilibrium(self):
        p=self.fixture();(self.root/'profile.csv').write_text('radius,cumulative_fraction\n0,0\n100,1\n')
        spec=self.root/'targets.json';spec.write_text(json.dumps(dict(components=[dict(galaxy=0,component=0,mass=1,cdf_csv='profile.csv')])))
        r=self.cli('initial',p,spec,'--out',self.root/'initout');self.assertEqual(r.returncode,0,r.stderr)
        report=json.loads((self.root/'initout/report.json').read_text());self.assertEqual(len(report['inputs']),3);self.assertEqual(report['status'],'INCOMPLETE')
    def test_batch_uses_manifest_relative_paths(self):
        self.fixture();manifest=self.root/'jobs.json';manifest.write_text(json.dumps(dict(jobs=[dict(args=['inspect','a.npz','--context','isolated'])])))
        r=self.cli('batch',manifest,'--out',self.root/'batch');self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads((self.root/'batch/batch.json').read_text())['jobs'][0]['status'],'INCOMPLETE')
    def test_reviewed_policy_failure_returns_nonzero(self):
        p=self.fixture(x=np.stack([self.x,self.x*2]));policy=self.root/'policy.json';policy.write_text(json.dumps(dict(reviewed=True,rationale='Synthetic regression fixture',limits={'isolated_height':.01})))
        r=self.cli('inspect',p,'--context','isolated','--policy',policy,'--out',self.root/'failout')
        self.assertEqual(r.returncode,2,r.stderr);self.assertEqual(json.loads((self.root/'failout/report.json').read_text())['status'],'FAIL')

if __name__=='__main__':unittest.main()
