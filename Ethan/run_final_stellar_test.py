"""Run the one planned stellar-resolution comparison and analyze its output."""
import argparse,hashlib,json,subprocess,sys,zipfile
from pathlib import Path
import numpy as np

EXPECTED='db34d50f26cc83e0a1149c9e2935d8f13ff620cc248b87d4fb004b6fdda41284'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path)
    p.add_argument('--zeno-root',default='/home/whisky/nbody-tools/zeno-reference')
    p.add_argument('--check-only',action='store_true')
    args=p.parse_args();here=Path(__file__).resolve().parent
    initial=here/'final_stellar_initial.npz'
    if hashlib.sha256(initial.read_bytes()).hexdigest()!=EXPECTED:p.error('Prepared initial-state checksum differs; extract the original package again')
    with np.load(initial,allow_pickle=False) as d:
        c=d['component'];m=d['m'];meta=json.loads(d['metadata'].item())
        if [int(np.sum(c==k)) for k in range(3)]!=[12288,4096,32768]:p.error('Wrong component counts')
        if len(np.unique(d['galaxy']))!=1 or d['x'].shape!=(1,49152,3):p.error('Expected the prepared isolated initial state')
        if not np.allclose([m[c==k].sum() for k in range(3)],[.1875,.0625,1.0000047683731454],rtol=1e-12):p.error('Component masses differ')
        if not np.allclose([meta['eps'],meta['dt'],meta['theta']],[.009375,1/512,.7],rtol=0,atol=1e-14) or meta['usequad'] is not True or meta['save_every']!=8:p.error('Force/integration settings differ')
        if not np.isfinite(d['x']).all() or not np.isfinite(d['v']).all() or not np.isfinite(m).all() or np.any(m<=0):p.error('Invalid phase-space state')
    print('Preflight passed: 12,288 disk + 4,096 bulge + 32,768 halo.',flush=True)
    print('Fixed: epsilon=0.009375, dt=1/512, theta=0.7, quadrupoles on; duration=1.25 (640 steps), save every 8.',flush=True)
    if args.check_only:return
    if args.out is None:p.error('--out is required unless using --check-only')
    out=args.out.resolve()
    if out.exists():p.error('Output folder exists; choose a fresh folder')
    # The existing project supplies Zeno bridge, WSL helpers and diagnostics.
    for filename in ['zeno_mice_run.py','zeno_snapshot_bridge.c']:
        if not (here/filename).exists():p.error(f'Place these package files in the existing Ethan folder; missing {filename}')
    command=[sys.executable,str(here/'mice_validation_run.py'),str(initial),'--mode','isolate',
             '--duration','1.25','--dt','0.001953125','--zeno-root',args.zeno_root,'--out',str(out)]
    subprocess.run(command,check=True)
    analysis=out/'analysis'
    subprocess.run([sys.executable,str(here/'analyze_isolated_galaxy.py'),'--trajectory',str(out/'run_trajectory.npz'),'--out',str(analysis)],check=True)
    archive=out/'analysis_to_share.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for file in sorted(analysis.iterdir()):
            if file.is_file():z.write(file,'analysis/'+file.name)
        for filename in ['isolation_summary.json','isolation_history.csv','config.json']:
            z.write(out/filename,filename)
    print('\nFinished. Send this small file:',archive,flush=True)
    print('The large trajectory remains in the run folder for animation work.',flush=True)


if __name__=='__main__':main()
