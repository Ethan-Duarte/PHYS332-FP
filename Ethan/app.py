"""N-Baddies V1.4: numerical studies, isolation screens, and 3D encounters."""
import json
import re
from datetime import datetime
from functools import partial
from pathlib import Path
from threading import Event, Lock

import gradio as gr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from nbody.animation import render
from nbody.encounters import parabolic_pair, relative_pair
from nbody.galaxies.generate import galaxy
from nbody.integrate import evolve
from nbody.io import load_trajectory, save_trajectory, load, save
from nbody.galaxies.sources.native import plummer
from nbody.diagnostics import inspect_isolation, encounter_report, structure
from nbody.solver_study import compare_trajectories
from nbody.viewer import write_viewer, viewer_iframe
from nbody.benchmarks import compare

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
RUNS = RESULTS / 'runs'
PRESET = ROOT / 'presets' / 'mice.json'
GALAXIES = RESULTS / 'galaxies'
RUNS.mkdir(parents=True, exist_ok=True)
GALAXIES.mkdir(parents=True, exist_ok=True)
CANCEL = Event()
RUN_LOCK = Lock()


def available_galaxies():
    return sorted((p.name for p in GALAXIES.glob('*.npz')), reverse=True)

def available_runs():
    return sorted((str(p.relative_to(RESULTS)) for p in RESULTS.rglob('*_trajectory.npz')), reverse=True)


def force_solver(name, eps, theta, cancel=None):
    if name == 'Barnes–Hut (Numba)':
        from nbody.gravity.fast_barnes_hut import acceleration
        return partial(acceleration, eps=eps, theta=theta, cancel=cancel)
    if name == 'Barnes–Hut (Python, slow)':
        from nbody.gravity.barnes_hut import acceleration
        return partial(acceleration, eps=eps, theta=theta, cancel=cancel)
    if name == 'Direct (NumPy)':
        from nbody.gravity.direct import acceleration
    else:
        from nbody.gravity.fast_direct import acceleration
    return partial(acceleration, eps=eps)


def separation_plot(times, positions, component, origins):
    centers = []
    for gid in (0, 1):
        mask = (component == 1) & (origins == gid)
        centers.append(positions[:, mask, :].mean(axis=1))
    sep = np.linalg.norm(centers[1] - centers[0], axis=1)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(times, sep, lw=2)
    ax.set(xlabel='Time (code units)', ylabel='Bulge-centroid separation',
           title='Sampled bulge-center separation')
    ax.grid(alpha=.2)
    fig.tight_layout()
    summary = (f'**{len(origins):,} live particles** · **{len(times)} frames** · '
               f't = {times[-1]:.3f} · minimum sampled bulge separation = {sep.min():.3f}')
    return fig, summary


def isolation_plot(times, positions, velocities, masses, component):
    from nbody.particles import Particles
    samples = [structure(Particles(x, v, masses, component, np.zeros(len(masses),dtype=int)))
               for x,v in zip(positions,velocities)]
    fig, ax = plt.subplots(figsize=(6,3))
    t = np.asarray(times)
    ax.plot(t, [s['r50']/samples[0]['r50'] for s in samples], label='Half-mass radius')
    if 'z_rms' in samples[0]:
        ax.plot(t, [s['z_rms']/samples[0]['z_rms'] for s in samples], label='Disk thickness')
    ax.axhline(1,ls=':',color='gray')
    ax.set(xlabel='Time (code units)',ylabel='Ratio to initial',title='Isolated structural evolution')
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout()
    return fig, f'**{len(masses):,} particles** · isolated evolution to t = {t[-1]:.3f}'


def visualizations(trajectory, folder, lim, show_halo, max_points, fps, export):
    t, x, v, m, comp, ids, meta = load_trajectory(trajectory)
    html = write_viewer(folder / 'viewer_3d.html', t, x, comp, ids,
                        limit=lim, show_halo=show_halo, max_points=int(max_points))
    viewer = viewer_iframe(html)
    curve, summary = (isolation_plot(t,x,v,m,comp) if meta.get('setup')=='Isolated galaxy'
                      else separation_plot(t, x, comp, ids))
    gif = mp4 = None
    # Optional 2D video for portable sharing. Interactive 3D is the primary viewer.
    if export in ('GIF', 'Both'):
        gif = folder / 'animation.gif'
        render(t, x, comp, ids, gif, fps=int(fps), lim=lim, view='xy', show_halo=show_halo)
    if export in ('MP4', 'Both'):
        mp4 = folder / 'animation.mp4'
        render(t, x, comp, ids, mp4, fps=int(fps), lim=lim, view='xy', show_halo=show_halo)
    return viewer, curve, str(gif) if gif else None, str(mp4) if mp4 else None, str(html), summary


FIELD_NAMES, FIELD_CONTROLS = [], []
def field(name, control):
    FIELD_NAMES.append(name)
    FIELD_CONTROLS.append(control)
    return control


def _status(message, running=True):
    return (gr.update(), gr.update(), gr.update(), gr.update(), gr.update(),
            gr.update(), gr.update(), message, gr.update(interactive=not running))


def cancel_active():
    CANCEL.set()
    return 'Cancellation requested. The current force evaluation must finish before the run stops.'


def simulate(*values, progress=gr.Progress()):
    p = dict(zip(FIELD_NAMES, values))
    if not RUN_LOCK.acquire(blocking=False):
        yield _status('Another simulation is already running; cancel it before starting a new one.', False)
        return
    CANCEL.clear()
    snapshots=[]
    try:
        isolated = p['setup']=='Isolated galaxy'
        if p['setup']=='Mice encounter' and p['pericenter']>=p['separation']:
            raise gr.Error('Nominal pericenter must be smaller than initial separation.')
        factor=1 if isolated else 2
        n=factor*(int(p['disk'])+int(p['bulge'])+int(p['halo']))
        if p['solver']=='Direct (NumPy)' and n>1500:
            raise gr.Error('NumPy direct uses N² temporary memory. Choose compiled direct or Barnes–Hut.')
        tag=re.sub(r'[^a-zA-Z0-9_-]+','_',p['name']).strip('_') or 'experiment'
        folder=RUNS/f'{tag}_{datetime.now():%Y%m%d_%H%M%S_%f}'
        folder.mkdir(parents=True,exist_ok=True)
        config=json.loads(PRESET.read_text())
        for key,section,parameter in (
            ('disk_mass','disk','mass'),('disk_scale','disk','scale'),
            ('disk_height','disk','height'),('bulge_mass','bulge','mass'),
            ('bulge_scale','bulge','scale'),('halo_mass','halo','mass'),
            ('halo_scale','halo','scale')):
            config[section][parameter]=float(p[key])
        config['seed']=int(p['seed'])
        config['particles']={k:int(p[k]) for k in ('disk','bulge','halo')}
        config['Q']=float(p['Q'])
        orientations=[[p['inc_a'],p['node_a']],[p['inc_b'],p['node_b']]]
        config['encounter']=dict(setup=p['setup'],separation=p['separation'],
                                 pericenter=p['pericenter'],orientations=orientations,
                                 relative_position=[p[k] for k in ('rx','ry','rz')],
                                 relative_velocity=[p[k] for k in ('vx','vy','vz')])
        backend=p['generator']
        saved_config=None
        yield _status('Constructing initial conditions…')
        if p['model']=='Plummer (analytic DF)':
            if not isolated: raise gr.Error('Plummer is currently offered for isolated tests only.')
            initial=plummer(sum(config['particles'].values()),
                            mass=sum(config[k]['mass'] for k in ('disk','bulge','halo')),
                            a=float(p['halo_scale']),seed=int(p['seed']))
            backend='analytic Plummer DF'
        else:
            if backend=='Saved snapshot':
                if not p['snapshot']: raise gr.Error('Choose a saved particle snapshot first.')
                g1, source=load(GALAXIES/Path(p['snapshot']).name)
                backend='saved: '+source.get('backend','unknown')
                saved_config=source.get('config')
            else:
                backend={'Native preview (not equilibrium)':'preview',
                         'AGAMA (equilibrium candidate)':'agama'}[backend]
                g1=galaxy(config,backend)
                if backend=='agama':
                    filename=f'agama_{datetime.now():%Y%m%d_%H%M%S}.npz'
                    save(GALAXIES/filename,g1,dict(backend='agama',config=config,
                         caveat='Candidate: requires isolated tests and resolution study'))
            if isolated:
                initial=g1
            else:
                if backend=='preview':
                    g2=galaxy({**config,'seed':config['seed']+1},'preview')
                else:
                    # Reuse the saved equilibrium realization (not independent noise).
                    g2=g1
                if p['setup']=='Custom encounter':
                    initial=relative_pair(g1,g2,config['encounter']['relative_position'],
                                         config['encounter']['relative_velocity'],orientations)
                else:
                    initial=parabolic_pair(g1,g2,p['separation'],p['pericenter'],orientations)
        force=force_solver(p['solver'],float(p['eps']),float(p['theta']),CANCEL)
        steps,save_every=int(p['steps']),int(p['save_every'])
        progress(0,desc='Integrating gravitational forces')
        yield _status(f'Integrating **{p["setup"]}** · {len(initial.m):,} particles · **{p["solver"]}**…')
        try:
            for time,particles in evolve(initial,force,float(p['dt']),steps,save_every,cancelled=CANCEL):
                snapshots.append((time,particles))
                progress(min(time/(steps*float(p['dt'])),.98),desc=f't = {time:.3f}')
                yield _status(f'{p["setup"]} · {len(initial.m):,} particles · '
                              f't={time:.3f}/{steps*float(p["dt"]):.3f} · '
                              f'{len(snapshots)} frames saved')
        except InterruptedError:
            CANCEL.set()
        if not snapshots:
            yield _status('Cancelled before any snapshot was available.',False)
            return
        metadata=dict(config=config,backend=backend,model=p['model'],setup=p['setup'],
                      solver=p['solver'],eps=float(p['eps']),theta=float(p['theta']),
                      dt=float(p['dt']),steps=steps,save_every=save_every,
                      cancelled=CANCEL.is_set(),
                      caveat='Native preview is NOT equilibrium; AGAMA output requires isolated validation')
        trajectory=folder/'run_trajectory.npz'
        save(folder/'initial_particles.npz',initial,metadata)
        save_trajectory(trajectory,snapshots,metadata)
        (folder/'config.json').write_text(json.dumps(metadata,indent=2))
        progress(.99,desc='Computing run diagnostics')
        if isolated:
            report,_=inspect_isolation(snapshots,float(p['eps']),folder,backend)
            if p['model']=='Disk–bulge–halo' and (backend in ('preview','agama') or saved_config):
                from nbody.diagnostics import profile_checks
                profile_checks(initial,saved_config or config,folder,
                               'preview' if 'preview' in backend else 'agama')
            scientific=f'Isolation **{report["verdict"]}** · '
            scientific+=f'Δr50 = {100*report["relative_changes"]["r50"]:+.2f}% · '
            scientific+=f'ΔE/|E0| = {report["relative_energy_change"]:+.3g}' if report['relative_energy_change'] is not None else ''
            scientific+='\n\n**This screening alone does not validate the galaxy.**'
            report_file='isolation_report.json'
        else:
            report=encounter_report(snapshots,float(p['eps']),folder)
            scientific=f'Minimum sampled bulge separation: {report["min_bulge_separation"]:.4f} at t={report["time_min_separation"]:.3f}.'
            report_file='encounter_report.json'
        yield _status('Preparing the interactive 3D viewer…')
        view,curve,gif,mp4,html,summary=visualizations(trajectory,folder,float(p['lim']),
                            bool(p['show_halo']),int(p['max_points']),int(p['fps']),
                            'None' if CANCEL.is_set() else p['export'])
        label='**Cancelled; partial trajectory saved.**' if CANCEL.is_set() else '**Simulation complete.**'
        progress(1,desc='Saved')
        yield (view,curve,gif,mp4,str(trajectory),str(folder/'config.json'),html,
               f'{label}\n\n{summary}\n\n{scientific}\n\n'
               f'Diagnostics: `{report_file}` and CSV in `results/runs/{folder.name}/`\n\n'
               f'Backend: **{backend}** · Saved to `results/runs/{folder.name}/`',
               gr.update(interactive=True))
    except Exception as exc:
        yield _status(f'**Run failed:** `{type(exc).__name__}: {exc}`. See launcher window for details.',False)
    finally:
        RUN_LOCK.release()
        CANCEL.clear()

def replay_saved(selected, lim, show_halo, max_points, fps, export, progress=gr.Progress()):
    if not selected:
        raise gr.Error('Select a saved trajectory.')
    path = (RESULTS / selected).resolve()
    if not path.is_relative_to(RESULTS.resolve()):
        raise gr.Error('Invalid trajectory selection.')
    folder = RUNS / ('replay_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    folder.mkdir(parents=True, exist_ok=True)
    progress(.2, desc='Rebuilding 3D frames')
    view, curve, gif, mp4, html, summary = visualizations(path, folder, lim, show_halo,
                                                           max_points, fps, export)
    return view, curve, gif, mp4, html, summary


def size_text(disk, bulge, halo, steps, solver, setup):
    n = (1 if setup=='Isolated galaxy' else 2)*(disk+bulge+halo)
    note = 'Pure-Python Barnes–Hut can be very slow; use the compiled version.' if 'Python' in solver else ''
    return f'**{n:,} live particles** · {int(steps):,} integration steps. {note}'


def run_benchmark(n, theta, eps):
    from nbody.galaxies.generate import galaxy
    config = json.loads(PRESET.read_text())
    d,b,h = int(n*.60/2), max(10,int(n*.08/2)), int(n*.32/2)
    config['particles'] = dict(disk=d,bulge=b,halo=h)
    a = galaxy(config,'preview')
    z = galaxy({**config,'seed':config['seed']+1},'preview')
    p = parabolic_pair(a,z,1.,.25,config['encounter']['orientations'])
    return compare(p.x,p.m,eps,theta), f'Benchmark: **{len(p.m)} particles**, same softened force law. Force errors are relative to compiled direct summation.'


def run_full_study(n,steps,theta,eps,progress=gr.Progress()):
    if not RUN_LOCK.acquire(blocking=False):
        raise gr.Error('An encounter is already running.')
    try:
        cfg=json.loads(PRESET.read_text())
        per=max(60,int(n)//2)
        cfg['particles']=dict(disk=int(per*.6),bulge=max(10,int(per*.08)),halo=int(per*.32))
        one=galaxy(cfg,'preview')
        two=galaxy({**cfg,'seed':cfg['seed']+1},'preview')
        initial=parabolic_pair(one,two,1.,.25,cfg['encounter']['orientations'])
        folder=RUNS/f'solver_comparison_{datetime.now():%Y%m%d_%H%M%S_%f}'
        results,plot=compare_trajectories(initial,folder,dt=.006,steps=int(steps),
                                          save_every=max(1,int(steps)//10),eps=float(eps),theta=float(theta))
        table=[[key,round(v['seconds'],4),round(100*v['force_p95_error'],3),
                round(v['rms_position_error_final'],7),round(100*v['energy_relative_change'],4)]
               for key,v in results.items()]
        return table,str(plot),str(folder/'solver_comparison.json'),(
            f'Compared **identical initial particles** with {len(initial.m)} live particles. '
            'Force errors and trajectories are relative to compiled direct summation. '
            '**These are preview galaxies, not validated equilibrium models.**')
    finally:
        RUN_LOCK.release()


with gr.Blocks(title='N-Baddies | 3D Galaxy Lab') as app:
    gr.Markdown('# N-Baddies V1.4 · Galaxy Dynamics Lab\nGenerate isolated systems, investigate a Barnes-inspired Mice encounter, compare force solvers, and replay results in 3D.')
    gr.Markdown('**Scientific status:** Native disk velocities remain approximate. AGAMA is an optional DF-based equilibrium candidate, not yet run-verified here. Isolation screening and convergence are required before claiming a physically validated Mice reconstruction.')
    with gr.Tab('Run a simulation'):
        with gr.Row():
            with gr.Column(scale=1, min_width=310):
                field('name', gr.Textbox(label='Run name', value='mice_test'))
                field('setup', gr.Radio(['Isolated galaxy','Mice encounter','Custom encounter'], value='Mice encounter', label='Experiment'))
                with gr.Accordion('Initial-condition model',open=True):
                    field('model', gr.Radio(['Disk–bulge–halo','Plummer (analytic DF)'],value='Disk–bulge–halo',label='Galaxy model (Plummer for isolated tests)'))
                    field('generator', gr.Radio(['Native preview (not equilibrium)','AGAMA (equilibrium candidate)','Saved snapshot'],value='Native preview (not equilibrium)',label='Disk generator'))
                    field('snapshot', gr.Dropdown(choices=available_galaxies(),value=None,label='Saved galaxy .npz (only used for Saved snapshot)'))
                    refresh_galaxies=gr.Button('Refresh saved galaxies',size='sm')
                    field('Q',gr.Slider(1,3,value=1.5,step=.1,label='Target Toomre Q near 1 disk scale (AGAMA only)'))
                with gr.Accordion('Galaxy particles and seed', open=True):
                    field('disk', gr.Slider(50, 5000, value=800, step=50, label='Disk particles / galaxy'))
                    field('bulge', gr.Slider(10, 500, value=100, step=10, label='Bulge particles / galaxy'))
                    field('halo', gr.Slider(50, 3000, value=350, step=50, label='Dark-matter particles / galaxy'))
                    field('seed', gr.Number(label='Seed', value=42, precision=0))
                with gr.Accordion('Galaxy structure · G = 1 code units', open=False):
                    field('disk_mass', gr.Number(label='Disk mass', value=.1875, minimum=.0001))
                    field('bulge_mass', gr.Number(label='Bulge mass', value=.0625, minimum=.0001))
                    field('halo_mass', gr.Number(label='Halo mass', value=1.0, minimum=.0001))
                    field('disk_scale', gr.Number(label='Disk radial scale', value=.0833333333, minimum=.0001))
                    field('disk_height', gr.Number(label='Disk thickness scale', value=.005, minimum=.0001))
                    field('bulge_scale', gr.Number(label='Bulge scale', value=.04168, minimum=.0001))
                    field('halo_scale', gr.Number(label='Halo scale', value=.1, minimum=.0001))
                with gr.Accordion('Encounter geometry', open=True):
                    with gr.Group(visible=True) as mice_group:
                        field('separation', gr.Slider(.3, 4, value=1, step=.05, label='Initial separation'))
                        field('pericenter', gr.Slider(.02, 1, value=.25, step=.01, label='Nominal pericenter'))
                    with gr.Group(visible=False) as custom_group:
                        gr.Markdown('**Galaxy B relative to A**, in code units. Set initial displacement and relative velocity; the system is placed at its barycenter.')
                        with gr.Row():
                            field('rx', gr.Number(value=1, label='Δx'))
                            field('ry', gr.Number(value=0, label='Δy'))
                            field('rz', gr.Number(value=0, label='Δz'))
                        with gr.Row():
                            field('vx', gr.Number(value=-1.3, label='Δvx'))
                            field('vy', gr.Number(value=.7, label='Δvy'))
                            field('vz', gr.Number(value=0, label='Δvz'))
                    with gr.Row():
                        field('inc_a', gr.Slider(-180,180,value=25,step=5,label='A inclination (°)'))
                        field('node_a', gr.Slider(-180,180,value=-30,step=5,label='A node (°)'))
                    with gr.Row():
                        field('inc_b', gr.Slider(-180,180,value=40,step=5,label='B inclination (°)'))
                        field('node_b', gr.Slider(-180,180,value=60,step=5,label='B node (°)'))
                with gr.Accordion('Gravity and integration', open=True):
                    field('solver', gr.Radio(['Fast direct (Numba)','Barnes–Hut (Numba)','Direct (NumPy)'], value='Fast direct (Numba)', label='Gravity solver'))
                    field('eps', gr.Number(value=.0125, minimum=.00001,label='Softening ε (Barnes reference = .0125)'))
                    field('theta', gr.Slider(.2,1,value=.5,step=.05,label='Barnes–Hut opening θ'))
                    field('dt', gr.Number(value=.006,minimum=.00001,label='Timestep Δt (Barnes reference = 1/512 ≈ .001953)'))
                    field('steps', gr.Slider(20,1200,value=240,step=20,label='Timesteps'))
                    field('save_every', gr.Slider(1,50,value=4,step=1,label='Save every N steps'))
                    reference_button = gr.Button('Apply Barnes (2004) softening and timestep',size='sm')
                with gr.Accordion('3D viewer and exports', open=True):
                    field('lim', gr.Slider(.5,5,value=1.5,step=.1,label='3D axis half-width'))
                    field('max_points', gr.Slider(300,6000,value=3000,step=100,label='Maximum displayed particles'))
                    field('show_halo', gr.Checkbox(value=False,label='Display dark matter (always simulated)'))
                    field('fps', gr.Slider(5,40,value=15,step=1,label='Export FPS'))
                    field('export', gr.Radio(['None','GIF','MP4','Both'],value='None',label='Optional 2D video export'))
                work = gr.Markdown('')
                with gr.Row():
                    run_button = gr.Button('Run simulation', variant='primary', scale=3)
                    stop_button = gr.Button('Cancel run', variant='stop', scale=1)
            with gr.Column(scale=2, min_width=530):
                status = gr.Markdown('Configure a run, then click **Run simulation**. You can rotate and zoom the 3D viewer with your mouse.')
                viewer = gr.HTML(label='Interactive 3D viewer · drag to rotate; use the controls below')
                curve = gr.Plot(label='Encounter separation / isolated structure')
                with gr.Row():
                    gif_file = gr.File(label='Optional GIF')
                    mp4_file = gr.File(label='Optional MP4')
                with gr.Row():
                    trajectory_file = gr.File(label='Full particle trajectory')
                    config_file = gr.File(label='Run configuration')
                    html_file = gr.File(label='Standalone interactive 3D viewer (.html)')
    with gr.Tab('Replay a saved run'):
        gr.Markdown('Change 3D camera bounds, displayed particle count, and halo visibility **without rerunning gravity**.')
        with gr.Row():
            saved = gr.Dropdown(choices=available_runs(), value=available_runs()[0] if available_runs() else None, label='Saved trajectory', scale=4)
            refresh = gr.Button('Refresh', scale=1)
        with gr.Row():
            replay_lim = gr.Slider(.5,5,value=1.5,step=.1,label='3D axis half-width')
            replay_halo = gr.Checkbox(False,label='Show halo')
            replay_max = gr.Slider(300,6000,value=3000,step=100,label='Max viewer particles')
        with gr.Row():
            replay_fps = gr.Slider(5,40,value=15,step=1,label='Export FPS')
            replay_export = gr.Radio(['None','GIF','MP4','Both'],value='None',label='Optional 2D video')
            replay_button = gr.Button('Open trajectory in 3D',variant='primary')
        replay_summary = gr.Markdown()
        replay_view = gr.HTML(label='Interactive 3D trajectory')
        replay_curve = gr.Plot(label='Encounter separation / isolated structure')
        with gr.Row():
            replay_html = gr.File(label='Download interactive 3D HTML')
            replay_gif = gr.File(label='Optional GIF')
            replay_mp4 = gr.File(label='Optional MP4')
    with gr.Tab('Compare gravitational solvers'):
        gr.Markdown('Compare accelerations at a **single initial Mice-style snapshot**. Fast direct summation is the reference. These timings exclude Numba compilation, and are not yet an end-to-end merger convergence test.')
        with gr.Row():
            bench_n = gr.Slider(100,2000,value=600,step=100,label='Total live particles')
            bench_theta = gr.Slider(.2,1,value=.5,step=.05,label='Barnes–Hut θ')
            bench_eps = gr.Number(value=.01,label='Softening ε')
            bench_button = gr.Button('Compare forces',variant='primary')
        bench_result = gr.Dataframe(headers=['Solver','Seconds','Median error %','95th percentile error %','Max absolute error'],interactive=False)
        bench_note = gr.Markdown()
    with gr.Tab('Full-trajectory solver study'):
        gr.Markdown('Compare **identical initial particles** using Fast direct, Barnes–Hut, and NumPy direct (small N). Reports include force errors, complete run times, energy drift, and projected galaxy-core trajectories.')
        with gr.Row():
            full_n=gr.Slider(120,1200,value=400,step=40,label='Total particles')
            full_steps=gr.Slider(10,150,value=30,step=10,label='Integration steps')
            full_theta=gr.Slider(.2,1,value=.5,step=.05,label='BH opening θ')
            full_eps=gr.Number(value=.0125,label='Softening ε')
        full_run=gr.Button('Run matched solver study',variant='primary')
        full_table=gr.Dataframe(headers=['Solver','Total seconds','P95 force error (%)','Final RMS position diff','Energy drift (%)'],interactive=False)
        full_plot=gr.Image(label='Separation curves and runtime')
        full_json=gr.File(label='Download complete comparison JSON')
        full_note=gr.Markdown()
    with gr.Tab('Guide and scientific status'):
        gr.Markdown('''
### Basic use
1. Start with **Mice encounter** and **Fast direct (Numba)**, then try **Barnes–Hut (Numba)**.
2. Click **Run simulation**, then drag the 3D plot to orbit the camera. Use the Play/Pause button and timeline beneath it.
3. To stop a long run, click **Cancel run**. Cancellation is checked between integration steps (and after each compiled force evaluation). A partial trajectory will be saved when available.
4. Use **Replay a saved run** to regenerate the interactive view, without re-integrating.
5. Use **Custom encounter** to specify the initial position and velocity of galaxy B relative to galaxy A. Both galaxies are centered around the barycenter.
6. All runs preserve numerical settings and full particle trajectories under `results/runs/`.

### Which solver?
- **Direct (NumPy):** exact softened all-pairs force within floating-point precision; uses O(N²) memory.
- **Fast direct (Numba):** same pairwise force law, compiled loops, O(N) extra memory.
- **Barnes–Hut (Numba):** approximate hierarchical forces with compiled traversal.
- **Barnes–Hut (Python, slow):** original reference code; use only for tiny diagnostic cases.

### Physics limitations
Native disk–bulge–halo velocities are heuristic. AGAMA can be selected explicitly if installed, but its code-unit normalization and equilibrium particle realizations need testing on your machine. The published Mice geometry is an approximate initialization, not a completed reconstruction.
This GUI is intended for exploring and benchmarking the numerical system; a publication-quality Mice reconstruction requires equilibrium initial conditions, convergence testing, and observational comparisons.
''')

    setup = FIELD_CONTROLS[FIELD_NAMES.index('setup')]
    setup.change(lambda mode: (gr.update(visible=mode=='Mice encounter'),
                               gr.update(visible=mode=='Custom encounter')),
                 inputs=setup, outputs=[mice_group,custom_group])
    refresh_galaxies.click(lambda:gr.update(choices=available_galaxies()),outputs=FIELD_CONTROLS[FIELD_NAMES.index('snapshot')])
    full_run.click(run_full_study,inputs=[full_n,full_steps,full_theta,full_eps],outputs=[full_table,full_plot,full_json,full_note])
    reference_button.click(lambda:(.0125,1/512), outputs=[FIELD_CONTROLS[FIELD_NAMES.index('eps')],FIELD_CONTROLS[FIELD_NAMES.index('dt')]])
    solver = FIELD_CONTROLS[FIELD_NAMES.index('solver')]
    disk = FIELD_CONTROLS[FIELD_NAMES.index('disk')]
    bulge = FIELD_CONTROLS[FIELD_NAMES.index('bulge')]
    halo = FIELD_CONTROLS[FIELD_NAMES.index('halo')]
    steps = FIELD_CONTROLS[FIELD_NAMES.index('steps')]
    for control in (disk,bulge,halo,steps,solver,setup):
        control.change(size_text, inputs=[disk,bulge,halo,steps,solver,setup], outputs=work)
    app.load(size_text, inputs=[disk,bulge,halo,steps,solver,setup], outputs=work)
    run_button.click(simulate, inputs=FIELD_CONTROLS,
                     outputs=[viewer,curve,gif_file,mp4_file,trajectory_file,config_file,
                              html_file,status,run_button], concurrency_limit=1)
    stop_button.click(cancel_active, outputs=status, queue=False)
    refresh.click(lambda: gr.update(choices=available_runs()), outputs=saved)
    bench_button.click(run_benchmark,inputs=[bench_n,bench_theta,bench_eps],outputs=[bench_result,bench_note])
    replay_button.click(replay_saved,
                        inputs=[saved,replay_lim,replay_halo,replay_max,replay_fps,replay_export],
                        outputs=[replay_view,replay_curve,replay_gif,replay_mp4,replay_html,replay_summary])


if __name__=='__main__':
    app.queue(default_concurrency_limit=2).launch(
        inbrowser=True, share=False, server_name='127.0.0.1',
        allowed_paths=[str(RESULTS)])
