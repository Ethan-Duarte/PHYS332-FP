"""N-Baddies local browser interface. Launch with: python app.py"""
import json
import re
from datetime import datetime
from functools import partial
from pathlib import Path

import gradio as gr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nbody.animation import render
from nbody.encounters import parabolic_pair
from nbody.galaxies.generate import galaxy
from nbody.integrate import evolve
from nbody.io import load_trajectory, save_trajectory

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RUNS = RESULTS / "runs"
PRESET = ROOT / "presets" / "mice.json"
RUNS.mkdir(parents=True, exist_ok=True)


def available_runs():
    return sorted((str(p.relative_to(RESULTS)) for p in RESULTS.rglob("*_trajectory.npz")), reverse=True)


def force_solver(name, eps, theta):
    if name == "Fast direct (Numba)":
        from nbody.gravity.fast_direct import acceleration
    elif name == "Barnes–Hut":
        from nbody.gravity.barnes_hut import acceleration
        return partial(acceleration, eps=eps, theta=theta)
    else:
        from nbody.gravity.direct import acceleration
    return partial(acceleration, eps=eps)


def separation_plot(times, positions, component, origins):
    # With low bulge resolution these are noisy sampled bulge centroids, not exact core orbits.
    centers = []
    for gid in (0, 1):
        mask = (component == 1) & (origins == gid)
        centers.append(positions[:, mask, :].mean(axis=1))
    sep = np.linalg.norm(centers[1] - centers[0], axis=1)
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.plot(times, sep, lw=2)
    ax.set(xlabel="Time (code units)", ylabel="Sampled bulge-centroid separation",
           title="Separation during encounter")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    return fig, f"**{len(origins):,} live particles** · **{len(times)} saved frames** · t = {times[-1]:.3f} code units · min sampled bulge separation = {sep.min():.3f}"


def render_saved(path, view, lim, fps, show_halo, export, folder):
    t, x, v, m, comp, ids, meta = load_trajectory(path)
    gif = mp4 = None
    if export in ("GIF", "Both"):
        gif = folder / "animation.gif"
        render(t, x, comp, ids, gif, fps=fps, lim=lim, view=view, show_halo=show_halo)
    if export in ("MP4", "Both"):
        mp4 = folder / "animation.mp4"
        render(t, x, comp, ids, mp4, fps=fps, lim=lim, view=view, show_halo=show_halo)
    chart, stats = separation_plot(t, x, comp, ids)
    return (str(gif) if gif else None, str(mp4) if mp4 else None,
            str(gif or mp4), chart, stats)


# The controls below are collected as a dict to keep the experiment function short.
FIELD_NAMES = []
FIELD_CONTROLS = []
def field(name, component):
    FIELD_NAMES.append(name)
    FIELD_CONTROLS.append(component)
    return component


def simulate(*values, progress=gr.Progress()):
    p = dict(zip(FIELD_NAMES, values))
    if p["pericenter"] >= p["separation"]:
        raise gr.Error("Pericenter must be smaller than initial separation.")

    tag = re.sub(r"[^a-zA-Z0-9_-]+", "_", p["name"]).strip("_") or "encounter"
    folder = RUNS / f"{tag}_{datetime.now():%Y%m%d_%H%M%S_%f}"
    folder.mkdir(parents=True, exist_ok=True)

    config = json.loads(PRESET.read_text())
    for key, section, parameter in (
        ("disk_mass", "disk", "mass"), ("disk_scale", "disk", "scale"),
        ("disk_height", "disk", "height"), ("bulge_mass", "bulge", "mass"),
        ("bulge_scale", "bulge", "scale"), ("halo_mass", "halo", "mass"),
        ("halo_scale", "halo", "scale"),
    ):
        config[section][parameter] = p[key]
    config["seed"] = int(p["seed"])
    config["particles"] = {k: int(p[k]) for k in ("disk", "bulge", "halo")}
    config["encounter"] = {
        "separation": p["separation"], "pericenter": p["pericenter"],
        "orientations": [[p["inc_a"], p["node_a"]], [p["inc_b"], p["node_b"]]],
    }
    backend = "agama" if p["backend"].startswith("AGAMA") else "preview"
    progress(0, desc="Generating the two galaxies")
    g1 = galaxy(config, backend)
    g2 = galaxy({**config, "seed": config["seed"] + 1}, backend)
    initial = parabolic_pair(g1, g2, p["separation"], p["pericenter"], config["encounter"]["orientations"])
    force = force_solver(p["solver"], p["eps"], p["theta"])
    steps = int(p["steps"])
    save_every = int(p["save_every"])

    snapshots = []
    for time, particles in evolve(initial, force, p["dt"], steps, save_every):
        snapshots.append((time, particles))
        progress(min(time / (steps * p["dt"]), 0.98), desc=f"Evolving particles: t = {time:.3f}")

    metadata = {
        "config": config, "backend": backend, "solver": p["solver"],
        "eps": p["eps"], "theta": p["theta"] if p["solver"] == "Barnes–Hut" else None,
        "dt": p["dt"], "steps": steps, "save_every": save_every,
        "caveat": "Native preview kinematics are not equilibrium validated" if backend == "preview"
                  else "AGAMA integration and isolated stability still require verification",
    }
    trajectory = folder / "run_trajectory.npz"
    save_trajectory(trajectory, snapshots, metadata)
    (folder / "config.json").write_text(json.dumps(metadata, indent=2))
    progress(0.98, desc="Rendering the animation")
    gif, mp4, downloadable, chart, summary = render_saved(
        trajectory, p["view"], p["lim"], int(p["fps"]), p["show_halo"],
        p["export"], folder)
    progress(1, desc="Saved")
    summary += f"\n\n**Saved to:** `results/runs/{folder.name}/`"
    return gif, mp4, str(trajectory), str(folder / "config.json"), chart, summary


def replay_saved(selected, view, lim, fps, show_halo, export, progress=gr.Progress()):
    if not selected:
        raise gr.Error("Choose a saved trajectory first.")
    path = (RESULTS / selected).resolve()
    if not path.is_relative_to(RESULTS.resolve()):
        raise gr.Error("Invalid trajectory selection.")
    folder = RUNS / ("replay_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    folder.mkdir(parents=True, exist_ok=True)
    progress(0.1, desc="Loading saved particles")
    gif, mp4, download, fig, summary = render_saved(path, view, lim, int(fps), show_halo, export, folder)
    return gif, mp4, download, fig, summary


with gr.Blocks(title="N-Baddies | Galaxy Simulator") as app:
    gr.Markdown("# N-Baddies · Galaxy Encounter Lab\nRun, animate, and replay gravitational N-body experiments without editing Python files.")
    gr.Markdown("**Scientific status:** Native disk–bulge–halo velocities are approximate. These animations are for development and exploration, **not yet validated reconstructions of the Mice galaxies**.")

    with gr.Tab("Run an encounter"):
        with gr.Row():
            with gr.Column(scale=1):
                field("name", gr.Textbox(label="Run name", value="mice_test"))
                field("backend", gr.Dropdown(label="Galaxy generator", choices=["Native preview (not equilibrium)", "AGAMA (experimental; requires separate install)"], value="Native preview (not equilibrium)"))
                with gr.Accordion("Particles and seed", open=True):
                    field("disk", gr.Slider(50, 3000, value=800, step=50, label="Disk particles / galaxy"))
                    field("bulge", gr.Slider(10, 500, value=100, step=10, label="Bulge particles / galaxy"))
                    field("halo", gr.Slider(50, 2000, value=350, step=50, label="Halo particles / galaxy"))
                    field("seed", gr.Number(label="Random seed", value=42, precision=0))
                with gr.Accordion("Galaxy structure (Barnes-inspired code units)", open=False):
                    field("disk_mass", gr.Number(label="Disk mass", value=0.1875, minimum=0.0001))
                    field("bulge_mass", gr.Number(label="Bulge mass", value=0.0625, minimum=0.0001))
                    field("halo_mass", gr.Number(label="Halo mass", value=1.0, minimum=0.0001))
                    field("disk_scale", gr.Number(label="Disk radial scale", value=0.0833333333, minimum=0.0001))
                    field("disk_height", gr.Number(label="Disk thickness scale", value=0.005, minimum=0.0001))
                    field("bulge_scale", gr.Number(label="Bulge scale radius", value=0.04168, minimum=0.0001))
                    field("halo_scale", gr.Number(label="Halo scale radius", value=0.1, minimum=0.0001))
                with gr.Accordion("Encounter geometry", open=False):
                    field("separation", gr.Slider(0.3, 4.0, value=1.0, step=0.05, label="Initial separation"))
                    field("pericenter", gr.Slider(0.02, 1.0, value=0.25, step=0.01, label="Nominal pericenter"))
                    field("inc_a", gr.Slider(-180, 180, value=25, step=5, label="Galaxy A inclination (°)"))
                    field("node_a", gr.Slider(-180, 180, value=-30, step=5, label="Galaxy A node angle (°)"))
                    field("inc_b", gr.Slider(-180, 180, value=40, step=5, label="Galaxy B inclination (°)"))
                    field("node_b", gr.Slider(-180, 180, value=60, step=5, label="Galaxy B node angle (°)"))
                with gr.Accordion("Integration and force solver", open=True):
                    field("solver", gr.Dropdown(label="Gravity solver", choices=["Fast direct (Numba)", "Direct (NumPy)", "Barnes–Hut"], value="Fast direct (Numba)"))
                    field("eps", gr.Number(label="Softening ε", value=0.01, minimum=0.00001))
                    field("theta", gr.Slider(0.2, 1.0, value=0.5, step=0.05, label="Barnes–Hut opening angle θ"))
                    field("dt", gr.Number(label="Timestep Δt", value=0.006, minimum=0.00001))
                    field("steps", gr.Slider(20, 1200, value=240, step=20, label="Integration steps"))
                    field("save_every", gr.Slider(1, 50, value=4, step=1, label="Save a frame every N steps"))
                with gr.Accordion("Animation settings", open=True):
                    field("view", gr.Radio(["both", "xy", "xz"], value="both", label="Camera view"))
                    field("lim", gr.Slider(0.5, 4.0, value=1.5, step=0.1, label="Plot half-width"))
                    field("fps", gr.Slider(5, 40, value=15, step=1, label="Animation FPS"))
                    field("show_halo", gr.Checkbox(value=False, label="Show dark-matter particles"))
                    field("export", gr.Radio(["GIF", "MP4", "Both"], value="GIF", label="Export format"))
                run_button = gr.Button("Run simulation and animate", variant="primary", size="lg")
            with gr.Column(scale=2):
                status = gr.Markdown("Set parameters on the left, then run the simulation. The Numba solver compiles on first use.")
                gif_display = gr.Image(label="Animated GIF", type="filepath", interactive=False)
                video_display = gr.Video(label="MP4 animation")
                separation_display = gr.Plot(label="Core separation")
                with gr.Row():
                    trajectory_file = gr.File(label="Download trajectory (.npz)", interactive=False)
                    config_file = gr.File(label="Download run settings (.json)", interactive=False)
    with gr.Tab("Replay a saved run"):
        gr.Markdown("Re-render a saved trajectory with a different camera, zoom, or frame rate. **No gravitational calculations are repeated.** The bundled V1.1 trajectory is available immediately.")
        with gr.Row():
            saved = gr.Dropdown(choices=available_runs(), label="Saved trajectory", value=available_runs()[0] if available_runs() else None, scale=4)
            refresh = gr.Button("Refresh list", scale=1)
        with gr.Row():
            replay_view = gr.Radio(["both", "xy", "xz"], value="both", label="Camera")
            replay_lim = gr.Slider(0.5, 4, value=1.5, step=0.1, label="Zoom / half-width")
            replay_fps = gr.Slider(5, 40, value=15, step=1, label="FPS")
            replay_halo = gr.Checkbox(value=False, label="Show halo")
            replay_format = gr.Radio(["GIF", "MP4", "Both"], value="GIF", label="Export")
        replay_button = gr.Button("Render saved trajectory", variant="primary")
        replay_summary = gr.Markdown()
        with gr.Row():
            replay_gif = gr.Image(label="GIF", type="filepath", interactive=False)
            replay_video = gr.Video(label="MP4")
        replay_plot = gr.Plot(label="Core separation")
        replay_download = gr.File(label="Download rendered animation", interactive=False)
    with gr.Tab("How this works"):
        gr.Markdown("""
### Start here
1. Select the **Run an encounter** tab and leave the default settings initially.
2. Choose **Fast direct (Numba)** for a few thousand particles. The first run compiles the force calculation.
3. Click **Run simulation and animate**. Watch the GIF, inspect the separation plot, and download the particle history.
4. Go to **Replay a saved run** to adjust the view and frame rate without rerunning the physics.
5. Every run saves a trajectory and configuration inside `results/runs/`.

### What the parameters mean
- **Disk/bulge/halo particles** are counts *per galaxy*. Halos are always live even if hidden.
- **Disk/bulge/halo masses and scales** change the intended galaxy model. Native preview velocities remain heuristic after these changes.
- **Initial separation/pericenter** define an approximate parabolic point-mass encounter. The actual pericenter of the live, extended galaxies can differ.
- **Opening angle θ** only applies to Barnes–Hut; smaller values usually give better force accuracy but are slower.
- **Softening ε** regularizes close encounters; **Δt** controls integration accuracy.
- **Save every N steps** determines how many movie frames are recorded, not how many steps are integrated.
- **Plot half-width** is in Barnes code-length units; change it in Replay if objects leave the frame.

### Research status
The native preview samples the requested spatial density profiles, but it does **not** construct a self-consistent disk–bulge–halo equilibrium. More particles and a compelling animation do not make its predictions scientifically validated. AGAMA is optional and its adapter still needs a local runtime and isolated-galaxy stability checks.
""")

    run_button.click(simulate, inputs=FIELD_CONTROLS, outputs=[gif_display, video_display, trajectory_file, config_file, separation_display, status])
    refresh.click(lambda: gr.update(choices=available_runs()), outputs=saved)
    replay_button.click(replay_saved, inputs=[saved, replay_view, replay_lim, replay_fps, replay_halo, replay_format], outputs=[replay_gif, replay_video, replay_download, replay_plot, replay_summary])


if __name__ == "__main__":
    app.launch(inbrowser=True, share=False, server_name="127.0.0.1", allowed_paths=[str(RESULTS)])
