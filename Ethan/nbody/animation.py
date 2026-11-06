"""Animated stellar distributions from stored simulation frames."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import numpy as np


def render(times, positions, components, galaxies, path, fps=15, lim=1.2,
           view='both', show_halo=False):
    projections = {'xy': (0, 1), 'xz': (0, 2)}
    panels = list(projections) if view == 'both' else [view]
    fig, axs = plt.subplots(1, len(panels), figsize=(5.8*len(panels), 5.2),
                            squeeze=False)
    axs = axs[0]
    colors = ['#5bd5fb', '#ffb66e']
    displayed = ((components <= 1) | (components == 3)) if not show_halo else np.ones(len(components), dtype=bool)
    artists = []
    for ax, panel in zip(axs, panels):
        dims = projections[panel]
        series = []
        for gid in (0, 1):
            mask = displayed & (galaxies == gid)
            disk = mask & ((components == 0) | (components == 3))
            bulge = mask & (components == 1)
            halo = mask & (components == 2)
            for subset, size, alpha in ((disk, 1.1, 0.72), (bulge, 2.5, 0.86),
                                        (halo, 0.15, 0.1)):
                if not subset.any():
                    continue
                plot = ax.scatter([], [], s=size, c=colors[gid], alpha=alpha,
                                  linewidths=0, rasterized=True)
                series.append((plot, subset))
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel=panel[0], ylabel=panel[1])
        ax.set_aspect('equal')
        ax.set_facecolor('#080e1d')
        ax.grid(alpha=.08)
        artists.append((dims, series))
    title = fig.suptitle('')
    note = ('Analytic Plummer DF: finite-N stability requires validation' if np.any(components==3)
            else 'Preview only: isolated disk stability and Mice reconstruction are not validated')
    fig.text(.5, .012, note,
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .95))

    def update(k):
        for dims, series in artists:
            for plot, mask in series:
                plot.set_offsets(positions[k, mask][:, dims])
        title.set_text(f'N-Baddies: stellar encounter  |  t = {times[k]:.3f} (code units)')
        return [plot for _, series in artists for plot, _ in series] + [title]

    movie = FuncAnimation(fig, update, frames=len(times), interval=1000/fps,
                          blit=False, repeat=True)
    writer = FFMpegWriter(fps=fps, bitrate=2400) if str(path).endswith('.mp4') else PillowWriter(fps=fps)
    movie.save(str(path), writer=writer, dpi=110)
    plt.close(fig)
