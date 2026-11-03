"""3D WebGL viewer for saved N-body trajectories. No additional physics runs."""
import numpy as np
import plotly.graph_objects as go

COLORS = ('#5bd5fb', '#ffb66e')


def figure3d(times, x, component, origin, limit=1.5, show_halo=False,
             max_points=3000):
    """Interactive rotation/zoom + time slider, with consistent particle sampling."""
    visible = (component <= 1) if not show_halo else np.ones(len(component), bool)
    indices = np.flatnonzero(visible)
    if len(indices) > max_points:
        indices = np.random.default_rng(2026).choice(indices, max_points, replace=False)
        indices.sort()
    groups = [indices[origin[indices] == g] for g in sorted(set(origin[indices]))]

    def traces(k):
        return [go.Scatter3d(
            x=x[k, ix, 0], y=x[k, ix, 1], z=x[k, ix, 2],
            mode='markers', name=f'Galaxy {gid+1}',
            marker=dict(size=2.3, opacity=0.7, color=COLORS[gid % len(COLORS)]),
            hoverinfo='skip') for gid, ix in enumerate(groups)]

    frames = [go.Frame(name=str(k), data=traces(k), traces=list(range(len(groups))))
              for k in range(len(times))]
    steps = [dict(method='animate', label=f'{float(t):.2f}',
                  args=[[str(k)], dict(mode='immediate', frame=dict(duration=0, redraw=True),
                                      transition=dict(duration=0))])
             for k, t in enumerate(times)]
    fig = go.Figure(data=traces(0), frames=frames)
    fig.update_layout(
        template='plotly_dark', title=f'3D encounter · {len(indices):,} displayed / {len(component):,} live particles',
        paper_bgcolor='#080e1d', plot_bgcolor='#080e1d',
        scene=dict(xaxis=dict(range=[-limit, limit], title='x'),
                   yaxis=dict(range=[-limit, limit], title='y'),
                   zaxis=dict(range=[-limit, limit], title='z'),
                   aspectmode='cube', bgcolor='#080e1d'),
        margin=dict(t=60, b=20, l=5, r=5), height=610,
        updatemenus=[dict(type='buttons', direction='left', x=0, y=0,
                          buttons=[dict(label='▶ Play', method='animate',
                                        args=[None, dict(frame=dict(duration=90, redraw=True),
                                                         transition=dict(duration=0), fromcurrent=True)]),
                                   dict(label='❚❚ Pause', method='animate',
                                        args=[[None], dict(mode='immediate',
                                                           frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, x=0.12, y=0, len=0.85, currentvalue=dict(prefix='t = '), steps=steps)],
        uirevision='camera')
    return fig
