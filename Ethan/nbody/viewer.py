"""Interactive 3D particle playback, independent of Gradio's Plotly renderer."""
import json
from html import escape
from pathlib import Path
from urllib.parse import quote

import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

COLORS = ('#5bd5fb', '#ffb66e')


def _sample(component, origin, show_halo, max_points):
    indices = np.flatnonzero((component <= 1) if not show_halo else np.ones(len(component), bool))
    if len(indices) > max_points:
        indices = np.sort(np.random.default_rng(2026).choice(indices, max_points, replace=False))
    return indices, [indices[origin[indices] == g] for g in sorted(set(origin[indices]))]


def figure3d(times, x, component, origin, limit=1.5, show_halo=False, max_points=3000):
    """Plotly Figure, retained for analysis/legacy export. Use write_viewer() for playback."""
    indices, groups = _sample(component, origin, show_halo, max_points)

    def traces(k):
        return [go.Scatter3d(x=x[k, ix, 0], y=x[k, ix, 1], z=x[k, ix, 2],
                             mode='markers', name=f'Galaxy {gid + 1}',
                             marker=dict(size=2.3, opacity=.7, color=COLORS[gid % 2]),
                             hoverinfo='skip') for gid, ix in enumerate(groups)]

    frames = [go.Frame(name=str(k), data=traces(k), traces=list(range(len(groups))))
              for k in range(len(times))]
    steps = [dict(method='animate', label=f'{float(t):.2f}',
                  args=[[str(k)], dict(mode='immediate', frame=dict(duration=0, redraw=True),
                                      transition=dict(duration=0))])
             for k, t in enumerate(times)]
    fig = go.Figure(data=traces(0), frames=frames)
    fig.update_layout(template='plotly_dark',
                      title=f'3D encounter · {len(indices):,} displayed / {len(component):,} live particles',
                      scene=dict(xaxis=dict(range=[-limit, limit]),
                                 yaxis=dict(range=[-limit, limit]),
                                 zaxis=dict(range=[-limit, limit]), aspectmode='cube'),
                      height=610, margin=dict(l=0, r=0, b=30, t=45),
                      updatemenus=[dict(type='buttons', x=0, y=0,
                                        buttons=[dict(label='Play', method='animate', args=[None])])],
                      sliders=[dict(steps=steps, active=0)])
    return fig


def write_viewer(path, times, x, component, origin, limit=1.5, show_halo=False,
                 max_points=3000, max_frames=180):
    """Write a self-contained HTML viewer with reliable JavaScript-driven playback.

    We don't rely on Plotly's animation frames inside Gradio's gr.Plot component:
    the browser controls directly restyle the three-dimensional traces.
    """
    indices, groups = _sample(component, origin, show_halo, int(max_points))
    sample_times = np.unique(np.linspace(0, len(times) - 1,
                                         min(len(times), max_frames)).round().astype(int))
    data = [[np.round(x[k, ix], 6).T.tolist() for ix in groups] for k in sample_times]
    time_data = np.round(times[sample_times], 6).tolist()
    payload = json.dumps(dict(frames=data, times=time_data), separators=(',', ':'))
    names = json.dumps([f'Galaxy {i + 1}' for i in range(len(groups))])
    colors = json.dumps([COLORS[i % len(COLORS)] for i in range(len(groups))])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    html = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>
*{box-sizing:border-box}body{margin:0;background:#080e1d;color:#e8efff;font:14px system-ui,Arial,sans-serif}
#plot{width:100%;height:615px}.controls{display:flex;align-items:center;gap:12px;padding:8px 18px 18px;flex-wrap:wrap}
button{background:#243351;color:white;border:1px solid #445677;border-radius:6px;padding:8px 16px;cursor:pointer}
button:hover{background:#354b76}input[type=range]{flex:1;min-width:180px;accent-color:#58c7ff}
small{color:#aabbd3}#clock{font-variant-numeric:tabular-nums;min-width:105px}
</style></head><body>
<div id="plot"></div>
<div class="controls"><button id="play">▶ Play</button><button id="prev">◀</button>
<button id="next">▶</button><input id="timeline" type="range" min="0" max="0" step="1" value="0">
<span id="clock"></span><small id="info"></small></div>
<script>__PLOTLY_SCRIPT__</script>
<script>
const payload=__PAYLOAD__;
const frames=payload.frames, times=payload.times;
const colors=__COLORS__, names=__NAMES__;
const plot=document.getElementById('plot'), slider=document.getElementById('timeline');
const play=document.getElementById('play'), clock=document.getElementById('clock');
slider.max=frames.length-1;
const limit=__LIMIT__;
const traceIds=frames[0].map((_,i)=>i);
const traceData=frames[0].map((coords,i)=>({type:'scatter3d', mode:'markers',
 x:coords[0], y:coords[1], z:coords[2], name:names[i], hoverinfo:'skip',
 marker:{size:2.3,opacity:.7,color:colors[i]}}));
const layout={paper_bgcolor:'#080e1d',plot_bgcolor:'#080e1d',font:{color:'#e8efff'},
 title:{text:'3D encounter · __PARTICLE_COUNT__ displayed / __LIVE_COUNT__ live particles'},
 scene:{xaxis:{range:[-limit,limit],title:'x'},yaxis:{range:[-limit,limit],title:'y'},
 zaxis:{range:[-limit,limit],title:'z'},aspectmode:'cube',bgcolor:'#080e1d'},
 margin:{l:0,r:0,b:0,t:50},uirevision:'keep-camera',showlegend:true};
let current=0,requested=0,drawing=false,timer=null;
function label(i){slider.value=i;clock.textContent='t = '+times[i].toFixed(3);
 document.getElementById('info').textContent=(i+1)+' / '+frames.length+' frames';}
async function seek(i){requested=Math.max(0,Math.min(frames.length-1,i));
 if(drawing)return;drawing=true;
 do{const k=requested, f=frames[k];
  await Plotly.restyle(plot,{x:f.map(g=>g[0]),y:f.map(g=>g[1]),z:f.map(g=>g[2])},traceIds);
  current=k;label(k);
 }while(current!==requested);
 drawing=false;
}
function pause(){if(timer){clearInterval(timer);timer=null;}play.textContent='▶ Play';}
async function toggle(){if(timer){pause();return;}if(current>=frames.length-1)await seek(0);
 play.textContent='❚❚ Pause';timer=setInterval(()=>{
  if(current>=frames.length-1){pause();return;}
  seek(Math.min(frames.length-1,requested+1));},100);}
slider.addEventListener('input',()=>{pause();seek(Number(slider.value));});
play.addEventListener('click',toggle);
document.getElementById('prev').addEventListener('click',()=>{pause();seek(current-1)});
document.getElementById('next').addEventListener('click',()=>{pause();seek(current+1)});
Plotly.newPlot(plot,traceData,layout,{responsive:true,displaylogo:false}).then(()=>label(0));
</script></body></html>'''
    for key, value in dict(__PLOTLY_SCRIPT__=get_plotlyjs(), __PAYLOAD__=payload,
                           __COLORS__=colors, __NAMES__=names,
                           __LIMIT__=repr(float(limit)),
                           __PARTICLE_COUNT__=str(len(indices)),
                           __LIVE_COUNT__=str(len(component))).items():
        html = html.replace(key, value)
    path.write_text(html, encoding='utf8')
    return path


def viewer_iframe(path):
    """Show the standalone HTML inside Gradio rather than a frame-stripping plot widget."""
    url = '/gradio_api/file=' + quote(Path(path).resolve().as_posix(), safe='/:')
    return (f'<iframe src="{escape(url, quote=True)}" style="width:100%;height:720px;'
            'border:0;border-radius:8px;background:#080e1d" loading="eager" '
            'title="Interactive 3D galaxy trajectory"></iframe>')
