import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rot(yaw, pitch, roll=0.0):
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll), np.sin(roll)
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1.]])
    ry = np.array([[cp, 0, sp], [0, 1., 0], [-sp, 0, cp]])
    rx = np.array([[1., 0, 0], [0, cr, -sr], [0, sr, cr]])
    # Camera azimuth changes the line of sight; elevation is above xy plane.
    right = np.array([-sy, cy, 0.])
    up = np.array([-sp * cy, -sp * sy, cp])
    sight = np.array([cp * cy, cp * sy, sp])
    return np.stack([right, up, sight])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", type=Path)
    ap.add_argument("--out", type=Path, default=Path("mice_view_search.png"))
    ap.add_argument('--azimuths', type=float, nargs='+', default=[0,60,120,180])
    ap.add_argument('--elevations', type=float, nargs='+', default=[0,20,40])
    ap.add_argument('--after', type=float, nargs='+', default=[.5,.75,.875,1.0])
    args = ap.parse_args()
    d = np.load(args.trajectory, allow_pickle=True)
    t, x, c, g = d["t"], d["x"], d["component"], d["galaxy"]

    centers = [np.median(x[:, (c == 1) & (g == gg)], axis=1) for gg in (0, 1)]
    peri = float(t[np.argmin(np.linalg.norm(centers[1] - centers[0], axis=1))])
    frame_ids = list(dict.fromkeys(int(np.argmin(abs(t - (peri + q)))) for q in args.after))
    views = [(np.deg2rad(y), np.deg2rad(p)) for p in args.elevations for y in args.azimuths]
    fig, axes = plt.subplots(len(frame_ids), len(views), figsize=(max(12,2*len(views)), 2.5*len(frame_ids)), squeeze=False,
                             facecolor="#08090e")
    colors = {0: "#62d5e8", 1: "#f465b5"}
    disk = c == 0
    for r, fi in enumerate(frame_ids):
        for j, (yaw, pitch) in enumerate(views):
            ax = axes[r, j]
            q = x[fi, disk] @ rot(yaw, pitch).T
            # Global frame: do not recenter each galaxy.
            for gg in (0, 1):
                mask = g[disk] == gg
                ax.scatter(q[mask, 0], q[mask, 1], s=0.18, alpha=.42,
                           color=colors[gg], linewidths=0, rasterized=True)
            ax.set_aspect("equal", adjustable="box")
            ax.set_facecolor("#08090e")
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
            for s in ax.spines.values(): s.set_visible(False)
            if j == 0: ax.set_ylabel(f"t={t[fi]:.3f}\nafter passage={t[fi]-peri:.3f}", color="white", fontsize=10)
            if r == 0: ax.set_title(f"az {np.rad2deg(yaw):.0f}\nelev {np.rad2deg(pitch):.0f}", color="#aab", fontsize=8)
    fig.suptitle("Global-frame camera search — collisionless Mice trial", color="white", fontsize=16)
    fig.text(.5, .015, "cyan/pink = disk particles from the two galaxies; no per-galaxy recentering", ha="center", color="#aab", fontsize=10)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved {args.out}")


if __name__ == "__main__": main()
