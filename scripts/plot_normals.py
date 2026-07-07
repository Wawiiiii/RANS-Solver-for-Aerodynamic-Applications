#!/usr/bin/env python3
"""
Visual check of the RansSolver face normals.

Reads the two files rans_solver_test writes:
  rans_geom_nodes.dat    "i j x y"           (mesh nodes)
  rans_geom_normals.dat  "fx fy nx ny length" (outward face normals)

Draws the grid and an arrow at each face center pointing along its outward unit
normal, scaled by the face length so wall cells (tiny) and outer cells (large)
both read cleanly. A full view plus a leading-edge zoom.

Usage:
    python plot_normals.py [--nodes rans_geom_nodes.dat]
                           [--normals rans_geom_normals.dat]
                           [--save out.png] [--no-show] [--scale 0.6]
"""

import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt


def load_nodes(path):
    ii, jj, xx, yy = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            ii.append(int(p[0])); jj.append(int(p[1]))
            xx.append(float(p[2])); yy.append(float(p[3]))
    ii = np.array(ii); jj = np.array(jj)
    nt = ii.max() + 1
    nr = jj.max() + 1
    x = np.full((nr, nt), np.nan); y = np.full((nr, nt), np.nan)
    x[jj, ii] = xx; y[jj, ii] = yy
    return x, y


def load_normals(path):
    data = np.loadtxt(path)  # columns: fx fy nx ny length
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]


def draw(ax, x, y, fx, fy, nx, ny, length, scale, title, xlim=None, ylim=None):
    nr, nt = x.shape
    for j in range(nr):
        ax.plot(x[j], y[j], "-", color="#9fb3c8", lw=0.4)
    for i in range(nt):
        ax.plot(x[:, i], y[:, i], "-", color="#9fb3c8", lw=0.4)
    ax.plot(x[0], y[0], "-", color="k", lw=1.2)  # airfoil surface

    # Arrow vector = unit normal * face length * scale.
    ax.quiver(fx, fy, nx * length * scale, ny * length * scale,
              angles="xy", scale_units="xy", scale=1.0,
              color="#c1121f", width=0.003)

    ax.set_aspect("equal")
    ax.set_title(title)
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)


def main():
    ap = argparse.ArgumentParser(description="Visualize RansSolver face normals.")
    ap.add_argument("--nodes", default="rans_geom_nodes.dat")
    ap.add_argument("--normals", default="rans_geom_normals.dat")
    ap.add_argument("--save", default=None)
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--scale", type=float, default=0.6,
                    help="arrow length as a fraction of face length")
    args = ap.parse_args()

    x, y = load_nodes(args.nodes)
    fx, fy, nx, ny, length = load_normals(args.normals)

    fig, ax = plt.subplots(1, 2, figsize=(15, 7))

    draw(ax[0], x, y, fx, fy, nx, ny, length, args.scale, "Face normals (full)")

    # Leading-edge zoom.
    xle = x[0].min()
    draw(ax[1], x, y, fx, fy, nx, ny, length, args.scale,
         "Leading-edge zoom",
         xlim=(xle - 0.05, xle + 0.35),
         ylim=(-0.2, 0.2))

    plt.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=110)
        print(f"Saved {args.save}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
