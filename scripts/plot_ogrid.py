#!/usr/bin/env python3
"""
Plot a structured O-grid around the airfoil.

Reads a mesh written by the test drivers in "i j x y" format (one node per line,
blank lines between radial layers), e.g. mesh_rans_256x128.dat or mesh_vj_256x256.dat.

Usage:
    python plot_ogrid.py <mesh.dat> [--save out.png] [--stride-i N] [--stride-j N]
                                    [--bl-layers K] [--no-show]

Examples:
    python plot_ogrid.py build/mesh_rans_256x128.dat
    python plot_ogrid.py build/mesh_rans_256x128.dat --save ogrid.png --no-show
"""

import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt


def load_mesh(path):
    """Load an 'i j x y' indexed mesh into (nr, nt) x/y arrays."""
    ii, jj, xx, yy = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            ii.append(int(parts[0]))
            jj.append(int(parts[1]))
            xx.append(float(parts[2]))
            yy.append(float(parts[3]))

    if not ii:
        raise ValueError(f"No mesh data found in {path}")

    ii = np.asarray(ii)
    jj = np.asarray(jj)
    nt = ii.max() + 1
    nr = jj.max() + 1

    x = np.full((nr, nt), np.nan)
    y = np.full((nr, nt), np.nan)
    x[jj, ii] = xx
    y[jj, ii] = yy

    if np.isnan(x).any():
        raise ValueError(f"{path} is not a complete structured grid (missing nodes)")

    return x, y, nt, nr


def plot_grid(ax, x, y, si, sj, jmax=None, color="#1f4e8c", lw=0.4):
    """Draw grid lines: radial layers (const j) and radial rays (const i)."""
    nr, nt = x.shape
    jhi = nr if jmax is None else min(jmax, nr)

    # Circumferential lines (constant j)
    for j in range(0, jhi, sj):
        ax.plot(x[j, :], y[j, :], "-", color=color, lw=lw)

    # Radial lines (constant i)
    for i in range(0, nt, si):
        ax.plot(x[:jhi, i], y[:jhi, i], "-", color=color, lw=lw)

    # Emphasize the airfoil surface (j = 0)
    ax.plot(x[0, :], y[0, :], "-", color="k", lw=1.3)


def main():
    ap = argparse.ArgumentParser(description="Plot an O-grid around the airfoil.")
    ap.add_argument("mesh", help="Mesh .dat file (i j x y format)")
    ap.add_argument("--save", default=None, help="Save figure to this path (PNG)")
    ap.add_argument("--stride-i", type=int, default=0,
                    help="Draw every Nth radial ray (0 = auto)")
    ap.add_argument("--stride-j", type=int, default=0,
                    help="Draw every Nth circumferential layer (0 = auto)")
    ap.add_argument("--bl-layers", type=int, default=40,
                    help="Number of near-wall layers shown in the zoom panels")
    ap.add_argument("--no-show", action="store_true",
                    help="Do not open an interactive window")
    args = ap.parse_args()

    x, y, nt, nr = load_mesh(args.mesh)

    # Auto strides so busy meshes stay legible.
    si = args.stride_i if args.stride_i > 0 else max(1, nt // 80)
    sj = args.stride_j if args.stride_j > 0 else max(1, nr // 50)
    bl = min(args.bl_layers, nr)

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f"{args.mesh}   (nt={nt}, nr={nr})", fontsize=12)

    # --- Panel 1: full O-mesh ---
    plot_grid(ax[0], x, y, si, sj)
    ax[0].set_title("Full O-mesh")
    ax[0].set_aspect("equal")
    span = max(np.ptp(x[0, :]), np.ptp(y[0, :]))
    r = 2.5 * max(span, 1.0)
    cx, cy = np.mean(x[0, :]), np.mean(y[0, :])
    ax[0].set_xlim(cx - r, cx + r)
    ax[0].set_ylim(cy - r, cy + r)

    # --- Panel 2: airfoil + boundary-layer region ---
    plot_grid(ax[1], x, y, max(1, nt // 160), 1, jmax=bl, lw=0.5)
    ax[1].set_title(f"Boundary-layer region (first {bl} layers)")
    ax[1].set_aspect("equal")
    ax[1].set_xlim(x[0, :].min() - 0.1, x[0, :].max() + 0.1)
    yr = 1.6 * np.abs(y[0, :]).max()
    ax[1].set_ylim(cy - yr, cy + yr)

    # --- Panel 3: leading-edge zoom ---
    plot_grid(ax[2], x, y, max(1, nt // 300), 1, jmax=min(bl, nr), lw=0.6)
    ax[2].set_title("Leading-edge clustering")
    ax[2].set_aspect("equal")
    xle = x[0, :].min()
    ax[2].set_xlim(xle - 0.01, xle + 0.06)
    ax[2].set_ylim(cy - 0.04, cy + 0.04)

    plt.tight_layout(rect=(0, 0, 1, 0.96))

    if args.save:
        fig.savefig(args.save, dpi=120)
        print(f"Saved figure to {args.save}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
