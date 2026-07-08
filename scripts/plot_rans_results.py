#!/usr/bin/env python3
"""
Plot RANS solver output files.

Expected solver files:
  flowfield: "# i j xc yc rho u v p mach"
  wall:      "# i x_wall y_wall p cp tau_x tau_y"

Examples:
  python scripts/plot_rans_results.py
  python scripts/plot_rans_results.py --flow rans_test_flowfield.dat --wall rans_test_wall.dat
  python scripts/plot_rans_results.py --flow flowfield.dat --wall wall.dat --out plots --no-show
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


def load_flowfield(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 9:
        raise ValueError(f"{path} must contain columns: i j xc yc rho u v p mach")

    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    ni = i.max() + 1
    nj = j.max() + 1

    fields = {
        "x": data[:, 2],
        "y": data[:, 3],
        "rho": data[:, 4],
        "u": data[:, 5],
        "v": data[:, 6],
        "p": data[:, 7],
        "mach": data[:, 8],
    }
    fields["speed"] = np.sqrt(fields["u"] ** 2 + fields["v"] ** 2)
    fields["temperature"] = fields["p"] / fields["rho"]

    gridded = {}
    for name, values in fields.items():
        grid = np.full((nj, ni), np.nan)
        grid[j, i] = values
        gridded[name] = grid

    return gridded, ni, nj


def load_wall(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 7:
        raise ValueError(f"{path} must contain columns: i x_wall y_wall p cp tau_x tau_y")

    wall = {
        "i": data[:, 0].astype(int),
        "x": data[:, 1],
        "y": data[:, 2],
        "p": data[:, 3],
        "cp": data[:, 4],
        "tau_x": data[:, 5],
        "tau_y": data[:, 6],
    }
    wall["tau_mag"] = np.sqrt(wall["tau_x"] ** 2 + wall["tau_y"] ** 2)
    wall["s"] = surface_distance(wall["x"], wall["y"])
    return wall


def surface_distance(x, y):
    s = np.zeros_like(x)
    if len(x) > 1:
        ds = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
        s[1:] = np.cumsum(ds)
    return s


def finite_limits(values):
    good = np.isfinite(values)
    if not np.any(good):
        return None
    return float(np.nanmin(values)), float(np.nanmax(values))


def plot_scalar(ax, flow, field_name, title, cmap):
    x = flow["x"]
    y = flow["y"]
    q = flow[field_name]

    levels = 50
    cf = ax.tricontourf(x.ravel(), y.ravel(), q.ravel(), levels=levels, cmap=cmap)
    ax.tricontour(x.ravel(), y.ravel(), q.ravel(), levels=12, colors="k", linewidths=0.25, alpha=0.35)
    ax.plot(x[0, :], y[0, :], color="black", lw=1.2)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return cf


def save_field_plot(flow, field_name, title, cmap, out_path, no_show):
    fig, ax = plt.subplots(figsize=(9, 6))
    cf = plot_scalar(ax, flow, field_name, title, cmap)
    fig.colorbar(cf, ax=ax, label=field_name)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if not no_show:
        plt.show()
    plt.close(fig)


def save_overview(flow, out_path, no_show):
    panels = [
        ("mach", "Mach number", "viridis"),
        ("p", "Pressure", "coolwarm"),
        ("rho", "Density", "magma"),
        ("speed", "Speed", "plasma"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for ax, (field, title, cmap) in zip(axes.ravel(), panels):
        cf = plot_scalar(ax, flow, field, title, cmap)
        fig.colorbar(cf, ax=ax, shrink=0.85)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if not no_show:
        plt.show()
    plt.close(fig)


def save_velocity_vectors(flow, out_path, no_show, max_vectors):
    x = flow["x"]
    y = flow["y"]
    u = flow["u"]
    v = flow["v"]
    speed = flow["speed"]

    nj, ni = x.shape
    stride_i = max(1, ni // max_vectors)
    stride_j = max(1, nj // max_vectors)

    fig, ax = plt.subplots(figsize=(9, 6))
    cf = ax.tricontourf(x.ravel(), y.ravel(), speed.ravel(), levels=50, cmap="plasma")
    ax.quiver(
        x[::stride_j, ::stride_i],
        y[::stride_j, ::stride_i],
        u[::stride_j, ::stride_i],
        v[::stride_j, ::stride_i],
        color="white",
        alpha=0.85,
        width=0.0025,
    )
    ax.plot(x[0, :], y[0, :], color="black", lw=1.2)
    ax.set_title("Velocity vectors over speed")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(cf, ax=ax, label="speed")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if not no_show:
        plt.show()
    plt.close(fig)


def save_wall_plots(wall, out_path, no_show):
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(wall["s"], wall["cp"], color="#1f77b4", lw=1.6)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Cp")
    axes[0].set_title("Wall pressure coefficient")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(wall["s"], wall["p"], color="#d62728", lw=1.6)
    axes[1].set_ylabel("p")
    axes[1].set_title("Wall pressure")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(wall["s"], wall["tau_x"], label="tau_x", color="#2ca02c", lw=1.4)
    axes[2].plot(wall["s"], wall["tau_y"], label="tau_y", color="#9467bd", lw=1.4)
    axes[2].plot(wall["s"], wall["tau_mag"], label="|tau|", color="black", lw=1.4)
    axes[2].set_xlabel("surface distance")
    axes[2].set_ylabel("wall shear")
    axes[2].set_title("Wall shear traction")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if not no_show:
        plt.show()
    plt.close(fig)


def save_surface_map(wall, out_path, no_show):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sc = ax.scatter(wall["x"], wall["y"], c=wall["cp"], s=18, cmap="coolwarm")
    ax.plot(wall["x"], wall["y"], color="black", lw=0.8, alpha=0.55)
    ax.set_title("Airfoil surface colored by Cp")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(sc, ax=ax, label="Cp")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if not no_show:
        plt.show()
    plt.close(fig)


def print_summary(flow, wall):
    print("\nFlowfield ranges:")
    for name in ["mach", "p", "rho", "speed", "temperature", "u", "v"]:
        lo_hi = finite_limits(flow[name])
        if lo_hi is not None:
            print(f"  {name:12s}: {lo_hi[0]: .6e} to {lo_hi[1]: .6e}")

    if wall is not None:
        print("\nWall ranges:")
        for name in ["cp", "p", "tau_x", "tau_y", "tau_mag"]:
            lo_hi = finite_limits(wall[name])
            if lo_hi is not None:
                print(f"  {name:12s}: {lo_hi[0]: .6e} to {lo_hi[1]: .6e}")


def main():
    ap = argparse.ArgumentParser(description="Plot RANS flowfield and wall output.")
    ap.add_argument("--flow", default="rans_flowfield.dat",
                    help="Flowfield file written by write_flowfield")
    ap.add_argument("--wall", default="rans_wall.dat",
                    help="Wall file written by write_wall_data")
    ap.add_argument("--out", default="rans_plots",
                    help="Directory where PNG plots are written")
    ap.add_argument("--no-show", action="store_true",
                    help="Save plots without opening interactive windows")
    ap.add_argument("--max-vectors", type=int, default=32,
                    help="Approximate max vectors per direction in velocity plot")
    args = ap.parse_args()

    flow_path = Path(args.flow)
    wall_path = Path(args.wall)
    out_dir = Path(args.out)

    if not flow_path.exists():
        raise FileNotFoundError(f"Flowfield file not found: {flow_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    flow, ni, nj = load_flowfield(flow_path)
    wall = load_wall(wall_path) if wall_path.exists() else None

    print(f"Loaded flowfield: {flow_path} ({ni} x {nj} cells)")
    if wall is not None:
        print(f"Loaded wall data: {wall_path} ({len(wall['i'])} faces)")
    else:
        print(f"Wall file not found: {wall_path} -- skipping wall plots")

    print_summary(flow, wall)

    save_overview(flow, out_dir / "overview.png", args.no_show)
    save_field_plot(flow, "mach", "Mach number", "viridis", out_dir / "mach.png", args.no_show)
    save_field_plot(flow, "p", "Pressure", "coolwarm", out_dir / "pressure.png", args.no_show)
    save_field_plot(flow, "rho", "Density", "magma", out_dir / "density.png", args.no_show)
    save_field_plot(flow, "temperature", "Temperature p/rho", "inferno", out_dir / "temperature.png", args.no_show)
    save_field_plot(flow, "speed", "Speed", "plasma", out_dir / "speed.png", args.no_show)
    save_velocity_vectors(flow, out_dir / "velocity_vectors.png", args.no_show, args.max_vectors)

    if wall is not None:
        save_wall_plots(wall, out_dir / "wall_cp_shear.png", args.no_show)
        save_surface_map(wall, out_dir / "surface_cp.png", args.no_show)

    return 0


if __name__ == "__main__":
    sys.exit(main())
