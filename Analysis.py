#!/usr/bin/env python3
"""
Analysis.py — Droplet breakup parameter sweep analysis

Parses all droplet_output_d<D>_sr<SR>_t<T>.csv files under results/,
groups by (depth, SR), filters timed-out trials, and produces:

  1. daughter_area_grid.png   — grid of histograms (rows=depth, cols=SR)
  2. daughter_area_means.png  — mean +/- std vs SR, one line per depth
  3. daughter_areas_by_depth_sr.csv — long-format CSV for downstream fitting

Usage:
    python3 Analysis.py
    python3 Analysis.py --pattern "results/d60_*/droplet_output_*.csv"
    python3 Analysis.py --bins 20 --grid-fig my_grid.png
"""

import argparse
import glob
import math
import re
from collections import defaultdict
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (must match droplet_sim.cpp)
# ---------------------------------------------------------------------------
D_DROPS_INIT = 1.0
D_OBSTS      = D_DROPS_INIT / 6.0
A0_INIT      = math.pi / 4.0 * D_DROPS_INIT**2   # initial droplet area


def exit_y_threshold(depth_cells: int) -> float:
    """Reproduces C++ formula: obsts_box_bottom - 2 * D_drops_init."""
    obsts_box_bottom = -(depth_cells * D_OBSTS) / 2.0
    return obsts_box_bottom - 2.0 * D_DROPS_INIT


def polygon_area(xs: List[float], ys: List[float]) -> float:
    """Shoelace formula."""
    n = len(xs)
    a = 0.0
    for i in range(n):
        j = (i + 1) % n
        a += xs[i] * ys[j] - xs[j] * ys[i]
    return abs(a) * 0.5


def parse_last_frame(csv_path: str) -> Tuple[int, List[float], float]:
    """
    Stream through csv_path once.
    Returns (last_frame_index, [droplet_areas], max_y_across_all_vertices).

    Skips comment lines (starting with #) and the obstacle header block.
    Robust to partial/truncated files.
    """
    last_frame = None
    last_rows: List[str] = []
    max_y = -1e30

    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            try:
                frame = int(float(parts[0]))
            except ValueError:
                continue

            if last_frame is None or frame > last_frame:
                last_frame = frame
                last_rows = [line]
            elif frame == last_frame:
                last_rows.append(line)

    if last_frame is None:
        return -1, [], -1e30

    areas: List[float] = []
    for line in last_rows:
        parts = line.split(",")
        if len(parts) < 4:
            continue
        try:
            nv = int(float(parts[2]))
        except ValueError:
            continue

        coords = parts[4:]
        if len(coords) < 2 * nv:
            continue

        xs, ys = [], []
        ok = True
        for i in range(nv):
            try:
                x = float(coords[2 * i])
                y = float(coords[2 * i + 1])
            except ValueError:
                ok = False
                break
            xs.append(x)
            ys.append(y)
            if y > max_y:
                max_y = y

        if not ok:
            continue
        a = polygon_area(xs, ys)
        if a > 0.0 and math.isfinite(a):
            areas.append(a)

    return last_frame, areas, max_y


def parse_filename(path: str):
    """
    Extract (depth, sr_float) from filename.
    Handles both 2-decimal (sr0p40) and 3-decimal (sr0p325) tags.
    Returns None if either field is missing.
    """
    d_match  = re.search(r"_d(\d+)_", path)
    sr_match = re.search(r"_sr(\d+p\d+)_", path)
    if not d_match or not sr_match:
        return None, None
    depth = int(d_match.group(1))
    sr    = float(sr_match.group(1).replace("p", "."))
    return depth, sr


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern",  default="results/d*/droplet_output_*.csv")
    ap.add_argument("--bins",     type=int, default=15)
    ap.add_argument("--grid-fig", default="daughter_area_grid.png")
    ap.add_argument("--means-fig",default="daughter_area_means.png")
    ap.add_argument("--out-csv",  default="daughter_areas_by_depth_sr.csv")
    args = ap.parse_args()

    files = sorted(glob.glob(args.pattern))
    if not files:
        raise SystemExit(f"No files matched: {args.pattern}")

    # Group by (depth, SR)
    groups = defaultdict(list)
    skipped = 0
    for path in files:
        depth, sr = parse_filename(path)
        if depth is None:
            skipped += 1
            continue
        groups[(depth, sr)].append(path)

    if skipped:
        print(f"Skipped {skipped} files with unparseable filenames.")

    depths = sorted({k[0] for k in groups})
    srs    = sorted({k[1] for k in groups})
    print(f"Depths: {depths}")
    print(f"SRs:    {srs}")

    # Collect normalized daughter areas per cell
    cell_data = {}   # (depth, sr) -> np.array of a/A0
    for (depth, sr), paths in groups.items():
        threshold = exit_y_threshold(depth)
        all_norm = []
        kept = dropped = 0
        for path in paths:
            _, areas, max_y = parse_last_frame(path)
            if not areas:
                dropped += 1
                continue
            if max_y > threshold:
                # Droplet never fully exited — timed out
                dropped += 1
                continue
            norm = np.array(areas) / A0_INIT
            # Keep only daughters (area < A0)
            norm = norm[norm < 1.0]
            all_norm.extend(norm.tolist())
            kept += 1

        arr = np.array(all_norm) if all_norm else np.array([])
        cell_data[(depth, sr)] = arr
        mu  = arr.mean() if arr.size else float("nan")
        std = arr.std()  if arr.size else float("nan")
        print(f"  d={depth:3d}  SR={sr:.3f}  kept={kept:4d}  dropped={dropped:3d}"
              f"  N={arr.size:6d}  mean={mu:.4f}  std={std:.4f}")

    # ------------------------------------------------------------------
    # Figure 1: grid of histograms
    # ------------------------------------------------------------------
    nr, nc = len(depths), len(srs)
    fig, axes = plt.subplots(
        nr, nc,
        figsize=(2.4 * nc + 0.8, 2.0 * nr + 0.8),
        sharex=True, sharey=False,
        squeeze=False,
    )
    bin_edges = np.linspace(0, 1.05, args.bins + 1)

    for ri, depth in enumerate(depths):
        for ci, sr in enumerate(srs):
            ax = axes[ri][ci]
            arr = cell_data.get((depth, sr), np.array([]))
            if arr.size:
                ax.hist(arr, bins=bin_edges, density=True, color="steelblue",
                        edgecolor="white", linewidth=0.3)
                mu = arr.mean()
                ax.axvline(mu, color="crimson", lw=1.0, ls="--")
                ax.text(0.97, 0.95, f"N={arr.size}\nμ={mu:.3f}",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=6)
            else:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="gray")

            if ri == 0:
                ax.set_title(f"SR={sr:.3f}", fontsize=8)
            if ci == 0:
                ax.set_ylabel(f"d={depth}", fontsize=7)
            ax.set_xlim(0, 1.05)
            ax.tick_params(labelsize=6)

    fig.supxlabel("Normalized daughter area  a/A₀", fontsize=9)
    fig.supylabel("Depth (D_obsts)", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.grid_fig, dpi=200)
    print(f"Saved {args.grid_fig}")

    # ------------------------------------------------------------------
    # Figure 2: mean ± std vs SR, one line per depth
    # ------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    cmap = plt.get_cmap("viridis", len(depths))

    for ri, depth in enumerate(depths):
        means, stds, sr_vals = [], [], []
        for sr in srs:
            arr = cell_data.get((depth, sr), np.array([]))
            if arr.size >= 5:
                means.append(arr.mean())
                stds.append(arr.std())
                sr_vals.append(sr)
        if sr_vals:
            ax2.errorbar(sr_vals, means, yerr=stds,
                         label=f"d={depth}", color=cmap(ri),
                         marker="o", ms=4, capsize=3, lw=1.2)

    ax2.set_xlabel("Spacing ratio SR", fontsize=11)
    ax2.set_ylabel("Mean normalized daughter area  ⟨a/A₀⟩", fontsize=11)
    ax2.set_ylim(0, 1.05)
    ax2.legend(fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(args.means_fig, dpi=200)
    print(f"Saved {args.means_fig}")

    # ------------------------------------------------------------------
    # Long-format CSV for downstream fitting
    # ------------------------------------------------------------------
    with open(args.out_csv, "w") as f:
        f.write("depth,spacing_ratio,normalized_area\n")
        for (depth, sr), arr in sorted(cell_data.items()):
            for val in arr:
                f.write(f"{depth},{sr:.4f},{val:.8f}\n")
    print(f"Saved {args.out_csv}")


if __name__ == "__main__":
    main()
