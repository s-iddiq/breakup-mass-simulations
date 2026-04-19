#!/usr/bin/env python3
"""
fit_beta_delta_zi.py — Zero-inflated Beta model for daughter droplet areas

Model:
    P(X = 1)  = w                    (intact pass-throughs; delta mass at 1)
    X | X < 1 ~ Beta(alpha, beta)    (broken fragment distribution)

The split is principled: intact droplets land at exactly x = 1 within
floating-point precision (area is conserved exactly when no breakup occurs).
The threshold eps=1e-3 is numerical tolerance only, not a physical choice.

Usage:
    python3 fit_beta_delta_zi.py
    python3 fit_beta_delta_zi.py --deep 50 --eps 1e-3
    python3 fit_beta_delta_zi.py --csv daughter_areas_by_depth_sr.csv

Input:
    daughter_areas_by_depth_sr.csv   (produced by Analysis.py)
    Columns: depth, spacing_ratio, normalized_area

Output:
    zi_beta_fits.png     — histograms with fitted ZI-Beta overlaid per SR
    zi_beta_params.png   — w(SR), alpha(SR), beta(SR) parameter plots
    zi_beta_params.csv   — table of fitted parameters
"""

import argparse
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Model fitting
# ---------------------------------------------------------------------------

def fit_zi_beta(x: np.ndarray, eps: float = 1e-3):
    """
    Fit zero-inflated Beta to normalized area data.

    Parameters
    ----------
    x   : array of a/A0 values (should be in (0, 1])
    eps : numerical tolerance for identifying intact droplets (x >= 1 - eps)

    Returns
    -------
    w     : float, fraction of intact pass-throughs
    alpha : float, Beta shape parameter (bulk fragment distribution)
    beta  : float, Beta shape parameter (bulk fragment distribution)
    bulk  : np.ndarray, the fragment data (x < 1 - eps)
    """
    intact_mask = x >= (1.0 - eps)
    w    = intact_mask.mean()
    bulk = x[~intact_mask]

    if bulk.size < 10:
        return w, float("nan"), float("nan"), bulk

    # Clip to open interval to avoid boundary issues in Beta MLE
    bulk_c = np.clip(bulk, 1e-6, 1.0 - 1e-6)
    try:
        alpha, beta_p, _, _ = stats.beta.fit(bulk_c, floc=0, fscale=1)
    except Exception:
        alpha, beta_p = float("nan"), float("nan")

    return w, alpha, beta_p, bulk


def mom_beta(x: np.ndarray):
    """Method-of-moments Beta fit. More robust than MLE for truncated data."""
    mu  = x.mean()
    var = x.var(ddof=1)
    if var <= 0 or mu <= 0 or mu >= 1:
        return float("nan"), float("nan")
    c = mu * (1.0 - mu) / var - 1.0
    if c <= 0:
        return float("nan"), float("nan")
    return mu * c, (1.0 - mu) * c


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",   default="daughter_areas_by_depth_sr.csv")
    ap.add_argument("--deep",  type=int, default=50,
                    help="Pool depths >= this value for steady-state fits.")
    ap.add_argument("--eps",   type=float, default=1e-3,
                    help="Numerical tolerance for intact-droplet identification.")
    ap.add_argument("--bins",  type=int, default=20)
    ap.add_argument("--hist-fig",   default="zi_beta_fits.png")
    ap.add_argument("--param-fig",  default="zi_beta_params.png")
    ap.add_argument("--param-csv",  default="zi_beta_params.csv")
    args = ap.parse_args()

    # Load data
    data = np.genfromtxt(args.csv, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    depths = data[:, 0]
    srs    = data[:, 1]
    areas  = data[:, 2]

    all_srs = sorted(set(srs.tolist()))

    # Collect per-SR steady-state data
    sr_data = {}
    for sr in all_srs:
        mask = (depths >= args.deep) & (srs == sr)
        arr  = areas[mask]
        arr  = arr[(arr > 0) & np.isfinite(arr)]
        if arr.size >= 20:
            sr_data[sr] = arr

    if not sr_data:
        raise SystemExit(f"No data with depth >= {args.deep}. "
                         f"Try --deep with a lower value.")

    # Fit ZI-Beta for each SR
    results = {}
    print(f"\n{'SR':>6}  {'N':>6}  {'w':>6}  {'alpha_MLE':>10}  "
          f"{'beta_MLE':>10}  {'alpha_MoM':>10}  {'beta_MoM':>10}")
    print("-" * 70)
    for sr, arr in sorted(sr_data.items()):
        w, al, be, bulk = fit_zi_beta(arr, eps=args.eps)
        am, bm = mom_beta(bulk) if bulk.size >= 10 else (float("nan"), float("nan"))
        results[sr] = dict(sr=sr, N=arr.size, w=w,
                           alpha_mle=al, beta_mle=be,
                           alpha_mom=am, beta_mom=bm,
                           bulk=bulk)
        print(f"{sr:6.3f}  {arr.size:6d}  {w:6.3f}  "
              f"{al:10.3f}  {be:10.3f}  {am:10.3f}  {bm:10.3f}")

    # ------------------------------------------------------------------
    # Figure 1: histogram + ZI-Beta overlay per SR
    # ------------------------------------------------------------------
    n = len(results)
    ncols = min(n, 4)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.5 * ncols, 3.0 * nrows),
                             squeeze=False)
    bin_edges = np.linspace(0, 1.0, args.bins + 1)
    bw = bin_edges[1] - bin_edges[0]
    xgrid = np.linspace(1e-4, 1.0 - 1e-4, 400)

    for idx, (sr, res) in enumerate(sorted(results.items())):
        ri, ci = divmod(idx, ncols)
        ax = axes[ri][ci]
        arr  = sr_data[sr]
        bulk = res["bulk"]
        w    = res["w"]

        # Histogram of all data
        ax.hist(arr, bins=bin_edges, density=True,
                color="steelblue", alpha=0.6, edgecolor="white", lw=0.3,
                label="data")

        # Delta mass at x=1 rendered as rectangle
        if w > 0 and bulk.size > 0:
            ax.bar(1.0 - bw / 2, w / bw, width=bw,
                   color="crimson", alpha=0.5, label=f"w={w:.3f}")

        # Beta PDF on bulk (MoM)
        am, bm = res["alpha_mom"], res["beta_mom"]
        if math.isfinite(am) and am > 0 and bm > 0:
            pdf = stats.beta.pdf(xgrid, am, bm) * (1.0 - w)
            ax.plot(xgrid, pdf, "k-", lw=1.5,
                    label=f"β({am:.2f},{bm:.2f})")

        ax.set_title(f"SR={sr:.3f}  N={arr.size}", fontsize=9)
        ax.set_xlabel("a/A₀", fontsize=8)
        ax.set_ylabel("density", fontsize=8)
        ax.legend(fontsize=7)
        ax.set_xlim(0, 1.05)

    # Hide unused subplots
    for idx in range(len(results), nrows * ncols):
        ri, ci = divmod(idx, ncols)
        axes[ri][ci].set_visible(False)

    fig.suptitle(f"Zero-inflated Beta fits  (depths ≥ {args.deep})", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.hist_fig, dpi=200)
    print(f"\nSaved {args.hist_fig}")

    # ------------------------------------------------------------------
    # Figure 2: parameter trends w(SR), alpha(SR), beta(SR)
    # ------------------------------------------------------------------
    sr_vals  = sorted(results.keys())
    ws       = [results[s]["w"]         for s in sr_vals]
    alphas   = [results[s]["alpha_mom"] for s in sr_vals]
    betas    = [results[s]["beta_mom"]  for s in sr_vals]

    fig2, axes2 = plt.subplots(1, 3, figsize=(11, 3.5))

    axes2[0].plot(sr_vals, ws, "o-", color="crimson")
    axes2[0].set_xlabel("SR"); axes2[0].set_ylabel("w  (pass-through fraction)")
    axes2[0].set_title("Intact pass-throughs vs SR")

    axes2[1].plot(sr_vals, alphas, "o-", color="steelblue")
    axes2[1].set_xlabel("SR"); axes2[1].set_ylabel("α  (Beta shape)")
    axes2[1].set_title("α(SR)  [MoM]")

    axes2[2].plot(sr_vals, betas, "o-", color="seagreen")
    axes2[2].set_xlabel("SR"); axes2[2].set_ylabel("β  (Beta shape)")
    axes2[2].set_title("β(SR)  [MoM]")

    for ax in axes2:
        ax.grid(True, alpha=0.3)

    fig2.tight_layout()
    fig2.savefig(args.param_fig, dpi=200)
    print(f"Saved {args.param_fig}")

    # ------------------------------------------------------------------
    # Parameter CSV
    # ------------------------------------------------------------------
    with open(args.param_csv, "w") as f:
        f.write("sr,N,w,alpha_mle,beta_mle,alpha_mom,beta_mom\n")
        for sr in sorted(results.keys()):
            r = results[sr]
            f.write(f"{r['sr']:.4f},{r['N']},{r['w']:.6f},"
                    f"{r['alpha_mle']:.6f},{r['beta_mle']:.6f},"
                    f"{r['alpha_mom']:.6f},{r['beta_mom']:.6f}\n")
    print(f"Saved {args.param_csv}")


if __name__ == "__main__":
    main()
