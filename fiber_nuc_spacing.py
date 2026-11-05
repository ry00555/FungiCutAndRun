#!/usr/bin/env python3
"""
Nucleosome spacing from Fiber-seq nucleosome calls.

Two modes:

  sample   — one sample: spacing tables + plots
      python fiber_nuc_spacing.py sample \
          --nucs SAMPLE.nucleosomes.bed.gz --sample SAMPLE --out OUTDIR

  compare  — overlay all samples that 'sample' mode has finished
      python fiber_nuc_spacing.py compare --in RESULTS_DIR --out SUMMARY_DIR

Input BED (made by FiberTools_Pt4_NucAnalysis.sh from `ft extract --nuc`),
one line per nucleosome, reference coordinates:
  chrom  start  end  fiber  nuc_len  strand  fiber_start  fiber_end

What is measured (all within single fibers, so no averaging across cells):
  * adjacent inter-dyad distance — dyad of nucleosome i to i+1 (= local repeat length)
  * linker length               — gap between neighbouring nucleosomes
  * dyad-distance ladder        — ALL dyad pairs on a fiber up to --max-lag bp.
    Regular arrays give peaks at 1x, 2x, 3x ... NRL; the slope of
    peak position vs peak number is the nucleosome repeat length (NRL),
    the single-molecule equivalent of an MNase ladder.
  * per-fiber regularity        — median spacing and its coefficient of variation (CV)
  * nucleosome density          — nucleosomes per kb of fiber
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLS = ["chrom", "start", "end", "fiber", "nuc_len", "strand", "fiber_start", "fiber_end"]


# ─────────────────────────────────────────────────────────────────────────────
def smooth(y, w=21):
    if w <= 1:
        return y.astype(float)
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")


def ladder_peaks(hist, min_first=100, window=60, max_peaks=5):
    """Local maxima of the smoothed all-pairs histogram, spaced >= window bp."""
    s = smooth(hist, 21)
    peaks = []
    for i in range(min_first, len(s) - window):
        lo, hi = max(0, i - window), i + window + 1
        if s[i] > 0 and s[i] == s[lo:hi].max():
            if not peaks or i - peaks[-1] > window:
                peaks.append(i)
        if len(peaks) >= max_peaks:
            break
    return peaks, s


def run_sample(args):
    os.makedirs(args.out, exist_ok=True)
    S = args.sample
    print(f"[{S}] reading {args.nucs}")
    df = pd.read_csv(args.nucs, sep="\t", header=None, names=COLS,
                     dtype={"chrom": str, "fiber": str, "strand": str})
    if args.max_nuc_len:
        df = df[df["nuc_len"] <= args.max_nuc_len]      # mononucleosomes only
    if df.empty:
        sys.exit(f"[{S}] no nucleosomes in {args.nucs}")

    df = df.sort_values(["fiber", "start"], kind="mergesort").reset_index(drop=True)
    dyad = ((df["start"].values + df["end"].values) // 2).astype(np.int64)
    fib = pd.factorize(df["fiber"])[0]
    start = df["start"].values
    end = df["end"].values

    # ── adjacent distances + linkers ─────────────────────────────────────────
    same = fib[1:] == fib[:-1]
    adj = (dyad[1:] - dyad[:-1])[same]
    linker = (start[1:] - end[:-1])[same]

    # ── all-pairs ladder ─────────────────────────────────────────────────────
    lag = args.max_lag
    ladder = np.zeros(lag + 1, dtype=np.int64)
    k = 1
    while True:
        if k >= len(dyad):
            break
        same_k = fib[k:] == fib[:-k]
        d = (dyad[k:] - dyad[:-k])[same_k]
        d = d[(d >= 0) & (d <= lag)]
        if d.size == 0:
            break
        ladder += np.bincount(d, minlength=lag + 1)[: lag + 1]
        k += 1

    peaks, ladder_s = ladder_peaks(ladder)
    nrl = np.nan
    if len(peaks) >= 2:
        n = np.arange(1, len(peaks) + 1)
        nrl = float(np.polyfit(n, peaks, 1)[0])

    # ── per-fiber metrics ────────────────────────────────────────────────────
    fiber_df = pd.DataFrame({"fib": fib, "dyad": dyad,
                             "fs": df["fiber_start"].values, "fe": df["fiber_end"].values,
                             "fiber": df["fiber"].values})
    rows = []
    for f, g in fiber_df.groupby("fib", sort=False):
        dy = g["dyad"].values
        span = max(int(g["fe"].iloc[0]) - int(g["fs"].iloc[0]), 1)
        d = np.diff(dy)
        d = d[(d >= 100) & (d <= 260)]          # one-nucleosome steps only (skip gaps from missed calls)
        med = float(np.median(d)) if d.size else np.nan
        cv = float(np.std(d) / np.mean(d)) if d.size >= 3 else np.nan
        rows.append((g["fiber"].iloc[0], len(dy), span, len(dy) / span * 1000, med, cv))
    pf = pd.DataFrame(rows, columns=["fiber", "n_nucs", "fiber_span_bp",
                                     "nucs_per_kb", "median_spacing", "spacing_cv"])
    pf.to_csv(os.path.join(args.out, f"{S}.per_fiber_spacing.tsv.gz"),
              sep="\t", index=False)

    # ── histograms ───────────────────────────────────────────────────────────
    dist = np.arange(lag + 1)
    adj_c = np.bincount(adj[(adj >= 0) & (adj <= lag)], minlength=lag + 1)[: lag + 1]
    pd.DataFrame({"distance": dist, "adjacent_dyad": adj_c, "all_pairs_dyad": ladder}) \
      .to_csv(os.path.join(args.out, f"{S}.dyad_distance_hist.tsv"), sep="\t", index=False)
    lk = linker[(linker >= 0) & (linker <= 500)]
    pd.DataFrame({"linker_bp": np.arange(501),
                  "count": np.bincount(lk, minlength=501)[:501]}) \
      .to_csv(os.path.join(args.out, f"{S}.linker_hist.tsv"), sep="\t", index=False)

    # ── summary ──────────────────────────────────────────────────────────────
    adj_in = adj[(adj >= 100) & (adj <= 500)]
    mode_adj = int(np.argmax(smooth(adj_c, 11)[100:501]) + 100) if adj_in.size else np.nan
    multi = pf[pf["n_nucs"] >= 5]
    summ = {
        "sample": S,
        "fibers": len(pf),
        "nucleosomes": len(df),
        "nucs_per_kb_median": round(float(pf["nucs_per_kb"].median()), 3),
        "nuc_length_median": float(np.median(df["nuc_len"])),
        "adjacent_dyad_median": float(np.median(adj_in)) if adj_in.size else np.nan,
        "adjacent_dyad_mode": mode_adj,
        "linker_median": float(np.median(lk)) if lk.size else np.nan,
        "NRL_ladder_fit": round(nrl, 1) if not np.isnan(nrl) else np.nan,
        "ladder_peaks": ",".join(map(str, peaks)),
        "fibers_with_5plus_nucs": len(multi),
        "per_fiber_spacing_cv_median": round(float(multi["spacing_cv"].median()), 3)
                                        if len(multi) else np.nan,
    }
    pd.DataFrame([summ]).to_csv(os.path.join(args.out, f"{S}.spacing_summary.tsv"),
                                sep="\t", index=False)
    print(f"[{S}] " + "  ".join(f"{k}={v}" for k, v in summ.items() if k != "sample"))

    # ── plots ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"{S} — nucleosome spacing ({summ['fibers']:,} fibers, "
                 f"{summ['nucleosomes']:,} nucleosomes)", fontweight="bold")

    a = ax[0, 0]
    a.plot(dist, smooth(adj_c, 5), color="#2166AC")
    a.axvline(summ["adjacent_dyad_mode"], ls="--", c="gray",
              label=f"mode {summ['adjacent_dyad_mode']} bp")
    a.set(xlim=(0, 600), xlabel="Adjacent dyad-to-dyad distance (bp)", ylabel="Count",
          title="Local repeat length (neighbouring nucleosomes)")
    a.legend()

    a = ax[0, 1]
    a.plot(np.arange(501), smooth(np.bincount(lk, minlength=501)[:501], 5), color="#1B7837")
    a.axvline(summ["linker_median"], ls="--", c="gray",
              label=f"median {summ['linker_median']:.0f} bp")
    a.set(xlim=(0, 300), xlabel="Linker length (bp)", ylabel="Count", title="Linker lengths")
    a.legend()

    a = ax[1, 0]
    a.plot(dist, ladder_s, color="#762A83")
    for p in peaks:
        a.axvline(p, ls=":", c="gray")
    a.set(xlim=(0, lag), xlabel="Dyad-to-dyad distance, all pairs on a fiber (bp)",
          ylabel="Count (smoothed)",
          title=f"Spacing ladder — NRL ≈ {summ['NRL_ladder_fit']} bp")

    a = ax[1, 1]
    m = multi.dropna(subset=["median_spacing", "spacing_cv"])
    if len(m):
        a.hexbin(m["median_spacing"], m["spacing_cv"], gridsize=40, mincnt=1,
                 cmap="viridis", bins="log")
    a.set(xlim=(100, 400), ylim=(0, 1), xlabel="Per-fiber median spacing (bp)",
          ylabel="Spacing CV (lower = more regular)",
          title="Per-fiber regularity (fibers with ≥5 nucleosomes)")

    for a in ax.flat:
        a.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    out = os.path.join(args.out, f"{S}.nuc_spacing.png")
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close()
    print(f"[{S}] wrote {out}")


# ─────────────────────────────────────────────────────────────────────────────
def run_compare(args):
    os.makedirs(args.out, exist_ok=True)
    summ = sorted(glob.glob(os.path.join(args.inp, "**", "*.spacing_summary.tsv"), recursive=True))
    if not summ:
        sys.exit(f"No *.spacing_summary.tsv found under {args.inp}")
    tab = pd.concat([pd.read_csv(f, sep="\t") for f in summ], ignore_index=True)
    tab.to_csv(os.path.join(args.out, "all_samples.spacing_summary.tsv"), sep="\t", index=False)
    print(tab.to_string(index=False))

    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    cmap = plt.get_cmap("tab10" if len(summ) <= 10 else "tab20")
    for i, f in enumerate(summ):
        S = os.path.basename(f).replace(".spacing_summary.tsv", "")
        d = os.path.dirname(f)
        h = pd.read_csv(os.path.join(d, f"{S}.dyad_distance_hist.tsv"), sep="\t")
        l = pd.read_csv(os.path.join(d, f"{S}.linker_hist.tsv"), sep="\t")
        c = cmap(i % cmap.N)
        ax[0].plot(h["distance"], smooth(h["adjacent_dyad"].values, 5) /
                   max(h["adjacent_dyad"].sum(), 1), color=c, label=S, lw=1.2)
        ax[1].plot(l["linker_bp"], smooth(l["count"].values, 5) /
                   max(l["count"].sum(), 1), color=c, lw=1.2)
        ax[2].plot(h["distance"], smooth(h["all_pairs_dyad"].values, 21) /
                   max(h["all_pairs_dyad"].sum(), 1), color=c, lw=1.2)
    ax[0].set(xlim=(0, 600), xlabel="Adjacent dyad distance (bp)", ylabel="Fraction",
              title="Local repeat length")
    ax[1].set(xlim=(0, 300), xlabel="Linker length (bp)", title="Linker lengths")
    ax[2].set(xlim=(0, 1000), xlabel="Dyad distance, all pairs (bp)", title="Spacing ladder")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    ax[0].legend(fontsize=7, loc="upper right")
    plt.tight_layout()
    out = os.path.join(args.out, "all_samples.nuc_spacing_compare.png")
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close()

    fig, a = plt.subplots(figsize=(max(6, len(tab) * 0.6), 5))
    a.bar(range(len(tab)), tab["NRL_ladder_fit"], color="#2166AC")
    a.set_xticks(range(len(tab)))
    a.set_xticklabels(tab["sample"], rotation=60, ha="right", fontsize=7)
    a.set(ylabel="NRL (bp, ladder fit)", title="Nucleosome repeat length per sample")
    a.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out, "all_samples.NRL.png"), dpi=170, bbox_inches="tight")
    plt.close()
    print(f"wrote comparison plots to {args.out}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("sample")
    s.add_argument("--nucs", required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--max-lag", type=int, default=1000)
    s.add_argument("--max-nuc-len", type=int, default=0,
                   help="ignore nucleosome calls longer than this, e.g. dinucleosomes (0 = keep all)")
    c = sub.add_parser("compare")
    c.add_argument("--in", dest="inp", required=True)
    c.add_argument("--out", required=True)
    args = p.parse_args()
    run_sample(args) if args.mode == "sample" else run_compare(args)


if __name__ == "__main__":
    main()
