#!/usr/bin/env python3
"""
generate_individual_bait_pools.py
──────────────────────────────────
Generates AF3 JSON pool inputs for screening the N. crassa proteome
against four individual bait proteins:
  - SET7  (NCU07496) — PRC2 catalytic subunit, H3K27me3 writer
  - SUZ12 (NCU05460) — PRC2 structural subunit
  - EED   (NCU05300) — PRC2 H3K27me3 reader subunit
  - EAF3  (NCU06787) — H3K36me3 chromodomain/MRG reader

Each bait is chain A in every pool. Candidate proteins fill chains
B, C, D... up to --max_aa amino acid budget.

One pool set is generated per bait → 4 output directories, each with:
  - <BAIT>_pool001.json ... <BAIT>_poolNNN.json
  - <BAIT>_pool_key.csv   (pool_name, n_candidates, total_aa,
                            candidate_chain_B, candidate_chain_C, ...)
  - <BAIT>_pool_key.xlsx  (same, Excel format)

Usage
-----
    python3 generate_individual_bait_pools.py \\
        --proteome_fasta  Ncrassa_proteome.fasta \\
        --bait_fasta      baits.fasta \\
        --outdir          individual_bait_pools \\
        --max_aa          5000 \\
        --max_candidates  4 \\
        --seeds           42 137

    # Or run one bait at a time:
    python3 generate_individual_bait_pools.py \\
        --proteome_fasta  Ncrassa_proteome.fasta \\
        --bait_fasta      baits.fasta \\
        --baits           SET7 EED \\
        --outdir          individual_bait_pools \\
        --max_aa          5000

Input FASTA formats
-------------------
--proteome_fasta : standard FASTA, one entry per protein
    Header can be:
      >NCU07496 ...
      >NCU07496|gene_name|...
      >gene_name|NCU07496|...
    The script extracts the first field as the protein ID.

--bait_fasta : FASTA with exactly one entry per bait.
    Headers must match bait names:
      >SET7
      >SUZ12
      >EED
      >EAF3
    (case-insensitive)

Alternatively, embed bait sequences directly via --bait_sequences
(see BAIT_SEQUENCES dict below).
"""

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path


# ── Hard-coded bait sequences (optional fallback) ─────────────────────
# Fill in if you don't want to supply a --bait_fasta file.
# Leave empty string to require --bait_fasta.
BAIT_SEQUENCES = {
    "SET7":  "",   # NCU07496
    "SUZ12": "",   # NCU05460
    "EED":   "",   # NCU05300
    "EAF3":  "",   # NCU06787
}

BAIT_NCUS = {
    "SET7":  "NCU07496",
    "SUZ12": "NCU05460",
    "EED":   "NCU05300",
    "EAF3":  "NCU06787",
}

CHAIN_LETTERS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


def parse_fasta(path):
    """Returns {header_id: sequence}. Uses first whitespace-delimited token as ID."""
    seqs = {}
    curr_id = None
    curr_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if curr_id:
                    seqs[curr_id] = "".join(curr_seq)
                raw = line[1:].split()[0]
                curr_id = raw.split("|")[0].strip()
                curr_seq = []
            elif curr_id:
                curr_seq.append(re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "",
                                       line.upper()))
    if curr_id:
        seqs[curr_id] = "".join(curr_seq)
    return seqs


def clean_seq(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def make_json(pool_name, bait_name, bait_seq, candidates, seeds):
    """
    Creates AF3 JSON with bait as chain A, candidates as B, C, D...
    candidates: list of (protein_id, sequence)
    """
    sequences = [
        {"protein": {"id": "A", "sequence": bait_seq}}
    ]
    for i, (_, seq) in enumerate(candidates):
        sequences.append({
            "protein": {
                "id": CHAIN_LETTERS[i + 1],  # B, C, D...
                "sequence": seq
            }
        })
    return {
        "name":        pool_name,
        "dialect":     "alphafold3",
        "version":     1,
        "modelSeeds":  seeds,
        "sequences":   sequences,
    }


def generate_pools(bait_name, bait_seq, proteome, max_aa, max_candidates, seeds,
                   outdir, exclude_ids=None):
    """
    Generates pool JSONs for one bait against the proteome.
    Returns list of pool key rows.
    """
    exclude_ids = exclude_ids or set()
    bait_aa = len(bait_seq)
    budget   = max_aa - bait_aa  # remaining aa after bait

    if budget <= 0:
        print(f"  WARNING: bait {bait_name} ({bait_aa} aa) already exceeds "
              f"max_aa={max_aa}. No candidates can fit.", file=sys.stderr)
        return []

    os.makedirs(outdir, exist_ok=True)

    # Sort candidates by aa length ascending so smaller proteins pack better
    candidates_all = [
        (pid, seq) for pid, seq in sorted(proteome.items(),
                                           key=lambda x: len(x[1]))
        if pid not in exclude_ids and len(seq) >= 20
    ]

    pool_rows = []
    pool_num  = 0
    i = 0

    while i < len(candidates_all):
        pool_num += 1
        pool_name = f"{bait_name.lower()}_pool{pool_num:03d}"

        current_candidates = []
        current_aa = 0

        j = i
        while j < len(candidates_all) and len(current_candidates) < max_candidates:
            pid, seq = candidates_all[j]
            if current_aa + len(seq) <= budget:
                current_candidates.append((pid, seq))
                current_aa += len(seq)
                j += 1
            else:
                # Try next smaller protein
                found = False
                for k in range(j + 1, min(j + 50, len(candidates_all))):
                    kpid, kseq = candidates_all[k]
                    if current_aa + len(kseq) <= budget:
                        current_candidates.append((kpid, kseq))
                        current_aa += len(kseq)
                        candidates_all.pop(k)
                        found = True
                        break
                if not found:
                    break  # nothing fits, close pool

        if not current_candidates:
            i = j
            continue

        # Write JSON
        pool_json = make_json(pool_name, bait_name, bait_seq,
                               current_candidates, seeds)
        json_path = os.path.join(outdir, f"{pool_name}.json")
        with open(json_path, "w") as f:
            json.dump(pool_json, f, indent=2)

        # Pool key row
        row = {
            "pool_name":    pool_name,
            "bait":         bait_name,
            "bait_ncu":     BAIT_NCUS.get(bait_name, ""),
            "n_candidates": len(current_candidates),
            "total_aa":     bait_aa + current_aa,
        }
        for ci, (pid, _) in enumerate(current_candidates):
            row[f"chain_{CHAIN_LETTERS[ci+1]}"] = pid
        pool_rows.append(row)

        i = j  # advance past consumed candidates

    print(f"  {bait_name}: {pool_num} pools, {len(candidates_all)} candidates "
          f"(bait={bait_aa} aa, budget={budget} aa/pool)")
    return pool_rows


def write_pool_key(rows, outdir, bait_name, max_candidates):
    if not rows:
        return

    # Determine all chain columns used
    chain_cols = [f"chain_{CHAIN_LETTERS[i+1]}" for i in range(max_candidates)]
    fieldnames = ["pool_name", "bait", "bait_ncu", "n_candidates",
                  "total_aa"] + chain_cols

    csv_path = os.path.join(outdir, f"{bait_name}_pool_key.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Pool key CSV: {csv_path}")

    # Excel version
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        for col in chain_cols:
            if col not in df.columns:
                df[col] = ""
        df = df[fieldnames]
        xlsx_path = os.path.join(outdir, f"{bait_name}_pool_key.xlsx")
        df.to_excel(xlsx_path, index=False)
        print(f"  Pool key XLSX: {xlsx_path}")
    except ImportError:
        print("  (pandas not available — skipping .xlsx)")

    # Also write scraper-compatible no-header CSV (lowercase pool names)
    scraper_path = os.path.join(outdir, f"{bait_name}_pool_key_scraper.csv")
    with open(scraper_path, "w", newline="") as f:
        writer = csv.writer(f)
        for row in rows:
            line = [row["pool_name"]]
            for col in chain_cols:
                line.append(row.get(col, ""))
            writer.writerow(line)
    print(f"  Scraper CSV:   {scraper_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Generate AF3 pool JSONs for 4-bait individual subunit screens"
    )
    ap.add_argument("--proteome_fasta", required=True,
                    help="FASTA of all N. crassa proteins to screen")
    ap.add_argument("--bait_fasta",     default=None,
                    help="FASTA containing bait sequences (headers = bait names)")
    ap.add_argument("--baits",          nargs="+",
                    default=["SET7","SUZ12","EED","EAF3"],
                    help="Which baits to generate pools for (default: all 4)")
    ap.add_argument("--outdir",         default="individual_bait_pools",
                    help="Root output directory (one subdir per bait)")
    ap.add_argument("--max_aa",         type=int, default=5000,
                    help="Max total aa per pool including bait (default 5000)")
    ap.add_argument("--max_candidates", type=int, default=4,
                    help="Max candidate chains per pool (default 4)")
    ap.add_argument("--seeds",          type=int, nargs="+", default=[42, 137],
                    help="AF3 model seeds (default: 42 137)")
    ap.add_argument("--exclude_baits",  action="store_true",
                    help="Exclude all 4 bait proteins from candidate pools")
    args = ap.parse_args()

    # Load proteome
    print("Loading proteome FASTA...")
    proteome = parse_fasta(args.proteome_fasta)
    proteome = {pid: clean_seq(seq) for pid, seq in proteome.items()
                if len(clean_seq(seq)) >= 20}
    print(f"  {len(proteome)} proteins loaded")

    # Load bait sequences
    bait_seqs = {}
    if args.bait_fasta:
        print("Loading bait FASTA...")
        raw_baits = parse_fasta(args.bait_fasta)
        # Match case-insensitively
        for bait in args.baits:
            for key, seq in raw_baits.items():
                if key.upper() == bait.upper():
                    bait_seqs[bait] = clean_seq(seq)
                    print(f"  {bait}: {len(bait_seqs[bait])} aa")
                    break

    # Fall back to hard-coded sequences
    for bait in args.baits:
        if bait not in bait_seqs:
            if BAIT_SEQUENCES.get(bait):
                bait_seqs[bait] = clean_seq(BAIT_SEQUENCES[bait])
                print(f"  {bait}: {len(bait_seqs[bait])} aa (hard-coded)")
            else:
                print(f"  ERROR: no sequence found for bait '{bait}'. "
                      f"Provide --bait_fasta or fill BAIT_SEQUENCES dict.",
                      file=sys.stderr)
                sys.exit(1)

    # Exclude bait NCU IDs from candidate pools
    exclude_ids = set()
    if args.exclude_baits:
        exclude_ids = set(BAIT_NCUS.values())
        print(f"Excluding bait proteins from candidates: {exclude_ids}")

    # Generate pools per bait
    os.makedirs(args.outdir, exist_ok=True)
    summary = {}

    for bait in args.baits:
        print(f"\nGenerating pools for {bait}...")
        bait_outdir = os.path.join(args.outdir, f"{bait}_pools")
        os.makedirs(bait_outdir, exist_ok=True)

        rows = generate_pools(
            bait_name      = bait,
            bait_seq       = bait_seqs[bait],
            proteome       = dict(proteome),   # fresh copy
            max_aa         = args.max_aa,
            max_candidates = args.max_candidates,
            seeds          = args.seeds,
            outdir         = bait_outdir,
            exclude_ids    = exclude_ids,
        )

        write_pool_key(rows, bait_outdir, bait, args.max_candidates)
        summary[bait] = len(rows)

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    total_pools = 0
    for bait, n in summary.items():
        print(f"  {bait:<8} {n:>5} pools → {args.outdir}/{bait}_pools/")
        total_pools += n
    print(f"  {'TOTAL':<8} {total_pools:>5} pool JSONs")
    print(f"\nNext steps:")
    print(f"  1. scp -r {args.outdir}/ ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/.../PooledPPI/")
    print(f"  2. Update PROTEIN, INPUT_DIR, OUTPUT_DIR, TOTAL in PRC2_AF3_Pools.sh for each bait")
    print(f"  3. sbatch PRC2_AF3_Pools.sh  (once per bait or modify to loop)")


if __name__ == "__main__":
    main()
