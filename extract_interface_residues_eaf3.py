#!/usr/bin/env python3
"""
extract_interface_residues_eaf3.py
────────────────────────────────────
Extracts inter-chain residue contacts between EAF3 (chain A) and
candidate chains (B, C, D, ...) using AF3's contact_probs matrix
from confidences.json.

Key differences from PRC2 version:
  - Bait = single chain A (EAF3), not a trimer
  - Candidates on chains B, C, D, E, F, G, H, I
  - Pool key has columns: bin_name, n_interactors, bin_total_aa, interactor_ids
  - Output folder structure: eaf3_poolXXX/eaf3_poolXXX_data.json
                             eaf3_poolXXX/seed-1_sample-0/confidences.json

Usage
-----
    python3 extract_interface_residues_eaf3.py \
        --hits_csv     EAF3_AF3_ipTM_scores.csv \
        --pool_key     EAF3_pool_key.xlsx \
        --output_dir   /lustre2/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI/EAF3_AF3_PooledJSON_output \
        --out          EAF3_interface_residues_table.csv \
        --contact_prob_thresh 0.5 \
        --window 5 \
        --sample seed-1_sample-0
"""

import argparse
import csv
import json
import os

BAIT_CHAIN   = 'A'
BAIT_NAME    = 'EAF3'
CAND_CHAINS  = ['B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']


def load_pool_key(pool_key_path):
    """
    Returns {pool_name: {chain_letter: protein_name}}.

    Supports both:
      - Excel (.xlsx): columns bin_name, n_interactors, bin_total_aa, interactor_ids
        interactor_ids is semicolon-separated list of NCU IDs
      - CSV: same column structure
    """
    pool_map = {}

    if pool_key_path.endswith('.xlsx'):
        try:
            import pandas as pd
            df = pd.read_excel(pool_key_path)
        except ImportError:
            raise ImportError("pandas + openpyxl required for .xlsx: pip install pandas openpyxl")

        for _, row in df.iterrows():
            pool = str(row['bin_name']).strip().lower()
            ids_raw = str(row.get('interactor_ids', '')).strip()
            ids = [x.strip() for x in ids_raw.split(';') if x.strip() and x != 'nan']
            chain_names = {}
            for i, cand_id in enumerate(ids):
                if i < len(CAND_CHAINS):
                    chain_names[CAND_CHAINS[i]] = cand_id
            pool_map[pool] = chain_names

    else:
        # CSV — try header detection
        with open(pool_key_path) as f:
            first = f.readline()
        has_header = 'bin_name' in first or 'pool' in first.lower()

        with open(pool_key_path) as f:
            reader = csv.reader(f)
            if has_header:
                next(reader)
            for row in reader:
                if not row:
                    continue
                pool = row[0].strip().lower()
                ids_raw = row[3] if len(row) > 3 else ''
                ids = [x.strip() for x in ids_raw.split(';') if x.strip()]
                chain_names = {}
                for i, cand_id in enumerate(ids):
                    if i < len(CAND_CHAINS):
                        chain_names[CAND_CHAINS[i]] = cand_id
                pool_map[pool] = chain_names

    print(f"Pool key loaded: {len(pool_map)} pools")
    return pool_map


def load_sequences_from_data_json(data_json_path):
    """Returns {chain_id: sequence} from AF3 data.json."""
    with open(data_json_path) as f:
        d = json.load(f)
    chain_seqs = {}
    for entry in d.get("sequences", []):
        prot = entry.get("protein", {})
        ids  = prot.get("id")
        seq  = prot.get("sequence")
        if ids is None or seq is None:
            continue
        if isinstance(ids, str):
            ids = [ids]
        for cid in ids:
            chain_seqs[cid] = seq
    return chain_seqs


def window_seq(full_seq, resnum, window=5):
    """1-indexed resnum. Returns sequence window with centre in [brackets]."""
    idx = resnum - 1
    if idx < 0 or idx >= len(full_seq):
        return ""
    lo = max(0, idx - window)
    hi = min(len(full_seq), idx + window + 1)
    chars = []
    for i in range(lo, hi):
        chars.append(f"[{full_seq[i]}]" if i == idx else full_seq[i])
    return "".join(chars)


def find_contacts(confidences_path, chain_seqs, chain_names,
                  contact_prob_thresh=0.5, window=5):
    """
    Finds contacts between EAF3 (chain A) and all candidate chains.
    chain_names: {chain_letter: protein_name}
    """
    with open(confidences_path) as f:
        conf = json.load(f)

    chain_ids    = conf["token_chain_ids"]
    res_ids      = conf["token_res_ids"]
    contact_probs = conf["contact_probs"]
    n = len(chain_ids)

    # Candidate chains = any chain that is not A
    candidate_chains = sorted(set(chain_ids) - {BAIT_CHAIN})

    rows = []
    seen = set()
    for i in range(n):
        if chain_ids[i] != BAIT_CHAIN:
            continue
        for j in range(n):
            cj = chain_ids[j]
            if cj not in candidate_chains:
                continue
            prob = contact_probs[i][j]
            if prob < contact_prob_thresh:
                continue
            tres  = res_ids[i]
            cres  = res_ids[j]
            key   = (BAIT_CHAIN, tres, cj, cres)
            if key in seen:
                continue
            seen.add(key)

            bait_seq  = chain_seqs.get(BAIT_CHAIN, "")
            cand_seq  = chain_seqs.get(cj, "")
            cand_name = chain_names.get(cj, f"chain_{cj}")

            rows.append({
                "subunit":            BAIT_NAME,
                "subunit_chain":      BAIT_CHAIN,
                "subunit_resnum":     tres,
                "subunit_resname":    bait_seq[tres - 1] if tres - 1 < len(bait_seq) else "?",
                "subunit_context":    window_seq(bait_seq, tres, window),
                "candidate_name":     cand_name,
                "candidate_chain":    cj,
                "candidate_resnum":   cres,
                "candidate_resname":  cand_seq[cres - 1] if cres - 1 < len(cand_seq) else "?",
                "candidate_context":  window_seq(cand_seq, cres, window),
                "contact_prob":       round(prob, 3),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits_csv",            required=True,
                    help="EAF3_AF3_ipTM_scores.csv — must have column 'pool_name'")
    ap.add_argument("--pool_key",            required=True,
                    help="EAF3_pool_key.xlsx or .csv")
    ap.add_argument("--output_dir",          required=True,
                    help="Path to EAF3_AF3_PooledJSON_output/")
    ap.add_argument("--out",                 default="EAF3_interface_residues_table.csv")
    ap.add_argument("--contact_prob_thresh", type=float, default=0.5)
    ap.add_argument("--window",              type=int,   default=5)
    ap.add_argument("--sample",              default="seed-1_sample-0")
    ap.add_argument("--all_pools",           action="store_true",
                    help="Process all pools, not just hits (ignores --hits_csv filter)")
    args = ap.parse_args()

    pool_map = load_pool_key(args.pool_key)

    # Get pools to process
    if args.all_pools:
        pools = set(pool_map.keys())
        print(f"Processing all {len(pools)} pools")
    else:
        pools = set()
        with open(args.hits_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                pool = row.get("pool_name", "").strip().lower()
                if pool:
                    pools.add(pool)
        print(f"Processing {len(pools)} hit pools from {args.hits_csv}")

    all_rows = []
    missing  = []
    errors   = []

    for pool in sorted(pools):
        # EAF3 output structure:
        # eaf3_poolXXX/
        #   eaf3_poolXXX_data.json
        #   seed-1_sample-0/
        #     confidences.json
        conf_path = os.path.join(args.output_dir, pool,
                                 args.sample, "confidences.json")
        data_path = os.path.join(args.output_dir, pool,
                                 f"{pool}_data.json")

        if not os.path.exists(conf_path) or not os.path.exists(data_path):
            missing.append(pool)
            continue

        try:
            chain_seqs  = load_sequences_from_data_json(data_path)
            chain_names = pool_map.get(pool, {})
            rows        = find_contacts(conf_path, chain_seqs, chain_names,
                                        args.contact_prob_thresh, args.window)
            for r in rows:
                r["pool_name"] = pool
            all_rows.extend(rows)
            if rows:
                print(f"  {pool}: {len(rows)} contact pairs")
        except Exception as e:
            errors.append((pool, str(e)))
            print(f"  {pool}: ERROR — {e}")

    if missing:
        print(f"\n{len(missing)} pools missing files (skipped):")
        for m in missing[:10]:
            print(f"  {m}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")

    if errors:
        print(f"\n{len(errors)} pools with errors:")
        for pool, err in errors[:5]:
            print(f"  {pool}: {err}")

    if not all_rows:
        print("No contacts found.")
        return

    fieldnames = ["pool_name", "subunit", "subunit_chain", "subunit_resnum",
                  "subunit_resname", "subunit_context", "candidate_name",
                  "candidate_chain", "candidate_resnum", "candidate_resname",
                  "candidate_context", "contact_prob"]

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nWrote {len(all_rows)} contact rows → {args.out}")
    print(f"Unique candidates with contacts: "
          f"{len(set(r['candidate_name'] for r in all_rows))}")


if __name__ == "__main__":
    main()
