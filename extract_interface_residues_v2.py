#!/usr/bin/env python3
"""
extract_interface_residues.py
──────────────────────────────
For each completed AF3 pool hit, finds inter-chain residue contacts between
the PRC2 trimer chains (A/B/C = SET7/SUZ12/EED) and the candidate chain(s)
(D, E, F...) using AF3's own contact_probs matrix from confidences.json,
then extracts a +-5 residue window of local sequence context around each
contacting residue on both sides of the interface.

Looks up candidate protein names from the pool key CSV so each row in the
output includes both the chain letter and the actual protein name (NCU ID
or gene name).

No biopython or .cif files needed -- uses confidences.json only.

Usage
-----
    python3 extract_interface_residues_v2.py \
        --hits_csv PRC2_AF3_all_hits.csv \
        --pool_key PRC2_pool_key_lower.csv \
        --output_dir /scratch/.../PRC2_AF3_PooledJSON_output \
        --out interface_residues_table.csv \
        --contact_prob_thresh 0.5 \
        --window 5
"""

import argparse
import csv
import json
import os

TRIMER_CHAINS = {'A': 'SET7', 'B': 'SUZ12', 'C': 'EED'}
CHAIN_LETTERS = ['D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']


def load_pool_key(pool_key_path):
    """Returns {pool_name: {chain_letter: protein_name}} from the pool key CSV.
    Expects columns: pool_name, chainA, chainB, chainC, chainD, chainE, chainF, chainG
    (no header row, matching PRC2_pool_key_lower.csv format)."""
    pool_map = {}
    with open(pool_key_path) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            pool = row[0].strip().lower()
            # Columns 0=pool_name, 1=A(SET7), 2=B(SUZ12), 3=C(EED), 4+=candidates D,E,F...
            chain_names = {}
            for i, letter in enumerate(CHAIN_LETTERS):
                col_idx = i + 4
                if col_idx < len(row) and row[col_idx].strip():
                    chain_names[letter] = row[col_idx].strip()
            pool_map[pool] = chain_names
    return pool_map


def load_sequences_from_data_json(data_json_path):
    """Returns {chain_id: full_one_letter_sequence} from the AF3 input data.json."""
    with open(data_json_path) as f:
        d = json.load(f)
    chain_seqs = {}
    for seq_entry in d.get("sequences", []):
        prot = seq_entry.get("protein", {})
        ids = prot.get("id")
        seq = prot.get("sequence")
        if ids is None or seq is None:
            continue
        if isinstance(ids, str):
            ids = [ids]
        for cid in ids:
            chain_seqs[cid] = seq
    return chain_seqs


def window_seq(full_seq, resnum, window=5):
    """resnum is 1-indexed. Returns string with center residue in brackets."""
    idx = resnum - 1
    if idx < 0 or idx >= len(full_seq):
        return ""
    lo = max(0, idx - window)
    hi = min(len(full_seq), idx + window + 1)
    chars = []
    for i in range(lo, hi):
        if i == idx:
            chars.append(f"[{full_seq[i]}]")
        else:
            chars.append(full_seq[i])
    return "".join(chars)


def find_contacts(confidences_path, chain_seqs, chain_names,
                  contact_prob_thresh=0.5, window=5):
    """chain_names: {chain_letter: protein_name} from pool key."""
    with open(confidences_path) as f:
        conf = json.load(f)

    chain_ids = conf["token_chain_ids"]
    res_ids = conf["token_res_ids"]
    contact_probs = conf["contact_probs"]
    n = len(chain_ids)

    candidate_chain_letters = sorted(set(chain_ids) - set(TRIMER_CHAINS))

    rows = []
    seen_pairs = set()
    for i in range(n):
        ci = chain_ids[i]
        if ci not in TRIMER_CHAINS:
            continue
        for j in range(n):
            cj = chain_ids[j]
            if cj not in candidate_chain_letters:
                continue
            prob = contact_probs[i][j]
            if prob < contact_prob_thresh:
                continue

            tres = res_ids[i]
            cres = res_ids[j]
            key = (ci, tres, cj, cres)
            if key in seen_pairs:
                continue
            seen_pairs.add(key)

            trimer_seq = chain_seqs.get(ci, "")
            cand_seq = chain_seqs.get(cj, "")
            candidate_name = chain_names.get(cj, "unknown")

            rows.append({
                "subunit": TRIMER_CHAINS[ci],
                "subunit_chain": ci,
                "subunit_resnum": tres,
                "subunit_resname": trimer_seq[tres - 1] if tres - 1 < len(trimer_seq) else "?",
                "subunit_context": window_seq(trimer_seq, tres, window),
                "candidate_name": candidate_name,
                "candidate_chain": cj,
                "candidate_resnum": cres,
                "candidate_resname": cand_seq[cres - 1] if cres - 1 < len(cand_seq) else "?",
                "candidate_context": window_seq(cand_seq, cres, window),
                "contact_prob": round(prob, 3),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits_csv", required=True)
    ap.add_argument("--pool_key", required=True,
                    help="PRC2_pool_key_lower.csv — no header, cols: pool_name,A,B,C,D,E,F,G")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--out", default="interface_residues_table.csv")
    ap.add_argument("--contact_prob_thresh", type=float, default=0.5)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--sample", default="seed-1_sample-0")
    args = ap.parse_args()

    pool_map = load_pool_key(args.pool_key)
    print(f"Loaded pool key: {len(pool_map)} pools")

    pools = set()
    with open(args.hits_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pool = row.get("pool_name")
            if pool:
                pools.add(pool.strip().lower())
    print(f"Found {len(pools)} unique hit pools to process")

    all_rows = []
    missing = []
    for pool in sorted(pools):
        conf_path = os.path.join(args.output_dir, pool, args.sample, "confidences.json")
        data_path = os.path.join(args.output_dir, pool, f"{pool}_data.json")

        if not os.path.exists(conf_path) or not os.path.exists(data_path):
            missing.append(pool)
            continue

        try:
            chain_seqs  = load_sequences_from_data_json(data_path)
            chain_names = pool_map.get(pool, {})
            rows = find_contacts(conf_path, chain_seqs, chain_names,
                                 args.contact_prob_thresh, args.window)
            for r in rows:
                r["pool_name"] = pool
            all_rows.extend(rows)
            print(f"{pool}: {len(rows)} contact residue pairs")
        except Exception as e:
            print(f"{pool}: ERROR - {e}")

    if missing:
        print(f"\n{len(missing)} pools missing confidences.json/data.json (skipped):")
        for m in missing[:10]:
            print(f"  {m}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")

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

    print(f"\nWrote {len(all_rows)} contact rows to {args.out}")


if __name__ == "__main__":
    main()
