#!/usr/bin/env python3
"""
AF3pool_scraper_v1.py
────────────────────
Extracts per-chain ipTM and PAE scores from pooled AlphaFold3 outputs.
Combines chain-pair matrix extraction (from scraper) with global iptm,
ranking_score, mean, and std across all 5 samples.

Usage:
    python3 AF3pool_scraper_v1.py \
        key_csv \
        base_dir \
        output_csv \
        --row-index 1
"""

import os
import json
import csv
import argparse
import glob
import statistics


def load_key_csv(key_csv_path):
    """Load key CSV mapping top-level folder names to chain names."""
    mapping = {}
    with open(key_csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) > 1:
                folder_name = row[0].strip()
                chains = [c.strip() for c in row[1:] if c.strip()]
                mapping[folder_name] = chains
    return mapping


def extract_row(json_path, key, row_index):
    """Extract a specific row (0-based) of a matrix from JSON by key."""
    with open(json_path) as f:
        data = json.load(f)
    matrix = data.get(key, [])
    if not isinstance(matrix, list) or not matrix:
        return []
    if 0 <= row_index < len(matrix):
        return matrix[row_index]
    return []


def main(key_csv, base_dir, output_csv, row_index):
    true_index = row_index - 1

    key_mapping = load_key_csv(key_csv)

    with open(output_csv, "w", newline='') as out_file:
        writer = csv.writer(out_file)

        # ── Header ────────────────────────────────────────────────────────────
        header = (
            ["pool_name", "chain_name"] +
            [f"iptm_sample-{i}" for i in range(5)] +
            [f"pae_sample-{i}" for i in range(5)] +
            ["global_iptm_best", "ranking_score_best",
             "mean_iptm", "std_iptm",
             "mean_ranking_score", "std_ranking_score"]
        )
        writer.writerow(header)

        for top_folder, chains in key_mapping.items():
            top_folder_path = os.path.join(base_dir, top_folder)
            if not os.path.isdir(top_folder_path):
                print(f"Warning: folder '{top_folder}' not found in {base_dir}")
                continue

            # ── Chain-pair matrices across 5 samples ──────────────────────────
            iptm_scores = []
            pae_scores  = []

            for i in range(5):
                pattern = os.path.join(top_folder_path, f"*sample-{i}",
                                       "summary_confidences.json")
                files = glob.glob(pattern)

                if not files:
                    print(f"Warning: no summary_confidences.json for sample-{i} in {top_folder}")
                    iptm_scores.append([])
                    pae_scores.append([])
                    continue

                json_path = files[0]
                iptm_scores.append(extract_row(json_path, "chain_pair_iptm",    true_index))
                pae_scores.append( extract_row(json_path, "chain_pair_pae_min", true_index))

            # ── Global scores from top-level summary_confidences.json ─────────
            top_json = os.path.join(top_folder_path,
                                    f"{top_folder}_summary_confidences.json")
            try:
                with open(top_json) as f:
                    top = json.load(f)
                global_iptm_best     = top.get('iptm', 'NA')
                ranking_score_best   = top.get('ranking_score', 'NA')
            except Exception:
                global_iptm_best, ranking_score_best = 'NA', 'NA'

            # ── Mean and std across all 5 samples ─────────────────────────────
            iptm_vals          = []
            ranking_score_vals = []

            for i in range(5):
                pattern = os.path.join(top_folder_path, f"*sample-{i}",
                                       "summary_confidences.json")
                files = glob.glob(pattern)
                if not files:
                    continue
                try:
                    with open(files[0]) as f:
                        data = json.load(f)
                    iptm = data.get('iptm', None)
                    rs   = data.get('ranking_score', None)
                    if iptm is not None:
                        iptm_vals.append(iptm)
                    if rs is not None:
                        ranking_score_vals.append(rs)
                except Exception:
                    continue

            mean_iptm          = round(statistics.mean(iptm_vals), 4)           if iptm_vals           else 'NA'
            std_iptm           = round(statistics.stdev(iptm_vals), 4)          if len(iptm_vals) > 1  else 'NA'
            mean_ranking_score = round(statistics.mean(ranking_score_vals), 4)  if ranking_score_vals  else 'NA'
            std_ranking_score  = round(statistics.stdev(ranking_score_vals), 4) if len(ranking_score_vals) > 1 else 'NA'

            # ── Write one row per chain ───────────────────────────────────────
            for idx, chain_name in enumerate(chains):
                row = [top_folder, chain_name]

                # chain-pair iptm across 5 samples
                for scores in iptm_scores:
                    row.append(scores[idx] if idx < len(scores) else "")

                # chain-pair pae across 5 samples
                for scores in pae_scores:
                    row.append(scores[idx] if idx < len(scores) else "")

                # global scores
                row += [global_iptm_best, ranking_score_best,
                        mean_iptm, std_iptm,
                        mean_ranking_score, std_ranking_score]

                writer.writerow(row)

    print(f"Output written to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-chain ipTM/PAE plus global scores from pooled AF3 outputs."
    )
    parser.add_argument("key_csv",    help="Path to scraper key CSV")
    parser.add_argument("base_dir",   help="Path to base directory containing output folders")
    parser.add_argument("output_csv", help="Path for output CSV")
    parser.add_argument("--row-index", type=int, default=1,
                        help="1-based row index for bait chain (default: 1 = chain A)")
    args = parser.parse_args()

    main(args.key_csv, args.base_dir, args.output_csv, args.row_index)
