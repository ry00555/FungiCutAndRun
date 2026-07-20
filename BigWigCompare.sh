#!/bin/bash
#SBATCH --job-name=BWCompare
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=../BWDeep.%j.out
#SBATCH --error=../BWDeep.%j.err


set -euo pipefail

# ---- Environment ----
module load deepTools

THREADS=12


# ---- Paths ----
BASEDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/ChIPSeq/RemappedBW"
META_CSV="${BASEDIR}/BW_Meta_deduplicated.csv"
AVG_DIR="${BASEDIR}/AveragedBW"
NORM_DIR="${BASEDIR}/Log2InputNormBW"

mkdir -p "$AVG_DIR" "$NORM_DIR"

echo "== Checking deepTools is available =="
which bigwigAverage
which bigwigCompare

echo "== Reading ${META_CSV} and running merge + input-normalize directly =="

python3 - "$BASEDIR" "$META_CSV" "$AVG_DIR" "$NORM_DIR" "$THREADS" <<'PYEOF'
import sys
import subprocess
import pandas as pd

basedir, meta_csv, avg_dir, norm_dir, threads = sys.argv[1:6]

df = pd.read_csv(meta_csv)

if "bigwig_found" in df.columns:
    df = df[df["bigwig_found"] == True].copy()

if "bigwig_filename" in df.columns:
    df["bigwig_path"] = df["bigwig_filename"].apply(lambda f: f"{basedir}/{f}")
elif "bigwig_path" in df.columns:
    pass  # already has usable paths
else:
    raise SystemExit("Meta CSV needs either a 'bigwig_filename' or 'bigwig_path' column")

# Sanity check: this script assumes BW_Meta_deduplicated.csv, which already has
# one row per physical replicate+mark+condition+factor (no repeats from a
# sample being listed in multiple source metadata sheets, no conflicting
# labels). If that's not true of whatever CSV got pointed at here, stop rather
# than silently double-counting or merging on a contradiction.
dup_check = df.groupby(["RunID", "mark", "Condition", "Factor"]).size()
if (dup_check > 1).any():
    raise SystemExit(
        "META_CSV is not deduplicated (found repeated RunID+mark+Condition+Factor rows) -- "
        "point this script at BW_Meta_deduplicated.csv, or dedupe first."
    )


def get_dominant_input_condition(chip_condition, ds_df):
    """Resolve a ChIP condition's Input group via bamControl -- handles
    per-timepoint-own-Input and shared-Input-across-timepoints patterns
    without assuming a fixed naming rule."""
    chip_rows = ds_df[(ds_df["Condition"] == chip_condition) & (ds_df["Factor"] != "Input")]
    if chip_rows.empty or chip_rows["bamControl"].isna().all():
        return None

    input_conditions = []
    for bc in chip_rows["bamControl"].dropna():
        match = ds_df[(ds_df["bamReads"] == bc) & (ds_df["Factor"] == "Input")]
        if len(match) == 1:
            input_conditions.append(match["Condition"].iloc[0])

    if not input_conditions:
        return None
    return pd.Series(input_conditions).value_counts().idxmax()


def run(cmd):
    print(f"  $ {cmd}", flush=True)
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"  !! command failed (exit {result.returncode}): {cmd}", flush=True)
    return result.returncode == 0


# Tracks are keyed by (mark, condition) only -- NOT by source dataset. Every
# physical replicate is already deduplicated in the meta sheet, so a strain
# shared across multiple original projects (e.g. tetRGFP, tetRGFPsuz12) still
# produces exactly one merged/normalized track, not one redundant copy per
# project it happened to be referenced from.
input_merge_cache = {}  # (mark, input_condition) -> merged Input path, so shared Inputs are only merged once
n_done = 0
n_skipped = 0
n_failed = 0

for mark, ds_df in df.groupby("mark"):
    conditions = ds_df.loc[ds_df["Factor"] != "Input", "Condition"].unique()

    for condition in conditions:
        chip_rows = ds_df[(ds_df["Condition"] == condition) & (ds_df["Factor"] != "Input")]
        if chip_rows.empty:
            continue

        input_condition = get_dominant_input_condition(condition, ds_df)
        if input_condition is None:
            print(f"SKIP [{mark}] {condition}: could not resolve Input via bamControl", flush=True)
            n_skipped += 1
            continue

        input_rows = ds_df[(ds_df["Condition"] == input_condition) & (ds_df["Factor"] == "Input")]
        if input_rows.empty:
            print(f"SKIP [{mark}] {condition}: resolved Input condition '{input_condition}' has no bigwigs", flush=True)
            n_skipped += 1
            continue

        print(f"\n[{mark}] {condition}  (Input: {input_condition})", flush=True)

        chip_bws = chip_rows["bigwig_path"].tolist()
        chip_tag = f"{mark}_{condition}_ChIP".replace(" ", "_")
        chip_merged = f"{avg_dir}/{chip_tag}.bw"
        chip_bw_args = " ".join(f'"{p}"' for p in chip_bws)

        print(f"  [merge ChIP] {len(chip_bws)} replicate(s) -> {chip_tag}.bw", flush=True)
        ok = run(f'bigwigAverage -b {chip_bw_args} -o "{chip_merged}" -p {threads}')
        if not ok:
            n_failed += 1
            continue

        # Reuse an already-merged Input bigwig for this mark/input_condition
        # rather than rebuilding it once per ChIP condition that shares it.
        input_key = (mark, input_condition)
        if input_key in input_merge_cache:
            input_merged = input_merge_cache[input_key]
            print(f"  [reusing merged Input] {input_key}", flush=True)
        else:
            input_bws = input_rows["bigwig_path"].tolist()
            input_tag = f"{mark}_{input_condition}_Input".replace(" ", "_")
            input_merged = f"{avg_dir}/{input_tag}.bw"
            input_bw_args = " ".join(f'"{p}"' for p in input_bws)

            print(f"  [merge Input] {len(input_bws)} replicate(s) -> {input_tag}.bw", flush=True)
            ok = run(f'bigwigAverage -b {input_bw_args} -o "{input_merged}" -p {threads}')
            if not ok:
                n_failed += 1
                continue
            input_merge_cache[input_key] = input_merged

        final_tag = f"{mark}_{condition}_InputNorm".replace(" ", "_")
        final_path = f"{norm_dir}/{final_tag}.bw"

        print(f"  [log2 compare] -> {final_tag}.bw", flush=True)
        ok = run(f'bigwigCompare -b1 "{chip_merged}" -b2 "{input_merged}" --operation log2 -o "{final_path}" -p {threads}')
        if not ok:
            n_failed += 1
            continue

        n_done += 1

print(f"\n== Summary ==")
print(f"Completed:        {n_done}")
print(f"Skipped:          {n_skipped} (no resolvable Input -- see SKIP lines above)")
print(f"Failed:           {n_failed} (deepTools command errored -- see !! lines above)")
print(f"Unique merged Input tracks built: {len(input_merge_cache)}")
PYEOF

echo "== Done =="
echo "Averaged (per-condition merged) bigwigs: ${AVG_DIR}"
echo "Log2 input-normalized bigwigs:            ${NORM_DIR}"
 
