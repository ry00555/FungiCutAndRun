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
META_CSV="${BASEDIR}/EAF3_MetaV4_withBigWigs.csv"
AVG_DIR="${BASEDIR}/AveragedBW"
NORM_DIR="${BASEDIR}/Log2InputNormBW"
RATIO_DIR="${BASEDIR}/InputNormRatioBW"

mkdir -p "$AVG_DIR" "$NORM_DIR" "$RATIO_DIR"

echo "== Checking deepTools is available =="
which bigwigAverage
which bigwigCompare

echo "== Reading ${META_CSV} and running merge + input-normalize directly =="

python3 - "$BASEDIR" "$META_CSV" "$AVG_DIR" "$NORM_DIR" "$RATIO_DIR" "$THREADS" <<'PYEOF'
import sys
import subprocess
import pandas as pd

basedir, meta_csv, avg_dir, norm_dir, ratio_dir, threads = sys.argv[1:7]

df = pd.read_csv(meta_csv)

if "bigwig_found" in df.columns:
   df = df[df["bigwig_found"] == True].copy()

if "bigwig_filename" in df.columns:
   df["bigwig_path"] = df["bigwig_filename"].apply(lambda f: f"{basedir}/{f}")
elif "bigwig_path" in df.columns:
   pass  # already has usable paths
else:
   raise SystemExit("Meta CSV needs either a 'bigwig_filename' or 'bigwig_path' column")

# Safety net: collapse any literal leftover duplicate rows (same RunID +
# Condition + Factor) rather than erroring or double-counting them. This is
# only safe for exact duplicates -- a genuine label CONFLICT (same RunID,
# different Condition/Factor) is a different, more serious problem and is
# checked for separately below.
before = len(df)
df = df.drop_duplicates(subset=["RunID", "Condition", "Factor"], keep="first")
if before != len(df):
   print(f"Collapsed {before - len(df)} exact-duplicate rows (same RunID+Condition+Factor).", flush=True)

conflict_runids = []
for runid, g in df.groupby("RunID"):
   if g["Condition"].nunique() > 1 or g["Factor"].nunique() > 1:
       conflict_runids.append(runid)
       print(f"CONFLICT: RunID {runid} has inconsistent Condition/Factor:", flush=True)
       print(g[["Condition", "Factor", "SampleID"]].to_string(index=False), flush=True)

if conflict_runids:
   raise SystemExit(
       "\nOne or more RunIDs have conflicting Condition/Factor labels (see CONFLICT "
       "block(s) above) -- fix the source CSV before rerunning. Refusing to guess."
   )


def get_dominant_input_condition(chip_condition, chip_factor, all_df):
   """Resolve a ChIP condition's Input group via bamControl -- handles
   per-timepoint-own-Input and shared-Input-across-timepoints patterns
   without assuming a fixed naming rule. Filters by BOTH Condition and
   Factor, since the same Condition string can have rows for more than one
   ChIP factor (e.g. tetRGFP_0hr has both H3K27me3 and H3K36me3 replicates) --
   without this, their bamControl values would get blended together."""
   chip_rows = all_df[(all_df["Condition"] == chip_condition) & (all_df["Factor"] == chip_factor)]
   if chip_rows.empty or chip_rows["bamControl"].isna().all():
       return None

   input_conditions = []
   for bc in chip_rows["bamControl"].dropna():
       match = all_df[(all_df["bamReads"] == bc) & (all_df["Factor"] == "Input")]
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


# Tracks are keyed by (Factor, Condition) -- built directly from the Factor
# column (H3K27me3 / H3K36me3 / etc.) rather than a separate 'mark' column.
# Merged Input tracks are cached by resolved input_condition ALONE (not by
# factor), since the physical merged Input bigwig is identical regardless of
# which ChIP factor references it -- avoids redundantly rebuilding the same
# Input merge once per factor that happens to share it.
input_merge_cache = {}  # input_condition -> merged Input path
n_done = 0
n_skipped = 0
n_failed = 0

chip_factors = df.loc[df["Factor"] != "Input", "Factor"].unique()

for factor in chip_factors:
   conditions = df.loc[df["Factor"] == factor, "Condition"].unique()

   for condition in conditions:
       chip_rows = df[(df["Condition"] == condition) & (df["Factor"] == factor)]
       if chip_rows.empty:
           continue

       input_condition = get_dominant_input_condition(condition, factor, df)
       if input_condition is None:
           print(f"SKIP [{factor}] {condition}: could not resolve Input via bamControl", flush=True)
           n_skipped += 1
           continue

       input_rows = df[(df["Condition"] == input_condition) & (df["Factor"] == "Input")]
       if input_rows.empty:
           print(f"SKIP [{factor}] {condition}: resolved Input condition '{input_condition}' has no bigwigs", flush=True)
           n_skipped += 1
           continue

       print(f"\n[{factor}] {condition}  (Input: {input_condition})", flush=True)

       chip_bws = chip_rows["bigwig_path"].tolist()
       chip_tag = f"{factor}_{condition}_ChIP".replace(" ", "_")
       chip_merged = f"{avg_dir}/{chip_tag}.bw"
       chip_bw_args = " ".join(f'"{p}"' for p in chip_bws)

       print(f"  [merge ChIP] {len(chip_bws)} replicate(s) -> {chip_tag}.bw", flush=True)
       ok = run(f'bigwigAverage -b {chip_bw_args} -o "{chip_merged}" -p {threads}')
       if not ok:
           n_failed += 1
           continue

       if input_condition in input_merge_cache:
           input_merged = input_merge_cache[input_condition]
           print(f"  [reusing merged Input] {input_condition}", flush=True)
       else:
           input_bws = input_rows["bigwig_path"].tolist()
           input_tag = f"{input_condition}_Input".replace(" ", "_")
           input_merged = f"{avg_dir}/{input_tag}.bw"
           input_bw_args = " ".join(f'"{p}"' for p in input_bws)

           print(f"  [merge Input] {len(input_bws)} replicate(s) -> {input_tag}.bw", flush=True)
           ok = run(f'bigwigAverage -b {input_bw_args} -o "{input_merged}" -p {threads}')
           if not ok:
               n_failed += 1
               continue
           input_merge_cache[input_condition] = input_merged

       final_tag = f"{factor}_{condition}_InputNorm".replace(" ", "_")
       final_path = f"{norm_dir}/{final_tag}.bw"

       print(f"  [log2 compare] -> {final_tag}.bw", flush=True)
       ok = run(f'bigwigCompare -b1 "{chip_merged}" -b2 "{input_merged}" --operation log2 -o "{final_path}" -p {threads}')
       if not ok:
           n_failed += 1
           continue

       ratio_tag = f"{factor}_{condition}_InputNormRatio".replace(" ", "_")
       ratio_path = f"{ratio_dir}/{ratio_tag}.bw"

       print(f"  [ratio compare] -> {ratio_tag}.bw", flush=True)
       ok = run(f'bigwigCompare -b1 "{chip_merged}" -b2 "{input_merged}" --operation ratio -o "{ratio_path}" -p {threads}')
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
echo "Ratio input-normalized bigwigs:           ${RATIO_DIR}"
