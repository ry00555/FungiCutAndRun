#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90gb
#SBATCH --time=6:00:00
#SBATCH --output=../MacsPeakCalling.%j.out
#SBATCH --error=../MacsPeakCalling.%j.err

cd $SLURM_SUBMIT_DIR
# Paths

BAMDIR="/scratch/ry00555/RNASeqPaper/SortedBamFiles"
META="${BAMDIR}/SortedBamFiles_meta_132to149.txt"
OUTDIR="/scratch/ry00555/RNASeqPaper/MACSPeaks"

# Empty outlist
> "$OUTLIST"

# Load modules
ml MACS3
ml BEDTools

# Gap threshold for merging peaks (adjust as needed, e.g., 100 bp)
MERGE_GAP=100

# Declare associative array to store peaks per ID
declare -A id_to_peaks

# ---------------------------
# Step 1: Collect peak files per ID
# ---------------------------
tail -n +2 "$META" | while IFS=$'\t' read -r ChIPBam BamIndex Strain Antibody Rep ID Input InputIndex MACS; do
    peakfile="${OUTDIR}/${ID}_peaks.broadPeak"  # Adjust if MACS3 naming is different
    if [[ -f "$peakfile" ]]; then
        id_to_peaks["$ID"]+="${peakfile} "
    fi
done

# ---------------------------
# Step 2: Create consensus peak sets
# ---------------------------
for ID in "${!id_to_peaks[@]}"; do
    files=(${id_to_peaks[$ID]})
    echo "➡️ Creating consensus peaks for $ID: ${files[*]}"

    merged_files=()

    # Sort & merge each replicate to account for shifted peaks
    for f in "${files[@]}"; do
        sorted_file="${f}.sorted"
        merged_file="${f}.merged"

        # Sort
        bedtools sort -i "$f" > "$sorted_file"

        # Merge nearby peaks within MERGE_GAP
        bedtools merge -i "$sorted_file" -d $MERGE_GAP > "$merged_file"

        merged_files+=("$merged_file")
    done

    # Intersect merged replicates to get consensus
    consensus="${OUTDIR}/${ID}_consensus_peaks.broadPeak"
    bedtools intersect -a "${merged_files[0]}" -b "${merged_files[@]:1}" > "$consensus"

    echo "   ✅ Consensus peaks written to $consensus"

    # Append to outlist
    echo "$consensus" >> "$OUTLIST"
done

echo "✅ All consensus peaks generated. List in $OUTLIST"
