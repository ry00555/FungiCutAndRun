#!/bin/bash
#SBATCH --job-name=MergeMacsPeaks
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30gb
#SBATCH --time=6:00:00
#SBATCH --output=../MergeMacsPeaks.%j.out
#SBATCH --error=../MergeMacsPeaks.%j.err

cd $SLURM_SUBMIT_DIR
\set -euo pipefail
set -x  # Trace commands for debugging

# ------------------------------
# User-configurable variables
# ------------------------------
META="/scratch/ry00555/RNASeqPaper/SortedBamFiles/SortedBamFiles_meta_132to149.txt"
OUTDIR="/scratch/ry00555/RNASeqPaper/MACSPeaks"
OUTLIST="${OUTDIR}/consensus_peak_list.txt"
MERGE_GAP=100  # distance in bp to merge nearby peaks

mkdir -p "$OUTDIR"
> "$OUTLIST"  # Clear previous outlist

declare -A id_to_peaks

# ------------------------------
# Build array of peak files per ID
# ------------------------------
tail -n +2 "$META" | while IFS=$'\t' read -r ChIPBam BamIndex Strain Antibody Rep ID Input InputIndex MACS; do
    peakfile="$OUTDIR/$MACS"

    if [[ -f "$peakfile" ]]; then
        id_to_peaks["$ID"]+="$peakfile "
    else
        echo "⚠️ Warning: peak file not found: $peakfile"
    fi
done

ml BEDTools
# ------------------------------
# Generate consensus peaks
# ------------------------------
for ID in "${!id_to_peaks[@]}"; do
    files=(${id_to_peaks[$ID]})
    if [ ${#files[@]} -eq 0 ]; then
        echo "⚠️ No peak files found for $ID, skipping..."
        continue
    fi

    echo "➡️ Creating consensus peaks for $ID: ${files[*]}"

    merged_files=()
    for f in "${files[@]}"; do
        sorted_file="${f}.sorted"
        merged_file="${f}.merged"

        # Sort
        bedtools sort -i "$f" > "$sorted_file"

        # Merge nearby peaks
        bedtools merge -i "$sorted_file" -d $MERGE_GAP > "$merged_file"

        merged_files+=("$merged_file")
    done

    # Intersect merged replicates to get consensus
    consensus="${OUTDIR}/${ID}_consensus_peaks.broadPeak"
    if [ ${#merged_files[@]} -eq 1 ]; then
        cp "${merged_files[0]}" "$consensus"
    else
        bedtools intersect -a "${merged_files[0]}" -b "${merged_files[@]:1}" > "$consensus"
    fi

    echo "   ✅ Consensus peaks written to $consensus"

    # Append to outlist
    echo "$consensus" >> "$OUTLIST"
done
