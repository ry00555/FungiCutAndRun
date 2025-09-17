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
# Paths

#!/bin/bash

# ------------------------------
# CONFIG
# ------------------------------
META="/scratch/ry00555/RNASeqPaper/SortedBamFiles/SortedBamFiles_meta_132to149.txt"
OUTDIR="/scratch/ry00555/RNASeqPaper/MACSPeaks"
OUTLIST="${OUTDIR}/consensus_peak_list.txt"
MERGE_GAP=100  # merge nearby peaks within 100 bp

> "$OUTLIST"  # clear previous outlist

declare -A id_to_peaks

# ------------------------------
# BUILD ID -> peak file mapping
# ------------------------------
tail -n +2 "$META" | while IFS=$'\t' read -r ChIPBam BamIndex Strain Antibody Rep ID Input InputIndex MACS; do
    # Only consider broadPeak files in MACS column
    if [[ -n "$MACS" && "$MACS" == *.broadPeak ]]; then
        peakfile="${OUTDIR}/${MACS}"
        if [[ -f "$peakfile" ]]; then
            # Append to array
            id_to_peaks["$ID"]+="$peakfile "
        else
            echo "⚠️ Warning: peak file not found: $peakfile"
        fi
    fi
done

ml BEDTools
# ------------------------------
# CREATE CONSENSUS PEAKS
# ------------------------------
for ID in "${!id_to_peaks[@]}"; do
    files=(${id_to_peaks[$ID]})
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

    consensus="${OUTDIR}/${ID}_consensus_peaks.broadPeak"
    if [ ${#merged_files[@]} -gt 1 ]; then
        bedtools intersect -a "${merged_files[0]}" -b "${merged_files[@]:1}" > "$consensus"
    else
        cp "${merged_files[0]}" "$consensus"
    fi

    echo "   ✅ Consensus peaks written to $consensus"
    echo "$consensus" >> "$OUTLIST"

done
