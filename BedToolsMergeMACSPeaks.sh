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

META="/scratch/ry00555/RNASeqPaper/SortedBamFiles/SortedBamFiles_meta_132to149.txt"
OUTDIR="/scratch/ry00555/RNASeqPaper/MACSPeaks"

ml BEDTools

# Gap threshold for merging peaks (adjust as needed, e.g., 100 bp)
MERGE_GAP=100

# Declare associative array to store peaks per ID
declare -A id_to_peaks

# ---------------------------
# Step 1: Collect peak files per ID
# ---------------------------
tail -n +2 "$META" | while IFS=$'\t' read -r ChIPBam BamIndex Strain Antibody Rep ID Input InputIndex MACS; do
    # Full path to the MACS broadPeak file
    peakfile="${OUTDIR}/${MACS}"  # Use MACS column for exact filename

    if [[ -f "$peakfile" ]]; then
        # Append to the ID's list
        if [[ -z "${id_to_peaks[$ID]}" ]]; then
            id_to_peaks["$ID"]="$peakfile"
        else
            id_to_peaks["$ID"]+=" $peakfile"
        fi
    else
        echo "⚠️ Warning: peak file not found: $peakfile"
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
