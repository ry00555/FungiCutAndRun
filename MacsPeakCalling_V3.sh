#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90gb
#SBATCH --time=8:00:00
#SBATCH --output=../MacsPeakCalling.%j.out
#SBATCH --error=../MacsPeakCalling.%j.err

cd $SLURM_SUBMIT_DIR
# Paths

BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged.csv"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"


OUTLIST="${OUTDIR}/MACS_peak_files.txt"
> "$OUTLIST"
ml MACS3

# Remove potential carriage returns (Mac Excel export issue)
dos2unix "$META" 2>/dev/null || true

# Skip header (tail -n +2), read CSV line by line
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller; do

    # Skip empty lines
    [[ -z "$RunID" ]] && continue

    base=$(basename "$bamReads" .bam)
    chip_path="${BAMDIR}/${bamReads}"
    input_path="${BAMDIR}/${bamControl:-}"
    index_path="${BAMDIR}/${BamIndex:-}"
    prefix="${OUTDIR}/${base}"

    echo "➡️ Processing: $base"

    # Check BAM + BAI existence
    if [[ ! -f "$chip_path" ]]; then
        echo "⚠️ Missing BAM: $chip_path"
        continue
    fi
    if [[ ! -f "$index_path" ]]; then
        echo "⚠️ Missing BAM index: $index_path"
        continue
    fi

    # Clean up partials
    rm -f "${prefix}_peaks.broadPeak" "${prefix}_peaks.xls" 2>/dev/null || true

    # --- Run MACS3 ---
    if [[ -n "$bamControl" && -f "$input_path" ]]; then
        echo "   Using control: $input_path"
        macs3 callpeak \
            -t "$chip_path" \
            -c "$input_path" \
            -f BAMPE \
            -n "$base" \
            --broad \
            --broad-cutoff 0.1 \
            -g 41037538 \
            --outdir "$OUTDIR" \
            --min-length 800 \
            --max-gap 500
    else
        echo "   No control found → running without input"
        macs3 callpeak \
            -t "$chip_path" \
            -f BAMPE \
            -n "$base" \
            --broad \
            --broad-cutoff 0.1 \
            -g 41037538 \
            --outdir "$OUTDIR" \
            --min-length 800 \
            --max-gap 500
    fi

    # Add to output list
    echo "${prefix}_peaks.broadPeak" >> "$OUTLIST"

done

echo "✅ Peak calling complete. Outputs listed in $OUTLIST"
