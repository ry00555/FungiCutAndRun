#!/bin/bash
#SBATCH --job-name=Macs3_Rerun
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=../Macs3_Rerun.%j.out
#SBATCH --error=../Macs3_Rerun.%j.err

# --- Paths ---
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"

ml MACS3

# Convert potential Windows line endings
dos2unix "$META" 2>/dev/null || true

# --- Targeted rerun list ---
declare -a MISSING_PEAKS=(
"134-12_cdp-6_H3K36me3_rep1_peaks.broadPeak"
"134-15_cdp-6_H3K36me3_rep2_peaks.broadPeak"
"147-35_H3K56R_H3K27me3_rep7_peaks.broadPeak"
"147-36_H3K56R_H3K36me3_rep7_peaks.broadPeak"
"147-33_H3K56R_Input_rep7_peaks.broadPeak"
"148-130_nst4_H3K4ac_rep1_peaks.broadPeak"
"132-29_rco-1_Input_rep1_peaks.broadPeak"
"131-55_rtt109_Input_rep3_peaks.broadPeak"
"149-32_set-7_H3K36me3_rep8_peaks.broadPeak"
"149-19_set-7_H3K9me3_rep1_peaks.broadPeak"
)

echo "🚀 Recalling MACS3 for ${#MISSING_PEAKS[@]} missing peak files..."

# Loop through META, matching only entries in the missing list
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$RunID" ]] && continue
    peakfile=$(basename "$Peaks")

    # Check if this row is one of the missing ones
    if printf '%s\n' "${MISSING_PEAKS[@]}" | grep -qx "$peakfile"; then
        echo "➡️  Recalling peaks for: $DesiredPeakName"

        chip_path="${BAMDIR}/${bamReads}"
        input_path="${BAMDIR}/${bamControl:-}"
        prefix="${OUTDIR}/${DesiredPeakName}"

        # Check BAM
        if [[ ! -f "$chip_path" ]]; then
            echo "⚠️ Missing BAM file: $chip_path"
            continue
        fi

        # Clean up old output (if any)
        rm -f "${prefix}_peaks."* 2>/dev/null || true

        # Run MACS3 (with or without control)
        if [[ -n "$bamControl" && -f "$input_path" ]]; then
            echo "   Using control: $input_path"
            macs3 callpeak \
                -t "$chip_path" \
                -c "$input_path" \
                -f BAMPE \
                -n "$DesiredPeakName" \
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
                -n "$DesiredPeakName" \
                --broad \
                --broad-cutoff 0.1 \
                -g 41037538 \
                --outdir "$OUTDIR" \
                --min-length 800 \
                --max-gap 500
        fi

        echo "✅ Finished: $DesiredPeakName"
    fi
done

echo "🎯 All targeted MACS3 re-runs complete."
