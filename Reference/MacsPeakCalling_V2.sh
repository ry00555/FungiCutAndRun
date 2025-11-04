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
META="${BAMDIR}/SortedBamFiles_meta_132to149.txt"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"

for bam in $(find . -name "*.bam"); do
    bai="${bam}.bai"
    if [ ! -f "$bai" ]; then
        echo "Missing index for: $bam"
    fi
done

# File to collect all output peaks
OUTLIST="${OUTDIR}/MACS_peak_files.txt"
> "$OUTLIST"
ml MACS3

# Loop through metadata (skip header)
tail -n +2 "$META" | while IFS=$'\t' read -r RunID	bamReads	BamIndex	SampleID	Factor	Tissue	Condition	Replicate	bamControl	bamInputIndex	ControlID	Peaks	PeakCaller; do

    # Build full paths
    chip_path="${BAMDIR}/${bamReads}"
    input_path=""
    [[ -n "$Input" ]] && input_path="${BAMDIR}/${bamControl}"

    base=$(basename "$bamReads" .bam)
    prefix="${OUTDIR}/${base}"
    echo "➡️ Processing: $base"
    # Expected output files for a complete run
       expected=("${prefix}_peaks.broadPeak" "${prefix}_peaks.xls")

       # Check if outputs already exist and are complete
       all_exist=true
       for f in "${expected[@]}"; do
           if [[ ! -s "$f" ]]; then
               all_exist=false
               break
           fi
       done

       if $all_exist; then
           echo "   ✅ Skipping (MACS3 output already complete)"
           echo "${prefix}_peaks.broadPeak" >> "$OUTLIST"
           continue
       else
           # Cleanup partials
           echo "   ⚠️ Detected missing/partial output → cleaning up"
           for f in "${expected[@]}"; do
               [[ -f "$f" ]] && rm -f "$f"
           done
       fi

    # ✅ Skip if already processed
    if [[ -f "$peakfile" ]]; then
        echo "   🔄 Peaks already exist: $peakfile → skipping"
        continue
    fi

    # Check BAM and index
    if [[ ! -f "$chip_path" ]]; then
        echo "⚠️ Missing ChIP BAM: $chip_path"
        continue
    fi
    if [[ ! -f "${BAMDIR}/${BamIndex}" ]]; then
        echo "⚠️ Missing BAM index for: $chip_path"
        continue
    fi

        # Run with input if available
        if [[ -n "$Input" && -f "$input_path" ]]; then
            echo "   Using input: $input_path"
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
            echo "   No input found → calling without control"
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

        # Collect output file
        echo "${OUTDIR}/${base}_peaks.broadPeak" >> "$OUTLIST"
        echo "${prefix}_peaks.broadPeak" >> "$OUTLIST"

    done

    echo "✅ Peak calling complete. Outputs listed in $OUTLIST"
