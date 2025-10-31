#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=../MacsPeakCalling.%j.out
#SBATCH --error=../MacsPeakCalling.%j.err

cd $SLURM_SUBMIT_DIR
#Set Paths

BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"


OUTLIST="${OUTDIR}/MACS_peak_files.txt"
> "$OUTLIST"
ml MACS3

 #Remove potential carriage returns (Mac Excel export issue)
 dos2unix "$META" 2>/dev/null || true
 #### run this in command line before starting to make sure bam files have reads and they were not truncated / are not missing EOF

 for bam in "$BAMDIR"/*.bam; do
     base=$(basename "$bam")
     mapped=$(samtools view -c -F 4 "$bam")
     echo -e "$base\t$mapped"
 done

 echo -e "BAM_File\tMapped_Reads\tEOF_OK"

 for bam in "$BAMDIR"/*.bam; do
     base=$(basename "$bam")

     # Check EOF
     EOF_OK="OK"
     if ! samtools quickcheck -v "$bam" &>/dev/null; then
         EOF_OK="MISSING_EOF"
     fi

     # Count mapped reads
     if [[ "$EOF_OK" == "OK" ]]; then
         mapped=$(samtools view -c -F 4 "$bam")
     else
         mapped="NA"
     fi

     echo -e "${base}\t${mapped}\t${EOF_OK}"
 done

 #--- STEP 1: Rename existing peak files ---
 echo "🔄 Checking for existing peak files to rename..."
 tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
     [[ -z "$RunID" ]] && continue

     old_base=$(basename "$bamReads" .bam)
     new_base="$DesiredPeakName"

     for ext in broadPeak gappedPeak xls; do
         old_file="${OUTDIR}/${old_base}_peaks.${ext}"
         new_file="${OUTDIR}/${new_base}_peaks.${ext}"
         if [[ -f "$old_file" && ! -f "$new_file" ]]; then
             echo "Renaming: $old_file → $new_file"
             mv "$old_file" "$new_file"
         fi
     done
 done

 echo "✅ Renaming step complete."

 #--- STEP 2: Run MACS3 where needed ---
 echo "🚀 Starting MACS3 peak calling..."
 tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
     [[ -z "$RunID" ]] && continue

     chip_path="${BAMDIR}/${bamReads}"
     input_path="${BAMDIR}/${bamControl:-}"
     index_path="${BAMDIR}/${BamIndex:-}"
     prefix="${OUTDIR}/${DesiredPeakName}"
     peakfile="${prefix}_peaks.broadPeak"

     echo "➡️ Processing: $DesiredPeakName"

     # Check for BAM + BAI
     if [[ ! -f "$chip_path" ]]; then
         echo "⚠️ Missing BAM: $chip_path"
         continue
     fi
     if [[ ! -f "$index_path" ]]; then
         echo "⚠️ Missing BAM index: $index_path"
         continue
     fi

     # Define expected outputs
     expected=(
         "${prefix}_peaks.broadPeak"
         "${prefix}_peaks.xls"
         "${prefix}_peaks.gappedPeak"
     )

     # Skip if all output files exist and broadPeak is >=1200 bytes
     all_exist=true
     for f in "${expected[@]}"; do
         if [[ ! -s "$f" ]]; then
             all_exist=false
             break
         fi
     done

     # Check broadPeak size
     if [[ -f "${prefix}_peaks.broadPeak" ]]; then
         bsize=$(stat -c%s "${prefix}_peaks.broadPeak" 2>/dev/null || stat -f%z "${prefix}_peaks.broadPeak")
         if (( bsize < 1200 )); then
             echo "⚠️ BroadPeak file too small (${bsize} bytes) → will rerun MACS3"
             all_exist=false
         fi
     fi

     if $all_exist; then
         echo "   ✅ Skipping (MACS3 output already complete)"
         echo "$peakfile" >> "$OUTLIST"
         continue
     fi

     echo "   ⚠️ Running MACS3 for: $DesiredPeakName"

     # Cleanup any partial output
     for f in "${expected[@]}"; do
         [[ -f "$f" ]] && rm -f "$f"
     done

     # --- Run MACS3 ---
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

     # Record output
     echo "$peakfile" >> "$OUTLIST"
 done

 echo "✅ STEP 2: MACS3 peak calling complete. Outputs listed in $OUTLIST"
