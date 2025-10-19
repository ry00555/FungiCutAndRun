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
set -euo pipefail
set -x  # debug trace

# === Paths ===
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"

MERGE_GAP=100  # bp gap allowed between merged peaks
SUMMARY="${OUTDIR}/replicate_vs_consensus.tsv"
OUTLIST="${OUTDIR}/consensus_peak_list.txt"
FRiP_summary="/scratch/ry00555/RNASeqPaper/Oct2025/FRiP_summary.tsv"

> "$SUMMARY"
> "$OUTLIST"
> "FRiP_summary"

ml BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

echo -e "Tissue\tReplicate\tRepPeakFile\tNumPeaks\tNumOverlap\tFracOverlap\tConsensusFile\tNumConsensusPeaks" > "$SUMMARY"
echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP" > "$FRiP_summary"

dos2unix "$META" 2>/dev/null || true
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    bam="${BAMDIR}/${SampleID}.bam"
    peaks="${MACSDIR}/${Peaks}"

    [[ ! -s "$bam" || ! -s "$peaks" ]] && continue

    total=$(samtools view -c -F 260 "$bam")
    inpeaks=$(bedtools intersect -u -a "$bam" -b "$peaks" | wc -l)
    frip=$(awk "BEGIN{printf \"%.3f\", ($total>0)?$inpeaks/$total:0}")

    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total}\t${inpeaks}\t${frip}" >> "$FRiP_summary"
done

echo "✅ FRiP summary saved to $FRiP_summary"

# ---------------------------------------------
# Step 1: Collect replicate peak files per Tissue
# ---------------------------------------------
while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$Tissue" ]] && continue
    peak="${MACSDIR}/${Peaks}"
    if [[ -s "$peak" ]]; then
        tissue_to_peaks["$Tissue"]+="$peak "
    else
        echo "⚠️ Skipping missing peak file: $peak"
    fi
done < <(awk -F',' 'NR>1' "$META")

# ---------------------------------------------
# Step 2: Loop through tissues to build consensus + overlaps
# ---------------------------------------------
for tissue in "${!tissue_to_peaks[@]}"; do
    files=(${tissue_to_peaks[$tissue]})
    n=${#files[@]}
    if (( n < 1 )); then
        echo "⚠️ No valid peak files for $tissue, skipping..."
        continue
    fi

    echo "🧬 Building consensus for $tissue ($n replicates)"
    merged_files=()
    for f in "${files[@]}"; do
        sorted="${MACSDIR}/$(basename "$f").sorted"
        merged="${MACSDIR}/$(basename "$f").merged"
        bedtools sort -i "$f" > "$sorted"
        bedtools merge -i "$sorted" -d $MERGE_GAP > "$merged"
        merged_files+=("$merged")
    done

    consensus="${CHIPR_DIR}/${tissue}_consensus_peaks.broadPeak"
    if (( n == 1 )); then
        cp "${merged_files[0]}" "$consensus"
    else
        bedtools intersect -a "${merged_files[0]}" -b "${merged_files[@]:1}" > "$consensus"
    fi

    num_consensus=$(wc -l < "$consensus")
    echo "✅ $tissue: consensus has $num_consensus peaks"
    echo "$consensus" >> "$OUTLIST"

    # ---------------------------------------------
    # Step 3: Compute overlap for each replicate
    # ---------------------------------------------
    for f in "${files[@]}"; do
        rep_name=$(basename "$f")
        overlap="${OUTDIR}/${rep_name%.bed}_vs_consensus.bed"
        bedtools intersect -u -a "$f" -b "$consensus" > "$overlap"

        num_peaks=$(wc -l < "$f")
        num_overlap=$(wc -l < "$overlap")
        frac=$(awk "BEGIN{printf \"%.3f\", ($num_peaks>0)?$num_overlap/$num_peaks:0}")

        echo -e "${tissue}\t${rep_name}\t${f}\t${num_peaks}\t${num_overlap}\t${frac}\t${consensus}\t${num_consensus}" >> "$SUMMARY"
    done
done

echo "🎉 Done!"
echo "📄 Summary written to: $SUMMARY"
echo "📁 Consensus peaks in: $OUTDIR"
