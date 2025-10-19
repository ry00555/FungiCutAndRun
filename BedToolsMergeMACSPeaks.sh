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

#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Path configuration
# -----------------------------
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"
MERGE_GAP=100

SUMMARY="${OUTDIR}/replicate_vs_consensus.tsv"
OUTLIST="${OUTDIR}/consensus_peak_list.txt"
FRIP_SUMMARY="${OUTDIR}/FRiP_summary.tsv"

# -----------------------------
# Setup
# -----------------------------
mkdir -p "$OUTDIR"

# Initialize output files
echo -e "Tissue\tReplicate\tRepPeakFile\tNumPeaks\tNumOverlap\tFracOverlap\tConsensusFile\tNumConsensusPeaks" > "$SUMMARY"
echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP" > "$FRIP_SUMMARY"
: > "$OUTLIST"

ml BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

# Normalize CSV line endings (just in case)
dos2unix "$META" 2>/dev/null || true

# -----------------------------
# Step 1: Calculate FRiP per replicate
# -----------------------------
echo "🧮 Calculating FRiP for each replicate..."

while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$SampleID" || -z "$Tissue" ]] && continue

    bam="${BAMDIR}/${bamReads}"
    peaks="${MACSDIR}/${Peaks}"

    if [[ ! -s "$bam" ]]; then
        echo "⚠️ Missing BAM file: $bam"
        continue
    fi
    if [[ ! -s "$peaks" ]]; then
        echo "⚠️ Missing peak file: $peaks"
        continue
    fi

    total_reads=$(samtools view -c -F 260 "$bam")
    reads_in_peaks=$(bedtools intersect -u -a "$bam" -b "$peaks" | wc -l)
    frip=$(awk "BEGIN{printf \"%.4f\", ($total_reads>0)?$reads_in_peaks/$total_reads:0}")

    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}" >> "$FRIP_SUMMARY"
done < <(tail -n +2 "$META")

echo "✅ FRiP summary saved to: $FRIP_SUMMARY"

# -----------------------------
# Step 2: Build associative array of peak files per Tissue
# -----------------------------
declare -A tissue_to_peaks

while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$Tissue" ]] && continue
    peak="${MACSDIR}/${Peaks}"

    if [[ -s "$peak" ]]; then
        tissue_to_peaks["$Tissue"]+="$peak "
    else
        echo "⚠️ Skipping missing peak file: $peak"
    fi
done < <(tail -n +2 "$META")

# -----------------------------
# Step 3: Optional consensus peak listing
# -----------------------------
echo "🧩 Collecting consensus peak sets per tissue..."
for tissue in "${!tissue_to_peaks[@]}"; do
    echo "$tissue ${tissue_to_peaks[$tissue]}" >> "$OUTLIST"
done

echo "✅ Consensus peak file list saved to: $OUTLIST"
echo "✅ Summary outputs:"
echo "   - FRiP summary: $FRIP_SUMMARY"
echo "   - Peak list: $OUTLIST"
echo "   - Replicate summary: $SUMMARY"

# ---------------------------------------------
# Step4: Loop through tissues to build consensus + overlaps
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
    # Step 5: Compute overlap for each replicate
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
