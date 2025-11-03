#!/bin/bash
#SBATCH --job-name=BAMtoFASTQ
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --output=../BAMtoFASTQ.%j.out
#SBATCH --error=../BAMtoFASTQ.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

set -euo pipefail

# ================================
# Directories and config
# ================================
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
#REMAPPED_BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/Remapped/SortedBamFiles"
#FASTQDIR="/scratch/ry00555/RNASeqPaper/Oct2025/Remapped/FASTQ_from_BAM"
BIGWIGDIR="/scratch/ry00555/RNASeqPaper/Oct2025/BigWigs"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025"
META="${OUTDIR}/BAM_File_Metadata_with_index_merged_V2.csv"
#mkdir -p "$REMAPPED_BAMDIR" "$FASTQDIR" "$REMAPPED_BAMDIR/tempReps" "$BIGWIGDIR"
dos2unix "$META" 2>/dev/null || true

# ================================
# Load modules
# ================================
ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0
ml deepTools

# ================================
# Extract FASTQ from BAMs with mapped reads
# ================================
#tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do

#    bam="${BAMDIR}/${bamReads}"

    # Output FASTQs
  #  fq1="$FASTQDIR/${DesiredPeakName}_R1.fastq"
  #  fq2="$FASTQDIR/${DesiredPeakName}_R2.fastq"

    # ✅ If FASTQs already exist AND are non-empty → skip
  #  if [[ -s "$fq1" && -s "$fq2" ]]; then
  #      echo "✅ FASTQs already exist for $DesiredPeakName → skipping."
  #      continue
  #  fi

    # ✅ Skip if BAM missing
  #  if [[ ! -f "$bam" ]]; then
  #      echo "⚠️ BAM not found: $bam"
  #      continue
  #  fi

  #  echo "Extracting ALL reads from $bam → $DesiredPeakName ..."

  #  samtools fastq \
  #      -1 "$fq1" \
  #      -2 "$fq2" \
  #      -0 /dev/null \
  #      -s /dev/null \
  #      -n \
#        "$bam"

#done


# ================================
# Remap FASTQs to genome
# ================================
# FILES="$FASTQDIR/*_R1.fastq"
# for f in $FILES; do
#     file=$(basename "$f")
#     sample=${file/_R1.fastq/}
#     read2="${f/_R1.fastq/_R2.fastq}"
#
#     bam="${REMAPPED_BAMDIR}/${sample}.bam"
#   #  QualityBam="${REMAPPED_BAMDIR}/${sample}_Q30.bam"
#     bigwig="${BIGWIGDIR}/${sample}"
#
#     echo "Mapping $sample..."
#     bwa mem -M -v 3 -t $THREADS $GENOME "$f" "$read2" \
#         | samtools view -bhSu - \
#         | samtools sort -@ $THREADS -T "$REMAPPED_BAMDIR/tempReps" -o "$bam"
#
#     samtools index "$bam"
#
#     echo "Generating BigWigs for $sample..."
#     bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 \
#         --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
# done
#
# echo "✅ Remapping and BigWig generation complete."
# echo "Skipped BAMs (no mapped reads) logged at $OUTDIR/skipped_no_mapped_reads.log"
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do

    bam="${BAMDIR}/${bamReads}"

    # Skip missing BAMs
    if [[ ! -f "$bam" ]]; then
        echo "❌ BAM not found: $bam"
        continue
    fi

    # Skip empty BAMs
    read_count=$(samtools view -c "$bam")
    if [[ "$read_count" -eq 0 ]]; then
        echo "⚠️ BAM is empty, skipping: $bam"
        continue
    fi

    # Output name
    bigwig="${BIGWIGDIR}/${DesiredPeakName}.bin_${BIN}.smooth_${SMOOTH}.bw"

    echo "→ Generating BigWig for sample: $DesiredPeakName"

    bamCoverage \
        -p $THREADS \
        -b "$bam" \
        -bs $BIN \
        --minMappingQuality 10 \
        --smoothLength $SMOOTH \
        --normalizeUsing BPM \
        -o "$bigwig"

done

echo "✅ BigWig generation complete (no remapping used)."
