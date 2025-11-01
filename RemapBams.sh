#!/bin/bash
#SBATCH --job-name=BAMtoFASTQ
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=5:00:00
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
REMAPPED_BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/Remapped/SortedBamFiles"
FASTQDIR="/scratch/ry00555/RNASeqPaper/Oct2025/Remapped/FASTQ_from_BAM"
BIGWIGDIR="/scratch/ry00555/RNASeqPaper/Oct2025/Remapped/BigWigs"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025"
META="${OUTDIR}/BAM_File_Metadata_with_index_merged_V2.csv"
mkdir -p "$REMAPPED_BAMDIR" "$FASTQDIR" "$REMAPPED_BAMDIR/tempReps" "$BIGWIGDIR"
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
# tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do
#     bam="${BAMDIR}/${bamReads}"
#     if [[ ! -f "$bam" ]]; then
#         echo "⚠️ BAM file not found: $bam"
#         continue
#     fi
#
#     mapped=$(samtools view -c -F 4 "$bam")
#     if [[ $mapped -gt 0 ]]; then
#         echo "Extracting FASTQ from $bam (SampleID: $DesiredPeakName, mapped reads: $mapped)..."
#         samtools fastq -1 "$FASTQDIR/${DesiredPeakName}_R1.fastq" \
#                        -2 "$FASTQDIR/${DesiredPeakName}_R2.fastq" \
#                        -0 /dev/null -s /dev/null -n "$bam"
#     else
#         echo "⚠️ $bam has no mapped reads, skipping."
#     fi
# done


# ================================
# Remap FASTQs to genome
# ================================
FILES="$FASTQDIR/*_R1.fastq"
for f in $FILES; do
    file=$(basename "$f")
    sample=${file/_R1.fastq/}
    read2="${f/_R1.fastq/_R2.fastq}"

    bam="${REMAPPED_BAMDIR}/${sample}.bam"
    QualityBam="${REMAPPED_BAMDIR}/${sample}_Q30.bam"
    bigwig="${BIGWIGDIR}/${sample}"

    echo "Mapping $sample..."
    bwa mem -M -v 3 -t $THREADS $GENOME "$f" "$read2" \
        | samtools view -bhSu - \
        | samtools sort -@ $THREADS -T "$REMAPPED_BAMDIR/tempReps" -o "$bam"

    samtools index "$bam"

    echo "Generating BigWigs for $sample..."
    bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 \
        --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
done

echo "✅ Remapping and BigWig generation complete."
echo "Skipped BAMs (no mapped reads) logged at $OUTDIR/skipped_no_mapped_reads.log"
