#!/bin/bash
#SBATCH --job-name=ft_run_all
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=ft_run_all.%j.out
#SBATCH --error=ft_run_all.%j.err

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

IN_DIR="/scratch/ry00555/ONTRun9_10Combined/InputBams"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
mkdir -p "$OUT_DIR"

THREADS=8

for BAM in "$IN_DIR"/*_merged.bam; do
    [ -f "$BAM" ] || continue
    SAMPLE=$(basename "$BAM" .bam)
    SAMPLE_DIR="$OUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"
    echo "Processing $SAMPLE ..."

    # 1) Add nucleosomes (consumes predicted BAM)
    NUCS_BAM="${SAMPLE_DIR}/${SAMPLE}.nucs.bam"
    ft add-nucleosomes "$BAM" "$NUCS_BAM" || { echo "add-nucleosomes failed for $SAMPLE"; continue; }

    # 3) Make decorated BED/track files
    ft track-decorators --bed12 "$SAMPLE_DIR/tracks" "$NUCS_BAM" --decorator || echo "track-decorators failed for $SAMPLE"

    # 4) Make pileup (per-base or per-feature aggregation)
    ft pileup --m6a --cpg --fiber-coverage "$NUCS_BAM" --out "$SAMPLE_DIR/${SAMPLE}.pileup.bedgraph" || echo "pileup failed for $SAMPLE"

ft qc --m6a-per-msp "$NUCS_BAM" --out "$SAMPLE_DIR/${SAMPLE}.txt" || echo "qc failed for $SAMPLE"

    # Index the final BAM
    samtools index "$NUCS_BAM"

    echo "Done $SAMPLE"
done
