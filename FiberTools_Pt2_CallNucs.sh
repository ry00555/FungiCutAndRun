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

    # 1) Predict m6A (adds MM/ML tags)
    PRED="${SAMPLE_DIR}/${SAMPLE}.predicted.bam"
    fibertools predict-m6a \
        --threads $THREADS \
        --ml 125 \
        --keep \
        "$BAM" \
        "$PRED" || { echo "predict-m6a failed for $SAMPLE"; continue; }

    # 2) Add nucleosomes (consumes predicted BAM)
    NUCS_BAM="${SAMPLE_DIR}/${SAMPLE}.nucs.bam"
    fibertools add-nucleosomes \
        --input "$PRED" \
        --output "$NUCS_BAM" || { echo "add-nucleosomes failed for $SAMPLE"; continue; }

    # 3) Make decorated BED/track files
    fibertools track-decorators \
        --input "$NUCS_BAM" \
        --outdir "$SAMPLE_DIR/tracks" || echo "track-decorators failed for $SAMPLE"

    # 4) Make pileup (per-base or per-feature aggregation)
    fibertools pileup \
        --input "$NUCS_BAM" \
        --out "$SAMPLE_DIR/${SAMPLE}.pileup.bedgraph" || echo "pileup failed for $SAMPLE"

    # Index the final BAM
    samtools index "$NUCS_BAM"

    echo "Done $SAMPLE"
done
