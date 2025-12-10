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
OUT_DIR="/lustre2/scratch/ry00555/ONTRun10/fibertools_results"
mkdir -p "$OUT_DIR"

# Loop over merged BAMs
for BAM in "$IN_DIR"/*_merged.bam; do
    [ -f "$BAM" ] || continue
    SAMPLE=$(basename "$BAM" .bam)
    SAMPLE_DIR="$OUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"
    echo "Processing $SAMPLE ..."

    # 1) Predict m6A (produces a new BAM with MM/ML tags or updates tags in-place)
    # replace the flags below with whatever flags you prefer after checking ft predict-m6a --help
    PRED="${SAMPLE_DIR}/${SAMPLE}.predicted.bam"
    ft predict-m6a --threads 8 --input "$BAM" --output "$PRED" || { echo "predict-m6a failed for $SAMPLE"; continue; }

    # 2) Add nucleosomes (consumes predicted BAM, writes nucleosome calls to BAM or BED)
    NUCS_BAM="${SAMPLE_DIR}/${SAMPLE}.nucs.bam"
    ft add-nucleosomes --input "$PRED" --output "$NUCS_BAM" || { echo "add-nucleosomes failed for $SAMPLE"; continue; }

    # 3) Make decorated BED/track files (bed of nucleosome regions, footprints, etc.)
    ft track-decorators --input "$NUCS_BAM" --outdir "$SAMPLE_DIR/tracks" || echo "track-decorators failed for $SAMPLE"

    # 4) Make pileup (per-base or per-feature aggregation) for genome browser (creates bedGraph-like output)
    ft pileup --input "$NUCS_BAM" --out "$SAMPLE_DIR/${SAMPLE}.pileup.bedgraph" || echo "pileup failed for $SAMPLE"

    # Index the BAMs
    samtools index "$NUCS_BAM"

    echo "Done $SAMPLE"
done

echo "All done."
