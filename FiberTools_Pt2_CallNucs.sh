#!/bin/bash
#SBATCH --job-name=ft_run_all
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
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
#ft add-nucleosomes "$BAM" "$NUCS_BAM" || { echo "add-nucleosomes failed for $SAMPLE"; continue; }
 #samtools index "$NUCS_BAM"

    # 3) Make decorated BED/track files
#ft track-decorators --bed12 "$SAMPLE_DIR/${SAMPLE}.nuctracks.bed" "$NUCS_BAM" --decorator "$SAMPLE_DIR/decorated_${SAMPLE}.nuctracks.bed" || echo "track-decorators failed for $SAMPLE"
#4) run qc
#ft qc --m6a-per-msp "$NUCS_BAM" "$SAMPLE_DIR/${SAMPLE}.txt" || echo "qc failed for $SAMPLE"

    # 5) Make pileup (per-base or per-feature aggregation)
ft pileup --m6a  --per-base  --fiber-coverage --out "$SAMPLE_DIR/${SAMPLE}.m6Apileup2.bedgraph" "$NUCS_BAM" || echo "pileup failed for $SAMPLE"
ft pileup --cpg --per-base --fiber-coverage --out "$SAMPLE_DIR/${SAMPLE}.5mcpileup2.bedgraph" "$NUCS_BAM" || echo "pileup failed for $SAMPLE"
ft pileup --fiber-coverage -per-base --out "$SAMPLE_DIR/${SAMPLE}.nucspileup2.bedgraph" "$NUCS_BAM" || echo "pileup failed for $SAMPLE"

#ft pileup --m6a --cpg --fiber-coverage --out "$SAMPLE_DIR/${SAMPLE}.totalinfo_pileup.bedgraph" "$NUCS_BAM" || echo "pileup failed for $SAMPLE"

#ft extract "$NUCS_BAM" --m6a "$SAMPLE_DIR"/"$SAMPLE"_m6a.bed --cpg "$SAMPLE_DIR"/"$SAMPLE"_cpg.bed --nuc "$SAMPLE_DIR"/"$SAMPLE"_nucleosome.bed --threads $THREADS || echo "extract failed for $SAMPLE"
#ft extract "$NUCS_BAM" --all "$SAMPLE_DIR"/"$SAMPLE"_totalinfo.bed --threads $THREADS || echo "extract failed for $SAMPLE"


    echo "Done $SAMPLE"
done
