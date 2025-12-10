#!/bin/bash
#SBATCH --job-name=ft_pt3
#SBATCH --partition=batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=ft_pt3.%j.out
#SBATCH --error=ft_pt3.%j.err


module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

WORKDIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
THREADS=8
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora.bed"
ml ucsc

for BAM in "$WORKDIR"/*_merged.nucs.bam; do
    [ -f "$BAM" ] || continue
    SAMPLE=$(basename "$BAM" _merged.nucs.bam)
    echo "Processing $SAMPLE ..."
    PILEUP="$WORKDIR/${SAMPLE}.pileup.bedgraph"
    BIGWIG="$WORKDIR/${SAMPLE}.pileup.bw"
    bedGraphToBigWig "$PILEUP" $GENOME "$BIGWIG" || echo "BigWig conversion failed for $SAMPLE"

TSS_MATRIX="$WORKDIR/${SAMPLE}.TSS.matrix.gz"
computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 1000 \
    -R "$TSS_BED" \
    -S "$BIGWIG" \
    -o "$TSS_MATRIX" \
    --skipZeros \
    --numberOfProcessors $THREADS

    plotProfile -m "$TSS_MATRIX" -out "$WORKDIR/${SAMPLE}.TSS_profile.png" \
        --perGroup --colors red

done
