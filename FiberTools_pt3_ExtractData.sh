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
mkdir -p "$OUT_DIR"
THREADS=8

for BAM in "$WORKDIR"/*_merged.nucs.bam; do
    [ -f "$BAM" ] || continue
    SAMPLE=$(basename "$BAM" _merged.nucs.bam)
    echo "Processing $SAMPLE ..."

    ft extract 
