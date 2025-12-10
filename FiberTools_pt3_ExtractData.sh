#!/bin/bash
#SBATCH --job-name=ft_pt3
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=ft_pt3.%j.out
#SBATCH --error=ft_pt3.%j.err

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

#!/bin/bash
#SBATCH --job-name=fibertools_TSSplots
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=TSSplots.out
#SBATCH --error=TSSplots.err

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

WORKDIR="/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
THREADS=8
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora.bed"

echo "Starting analysis..."
echo "WORKDIR: $WORKDIR"
echo "GENOME: $GENOME"
echo "TSS BED: $TSS_BED"
echo ""

### ─────────────────────────────────────────────
### FIND ALL BAMs (subdirectories included)
### ─────────────────────────────────────────────
find "$WORKDIR" -maxdepth 2 -type f -name "*_merged.nucs.bam" | while read -r BAM; do
    SAMPLE=$(basename "$BAM" _merged.nucs.bam)
    SAMPLE_DIR=$(dirname "$BAM")

    echo "────────────────────────────────────────"
    echo " Processing sample: $SAMPLE"
    echo " Directory: $SAMPLE_DIR"
    echo "────────────────────────────────────────"

    ### ─────────────────────────────────────────
    ### 1. Identify existing pileups
    ###    (nucs and m6A and 5mC)
    ### ─────────────────────────────────────────
    NUCPU="${SAMPLE_DIR}/${SAMPLE}.nucspileup.bedgraph"
    M6APU="${SAMPLE_DIR}/${SAMPLE}.m6Apileup.bedgraph"
    MCPPU="${SAMPLE_DIR}/${SAMPLE}.5mcpileup.bedgraph"

    ### ─────────────────────────────────────────
    ### 2. Convert BEDGRAPH → BigWig (if exists)
    ### ─────────────────────────────────────────
    for PU in "$NUCPU" "$M6APU" "$MCPPU"; do
        if [ -f "$PU" ]; then
            BW="${PU%.bedgraph}.bw"
            echo "Converting to bigwig: $PU → $BW"
            bedGraphToBigWig "$PU" "$GENOME" "$BW" || echo "BigWig conversion failed for $PU"
        else
            echo "Pileup not found: $PU (skipping)"
        fi
    done

    ### ─────────────────────────────────────────
    ### 3. TSS computeMatrix + plotProfile
    ###    (run for each available bigwig)
    ### ─────────────────────────────────────────
    for BW in "${SAMPLE_DIR}/${SAMPLE}.nucspileup.bw" \
              "${SAMPLE_DIR}/${SAMPLE}.m6Apileup.bw" \
              "${SAMPLE_DIR}/${SAMPLE}.5mcpileup.bw"; do

        if [ -f "$BW" ]; then
            FEATURE=$(basename "$BW" .bw | sed 's/^.*\.//')   # nucspileup / m6Apileup / 5mcpileup

            MATRIX="${SAMPLE_DIR}/${SAMPLE}.${FEATURE}.TSS.matrix.gz"
            PNG="${SAMPLE_DIR}/${SAMPLE}.${FEATURE}.TSS_profile.png"

            echo "Running computeMatrix for $BW ..."
            computeMatrix reference-point \
                --referencePoint TSS \
                -b 2000 -a 1000 \
                -R "$TSS_BED" \
                -S "$BW" \
                -o "$MATRIX" \
                --skipZeros \
                --numberOfProcessors $THREADS

            echo "Plotting profile: $PNG"
            plotProfile -m "$MATRIX" \
                -out "$PNG" \
                --perGroup \
                --plotTitle "${SAMPLE} ${FEATURE} TSS Profile"
        fi
    done

    echo "Done with sample: $SAMPLE"
    echo ""
done

echo "All samples completed!"
