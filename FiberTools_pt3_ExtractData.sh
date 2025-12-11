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

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

WORKDIR="/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
THREADS=8
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora.bed"

# iterate through all BAMs inside subdirectories
for BAM in "$WORKDIR"/*/*.nucs.bam; do
    DIR=$(dirname "$BAM")
    FILE=$(basename "$BAM")

    SAMPLE=${FILE%.nucs.bam}

    echo "────────────────────────────────────────"
    echo " Processing: $SAMPLE"
    echo " Directory : $DIR"
    echo "────────────────────────────────────────"

    # ----- Expected pileup files -----
    NUC_BED="$DIR/${SAMPLE}.nucspileup.bedgraph"
    M6A_BED="$DIR/${SAMPLE}.m6Apileup.bedgraph"
    CPG_BED="$DIR/${SAMPLE}.5mcpileup.bedgraph"

    # ---------- BigWig Conversion ----------
    for BEDGRAPH in "$NUC_BED" "$M6A_BED" "$CPG_BED"; do
        if [ -f "$BEDGRAPH" ]; then
            BW="${BEDGRAPH%.bedgraph}.bw"
            if [ ! -f "$BW" ]; then
                echo "Converting $(basename $BEDGRAPH) → $(basename $BW)"
                bedGraphToBigWig "$BEDGRAPH" "$GENOME" "$BW"
            fi
        else
            echo "Pileup not found: $BEDGRAPH (skipping)"
        fi
    done

    NUC_BW="${NUC_BED%.bedgraph}.bw"
    M6A_BW="${M6A_BED%.bedgraph}.bw"
    CPG_BW="${CPG_BED%.bedgraph}.bw"

    # ---------- Generate TSS matrices ----------
    # 1. nucleosome-only
    if [ -f "$NUC_BW" ]; then
        MATRIX="$DIR/${SAMPLE}.nuc.TSS.matrix.gz"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 1000 \
            -R "$TSS_BED" \
            -S "$NUC_BW" \
            -o "$MATRIX" \
            --skipZeros \
            --numberOfProcessors "$THREADS"

        plotProfile \
            -m "$MATRIX" \
            -out "$DIR/${SAMPLE}.nuc.TSS_profile.png" \
            --plotTitle "${SAMPLE} Nucleosomes" --colorMap 'Reds'
    fi

    # 2. m6A-only
    if [ -f "$M6A_BW" ]; then
        MATRIX="$DIR/${SAMPLE}.m6A.TSS.matrix.gz"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 1000 \
            -R "$TSS_BED" \
            -S "$M6A_BW" \
            -o "$MATRIX" \
            --skipZeros \
            --numberOfProcessors "$THREADS"

        plotProfile \
            -m "$MATRIX" \
            -out "$DIR/${SAMPLE}.m6A.TSS_profile.png" \
            --plotTitle "${SAMPLE} m6A" --colorMap 'Greens'
    fi

    # 3. 5mC-only
    if [ -f "$CPG_BW" ]; then
        MATRIX="$DIR/${SAMPLE}.5mC.TSS.matrix.gz"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 1000 \
            -R "$TSS_BED" \
            -S "$CPG_BW" \
            -o "$MATRIX" \
            --skipZeros \
            --numberOfProcessors "$THREADS"

        plotProfile \
            -m "$MATRIX" \
            -out "$DIR/${SAMPLE}.5mC.TSS_profile.png" \
            --plotTitle "${SAMPLE} 5mC" --colorMap 'Blues'


          fi

done
