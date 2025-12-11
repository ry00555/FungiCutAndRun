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


# ======================================
# Function: sort bedGraph properly
# ======================================
sort_bedgraph() {
    IN="$1"
    if [ -f "$IN" ]; then
        echo "Sorting $(basename "$IN")"
        TMP="${IN}.tmp"
        sort -k1,1 -k2,2n "$IN" > "$TMP"
        mv "$TMP" "$IN"
    fi
}

# ======================================
# Iterate over all BAMs in subdirectories
# ======================================
for BAM in "$WORKDIR"/*/*.nucs.bam; do
    [ -f "$BAM" ] || continue

    DIR=$(dirname "$BAM")
    FILE=$(basename "$BAM")
    SAMPLE=${FILE%.nucs.bam}

    echo "────────────────────────────────────────"
    echo " Processing sample: $SAMPLE"
    echo " Directory        : $DIR"
    echo "────────────────────────────────────────"

    # ----- Expected pileup files -----
    NUC_BED="$DIR/${SAMPLE}.nucspileup.bedgraph"
    M6A_BED="$DIR/${SAMPLE}.m6Apileup.bedgraph"
    CPG_BED="$DIR/${SAMPLE}.5mcpileup.bedgraph"

    # ---------- Sort and convert bedGraphs ----------
    for BEDGRAPH in "$NUC_BED" "$M6A_BED" "$CPG_BED"; do
        if [ -f "$BEDGRAPH" ]; then
            sort_bedgraph "$BEDGRAPH"
            BW="${BEDGRAPH%.bedgraph}.bw"
            if [ ! -f "$BW" ]; then
                echo "Converting $(basename $BEDGRAPH) → $(basename $BW)"
                bedGraphToBigWig "$BEDGRAPH" "$GENOME" "$BW"
            fi
        else
            echo "Pileup not found: $BEDGRAPH (skipping)"
        fi
    done

    # ---------- BigWig paths ----------
    NUC_BW="${NUC_BED%.bedgraph}.bw"
    M6A_BW="${M6A_BED%.bedgraph}.bw"
    CPG_BW="${CPG_BED%.bedgraph}.bw"

    # ---------- Generate TSS matrices & plots ----------
    for FEATURE in nuc m6A 5mC; do
        BW_VAR="${FEATURE^^}_BW"      # NUC_BW, M6A_BW, CPG_BW
        MATRIX="$DIR/${SAMPLE}.${FEATURE}.TSS.matrix.gz"
        PLOT="$DIR/${SAMPLE}.${FEATURE}.TSS_profile.png"
        COLOR="Reds"
        [ "$FEATURE" == "m6A" ] && COLOR="Greens"
        [ "$FEATURE" == "5mC" ] && COLOR="Blues"

        BW="${!BW_VAR}"
        if [ -f "$BW" ]; then
            echo "Generating TSS matrix and plot for $FEATURE"
            computeMatrix reference-point \
                --referencePoint TSS \
                -b 2000 -a 1000 \
                -R "$TSS_BED" \
                -S "$BW" \
                -o "$MATRIX" \
                --skipZeros \
                --numberOfProcessors "$THREADS"

            plotProfile \
                -m "$MATRIX" \
                -out "$PLOT" \
                --plotTitle "${SAMPLE} ${FEATURE}" \
                --colorMap "$COLOR"
        fi
    done

    # ---------- Optional: combined TSS plot ----------
    COMBINED_MATRIX="$DIR/${SAMPLE}.combined.TSS.matrix.gz"
    COMBINED_PLOT="$DIR/${SAMPLE}.combined.TSS_profile.png"
    EXISTING_BWS=()
    [ -f "$NUC_BW" ] && EXISTING_BWS+=("$NUC_BW")
    [ -f "$M6A_BW" ] && EXISTING_BWS+=("$M6A_BW")
    [ -f "$CPG_BW" ] && EXISTING_BWS+=("$CPG_BW")

    if [ "${#EXISTING_BWS[@]}" -gt 1 ]; then
        echo "Generating combined TSS matrix and plot"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 1000 \
            -R "$TSS_BED" \
            -S "${EXISTING_BWS[@]}" \
            -o "$COMBINED_MATRIX" \
            --skipZeros \
            --numberOfProcessors "$THREADS"

        plotProfile \
            -m "$COMBINED_MATRIX" \
            -out "$COMBINED_PLOT" \
            --plotTitle "${SAMPLE} Combined" \
            --perGroup \
            --colors "red green blue"
    fi

done

echo "All done!"
