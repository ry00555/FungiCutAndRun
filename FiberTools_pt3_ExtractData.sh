#!/bin/bash
#SBATCH --job-name=FiberTools_Pt3_Extract
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=85G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_Pt3_Extract.%j.out
#SBATCH --error=FiberTools_Pt3_Extract.%j.err

# =============================================================================
# FiberTools Part 3 — Convert pileups to BigWig, generate TSS heatmaps
# =============================================================================
# Loops over ONTRun9 and ONTRun10 results separately.
# For each sample:
#   1. Sort bedgraphs
#   2. Convert bedgraph → BigWig (bedGraphToBigWig)
#   3. computeMatrix + plotHeatmap around TSS for m6A, 5mC, nucleosomes
#
# Input layout (from Part 2):
#   fibertools_results/
#     ONTRun9/<sample>/  *.nucspileup.bedgraph  *.m6Apileup.bedgraph  *.5mcpileup.bedgraph
#     ONTRun10/<sample>/ (same)
#
# Output (added to each sample dir):
#   *.nuc.bw  *.m6A.bw  *.5mC.bw
#   *.nuc.TSS_profile.png
#   *.m6A.TSS_profile.png
#   *.5mC.TSS_profile.png
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

WORKDIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora.bed"

echo "============================================================"
echo " FiberTools Pt3: BigWig conversion + TSS heatmaps"
echo " Working dir: $WORKDIR"
echo "============================================================"

# ── Helper: sort a bedgraph in place ────────────────────────────
sort_bedgraph() {
    local IN="$1"
    if [ -f "$IN" ]; then
        local TMP="${IN}.tmp"
        sort -k1,1 -k2,2n "$IN" > "$TMP" && mv "$TMP" "$IN"
    fi
}

# ── Helper: bedgraph → BigWig ────────────────────────────────────
# Clips any coordinates that exceed chrom sizes (common with ONT data),
# keeps only 4 columns, then converts.
make_bigwig() {
    local BG="$1"
    local BW="$2"
    local GENOME="$3"

    if [ -f "$BW" ]; then
        echo "      ⏭️  $(basename $BW) already exists"
        return 0
    fi

    if [ ! -f "$BG" ]; then
        echo "      ⚠️  bedgraph not found: $(basename $BG) — skipping"
        return 0
    fi

    echo "      Converting $(basename $BG) → $(basename $BW)"

    local BG4="${BG%.bedgraph}.4col.bedgraph"

    # Sort, keep 4 cols, clip coordinates to chrom sizes
    sort_bedgraph "$BG"
    cut -f1-4 "$BG" \
        | awk -v genome="$GENOME" '
            BEGIN {
                while ((getline line < genome) > 0) {
                    split(line, a, "\t"); chromlen[a[1]] = a[2]
                }
            }
            {
                if ($2 < 0) $2 = 0
                if (chromlen[$1] && $3 > chromlen[$1]) $3 = chromlen[$1]
                if ($2 < $3) print $1"\t"$2"\t"$3"\t"$4
            }
        ' > "$BG4"

    bedGraphToBigWig "$BG4" "$GENOME" "$BW" \
        && rm -f "$BG4" \
        || { echo "      ❌ bedGraphToBigWig failed for $(basename $BG)"; rm -f "$BG4"; return 1; }
}

# ── Main loop ────────────────────────────────────────────────────
for RUN in ONTRun9 ONTRun10; do

    RUN_DIR="$WORKDIR/$RUN"
    if [ ! -d "$RUN_DIR" ]; then
        echo "⚠️  $RUN_DIR not found — skipping $RUN"
        continue
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Processing $RUN"
    echo "╚══════════════════════════════════════════════════════╝"

    for NUCS_BAM in "$RUN_DIR"/*/*.nucs.bam; do
        [ -f "$NUCS_BAM" ] || continue

        DIR=$(dirname "$NUCS_BAM")
        FILE=$(basename "$NUCS_BAM")
        SAMPLE="${FILE%.nucs.bam}"

        echo ""
        echo "  ── $SAMPLE ──"

        # Expected bedgraph paths (matching Part 2 output names)
        NUC_BG="$DIR/${SAMPLE}.nucspileup.bedgraph"
        M6A_BG="$DIR/${SAMPLE}.m6Apileup.bedgraph"
        CPG_BG="$DIR/${SAMPLE}.5mcpileup.bedgraph"

        # BigWig paths
        NUC_BW="${NUC_BG%.bedgraph}.bw"
        M6A_BW="${M6A_BG%.bedgraph}.bw"
        CPG_BW="${CPG_BG%.bedgraph}.bw"

        # ── Convert bedgraphs → BigWigs ──────────────────────────
        echo "    [1/2] BigWig conversion"
        make_bigwig "$NUC_BG" "$NUC_BW" "$GENOME"
        make_bigwig "$M6A_BG" "$M6A_BW" "$GENOME"
        make_bigwig "$CPG_BG" "$CPG_BW" "$GENOME"

        # ── computeMatrix + plotHeatmap ──────────────────────────
        echo "    [2/2] TSS matrices + heatmaps"

        # 1. Nucleosome
        if [ -f "$NUC_BW" ]; then
            NUC_MATRIX="$DIR/${SAMPLE}.nuc.TSS.gz"
            if [ ! -f "$NUC_MATRIX" ]; then
                computeMatrix reference-point \
                    --referencePoint TSS \
                    -b 2000 -a 1000 \
                    -R "$TSS_BED" \
                    -S "$NUC_BW" \
                    -o "$NUC_MATRIX" \
                    --outFileNameMatrix "$DIR/${SAMPLE}.nuc.TSS.tab" \
                    --skipZeros -p 12
            fi
            plotHeatmap \
                -m "$NUC_MATRIX" \
                -out "$DIR/${SAMPLE}.nuc.TSS_profile.png" \
                --plotTitle "${SAMPLE} Nucleosomes" \
                --colorMap 'Reds' \
                --whatToShow 'heatmap and colorbar'
            echo "      ✅  nucleosome heatmap done"
        fi

        # 2. m6A
        if [ -f "$M6A_BW" ]; then
            M6A_MATRIX="$DIR/${SAMPLE}.m6A.TSS.gz"
            if [ ! -f "$M6A_MATRIX" ]; then
                computeMatrix reference-point \
                    --referencePoint TSS \
                    -b 2000 -a 1000 \
                    -R "$TSS_BED" \
                    -S "$M6A_BW" \
                    -o "$M6A_MATRIX" \
                    --outFileNameMatrix "$DIR/${SAMPLE}.m6A.TSS.tab" \
                    --skipZeros -p 12
            fi
            plotHeatmap \
                -m "$M6A_MATRIX" \
                -out "$DIR/${SAMPLE}.m6A.TSS_profile.png" \
                --plotTitle "${SAMPLE} m6A" \
                --colorMap 'Greens' \
                --whatToShow 'heatmap and colorbar'
            echo "      ✅  m6A heatmap done"
        fi

        # 3. 5mC  (note: uses $CPG_MATRIX not $NUC_MATRIX — bug fixed)
        if [ -f "$CPG_BW" ]; then
            CPG_MATRIX="$DIR/${SAMPLE}.5mC.TSS.gz"
            if [ ! -f "$CPG_MATRIX" ]; then
                computeMatrix reference-point \
                    --referencePoint TSS \
                    -b 2000 -a 1000 \
                    -R "$TSS_BED" \
                    -S "$CPG_BW" \
                    -o "$CPG_MATRIX" \
                    --outFileNameMatrix "$DIR/${SAMPLE}.5mC.TSS.tab" \
                    --skipZeros -p 12
            fi
            plotHeatmap \
                -m "$CPG_MATRIX" \
                -out "$DIR/${SAMPLE}.5mC.TSS_profile.png" \
                --plotTitle "${SAMPLE} 5mC" \
                --colorMap 'Blues' \
                --whatToShow 'heatmap and colorbar'
            echo "      ✅  5mC heatmap done"
        fi
        if [ ! -f "${DIR}/${SAMPLE}_nucleosome.bed" ]; then
            ft extract "$NUCS_BAM" \
                --nuc "${DIR}/${SAMPLE}_nucleosome.bed" \
                --m6a "${DIR}/${SAMPLE}_m6a.bed" \
                --cpg "${DIR}/${SAMPLE}_cpg.bed" \
                --threads 8
        fi
    done
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " Part 3 complete. Output files per sample:"
echo "============================================================"
for RUN in ONTRun9 ONTRun10; do
    echo ""
    echo "  $RUN:"
    find "$WORKDIR/$RUN" -name "*.TSS_profile.png" | sort | while read -r PNG; do
        echo "    $(basename $PNG)"
    done
done

echo ""
echo "All done!"
