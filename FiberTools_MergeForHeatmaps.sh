#!/usr/bin/env bash
#SBATCH --job-name=FiberTools_MergeHeatmaps
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_MergeHeatmaps.%j.out
#SBATCH --error=FiberTools_MergeHeatmaps.%j.err

# =============================================================================
# Merge nucs.bam files across replicates + both runs for heatmap generation
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

FT_RESULTS="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedForHeatmaps"
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora.bed"

mkdir -p "$OUT_DIR"

# =============================================================================
# Group definitions — one line per group
# Format: GROUP_NAME|barcode1 barcode2 barcode3
# Barcodes matched against sample dir names from Part 2 (both runs)
# =============================================================================
GROUPS=(
    "WT_Eddie|barcode01 barcode02 barcode13"
    "WT_Rochelle|barcode07 barcode08 barcode18"
    "cac-1|barcode03 barcode04 barcode14"
    "cac-2|barcode05 barcode06 barcode15"
    "rtt109|barcode09 barcode10 barcode19"
    "rtt109FLAG|barcode11 barcode12 barcode20"
    "gDNA_Hia5_ctrl|barcode16 barcode21"
    "WizardgDNA_HMW|barcode17 barcode22"
)

# =============================================================================
# Helper: bedgraph → BigWig
# =============================================================================
make_bigwig() {
    local BG="$1"
    local BW="$2"

    [ -f "$BW" ] && { echo "      ⏭️  $(basename $BW) exists"; return 0; }
    [ -f "$BG" ] || { echo "      ⚠️  $(basename $BG) not found — skipping"; return 0; }

    local BG4="${BG%.bedgraph}.4col.bedgraph"
    sort -k1,1 -k2,2n "$BG" \
        | cut -f1-4 \
        | awk -v genome="$GENOME" '
            BEGIN { while ((getline line < genome) > 0) {
                split(line,a,"\t"); chromlen[a[1]]=a[2] } }
            { if ($2<0) $2=0
              if (chromlen[$1] && $3>chromlen[$1]) $3=chromlen[$1]
              if ($2<$3) print $1"\t"$2"\t"$3"\t"$4 }
        ' > "$BG4"

    bedGraphToBigWig "$BG4" "$GENOME" "$BW" && rm -f "$BG4" \
        || { echo "      ❌ bedGraphToBigWig failed"; rm -f "$BG4"; return 1; }
    echo "      ✅  $(basename $BW)"
}

# =============================================================================
# Helper: computeMatrix + plotHeatmap
# =============================================================================
make_heatmap() {
    local BW="$1"
    local LABEL="$2"
    local SUFFIX="$3"
    local CMAP="$4"
    local DIR="$5"

    [ -f "$BW" ] || { echo "      ⚠️  $(basename $BW) not found — skipping heatmap"; return 0; }

    local MATRIX="$DIR/${LABEL}.${SUFFIX}.TSS.gz"
    local PNG="$DIR/${LABEL}.${SUFFIX}.TSS_profile.png"

    if [ ! -f "$MATRIX" ]; then
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 1000 \
            -R "$TSS_BED" \
            -S "$BW" \
            -o "$MATRIX" \
            --outFileNameMatrix "$DIR/${LABEL}.${SUFFIX}.TSS.tab" \
            --skipZeros -p 12
    fi

    plotHeatmap \
        -m "$MATRIX" \
        -out "$PNG" \
        --plotTitle "${LABEL} ${SUFFIX}" \
        --colorMap "$CMAP" \
        --whatToShow 'heatmap and colorbar'

    echo "      ✅  $(basename $PNG)"
}

# =============================================================================
# Main loop
# =============================================================================
echo "============================================================"
echo " Merging replicates + runs for heatmap generation"
echo " Output: $OUT_DIR"
echo "============================================================"

for ENTRY in "${GROUPS[@]}"; do

    GROUP="${ENTRY%%|*}"
    BARCODES="${ENTRY##*|}"

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Group: $GROUP  (barcodes: $BARCODES)"
    echo "╚══════════════════════════════════════════════════════╝"

    GROUP_DIR="$OUT_DIR/${GROUP}"
    MERGED_BAM="$OUT_DIR/${GROUP}_merged.nucs.bam"
    mkdir -p "$GROUP_DIR"

    # ── Collect matching nucs.bam files from both runs ───────────
    BAMS_TO_MERGE=()
    for RUN in ONTRun9 ONTRun10; do
        for BARCODE in $BARCODES; do
            for SAMPLE_DIR in "$FT_RESULTS/$RUN"/*${BARCODE}*/; do
                [ -d "$SAMPLE_DIR" ] || continue
                NUCS_BAM=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*.nucs.bam" | head -1)
                if [ -n "$NUCS_BAM" ] && [ -f "$NUCS_BAM" ]; then
                    echo "  Found: $(basename $NUCS_BAM)"
                    BAMS_TO_MERGE+=("$NUCS_BAM")
                else
                    echo "  ⚠️  No nucs.bam in $SAMPLE_DIR"
                fi
            done
        done
    done

    if [ ${#BAMS_TO_MERGE[@]} -eq 0 ]; then
        echo "  ❌  No nucs.bam files found for $GROUP — skipping"
        continue
    fi

    echo "  Total: ${#BAMS_TO_MERGE[@]} nucs.bam(s) to merge"

    # ── Merge ────────────────────────────────────────────────────
    if [ ! -f "$MERGED_BAM" ]; then
        if [ ${#BAMS_TO_MERGE[@]} -eq 1 ]; then
            cp "${BAMS_TO_MERGE[0]}" "$MERGED_BAM"
        else
            samtools merge -f -@ 8 "$MERGED_BAM" "${BAMS_TO_MERGE[@]}"
        fi
        samtools index "$MERGED_BAM"
        echo "  ✅  Merged: $(basename $MERGED_BAM)"
    else
        echo "  ⏭️  Merged BAM already exists"
    fi

    # ── Pileups ──────────────────────────────────────────────────
    echo "  Running pileups..."
    M6A_BG="$GROUP_DIR/${GROUP}.m6Apileup.bedgraph"
    CPG_BG="$GROUP_DIR/${GROUP}.5mcpileup.bedgraph"
    NUC_BG="$GROUP_DIR/${GROUP}.nucspileup.bedgraph"

    [ -f "$M6A_BG" ] || ft pileup --m6a --per-base --fiber-coverage \
        --out "$M6A_BG" "$MERGED_BAM" || echo "  ❌ m6A pileup failed"

    [ -f "$CPG_BG" ] || ft pileup --cpg --per-base --fiber-coverage \
        --out "$CPG_BG" "$MERGED_BAM" || echo "  ❌ 5mC pileup failed"

    [ -f "$NUC_BG" ] || ft pileup --per-base --fiber-coverage \
        --out "$NUC_BG" "$MERGED_BAM" || echo "  ❌ nuc pileup failed"

    # ── BigWigs ──────────────────────────────────────────────────
    echo "  Converting to BigWig..."
    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    make_bigwig "$NUC_BG" "$NUC_BW"
    make_bigwig "$M6A_BG" "$M6A_BW"
    make_bigwig "$CPG_BG" "$CPG_BW"

    # ── Heatmaps ─────────────────────────────────────────────────
    echo "  Generating heatmaps..."
    make_heatmap "$NUC_BW" "$GROUP" "nuc" "Reds"   "$GROUP_DIR"
    make_heatmap "$M6A_BW" "$GROUP" "m6A" "Greens" "$GROUP_DIR"
    make_heatmap "$CPG_BW" "$GROUP" "5mC" "Blues"  "$GROUP_DIR"

done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " All done. Heatmaps generated:"
echo "============================================================"
find "$OUT_DIR" -name "*.TSS_profile.png" | sort | while read -r PNG; do
    echo "  $(basename $PNG)"
done
