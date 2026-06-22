#!/bin/bash
#SBATCH --job-name=FiberTools_Pt2_CallNucs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_Pt2_CallNucs.%j.out
#SBATCH --error=FiberTools_Pt2_CallNucs.%j.err

# =============================================================================
# FiberTools Part 2 — Call nucleosomes and generate pileups
# =============================================================================
# Loops over MergedBams/ONTRun9 and MergedBams/ONTRun10 separately.
# For each merged BAM:
#   1. ft add-nucleosomes  → *.nucs.bam
#   2. ft pileup (m6A, 5mC, nucleosome coverage)
#
# Output layout:
#   fibertools_results/
#     ONTRun9/
#       ONTRun9_WT_Eddie_barcode01_merged/
#         *.nucs.bam
#         *.m6Apileup.bedgraph
#         *.5mcpileup.bedgraph
#         *.nucspileup.bedgraph
#     ONTRun10/
#       (same structure)
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

MERGED_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedBams"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
THREADS=8

mkdir -p "$OUT_DIR"

echo "============================================================"
echo " FiberTools Pt2: Call nucleosomes + pileups"
echo " Input : $MERGED_DIR"
echo " Output: $OUT_DIR"
echo "============================================================"

for RUN in ONTRun9 ONTRun10; do

    RUN_IN="$MERGED_DIR/$RUN"
    RUN_OUT="$OUT_DIR/$RUN"
    mkdir -p "$RUN_OUT"

    if [ ! -d "$RUN_IN" ]; then
        echo "⚠️  $RUN_IN not found — skipping $RUN"
        continue
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Processing $RUN"
    echo "╚══════════════════════════════════════════════════════╝"

    for BAM in "$RUN_IN"/*_merged.bam; do
        [ -f "$BAM" ] || continue

        SAMPLE=$(basename "$BAM" .bam)          # e.g. ONTRun9_WT_Eddie_barcode01_merged
        SAMPLE_DIR="$RUN_OUT/$SAMPLE"
        mkdir -p "$SAMPLE_DIR"

        NUCS_BAM="$SAMPLE_DIR/${SAMPLE}.nucs.bam"

        echo ""
        echo "  ── $SAMPLE ──"

        # ── Step 1: add-nucleosomes ──────────────────────────────
        if [ ! -f "$NUCS_BAM" ]; then
            echo "    ft add-nucleosomes..."
            ft add-nucleosomes "$BAM" "$NUCS_BAM" \
                || { echo "    ❌ add-nucleosomes failed for $SAMPLE"; continue; }
            samtools index "$NUCS_BAM"
        else
            echo "    ⏭️  nucs.bam already exists — skipping add-nucleosomes"
        fi

        # ── Step 2: pileups ─────────────────────────────────────
        M6A_BG="$SAMPLE_DIR/${SAMPLE}.m6Apileup.bedgraph"
        CPG_BG="$SAMPLE_DIR/${SAMPLE}.5mcpileup.bedgraph"
        NUC_BG="$SAMPLE_DIR/${SAMPLE}.nucspileup.bedgraph"

        if [ ! -f "$M6A_BG" ]; then
            echo "    ft pileup — m6A..."
            ft pileup --m6a --per-base --fiber-coverage \
                --out "$M6A_BG" "$NUCS_BAM" \
                || echo "    ❌ m6A pileup failed for $SAMPLE"
        else
            echo "    ⏭️  m6A pileup exists — skipping"
        fi

        if [ ! -f "$CPG_BG" ]; then
            echo "    ft pileup — 5mC..."
            ft pileup --cpg --per-base --fiber-coverage \
                --out "$CPG_BG" "$NUCS_BAM" \
                || echo "    ❌ 5mC pileup failed for $SAMPLE"
        else
            echo "    ⏭️  5mC pileup exists — skipping"
        fi

        if [ ! -f "$NUC_BG" ]; then
            echo "    ft pileup — nucleosome coverage..."
            ft pileup --per-base --fiber-coverage \
                --out "$NUC_BG" "$NUCS_BAM" \
                || echo "    ❌ nuc pileup failed for $SAMPLE"
        else
            echo "    ⏭️  nuc pileup exists — skipping"
        fi

        echo "    ✅  $SAMPLE done"

    done
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " All done. Results:"
echo "============================================================"
for RUN in ONTRun9 ONTRun10; do
    echo ""
    echo "  $RUN:"
    find "$OUT_DIR/$RUN" -name "*.nucs.bam" | sort | while read -r BAM; do
        SIZE=$(du -h "$BAM" | cut -f1)
        echo "    $SIZE  $(basename $BAM)"
    done
done

echo ""
echo "Next: run FiberTools_Pt3_ExtractData.sh"
echo "  Point WORKDIR to: $OUT_DIR"
