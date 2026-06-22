#!/bin/bash
#SBATCH --job-name=FiberTools_Pt1_MergeBams
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=FiberTools_Pt1_MergeBams.%j.out
#SBATCH --error=FiberTools_Pt1_MergeBams.%j.err

# =============================================================================
# FiberTools Part 1 — Merge BAMs across ONTRun9 + ONTRun10
# =============================================================================
# For each barcode (01–22), merges all split BAMs from both runs into a single
# named merged BAM using sample/strain info from the metadata CSV.
#
# Input layout:
#   InputBams/
#     ONTRun9/barcode01/*.bam
#     ONTRun10/barcode01/*.bam
#     (same barcodes in both runs = same biological sample)
#
# Output layout:
#   MergedBams/
#     WT_barcode01_merged.bam
#     WT_barcode02_merged.bam
#     cac-1_barcode03_merged.bam
#     ...
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

BASE_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/InputBams"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedBams"
SORTED_TMP="${BASE_DIR}/sorted_q10_tmp"

mkdir -p "$OUT_DIR"
mkdir -p "$SORTED_TMP"

# =============================================================================
# Barcode → sample name mapping
# Derived from FiberMeta.csv (ONT9 metadata; same barcodes used for ONT10)
# Format: [barcode]="SampleName"
# Barcodes 23 and 24 are empty in the metadata — skipped.
# gDNA controls (16, 17, 21, 22) and HMW prep included — remove if unwanted.
# =============================================================================
declare -A SAMPLE_MAP=(
    [barcode01]="WT_Eddie_barcode01"
    [barcode02]="WT_Eddie_barcode02"
    [barcode03]="cac-1_Eddie_barcode03"
    [barcode04]="cac-1_Eddie_barcode04"
    [barcode05]="cac-2_Eddie_barcode05"
    [barcode06]="cac-2_Eddie_barcode06"
    [barcode07]="WT_Rochelle_barcode07"
    [barcode08]="WT_Rochelle_barcode08"
    [barcode09]="rtt109_Rochelle_barcode09"
    [barcode10]="rtt109_Rochelle_barcode10"
    [barcode11]="rtt109FLAG_Rochelle_barcode11"
    [barcode12]="rtt109FLAG_Rochelle_barcode12"
    [barcode13]="WT_Eddie_barcode13"
    [barcode14]="cac-1_Eddie_barcode14"
    [barcode15]="cac-2_Eddie_barcode15"
    [barcode16]="gDNA_Hia5_ctrl_Eddie_barcode16"
    [barcode17]="WizardgDNA_HMW_Eddie_barcode17"
    [barcode18]="WT_Rochelle_barcode18"
    [barcode19]="rtt109_Rochelle_barcode19"
    [barcode20]="rtt109FLAG_Rochelle_barcode20"
    [barcode21]="gDNA_Hia5_ctrl_Rochelle_barcode21"
    [barcode22]="WizardgDNA_HMW_Rochelle_barcode22"
)

echo "============================================================"
echo " FiberTools Pt1: Merge BAMs — ONTRun9 + ONTRun10 combined"
echo " Output dir: $OUT_DIR"
echo "============================================================"

for BARCODE in "${!SAMPLE_MAP[@]}"; do

    SAMPLE="${SAMPLE_MAP[$BARCODE]}"
    FINAL_BAM="$OUT_DIR/${SAMPLE}_merged.bam"

    # Skip if already merged
    if [ -f "$FINAL_BAM" ] && [ -f "${FINAL_BAM}.bai" ]; then
        echo "⏭️   Skipping $SAMPLE — merged BAM already exists"
        continue
    fi

    echo ""
    echo "────────────────────────────────────────────────────────"
    echo " Barcode : $BARCODE"
    echo " Sample  : $SAMPLE"
    echo "────────────────────────────────────────────────────────"

    SORTED_BAMS=()

    # ── Collect BAMs from both runs ──────────────────────────────
    for RUN in ONTRun9 ONTRun10; do
        SRC_DIR="$BASE_DIR/$RUN/$BARCODE"

        if [ ! -d "$SRC_DIR" ]; then
            echo "  ⚠️  $RUN/$BARCODE not found — skipping this run"
            continue
        fi

        BAM_COUNT=$(find "$SRC_DIR" -maxdepth 1 -name "*.bam" ! -name "*.bam.bai" | wc -l)
        echo "  $RUN: found $BAM_COUNT BAM file(s)"

        for BAM in "$SRC_DIR"/*.bam; do
            [ -f "$BAM" ] || continue
            [ "${BAM}" != "${BAM%.bai}" ] && continue   # skip .bai files

            BAM_BASE=$(basename "$BAM" .bam)
            SORTED_BAM="$SORTED_TMP/${SAMPLE}_${RUN}_${BAM_BASE}_q10_sorted.bam"

            if [ ! -f "$SORTED_BAM" ]; then
                echo "    Filtering + sorting: $(basename $BAM)"
                samtools view -b -q 10 "$BAM" \
                    | samtools sort -@ 4 -o "$SORTED_BAM"
            else
                echo "    Already sorted: $(basename $SORTED_BAM)"
            fi

            SORTED_BAMS+=("$SORTED_BAM")
        done
    done

    # ── Merge all sorted BAMs for this barcode ───────────────────
    if [ ${#SORTED_BAMS[@]} -eq 0 ]; then
        echo "  ❌  No BAMs found for $BARCODE in either run — skipping"
        continue
    fi

    echo "  Merging ${#SORTED_BAMS[@]} sorted BAM(s) → $(basename $FINAL_BAM)"

    if [ ${#SORTED_BAMS[@]} -eq 1 ]; then
        # Only one BAM — just copy rather than merge
        cp "${SORTED_BAMS[0]}" "$FINAL_BAM"
    else
        samtools merge -f -@ 4 "$FINAL_BAM" "${SORTED_BAMS[@]}"
    fi

    samtools index "$FINAL_BAM"
    echo "  ✅  Done: $FINAL_BAM"

done

# =============================================================================
# Clean up temp sorted BAMs (comment out if you want to keep them)
# =============================================================================
echo ""
echo "Cleaning up temp sorted BAMs..."
rm -rf "$SORTED_TMP"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " Merge complete. Final BAMs in: $OUT_DIR"
echo "============================================================"
find "$OUT_DIR" -name "*_merged.bam" | sort | while read -r BAM; do
    READS=$(samtools view -c "$BAM" 2>/dev/null || echo "?")
    SIZE=$(du -h "$BAM" | cut -f1)
    echo "  $SIZE  $(basename $BAM)  ($READS reads)"
done

echo ""
echo "Next step: update IN_DIR in FiberTools_Pt2_CallNucs.sh to:"
echo "  $OUT_DIR"
