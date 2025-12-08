#!/bin/bash
#SBATCH --job-name=FiberTools_Pt1_MergeBams
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=../FiberTools_Pt1_MergeBams.%j.out
#SBATCH --error=../FiberTools_Pt1_MergeBams.%j.err

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools
BASE_DIR="/lustre2/scratch/ry00555/ONTRun10/bam_pass"
OUT_DIR="$BASE_DIR/merged_bams_passed"
META="/lustre2/scratch/ry00555/ONTRun9/ONTRun9.txt"

mkdir -p "$OUT_DIR"

# ---------------------------
# Step 0: Rename barcode folders based on sample sheet
# ---------------------------
echo "🔄 Renaming barcode folders based on sample sheet..."
tail -n +2 "$META" | while IFS=$'\t' read -r SAMPLE_ID EXPERIMENT_ID INDEX_KIT INDEX OWNER STRAIN EXP_TYPE ANTIBODY GDNA_KIT ELUTION_CONC ELUTION_VOL TOTAL_DNA NOTES FINAL_CONC FINAL_ELUTION UL_200NG H2O_TO_ADD; do
    # Extract the barcode from SAMPLE_ID (assumes format ONT9_barcodeXX_...)
    BARCODE=$(echo "$SAMPLE_ID" | cut -d'_' -f2)
    if [ -d "$BASE_DIR/$BARCODE" ]; then
        echo "Renaming $BARCODE --> $SAMPLE_ID"
        mv "$BASE_DIR/$BARCODE" "$BASE_DIR/$SAMPLE_ID"
    fi
done

# ---------------------------
# Step 1-2: Filter, sort, and merge BAMs
# ---------------------------
for BARCODE_DIR in "$BASE_DIR"/barcode*/; do
    BARCODE_NAME=$(basename "$BARCODE_DIR")
    SORTED_DIR="${BARCODE_DIR}sorted_q9"
    mkdir -p "$SORTED_DIR"

    # Filter + Sort each BAM
    SORTED_BAMS=()
    for BAM in "$BARCODE_DIR"/*.bam; do
        [ -f "$BAM" ] || continue  # skip if no BAMs
        BAM_BASENAME=$(basename "$BAM" .bam)
        SORTED_BAM="$SORTED_DIR/${BAM_BASENAME}_q9_sorted.bam"
        samtools view -b -q 9 "$BAM" | samtools sort -@ 4 -o "$SORTED_BAM"
        SORTED_BAMS+=("$SORTED_BAM")
    done

    # Merge all sorted BAMs at once
    FINAL_BAM="$OUT_DIR/${BARCODE_NAME}_merged.bam"
    if [ ${#SORTED_BAMS[@]} -gt 0 ]; then
        samtools merge "$FINAL_BAM" "${SORTED_BAMS[@]}"
        samtools index "$FINAL_BAM"
        echo "✅ Merged and indexed: $FINAL_BAM"
    else
        echo "⚠️ No BAM files found for $BARCODE_NAME, skipping..."
    fi
done

echo "🎉 All barcodes processed! Merged BAMs at: $OUT_DIR"
