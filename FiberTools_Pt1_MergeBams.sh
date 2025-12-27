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
BASE_DIR="/project/zallab/SequencingArchive/Oxford_Nanopore_Runs/ONTRun11/bam_pass"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun11/merged_bams_passed"
META="/lustre2/scratch/ry00555/ONTRun11/ONTRun11.txt"

mkdir -p "$OUT_DIR"

# Step 1-2: Filter, sort, and merge BAMs
for BARCODE_DIR in "$BASE_DIR"/barcode*/; do
    BARCODE_NAME=$(basename "$BARCODE_DIR")
    SORTED_DIR="${OUT_DIR}sorted_q10"
    mkdir -p "$SORTED_DIR"

    # Filter + Sort each BAM
    SORTED_BAMS=()
    for BAM in "$BARCODE_DIR"/*.bam; do
        [ -f "$BAM" ] || continue  # skip if no BAMs
        BAM_BASENAME=$(basename "$BAM" .bam)
        SORTED_BAM="$SORTED_DIR/${BAM_BASENAME}_q10_sorted.bam"
        samtools view -b -q 9 "$BAM" | samtools sort -@ 4 -o "$SORTED_BAM"
        SORTED_BAMS+=("$SORTED_BAM")
    done

    # Merge all sorted BAMs at once
    if [ ${#SORTED_BAMS[@]} -gt 0 ]; then
        # Find the sample_id corresponding to this barcode
        SAMPLE_ID=$(awk -v barcode="$BARCODE_NAME" -F'\t' '$1 ~ barcode {print $1}' "$META")
        [ -z "$SAMPLE_ID" ] && SAMPLE_ID="$BARCODE_NAME"  # fallback if not found

        FINAL_BAM="$OUT_DIR/${SAMPLE_ID}_merged.bam"
        samtools merge "$FINAL_BAM" "${SORTED_BAMS[@]}"
        samtools index "$FINAL_BAM"
        echo "✅ Merged and indexed: $FINAL_BAM"
    else
        echo "⚠️ No BAM files found for $BARCODE_NAME, skipping..."
    fi
done

echo "🎉 All barcodes processed! Merged BAMs at: $OUT_DIR"
