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
BASE_DIR="/lustre2/scratch/ry00555/ONTRun9/bam_pass"
OUT_DIR="$BASE_DIR/merged_bams_passed"
mkdir -p "$OUT_DIR"

MAX_FILES_PER_BATCH=100  # adjust as needed

for BARCODE_DIR in "$BASE_DIR"/barcode*/; do
    BARCODE_NAME=$(basename "$BARCODE_DIR")
    SORTED_DIR="${BARCODE_DIR}sorted_q9"
    mkdir -p "$SORTED_DIR"

    # Step 1: Filter + Sort
    for BAM in "$BARCODE_DIR"/*.bam; do
        BAM_BASENAME=$(basename "$BAM" .bam)
        samtools view -b -q 9 "$BAM" | samtools sort -o "$SORTED_DIR/${BAM_BASENAME}_q9_sorted.bam"
    done

    # Step 2: Batch merge
    TMP_MERGED_LIST=()
    BATCH_INDEX=0
    FILE_COUNT=0
    BATCH_FILES=()

    for FILE in "$SORTED_DIR"/*.bam; do
        BATCH_FILES+=("$FILE")
        FILE_COUNT=$((FILE_COUNT + 1))

        if [[ $FILE_COUNT -eq $MAX_FILES_PER_BATCH ]]; then
            BATCH_OUT="$SORTED_DIR/batch_${BATCH_INDEX}.bam"
            samtools merge "$BATCH_OUT" "${BATCH_FILES[@]}"
            TMP_MERGED_LIST+=("$BATCH_OUT")
            BATCH_INDEX=$((BATCH_INDEX + 1))
            FILE_COUNT=0
            BATCH_FILES=()
        fi
    done

    # Merge remaining files (last batch)
    if [[ ${#BATCH_FILES[@]} -gt 0 ]]; then
        BATCH_OUT="$SORTED_DIR/batch_${BATCH_INDEX}.bam"
        samtools merge "$BATCH_OUT" "${BATCH_FILES[@]}"
        TMP_MERGED_LIST+=("$BATCH_OUT")
    fi

    # Step 3: Merge all batches into final BAM
    samtools merge "$OUT_DIR/${BARCODE_NAME}_merged.bam" "${TMP_MERGED_LIST[@]}"
    samtools index "$OUT_DIR/${BARCODE_NAME}_merged.bam"
done

echo "✅ Done! Batches merged and indexed at: $OUT_DIR"
