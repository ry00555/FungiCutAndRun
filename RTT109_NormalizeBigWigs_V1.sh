#!/bin/bash
#SBATCH --job-name=RTT109NormalizeBigWigs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=../RTT109NormalizeBigWigs.%j.out
#SBATCH --error=../RTT109NormalizeBigWigs.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

OUTDIR="/scratch/ry00555/RTT109PaperFigures"

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

OUTDIR="/scratch/ry00555/RTT109PaperFigures"
BW_DIR="${OUTDIR}/BigWigs"
OUT_NORM="${BW_DIR}/NormalizedBigWigs"
MERGED_DIR="${OUT_NORM}/MergedBigWigs"
ml deepTools

mkdir -p "$OUT_NORM"

echo "🔍 Processing meta file: $PAIRFILE"
echo "   BigWigs directory: $BW_DIR"
echo

# Read header line
read -r header < "$PAIRFILE"
echo "Header columns: $header"

# Initialize arrays
missing_files=()
processed=()

# Skip header and loop through rows
tail -n +2 "$PAIRFILE" | while IFS=$'\t' read -r RunSample ID Strain Antibody Rep BigWig Input; do
    chip_path="${BW_DIR}/${BigWig}"
    input_path="${BW_DIR}/${Input}"

    # Skip if ChIP file is missing
#     if [[ ! -f "$chip_path" ]]; then
#         echo "⚠️ Missing ChIP file: $chip_path"
#         missing_files+=("$BigWig")
#         continue
#     fi
#
#     outbase="${OUT_NORM}/${ID}_R${Rep}"
#
#     # If Input exists, normalize
#     if [[ -n "$Input" && -f "$input_path" ]]; then
#         echo "Normalizing $BigWig vs $Input → ${outbase}_foldchange.bw"
#         bigwigCompare \
#           -b1 "$chip_path" \
#           -b2 "$input_path" \
#           --operation ratio \
#           --pseudocount 1 \
#           --binSize 25 \
#           -o "${outbase}_foldchange.bw" \
#           --skipZeroOverZero \
#           --verbose
#     else
#         # No Input → copy ChIP as normalized
#         echo "⚠️ No Input for $BigWig → copying as ${outbase}.bw"
#         cp "$chip_path" "${outbase}.bw"
#     fi
#
#     processed+=("$ID R$Rep")
# done
#
# # Summary
# echo
# echo "========== PROCESSING SUMMARY =========="
# echo "✅ Processed files: ${#processed[@]}"
# for f in "${processed[@]}"; do
#     echo "   $f"
# done
#
# echo
# echo "⚠️ Missing ChIP files: ${#missing_files[@]}"
# for m in "${missing_files[@]}"; do
#     echo "   $m"
# done

ml ucsc  # load UCSC tools

# Find all bigwig files
bw_files=($(ls "$OUT_NORM"/*.bw 2>/dev/null))
echo "Found ${#bw_files[@]} .bw files"

# Loop over unique IDs
for id in $(for f in "${bw_files[@]}"; do
        basename "$f" .bw | sed -E 's/_R[0-9]+(_foldchange)?$//'
      done | sort | uniq); do

    echo "Processing ID: $id"

    # Gather all replicates for this ID
    rep_files=($(ls "$OUT_NORM/${id}"_R*.bw "$OUT_NORM/${id}"_R*_foldchange.bw 2>/dev/null))
    echo "Replicates found: ${rep_files[@]}"

    if [ ${#rep_files[@]} -eq 0 ]; then
        echo "No files found for $id, skipping..."
        continue
    fi

    # Print command before running
    echo "Running multiBigwigSummary for $id..."
    echo "multiBigwigSummary bins -b ${rep_files[@]} --binSize 25 --outRawCounts $MERGED_DIR/${id}_summary.tab -o $MERGED_DIR/${id}_summary.npz"

    # Uncomment when ready
    # multiBigwigSummary bins \
    #     -b "${rep_files[@]}" \
    #     --binSize 25 \
    #     --outRawCounts "$MERGED_DIR/${id}_summary.tab" \
    #     -o "$MERGED_DIR/${id}_summary.npz"

    # Check if .tab exists before converting
    if [ -f "$MERGED_DIR/${id}_summary.tab" ]; then
        echo "Converting to BEDGraph..."
        awk 'NR>1 {print $1 "\t" $2 "\t" $3 "\t" $4}' "$MERGED_DIR/${id}_summary.tab" \
            | sort -k1,1 -k2,2n \
            > "$MERGED_DIR/${id}_summary.bedGraph"

        echo "Creating bigWig..."
        bedGraphToBigWig "$MERGED_DIR/${id}_summary.bedGraph" \
            "/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt" \
            "$MERGED_DIR/${id}_summary.bw"
    else
        echo "Summary tab not found for $id, skipping conversion"
    fi

done

echo "Done!"

# multiBigwigSummary BED-file \
#   --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K27me3*.bw \
#   --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
#   -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_K27genes_signal_matrix.npz \
#   --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_K27genesonly_signal_matrix.tab
#
#   multiBigwigSummary BED-file \
#     --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
#     --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
#     -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3ChIP_K27genes_signal_matrix.npz \
#     --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigsH_3K36me3ChIP__K27genesonly_signal_matrix.tab
#
# multiBigwigSummary BED-file \
#       --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
#       --BED "/scratch/ry00555/GeneList_BedFiles/NonK27genes.bed" \
#       -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NonK27genes_signal_matrix.npz \
#       --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_NonK27genesonly_signal_matrix.tab
#
#       multiBigwigSummary BED-file \
#         --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K27me3*.bw \
#         --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
#         -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_K27regions_signal_matrix.npz \
#         --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_H3K27me3regions_signal_matrix.tab
#
#
#         multiBigwigSummary BED-file \
#           --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
#           --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
#           -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3_K27regions_signal_matrix.npz \
#           --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3_NormalizedBigWigs_H3K27me3regions_signal_matrix.tab
