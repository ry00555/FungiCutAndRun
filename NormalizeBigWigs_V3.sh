#!/bin/bash
#SBATCH --job-name=NormalizeBigWigs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=../NormalizeBigWigs.%j.out
#SBATCH --error=../NormalizeBigWigs.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

OUTDIR="/scratch/ry00555/RNASeqPaper"

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

#mkdir -p "${OUTDIR}/NormalizedBigWigs/Run135toRun150"

#!/bin/bash

PAIRFILE="${OUTDIR}/chip_input_pairs_Run135_Run150_H3K27me3only.txt"
BW_DIR="${OUTDIR}/CandidateBigWigs"
OUT_NORM="${OUTDIR}/NormalizedBigWigs/Run135toRun150"

ml deepTools

mkdir -p "$OUT_NORM"

echo "🔍 Checking pairs listed in: $PAIRFILE"
echo "   against files in: $BW_DIR"
echo

missing_pairs=()
processed_pairs=()

while IFS=$'\t' read -r chip_bw input_bw; do
    chip_path="${BW_DIR}/${chip_bw}"
    input_path="${BW_DIR}/${input_bw}"

    # Skip if ChIP file is missing
    if [[ ! -f "$chip_path" ]]; then
        echo "⚠️ Missing ChIP file: $chip_path"
        missing_pairs+=("$chip_bw")
        continue
    fi

    outname=$(basename "$chip_bw" .bw)_norm
    if [[ -n "$input_bw" && -f "$input_path" ]]; then
        # Input exists → do full normalization
        outname="${outname}_foldchange.bw"
        echo "Normalizing $chip_bw against $input_bw → $outname"
        bigwigCompare \
          -b1 "$chip_path" \
          -b2 "$input_path" \
          --operation ratio \
          --pseudocount 1 \
          --binSize 25 \
          -o "${OUT_NORM}/${outname}" \
          --skipZeroOverZero \
          --verbose

    else
        # No Input → just copy or rename the ChIP file
        outname="${outname}.bw"
        echo "⚠️ No Input for $chip_bw → copying as $outname"
        cp "$chip_path" "${OUT_NORM}/${outname}"
    fi

    processed_pairs+=("$chip_bw vs ${input_bw:-none}")
done < "$PAIRFILE"

echo
echo "========== PROCESSING SUMMARY =========="
echo "✅ Processed pairs/files: ${#processed_pairs[@]}"
for p in "${processed_pairs[@]}"; do
    echo "   $p"
done

echo
echo "⚠️ Missing ChIP files: ${#missing_pairs[@]}"
for m in "${missing_pairs[@]}"; do
    echo "   $m"
done

multiBigwigSummary BED-file \
  --bwfiles ${OUTDIR}/NormalizedBigWigs/Run135toRun150/*.bw \
  --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
  -out ${OUTDIR}/NormalizedBigWigs/Run135toRun150/K27genes_signal_matrix.npz \
  --outRawCounts ${OUTDIR}/NormalizedBigWigs/Run135toRun150/NormalizedBigWigs_K27genesonly_signal_matrix.tab

multiBigwigSummary BED-file \
      --bwfiles ${OUTDIR}/NormalizedBigWigs/Run135toRun150/*.bw \
      --BED "/scratch/ry00555/GeneList_BedFiles/NonK27genes.bed" \
      -out ${OUTDIR}/NormalizedBigWigs/Run135toRun150/NonK27genes_signal_matrix.npz \
      --outRawCounts ${OUTDIR}/NormalizedBigWigs/Run135toRun150/NormalizedBigWigs_NonK27genesonly_signal_matrix.tab

      multiBigwigSummary BED-file \
        --bwfiles ${OUTDIR}/NormalizedBigWigs/Run135toRun150/*.bw \
        --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
        -out ${OUTDIR}/NormalizedBigWigs/Run135toRun150/K27regions_signal_matrix.npz \
        --outRawCounts ${OUTDIR}/NormalizedBigWigs/Run135toRun150/NormalizedBigWigs_H3K27me3regions_signal_matrix.tab
