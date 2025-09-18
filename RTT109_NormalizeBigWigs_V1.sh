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

PAIRFILE="${OUTDIR}/RTT109_meta_biwigs.txt"
BW_DIR="${OUTDIR}/BigWigs"
OUT_NORM="${OUTDIR}/BigWigs/NormalizedBigWigs"

ml deepTools

mkdir -p "$OUT_NORM"

echo "🔍 Checking IDs listed in: $PAIRFILE"
echo "   against files in: $BW_DIR"
echo

missing_pairs=()
processed_ids=()

declare -A chip_reps
declare -A input_map

# Read file: ID, BigWig, Input
tail -n +2 "$PAIRFILE" | while IFS=$'\t' read -r ID BigWig Input; do
    chip_path="${BW_DIR}/${BigWig}"
    input_path="${BW_DIR}/${Input}"

    # Track missing ChIP files
    if [[ ! -f "$chip_path" ]]; then
        echo "⚠️ Missing ChIP file: $chip_path"
        missing_pairs+=("$BigWig")
        continue
    fi

    chip_reps[$ID]+=" $chip_path"
    if [[ -n "$Input" && -f "$input_path" ]]; then
        input_map[$ID]="$input_path"
    fi
done

# Now process each ID
for ID in "${!chip_reps[@]}"; do
    chips="${chip_reps[$ID]}"
    input="${input_map[$ID]}"
    outbase="${OUT_NORM}/${ID}"

    # If multiple replicates: average them
    if [[ $(echo "$chips" | wc -w) -gt 1 ]]; then
        echo "👉 Averaging replicates for $ID"
        bigwigMerge $chips "${outbase}_merged.bedGraph"
        awk '{print $1, $2, $3, $4/NR}' "${outbase}_merged.bedGraph" > "${outbase}_avg.bedGraph"
        bedGraphToBigWig "${outbase}_avg.bedGraph" chrom.sizes "${outbase}_ChIP_avg.bw"
        chip_final="${outbase}_ChIP_avg.bw"
        rm "${outbase}_merged.bedGraph" "${outbase}_avg.bedGraph"
    else
        chip_final="$chips"
    fi

    # Normalize against Input if present
    if [[ -n "$input" ]]; then
        echo "Normalizing $ID: $(basename "$chip_final") vs $(basename "$input")"
        bigwigCompare \
          -b1 "$chip_final" \
          -b2 "$input" \
          --operation ratio \
          --pseudocount 1 \
          --binSize 25 \
          -o "${outbase}_foldchange.bw" \
          --skipZeroOverZero \
          --verbose
    else
        echo "⚠️ No Input for $ID → keeping ChIP only"
        cp "$chip_final" "${outbase}_noInput.bw"
    fi

    processed_ids+=("$ID")
done

echo
echo "========== PROCESSING SUMMARY =========="
echo "✅ Processed IDs: ${#processed_ids[@]}"
for p in "${processed_ids[@]}"; do
    echo "   $p"
done

echo
echo "⚠️ Missing ChIP files: ${#missing_pairs[@]}"
for m in "${missing_pairs[@]}"; do
    echo "   $m"
done


multiBigwigSummary BED-file \
  --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K27me3*.bw \
  --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
  -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_K27genes_signal_matrix.npz \
  --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_K27genesonly_signal_matrix.tab

  multiBigwigSummary BED-file \
    --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
    --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
    -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3ChIP_K27genes_signal_matrix.npz \
    --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigsH_3K36me3ChIP__K27genesonly_signal_matrix.tab

multiBigwigSummary BED-file \
      --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
      --BED "/scratch/ry00555/GeneList_BedFiles/NonK27genes.bed" \
      -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NonK27genes_signal_matrix.npz \
      --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_NonK27genesonly_signal_matrix.tab

      multiBigwigSummary BED-file \
        --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K27me3*.bw \
        --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
        -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_K27regions_signal_matrix.npz \
        --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_NormalizedBigWigs_H3K27me3regions_signal_matrix.tab


        multiBigwigSummary BED-file \
          --bwfiles ${OUTDIR}/BigWigs/NormalizedBigWigs/*H3K36me3*.bw \
          --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
          -out ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3_K27regions_signal_matrix.npz \
          --outRawCounts ${OUTDIR}/BigWigs/NormalizedBigWigs/RTT109_H3K36me3_NormalizedBigWigs_H3K27me3regions_signal_matrix.tab
