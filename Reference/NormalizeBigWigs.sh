#!/bin/bash
#SBATCH --job-name=NormalizeBigWigs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=../NormalizeBigWigs.%j.out
#SBATCH --error=../NormalizeBigWigs.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

OUTDIR="/scratch/ry00555/RNASeqPaper"

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

mkdir -p "${OUTDIR}/NormalizedBigWigs"

PAIRFILE="${OUTDIR}/chip_input_pairs.txt"

ml deepTools
while IFS=$'\t' read -r chip_bw input_bw; do
    chip_path="${OUTDIR}/CandidateBigWigs/$chip_bw"
    input_path="${OUTDIR}/CandidateBigWigs/$input_bw"

    if [[ -f "$chip_path" && -f "$input_path" ]]; then
        outname=$(basename "$chip_bw" .bw)_norm_foldchange.bw
        echo "Normalizing $chip_bw against $input_bw → $outname"

        bigwigCompare \
          -b1 "$chip_path" \
          -b2 "$input_path" \
          --operation ratio \
          --pseudocount 1 \
          --smoothLength 150 \
          --binSize 25 \
          -o "${OUTDIR}/NormalizedBigWigs/$outname" \
          --skipZeroOverZero \
          --verbose
    else
        echo "⚠️ Missing file(s): $chip_path or $input_path" >&2
    fi
done < "$PAIRFILE"

multiBigwigSummary BED-file \
  --bwfiles ${OUTDIR}/NormalizedBigWigs/*.bw \
  --BED "/scratch/ry00555/heatmapPRC2genes.bed" \
  -out ${OUTDIR}/region_signal_matrix.npz \
  --outRawCounts ${OUTDIR}/NormalizedBigWigs_K27genesonly_signal_matrix.tab

  multiBigwigSummary BED-file \
    --bwfiles *.bw \
    --BED "/scratch/ry00555/heatmapPRC2genes.bed" \
    -out L2FC_region_signal_matrix.npz \
    --outRawCounts NormalizedBigWigs_K27genesonly_signal_matrix.tab
