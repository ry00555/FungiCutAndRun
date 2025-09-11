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

PAIRFILE="${OUTDIR}/chip_input_pairs_Run135_Run150_H3K27me3only.txt"

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
          -o "${OUTDIR}/NormalizedBigWigs/Run135toRun150/$outname" \
          --skipZeroOverZero \
          --verbose
    else
        echo "⚠️ Missing file(s): $chip_path or $input_path" >&2
    fi
done < "$PAIRFILE"

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
