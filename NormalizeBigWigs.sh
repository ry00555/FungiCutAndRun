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

source config.txt
OUTDIR="/scratch/ry00555/RNASeqPaper/"

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

mkdir -p "${OUTDIR}/NormalizedBigWigs"


ml deepTools
# Loop through all H3K27me3 bigwigs
# Loop through all H3K27me3 BigWigs in the input directory
for chip_bw in ${OUTDIR}/CandidateBigWigs/*H3K27me3*.bw; do
  # Attempt to infer the matching Input BigWig
  base=$(basename "$chip_bw" | sed 's/H3K27me3.*/Input/')
  input_bw=$(find ${OUTDIR}/CandidateBigWigs -name "*${base}*.bw" | head -n 1)

  if [[ -f "$input_bw" ]]; then
    outname=$(basename "$chip_bw" .bw)_norm_foldchange.bw
    echo "Normalizing $chip_bw against $input_bw → $outname"

    bigwigCompare \
      -b1 "$chip_bw" \
      -b2 "$input_bw" \
      --operation ratio \
      --pseudocount 1 \
      --smoothLength 150 \
      --binSize 25 \
      -o "${OUTDIR}/NormalizedBigWigs/$outname" \
      --skipZeroOverZero \
      --verbose
  else
    echo "⚠️ No matching input found for $chip_bw" >&2
  fi
done
