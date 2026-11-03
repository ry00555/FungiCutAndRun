#!/bin/bash
# Submits one MapData.sh job per sample, so all samples map at the same time.
# Run with bash (NOT sbatch), from the folder containing MapData.sh and config.txt:
#
#   bash submitMappingJob.sh                    # map EVERY sample in fastqPath
#   bash submitMappingJob.sh sampleList.txt     # map only the samples in the list
#
# A sample name is the FastQ file name without _R1_001.fastq.gz,
# e.g. 158-1_ChIP_WT_H3K27me3_Rep1_S1_L004

fastqPath="/scratch/ry00555/Run158/FASTQ"
outdir="/scratch/ry00555/Run158"

[[ -f MapData.sh && -f config.txt ]] || { echo "ERROR: run this from the folder containing MapData.sh and config.txt"; exit 1; }
[[ -d "$fastqPath" ]] || { echo "ERROR: FastQ folder not found: $fastqPath"; exit 1; }

# Same output folders as MapCutAndRun.sh
mkdir -p "${outdir}/logs" "${outdir}/TrimmedReads" "${outdir}/SortedBamFiles" \
         "${outdir}/BigWigs" "${outdir}/MACSPeaks"

# Sample list: from the file given, or every R1 file in fastqPath
if [[ -n "$1" ]]; then
  samples=$(tr -d '\r' < "$1" | grep -v '^#' | grep -v '^[[:space:]]*$')
else
  samples=$(ls "$fastqPath" | grep '_R1_001\.fastq\.gz$' | sed 's/_R1_001\.fastq\.gz$//')
fi
[[ -n "$samples" ]] || { echo "ERROR: no samples found (looking for *_R1_001.fastq.gz in $fastqPath)"; exit 1; }

n=0
for line in $samples; do
  if [[ ! -s "${fastqPath}/${line}_R1_001.fastq.gz" || ! -s "${fastqPath}/${line}_R2_001.fastq.gz" ]]; then
    echo "SKIPPED $line: R1 or R2 file missing in ${fastqPath}"
    continue
  fi

  sbatch \
    --job-name="${line}" \
    --output="${outdir}/logs/%x.out" \
    --error="${outdir}/logs/%x.err" \
    --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" \
    MapData.sh
  echo "$line mapping job submitted"
  n=$((n+1))
  sleep 1
done

echo "Submitted $n jobs. Check with: squeue -u $USER"
