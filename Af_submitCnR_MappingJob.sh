#!/bin/bash
# Submits one CUT&RUN mapping job (MapCutRun.sh) per sample in the list.
# Usage (from the folder with MapCutRun.sh and config.txt):
#   bash submitCutRunJob.sh sampleList_CutRun.txt

if [[ -z "$1" ]]; then
  echo "
  ################
  ERROR: You must include the path to your sample list.
  eg. bash submitCutRunJob.sh sampleList_CutRun.txt
  ################
  "
  exit 1
fi

[[ -f MapCutRun.sh && -f config.txt ]] || { echo "ERROR: run this from the folder containing Map_AfCutnRun.sh and config.txt"; exit 1; }

fastqPath="/scratch/ry00555/Afum/CEA10/FastQ"   # <name>_R1.fastq.gz or <name>_R1_001.fastq.gz
outdir="/scratch/ry00555/Afum/CEA10"            # results go to $outdir/CutRun

mkdir -p "${outdir}/CutNRun/logs"

while read -r line || [[ -n "$line" ]]; do
  line="${line//$'\r'/}"            # strip Windows/Excel line endings
  [[ -z "$line" || "$line" == \#* ]] && continue

  if [[ ! -s "${fastqPath}/${line}_R1.fastq.gz" && ! -s "${fastqPath}/${line}_R1_001.fastq.gz" ]]; then
    echo "SKIPPED $line: no R1 file in ${fastqPath}"
    continue
  fi

  sbatch \
    --job-name="CR_${line}" \
    --output="${outdir}/CutNRun/logs/%x.out" \
    --error="${outdir}/CutNRun/logs/%x.err" \
    --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" \
    Map_AfCutnRun.sh
  echo "$line CUT&RUN job submitted"
  sleep 1
done < "$1"
