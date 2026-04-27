#!/bin/bash
#SBATCH --job-name=NCBI_FASTQ_pull
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24gb
#SBATCH --time=36:00:00
#SBATCH --output=../NCBI_FASTQ_pull.%j.out
#SBATCH --error=../NCBI_FASTQ_pull.%j.err
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --mail-type=ALL

cd "$SLURM_SUBMIT_DIR"

OUTDIR="/home/ry00555/RTT109Paper/GeoAccessions/RNASeq"
SRA_DIR="${OUTDIR}/SRA"
FASTQ_DIR="${OUTDIR}/FASTQ"
FAILED_LOG="${OUTDIR}/failed_srr_ids.txt"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

mkdir -p "$SRA_DIR" "$FASTQ_DIR"
: > "$FAILED_LOG"

ml purge
ml SRA-Toolkit
ml pigz 2>/dev/null || true

readarray -t SRR_IDS <<'EOF'
SRR8269709
SRR8269844
SRR8269736
SRR8269785
SRR8269838
SRR8269750
SRR8269689
SRR8269850
SRR8269692
SRR8269834
SRR8269815
SRR8269722
SRR8269775
SRR8269810
SRR7970564
SRR7970565
SRR7970613
SRR7970719
SRR7970723
SRR7970724
SRR7970603
SRR7970606
SRR7970610
SRR7970614
SRR7970615
SRR7970706
SRR9027639
SRR9027640
SRR9027641
SRR10916318
SRR10916319
SRR10916320
SRR10956838
SRR10956839
SRR10956840
SRR10916324
SRR10916325
SRR10916326
SRR10916163
SRR10916164
SRR10916165
SRR10958908
SRR10958909
SRR10958910
SRR9027682
SRR9027762
SRR9027766
EOF

for SRR_ID in "${SRR_IDS[@]}"; do
  echo "========== ${SRR_ID} =========="

  if [[ -s "${FASTQ_DIR}/${SRR_ID}_1.fastq.gz" || -s "${FASTQ_DIR}/${SRR_ID}.fastq.gz" ]]; then
    echo "FASTQ already exists for ${SRR_ID}, skipping."
    continue
  fi

  if ! prefetch -O "$SRA_DIR" "$SRR_ID"; then
    echo -e "${SRR_ID}\tprefetch_failed" >> "$FAILED_LOG"
    continue
  fi

  SRA_FILE="${SRA_DIR}/${SRR_ID}/${SRR_ID}.sra"
  if [[ ! -s "$SRA_FILE" ]]; then
    echo -e "${SRR_ID}\tmissing_sra_file" >> "$FAILED_LOG"
    continue
  fi

  if command -v fasterq-dump >/dev/null 2>&1; then
    if ! fasterq-dump "$SRA_FILE" --split-files --threads "$THREADS" --outdir "$FASTQ_DIR"; then
      echo -e "${SRR_ID}\tfasterq_dump_failed" >> "$FAILED_LOG"
      continue
    fi

    shopt -s nullglob
    fq_files=( "$FASTQ_DIR/${SRR_ID}"*.fastq )
    if ((${#fq_files[@]} == 0)); then
      echo -e "${SRR_ID}\tno_fastq_created" >> "$FAILED_LOG"
      shopt -u nullglob
      continue
    fi

    if command -v pigz >/dev/null 2>&1; then
      pigz -p "$THREADS" "${fq_files[@]}"
    else
      gzip "${fq_files[@]}"
    fi
    shopt -u nullglob
  else
    if ! fastq-dump --split-files --gzip "$SRA_FILE" -O "$FASTQ_DIR"; then
      echo -e "${SRR_ID}\tfastq_dump_failed" >> "$FAILED_LOG"
      continue
    fi
  fi
done

echo "Done. FASTQs in: $FASTQ_DIR"
if [[ -s "$FAILED_LOG" ]]; then
  echo "Some SRRs failed. Check: $FAILED_LOG"
else
  echo "All SRRs completed successfully."
fi
