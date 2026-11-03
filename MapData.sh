#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=../Run158/logs/%x.out
#SBATCH --error=../Run158/logs/%x.err

## Variables passed in from submitMappingJob.sh:
##   accession = FastQ name without _R1_001.fastq.gz (e.g. 158-1_ChIP_WT_H3K27me3_Rep1_S1_L004)
##   fastqPath = /scratch/ry00555/Run158/FASTQ
##   outdir    = /scratch/ry00555/Run158
## From config.txt: THREADS, GENOME, BIN, SMOOTH.
## (config.txts FASTQ and OutputFolderName lines are NOT used; the paths come from submitMappingJob.sh.)
##
## Output (same folders and names as MapCutAndRun.sh):
##   TrimmedReads/<name>_val_1.fq.gz, _val_2.fq.gz + FastQC
##   SortedBamFiles/<name>.bam (+ .bai)
##   BigWigs/<name>.bin_<BIN>.smooth_<SMOOTH>Bulk.bw
##   MACSPeaks/<name>_peaks.broadPeak etc.
## where <name> is the sample name with the _S##_L00# ending removed.

cd $SLURM_SUBMIT_DIR
source config.txt
set -euo pipefail

################### input files ###################
# Exact file names (no wildcards), so one sample can never pick up another's reads
read1="${fastqPath}/${accession}_R1_001.fastq.gz"
read2="${fastqPath}/${accession}_R2_001.fastq.gz"
for f in "$read1" "$read2"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# Short name for outputs: drop the _S##_L00# ending (and any extra underscores before it)
name=$(echo "$accession" | sed -E 's/_+S[0-9]+_L[0-9]{3}$//')

################### output folders ###################
trimmed="${outdir}/TrimmedReads"
BAMDIR="${outdir}/SortedBamFiles"
BWDIR="${outdir}/BigWigs"
PEAKDIR="${outdir}/MACSPeaks"
tmp="${outdir}/tempFile/${name}"
mkdir -p "$trimmed" "$BAMDIR" "$BWDIR" "$PEAKDIR" "$tmp"

bam="${BAMDIR}/${name}.bam"
bigwig="${BWDIR}/${name}"

################### trimming ###################
# remove adaptors, trim low-quality bases (Phred 20), keep reads >= 25 bp
module load Trim_Galore/0.6.10-GCCcore-12.3.0
trim_galore --illumina --fastqc --paired --length 25 \
  --basename "${name}" --gzip -o "$trimmed" "$read1" "$read2"

################### mapping ###################
module load SAMtools/1.21-GCC-13.3.0
module load BWA/0.7.18-GCCcore-13.3.0

bwa mem -M -v 3 -t $THREADS "$GENOME" \
    "${trimmed}/${name}_val_1.fq.gz" "${trimmed}/${name}_val_2.fq.gz" \
  | samtools view -bhSu - \
  | samtools sort -@ $THREADS -T "${tmp}/${name}" -o "$bam" -
samtools index "$bam"
rm -r "$tmp"

################### bigWig ###################
module load deepTools/3.5.5-gfbf-2023a
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 20 \
  --smoothLength $SMOOTH -of bigwig -b "$bam" \
  -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

################### peaks ###################
module load MACS3/3.0.1-gfbf-2023a
# No control (same as MapCutAndRun.sh). -g 41037538 = Neurospora crassa genome size.
macs3 callpeak -t "$bam" -f BAMPE -n "${name}" --broad -g 41037538 \
  --broad-cutoff 0.1 --outdir "$PEAKDIR" --min-length 800 --max-gap 500

echo "${name} done"
