#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR
source config.txt

# Variables passed in from submitMappingJob.sh:
#   accession = SampleID (e.g. CEA10_H3K4me2_ChIP_1_5glucose_Rep1)
#   fastqPath = folder holding <SampleID>_R1.fastq.gz and <SampleID>_R2.fastq.gz
#   outdir    = output folder

set -euo pipefail

################### input files ###################
# Exact names (no wildcard), so one sample can never pick up another's reads
read1="${fastqPath}/${accession}_R1.fastq.gz"
read2="${fastqPath}/${accession}_R2.fastq.gz"

for f in "$read1" "$read2"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

################### output folders ###################
trimmed="${outdir}/TrimmedFastQs/${accession}"
tmp="${outdir}/tempFile/${accession}"
bamdir="${outdir}/bamFiles"
bwDir="${outdir}/bigWig"
PeakDir="${outdir}/Peaks/${accession}"
mkdir -p "$trimmed" "$tmp" "$bamdir" "$bwDir" "$PeakDir"

bam="${bamdir}/${accession}.bam"
bigwig="${bwDir}/${accession}"

################### trimming ###################
# remove adaptors, trim low-quality bases (Phred $TRIMQ), keep reads >= 25 bp
module load Trim_Galore/0.6.10-GCCcore-12.3.0
trim_galore --illumina --fastqc --paired --length 25 --quality $TRIMQ --cores 2 \
  --basename "${accession}" --gzip -o "$trimmed" "$read1" "$read2"

################### mapping ###################
module load SAMtools/1.21-GCC-13.3.0
module load BWA/0.7.18-GCCcore-13.3.0

bwa mem -M -v 3 -t $THREADS "$GENOME" \
    "${trimmed}/${accession}_val_1.fq.gz" "${trimmed}/${accession}_val_2.fq.gz" \
  | samtools view -bhSu - \
  | samtools sort -@ $THREADS -T "${tmp}/${accession}" -o "$bam" -
samtools index "$bam"
rm -r "$tmp"

################### bigWigs ###################
module load deepTools/3.5.5-gfbf-2023a
# all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality $MAPQ \
  --smoothLength $SMOOTH -of bigwig -b "$bam" \
  -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

# mononucleosome-sized fragments (bin 1, smooth 25 -> named to match)
bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --minMappingQuality $MAPQ \
  --smoothLength 25 -of bigwig -b "$bam" \
  -o "${bigwig}.bin_1.smooth_25_MNase.bw"

################### peaks ###################
module load MACS3/3.0.1-gfbf-2023a
# --nolambda: call peaks without a control
macs3 callpeak -t "$bam" -f BAMPE -n "${accession}" --broad -g $GSIZE \
  --broad-cutoff 0.1 --outdir "$PeakDir" --min-length 800 --max-gap 500 --nolambda
