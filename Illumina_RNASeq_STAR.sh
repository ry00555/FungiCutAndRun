#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=36:00:00
#SBATCH --output=../%x.out
#SBATCH --error=../%x.err

cd "$SLURM_SUBMIT_DIR"

THREADS=12

fastqPath="/lustre2/scratch/ry00555/Run155/2026_Run155_FastQ/RNA_FASTQ"
outdir="/lustre2/scratch/ry00555/Run155"

mkdir -p \
  "${outdir}/TrimmedFastQs" \
  "${outdir}/SortedBamFiles" \
  "${outdir}/counts" \
  "${outdir}/BigWigs" \
  "../MappingOutput/logs"

module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load STAR/2.7.11b-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load Subread/2.0.6-GCC-12.3.0
module load deepTools/3.5.5-gfbf-2023a

for read1 in "${fastqPath}"/*_R1_001.fastq.gz; do

  read2="${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  accession=$(basename "$read1" _R1_001.fastq.gz)

  echo "Processing $accession"

  if [ ! -f "$read2" ]; then
    echo "Missing R2 for $accession:"
    echo "  Expected: $read2"
    continue
  fi

  trimmed="${outdir}/TrimmedFastQs/${accession}"
  bamdir="${outdir}/SortedBamFiles/${accession}"
  bwDir="${outdir}/BigWigs/${accession}"

  mkdir -p "$trimmed" "$bamdir" "$bwDir"

  bam_prefix="${bamdir}/${accession}_"
  bam_out="${bam_prefix}Aligned.sortedByCoord.out.bam"
  bw_out="${bwDir}/${accession}.bw"
trim_galore --illumina --fastqc --paired --length 25 --basename "$accession" --gzip -o "$trimmed" "$read1" "$read2"

  STAR --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
    --outFileNamePrefix "$bam_prefix" \
    --readFilesIn "$trimmed/${accession}_val_1.fq.gz" "$trimmed/${accession}_val_2.fq.gz" \
    --readFilesCommand zcat \
    --alignIntronMax 10000 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingBinsN 100 \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --limitBAMsortRAM 19990000000

  samtools index -@ "$THREADS" "$bam_out"

  bamCoverage \
    -p "$THREADS" \
    -bs 50 \
    --normalizeUsing BPM \
    -of bigwig \
    -b "$bam_out" \
    -o "$bw_out"

done

echo "Running combined featureCounts..."

featureCounts -T "$THREADS" \
  -t CDS \
  -g gene_name \
  -s 2 \
  -p \
  --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o "${outdir}/counts/Run155_gene_counts.txt" \
  "${outdir}"/SortedBamFiles/*/*Aligned.sortedByCoord.out.bam

echo "Done."
echo "Combined counts:"
echo "${outdir}/counts/Run155_gene_counts_clean_names.txt"
