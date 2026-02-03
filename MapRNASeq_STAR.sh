#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=8:00:00
#SBATCH --output=../MappingOutput/logs/%x.out
#SBATCH --error=../MappingOutput/logs/%x.err

cd "$SLURM_SUBMIT_DIR"

# Use the CPUs you requested
THREADS=12

fastqPath="/lustre2/scratch/ry00555/RNASeqBamCoverage/Eaf3/FASTQ"
outdir="/lustre2/scratch/ry00555/RNASeqBamCoverage/Eaf3"

# (Optional) load modules once
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load STAR/2.7.11b-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load Subread/2.0.6-GCC-12.3.0
module load deepTools/3.5.5-gfbf-2023a

while read -r accession; do
  [ -z "$accession" ] && continue
  echo "Processing $accession"

  # Input FASTQs (your naming matches this perfectly)
  read1="${fastqPath}/${accession}_1.fastq.gz"
  read2="${fastqPath}/${accession}_2.fastq.gz"

  # Quick sanity check
  if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
    echo "Missing FASTQ(s) for $accession:"
    echo "  $read1"
    echo "  $read2"
    continue
  fi

  # Output directories (per accession)
  trimmed="${outdir}/TrimmedFastQs/${accession}"
  bamdir="${outdir}/bamFiles/${accession}"
  countsdir="${outdir}/counts/${accession}"
  bwDir="${outdir}/bigWig/${accession}"

  mkdir -p "$trimmed" "$bamdir" "$countsdir" "$bwDir"

  # Output file prefixes
  bam_prefix="${bamdir}/${accession}_"
  bam_out="${bam_prefix}Aligned.sortedByCoord.out.bam"
  counts_out="${countsdir}/${accession}_counts.txt"
  bw_out="${bwDir}/${accession}.bw"

  # Trim
  trim_galore --illumina --fastqc --paired --length 25 \
    --basename "${accession}" --gzip -o "$trimmed" \
    "$read1" "$read2"

  # Map with STAR
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

  # Index BAM
  samtools index -@ "$THREADS" "$bam_out"

  # Quantify
  featureCounts -T "$THREADS" \
    -t CDS -g gene_name -s 2 -p --primary \
    -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
    -o "$counts_out" \
    "$bam_out"

  # BigWig
  bamCoverage -p "$THREADS" -bs 50 --normalizeUsing BPM -of bigwig \
    -b "$bam_out" -o "$bw_out"

done < /lustre2/scratch/ry00555/RNASeqBamCoverage/Eaf3/FASTQ/accessions.txt
