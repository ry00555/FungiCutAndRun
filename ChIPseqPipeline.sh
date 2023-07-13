#!/bin/bash
#SBATCH --job-name=Run132ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun132.%j.out
#SBATCH --error=../MapCutAndRun132.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

OUTDIR=/scratch/ry00555/OutputRun132


# Process reads using trimGalore
ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz"

mkdir "${OUTDIR}/SortedBamFiles"
mkdir "${OUTDIR}/BigWigs"
mkdir "${OUTDIR}/Peaks"
mkdir "$OUTDIR/HomerTagDirectories"
mkdir "$OUTDIR/TdfFiles"

# Iterate over the files
for f in $FILES
do
  file=${f##*/}
  name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}

  read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
  bam="${OUTDIR}/SortedBamFiles/${name}.bam"
  bigwig="${OUTDIR}/BigWigs/${name}"

  ml SAMtools/1.9-GCC-8.3.0
  ml BWA/0.7.17-GCC-8.3.0

  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
  samtools index "$bam"

  ml deepTools/3.3.1-intel-2019b-Python-3.7.4

  bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
  bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
done

ml Homer/4.11-foss-2019b SAMtools/1.16.1-GCC-11.3.0

# Set common variables
OUTDIR="${OUTDIR}/Peaks"
TAGDIR="$OUTDIR/HomerTagDirectories"
BAMDIR="/${OUTDIR}/SortedBamFiles"

# Generate list of sample IDs based on a naming pattern
common_id="132"
start_num=17
end_num=95
samples=()
for ((num=start_num; num<=end_num; num++)); do
  sample_id="${common_id}-${num}"
  samples+=("$sample_id")
done

# Iterate over the sample IDs
for sample in "${samples[@]}"; do
  # Construct file paths
  bam="${BAMDIR}/${sample}_S${sample}_L001_R1_001_val_1.fq.gz.bam"
  input_bam="${BAMDIR}/${sample}_input_S${sample}_L001_R1_001_val_1.fq.gz.bam"

  # Make tag directory
  makeTagDirectory "${TAGDIR}/${sample}" "$bam"

  # Call peaks
  findPeaks "${TAGDIR}/${sample}" -style histone -region -size 150 -minDist 530 -o "${OUTDIR}/${sample}_peaks.txt" -i "$input_bam"
done

# Rest of the script...
#ml deepTools/3.5.1-intel-2020b-Python-3.8.6

# Directory containing the files
#file_dir="/path/to/files/directory"

# Iterate over each file in the directory
#for file_path in "${bigwig}"/*.bw; do
  # Get the base name of the file
#  file_name=$(basename "${file_path}")

  # Remove the file extension to create the sample ID
#  sample_id="${file_name%.*}"

  # Compute matrix
#  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R heatmapPRC2genes.bed --skipZeros -o "matrix_${sample_id}.gz"

  # Plot heatmap
#  plotHeatmap -m "matrix_${sample_id}.gz" -out "cacheatmap_H3K4me2_${sample_id}_hclust.png" --samplesLabel "${samples[*]}" --hclust 1 --colorMap Reds
#done
