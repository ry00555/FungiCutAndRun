#!/bin/bash
#SBATCH --job-name=ChIPSeqPipeline
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../ChIPSeqPipeline132.%j.out
#SBATCH --error=../ChIPSeqPipeline132.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

#source config.txt

OUTDIR=/scratch/ry00555/OutputRun132
ml SAMtools/1.16.1-GCC-10.2.0
ml BWA/0.7.17-GCC-10.2.0
ml Homer/4.11-foss-2020b
ml deepTools/3.5.1-foss-2020b-Python-3.8.6

# Process reads using trimGalore
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
#trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

#FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz"


#mkdir "${OUTDIR}/SortedBamFiles"
#mkdir "${OUTDIR}/BigWigs"
#mkdir "${OUTDIR}/Peaks"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"
#mkdir "${OUTDIR}/NormalizedBigWigs"

# Iterate over the files
#do
# file=${f##*/}
#  name=${file/%_S[1-99]*_R1_001_val_1.fq.gz/}

#  read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
#  bam="${OUTDIR}/SortedBamFiles/${name}.bam"
#  bigwig="${OUTDIR}/BigWigs/${name}"



  #bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
  #samtools index "$bam"

#  ml deepTools/3.3.1-intel-2019b-Python-3.7.4
#  conda install -c bioconda deeptools

#  bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#  bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
#done

# Set common variables
PEAKDIR="${OUTDIR}/Peaks"
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"

# Iterate over each BAM file in the directory
#for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
#  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
#  sample_id="${sample_id%%_Rep_1*}"

  # Make tag directory
#  makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"

  # Call peaks
#  findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${TAGDIR}/${sample_id}_peaks.txt"

  # Normalize to mitochondrial DNA which has no nucleosomes
  # Calculate read counts using samtools idxstats
  # mt_read_count=$(samtools idxstats "${bam_file}" | awk '$1=="MT"{print $3}')
  # reference_read_count=$(samtools idxstats "${bam_file}" | awk 'BEGIN{total=0}{if($1!="MT"){total+=$3}}END{print total}')

   # Calculate scaling factor
   #scaling_factor=$(bc <<< "scale=4; ${mt_read_count} / ${reference_read_count}")

   # Normalize the ChIP-seq signal using bamCoverage with the scaling factor
#   bamCoverage --scaleFactor $scaling_factor -of bigwig -b "${bam_file}".bam -o "${OUTDIR}/NormalizedBigWigs/${sample_id}_normalized.bw"


   #This is for unnormalized to mtDNA
#   bamCoverage -b "${bam_file}" -o "${OUTDIR}/BigWigs/${sample_id}.bw"
 #done



# Iterate over each file in the directory
for file_path in "${OUTDIR}/BigWigs"/*.bw; do
  # Get the base name of the file
  BW_name=$(basename "${file_path}")

  # Remove the file extension to create the sample ID
  BW_id="${BW_name%.*}"

  # Compute matrix
  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${OUTDIR}/BigWigs/${file_path}"*.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_${BW_id}.gz"

#For normalized to mtDNA note as of current there's only one variable for the bw file_path as we don't know if normalization will work yet
  #computeMatrix reference-point --referencePoint TSS -S "${file_path}" -R "/scratch/ry00555/neurospora.bed" -a 1500 -b 1500 --skipZeros -o "${OUTDIR}/Matrices/matrix_normalized_${BW_id}.gz"
  #plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_normalized_${BW_id}.gz" --outFileName "${BW_id}_normalized_hclust.png" \
  #            --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds

  # Plot heatmap
plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}.gz" --outFileName "${BW_id}_hclust.png" \
            --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds
done
