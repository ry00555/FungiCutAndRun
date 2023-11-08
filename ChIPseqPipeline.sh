#!/bin/bash
#SBATCH --job-name=ChIPSeqPipeline
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=../ChIPSeqPipeline131.%j.out
#SBATCH --error=../ChIPSeqPipeline131.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Set output directory specific for each sequencing experiment
OUTDIR=/scratch/ry00555/OutputRun131


#Load all the modules that are needed for the entire pipeline
ml SAMtools/1.16.1-GCC-10.2.0
ml BWA/0.7.17-GCC-10.2.0
ml Homer/4.11-foss-2020b
ml deepTools/3.5.1-foss-2020b-Python-3.8.6
# ml Perl/5.30.0-GCCcore-8.3.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
ml BEDTools/2.30.0-GCC-10.2.0

#This can be changed but make sure that there are appropriate paths and output directories based on the analysis specifically when peak calling and normalizing

#mkdir -p "${OUTDIR}/SortedBamFiles"
#mkdir -p "${OUTDIR}/BigWigs"
#mkdir -p "${OUTDIR}/Peaks"
#mkdir -p "$OUTDIR/HomerTagDirectories"
#mkdir -p "$OUTDIR/TdfFiles"
#mkdir -p "${OUTDIR}/Heatmaps"
#mkdir -p "${OUTDIR}/Matrices"
#mkdir -p "${OUTDIR}/NormalizedBigWigs"
#mkdir -p "${OUTDIR}/Beds"
mkdir -p "${OUTDIR}/Counts"


# Process reads using trimGalore
#trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz"

# Iterate over the files
#for f in $FILES
#do
 #file=${f##*/}
  #name=${file/%_S[1-99]*_R1_001_val_1.fq.gz/}

  #read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
#  bam="${OUTDIR}/SortedBamFiles/${name}.bam"
  #bigwig="${OUTDIR}/BigWigs/${name}"

#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
  #samtools index "$bam"
  #bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#  bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
#done

# Set common variables - NOTE this is my organization method, do what works best for you, but if you change these variables then you might need to change the paths and variables below
#PEAKDIR="${OUTDIR}/Peaks"
#TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
#BEDDIR="${OUTDIR}/Beds"

#Create mitochondrial normalization source files

#genomeFile="/scratch/ry00555/GCA_000182925.2_NC12_genomic.fna"
#samtools faidx $genomeFile
#Genome=$(basename $genomeFile .fna)

#create txt file with Chromosome Name and Length
#cut -f1,2 GCA_000182925.2_NC12_genomic.fna.fai > GCA_000182925.2_NC12_genomic.fna.sizes
#cut -f1,2 ${Genome}.fna.fai > ${Genome}.sizes


#Create bed file with Chromosome Name, Start Coordinate, End Coordinate
#perl -p -i -e 's/([A-Z]{2}[0-9]{6}\.1)\t(\d*)\t(\d*)\t(\d*)\t(\d*)/$1\t1\t$2/i' GCA_000182925.2_NC12_genomic.fna.fai

#convert fai file to chromsomal bed file
#cat GCA_000182925.2_NC12_genomic.fna.fai  > GCA_000182925.2_NC12_genomic.bed
#cut -f1,2 ${Genome}.fna.fai > ${Genome}.sizes

#create genomic file with desired window sizes
#windowSize=1000
#bedtools makewindows -g ${Genome}.sizes -w ${windowSize} > ${Genome}_1kbWindows.bed
#bedtools makewindows -g GCA_000182925.2_NC12_genomic.fna.sizes -w 1000 > GCA_000182925.2_NC12_genomic_1kbWindows.bed


# Iterate over each BAM file in the directory
for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
  sample_id="${sample_id%%_Rep_1*}"
 outputFile="${BEDDIR}/${sample_id}_normalized.bed"
  unnormalizedBigWig="${OUTDIR}/BigWigs/${sample_id}.bw"
  normalizedBigWig="${OUTDIR}/NormalizedBigWigs/${sample_id}_normalized.bw"


  # Make tag directory
 makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"

  # Call peaks

 findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${TAGDIR}/${sample_id}_peaks.txt"

 # Calculate coverage for 1kb windows
bedtools coverage -a "/scratch/ry00555/${Genome}_1kbWindows.bed" -b "${bam_file}" > "${BEDDIR}/${sample_id}_coverage.bed"

  # Check if the sample is an Input sample (for mitochondrial normalization)
  if [[ "${sample_id}" == *"Input"* ]]; then
    # Calculate the median coverage for mitochondrial chromosome from input file against the mitochondrial genome
    bedtools map -a "/scratch/ry00555/${Genome}_1kbWindows.bed" -b "${BEDDIR}/${sample_id}_coverage.bed" -c 4 -o median > "${BEDDIR}/${sample_id}_coverage.medianChromosomes.out"

    # Copy the unnormalized bigwig file (Input samples remain unnormalized)
    cp "${bam_file}" "${unnormalizedBigWig}"
  else
    # Calculate the median coverage for mitochondrial chromosome from ChIP-seq sample against the mitochondrial genome
    bedtools map -a "/scratch/ry00555/${Genome}_1kbWindows.bed" -b "${BEDDIR}/${sample_id}_coverage.bed" -c 4 -o median > "${BEDDIR}/${sample_id}_coverage.medianChromosomes.out"

    # Calculate read counts using samtools idxstats for the mitochondrial genome

  mt_read_count=$(samtools idxstats "${bam_file}" | awk -v mito_chr="KI" '$1 ~ mito_chr { total += $3 } END { print total }')

    # Calculate read counts using samtools idxstats for the non-mitochondrial genome
    reference_read_count=$(samtools idxstats "${bam_file}" | awk '$1 !~ /KI/ { total += $3 } END { print total }')

    # Append the read counts to a file
echo "${sample_id},${mt_read_count},${reference_read_count}" >> "${OUTDIR}/Counts/read_counts.csv"
fi
    Check if the read counts are valid (non-empty and numeric)
   if ! [[ "$mt_read_count" =~ ^[0-9]+$ && "$reference_read_count" =~ ^[0-9]+$ ]]; then
     echo "Error: Invalid read counts for ${bam_file}"
     continue
   fi

    Check if the read counts are greater than 0 to avoid division by zero
   if (( mt_read_count > 0 && reference_read_count > 0 )); then
      Calculate the scaling factor
     scaling_factor=$(awk "BEGIN {printf \"%.4f\", ${mt_read_count} / ${reference_read_count}}")

      Normalize the ChIP-seq signal using bamCoverage with the scaling factor and create normalized bigwig and bedgraph files
     bamCoverage --scaleFactor "${scaling_factor}" -of bedgraph -b "${bam_file}" -o "${OUTDIR}/NormalizedBigWigs/${sample_id}_normalized.bedgraph"
     bamCoverage --scaleFactor "${scaling_factor}" -of bigwig -b "${bam_file}" -o "${normalizedBigWig}"
   else
     echo "Error: Invalid read counts for ${bam_file}"
     continue
   fi

  Create unnormalized bigwig file for mtDNA
 bamCoverage -of bigwig -b "${bam_file}" -o "${unnormalizedBigWig}"

done


# Iterate over each file in the directory
for file_path in "${OUTDIR}/NormalizedBigWigs"/*.bw; do
  # Get the base name of the file
  BW_name=$(basename "${file_path}" .bw)

  # Remove the file extension to create the sample ID
  BW_id="${BW_name%.*}"
  # Replace special characters with underscores in the sample ID
 BW_id=${BW_id//[^a-zA-Z0-9]/_}

 # Limit the length of the sample ID to avoid long filenames
 BW_id=${BW_id:0:50}
 # Compute matrix for the reference-point TSS
 computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_${BW_id}.gz"

  # Compute matrix for the reference-point TSS with specific regions (e.g., PRC2 genes)
  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz"

  # Preprocess the matrix files to replace nan values with zeros
  zcat "${OUTDIR}/Matrices/matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/matrix_${BW_id}_processed.gz"
  zcat "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed.gz"

  # Plot heatmaps for the reference-point TSS
 plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds

  # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
 plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds
done
