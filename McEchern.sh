#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=GATK.%j.out
#SBATCH --error=GATK.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make output directory
OUTDIR="/scratch/ry00555/McEachern/"


# mkdir TrimmedReads
# mkdir BigWigs
# mkdir SortedBamFiles


#Load these modules that are compatible with GATK version 4.3
#module load BEDOPS/2.4.39-foss-2019b

#ml R/3.6.2-foss-2019b
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > klactis.fasta
#
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > $WorkDir/Genome/klactis.fasta
#
# module load SAMtools
# samtools faidx $WorkDir/Genome/klactis.fasta
# samtools faidx klactis.fasta
#
# module load Bowtie2
#
# faidx --transform bed klactis.fasta > klactis.bed
#
 #ml picard
#
# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
#       -R GCF_000002515.2_ASM251v1_genomic.fna \
#       -O GCF_000002515.2_ASM251v1_genomic.dict
#
#       java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
#       -I /scratch/ry00555/McEachern/Genome/Klactisgenes.bed \
#       -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#       -SD /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
#       -O /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.interval_list
#
#ml GATK
#
# #WGS uses 1000 bp bins
#       gatk PreprocessIntervals \
#       -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#       -L /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.interval_list \
#       --interval-merging-rule OVERLAPPING_ONLY \
#       --bin-length 1000 \
#       --padding 0 \
#       -O /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list
#
#       gatk AnnotateIntervals \
#       -R GCF_000002515.2_ASM251v1_genomic.fna \
#       -L klactis_preprocessed1000_intervals.interval_list \
#        --interval-merging-rule OVERLAPPING_ONLY \
#        -O GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv

#4/15/24 Potentially M1-7 are Kluveromyces marxianus so prep genomes just like klactis Reference: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001417885.1/
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > klactis.fasta
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_genomic.fna.gz | gunzip -c > Kluyveromycesmarxianus.fna
# module load SAMtools
# samtools faidx Kluyveromycesmarxianus.fna
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_genomic.gtf.gz | gunzip -c > Kluyveromycesmarxianus.gtf
# module load BEDOPS/2.4.39-foss-2019b
# gtf2bed < Kluyveromycesmarxianus.gtf > Kluyveromycesmarxianus.bed --max-mem 10G
# awk '$8 == "gene"' Kluyveromycesmarxianus.bed > Kluyveromycesmarxianus_genes.bed

# ml picard
# #
# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
# -R Kluyveromycesmarxianus.fna \
# -O Kluyveromycesmarxianus.dict
# #
# java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
# -I Kluyveromycesmarxianus_genes.bed \
# -R Kluyveromycesmarxianus.fna \
# -SD Kluyveromycesmarxianus.dict \
# -O Kluyveromycesmarxianus.interval_list
# #
# ml GATK
# #
# #WGS uses 1000 bp bins
# gatk PreprocessIntervals \
# -R Kluyveromycesmarxianus.fna \
# -L Kluyveromycesmarxianus.interval_list \
# --interval-merging-rule OVERLAPPING_ONLY \
# --bin-length 1000 \
# --padding 0 \
# -O Kluyveromycesmarxianus_preprocessed1000_intervals.interval_list
##########################################
#Now handle the samples
# ml Trim_Galore
# trim_galore --paired --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/113*fastq\.gz
#        # #
#      #
#       # mkdir "${OUTDIR}/SortedBamFiles"
#       # mkdir "${OUTDIR}/BigWigs"
#       # mkdir "${OUTDIR}/Peaks"
#      #mkdir "$OUTDIR/HomerTagDirectories"
#      #mkdir "$OUTDIR/TdfFiles"
#
#
#  FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1.fq.gz" # Don't forget the *
# #
# # # Iterate over the files
#  for f in $FILES
#  do
#      file=${f##*/}
#      name=${file/%_S[1-12]*R1_001_val_1.fq.gz/}
#      read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R1_001_val_2\.fq\.gz/g')
#      bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
# #bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
# ml SAMtools
# ml BWA
#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
# ml deepTools
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#  done

##On April 15th, ZL told me to rerun the analuysis on all of the Run 113 gDNA files and they were all single end so merge and cat then trim
#command line
#  cat 113-1-* > 113-1-gDNA-CBS2359_merged.fastq.gz for all 113-22 samples

#now trim the single end reads and make bam files and place them in the TrimmedReads directory. Use Run113 as the fastq directory so it doesn't redo 1 and 12

#ml Trim_Galore
#trim_galore  --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads 113*merged.fastq.gz
#trim_galore  --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/Run113/113*fastq\.gz
FILES="${OUTDIR}/TrimmedReads/113*merged.fq.gz" # Don't forget the *

#need to merge the fastq files after trimming	#need to merge the fastq files before trimming so do this cat 113-1* > ../113-1-gDNA-CBS2359_merged.fq.gz in Run113 directory don't do it again for samples 1 and 12

for f in $FILES
do
file=${f##*/}
name=${file/%_S[1-12]*merged.fq.gz/}
bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
ml SAMtools
ml BWA
bwa mem -M -v 3 -t $THREADS $GENOME $f | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
ml deepTools
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
done

# # Set the directory containing the sorted BAM files
# SORTED_BAM_DIR="/scratch/ry00555/McEachern/SortedBamFiles/"
# #
# # # Iterate over all BAM files in the directory
# ml picard
# for bam_file in $SORTED_BAM_DIR/113*.bam
# do
# #     # Get the base name of the BAM file
#      base_name=$(basename "$bam_file" .bam)
# #
# #     # Define the output file path
#      output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
# #
# #     # Run Picard to add or replace read groups
#      java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#      -I "$bam_file" \
#      -O "$output_file" \
#      -RGID 2 \
#      -RGLB lib1 \
#      -RGPL illumina \
#      -RGPU S34 \
#      -RGSM "${base_name%.*}"
#  done
# #
# # #mkdir CountTSVs
#  ml GATK
#   for bam_file in $SORTED_BAM_DIR/113*_output.bam
#   do
# # #   # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" _output.bam)
# # #   # Define the output file path
#     input_file="${SORTED_BAM_DIR}/${base_name}"
#  samtools index "$input_file"
# # #
#   gatk CollectReadCounts \
#   -I "$input_file" \
#   -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#   -L /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list \
#   --interval-merging-rule OVERLAPPING_ONLY \
#   -O /scratch/ry00555/McEachern/CountTSVs/$base_name.counts.tsv
#
#  done
#
# CountTSVsDIR="/scratch/ry00555/McEachern/CountTSVs/"
#
# # gatk CreateReadCountPanelOfNormals \
# # -I ${CountTSVsDIR}/113-1-gDNA-CBS2359_merged.counts.tsv \
# # -I ${CountTSVsDIR}/113-12-gDNA-7B520_merged.counts.tsv  \
# # --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# # -O ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5
# #
#  for count_files in $CountTSVsDIR/113*tsv
#  do
# #
# #   # Get the base name of the counts file
#       base_name=$(basename "$count_files" .counts.tsv)
# #  #   # Define the output file path
#    input_file="${CountTSVsDIR}/${base_name}"
#  gatk DenoiseReadCounts \
#  -I "$input_file" \
#  --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
#  --count-panel-of-normals ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5 \
#  --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
#  --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
# #
# done
# #
#
# ml R/3.6.2-foss-2019b
# ml GATK/4.3.0.0-GCCcore-8.3.0-Java-11
# for copy_ratios in ${OUTDIR}/CopyRatios/
# do
# # # Get the base name of the counts file
#     base_name=$(basename "$copy_ratios")
#  #   # Define the output file path
# # input_file="${OUTDIR}/CopyRatios/${base_name}"
# #
#  gatk PlotDenoisedCopyRatios \
#  --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
#  --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv  \
#  --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
#  --point-size-copy-ratio 1 \
#  --output-prefix ${base_name} \
#  --output ${OUTDIR}/PlotDenoisedCopyRatios
#  done
#
#  for copy_ratios in ${OUTDIR}/CopyRatios/*.denoisedCR.tsv
#  do
#
#  base_name=$(basename "$copy_ratios" .denoisedCR.tsv)
#
#  gatk ModelSegments \
#  --denoised-copy-ratios CopyRatios/${base_name}.denoisedCR.tsv \
#  --output-prefix ${base_name} \
#  -O ${OUTDIR}/ModelSegments
#  gatk PlotModeledSegments \
#  --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv \
#  --segments ModelSegments/${base_name}.modelFinal.seg \
#  --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
#         --point-size-copy-ratio 1 \
#         --output-prefix ${base_name} \
# -O ${OUTDIR}/PlotModelSegments
#  done
#
#
#  #Do the script for the Km samples
#  for bam_file in $SORTED_BAM_DIR/*_output.bam
# do
#     # Check if the file name contains M1-M7
#     if [[ "$bam_file" == *M1* || "$bam_file" == *M2* || "$bam_file" == *M3* || "$bam_file" == *M4* || "$bam_file" == *M5* || "$bam_file" == *M6* || "$bam_file" == *M7* ]]; then
#         # Get the base name of the BAM file
#         base_name=$(basename "$bam_file" _output.bam)
#         # Define the output file path
#         input_file="${SORTED_BAM_DIR}/${base_name}"
#         samtools index "$input_file"
#         gatk CollectReadCounts \
#             -I "$input_file" \
#             -R /scratch/ry00555/McEachern/Genome/Kluyveromycesmarxianus.fna \
#             -L /scratch/ry00555/McEachern/Genome/ Kluyveromycesmarxianus_preprocessed1000_intervals.interval_list \
#             --interval-merging-rule OVERLAPPING_ONLY \
#             -O /scratch/ry00555/McEachern/CountTSVs/$base_name.counts.tsv
#     fi
# done
