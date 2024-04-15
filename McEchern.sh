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

#ml R/3.6.2-foss-2019b
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > klactis.fasta
#
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > $WorkDir/Genome/klactis.fasta
#
# module load SAMtools
# samtools faidx $WorkDir/Genome/klactis.fasta
# samtools faidx klactis.fasta
#module load BEDOPS/2.4.39-foss-2019b
#gtf2bed < klactis.gtf > klactis.bed
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

#ml Trim_Galore
#trim_galore --paired --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/113*fastq\.gz
     # #
   #
    # mkdir "${OUTDIR}/SortedBamFiles"
    # mkdir "${OUTDIR}/BigWigs"
    # mkdir "${OUTDIR}/Peaks"
   #mkdir "$OUTDIR/HomerTagDirectories"
   #mkdir "$OUTDIR/TdfFiles"


# FILES="${OUTDIR}/TrimmedReads/113*_S11_L002_R1_001_trimmed.fq.gz" # Don't forget the *
#
# # Iterate over the files
#  for f in $FILES
#  do

# for Trimmedreads in ${OUTDIR}/TrimmedReads/113*_S11_L002_R1_001_trimmed.fq.gz
# do
#     # Get the base name of the BAM file
#     base_name=$(basename "$Trimmedreads" _S11_L002_R1_001_trimmed.fq.gz)
#Example 113-11-gDNA-CBS2Asm-Int3A_S11_L002_R1_001_trimmed.fq.gz


#      file=${f##*/}
#      name=${file/%_S[1-12]*R1_001_val_1.fq.gz/}
#      bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
#bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
# ml SAMtools
# ml BWA
#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
# ml deepTools
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#  done

#Command lines for the two 113 wildtype samples
# bam="/scratch/ry00555/McEachern/SortedBamFiles/113-1-gDNA-CBS2359_merged.bam"
# bigwig="/scratch/ry00555/McEachern/BigWigs/113-1-gDNA-CBS2359_merged"
# bam2="/scratch/ry00555/McEachern/SortedBamFiles/113-12-gDNA-7B520_merged.bam"
# bigwig2="/scratch/ry00555/McEachern/BigWigs/113-12-gDNA-7B520_merged"
#
#ml SAMtools
# ml BWA
#  bwa mem -M -v 3 -t $THREADS $GENOME ${OUTDIR}/TrimmedReads/113-1-gDNA-CBS2359_merged_trimmed.fq.gz | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
#  bwa mem -M -v 3 -t $THREADS $GENOME ${OUTDIR}/TrimmedReads/113-12-gDNA-7B520_merged_trimmed.fq.gz | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam2" -
#
# samtools index "$bam"
# samtools index "$bam2"
#
# ml deepTools
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam2" -o "${bigwig2}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"


#Redo all the Run 113 samples
#ml Trim_Galore
#trim_galore  --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/Run113/113*fastq\.gz
# FILES="${OUTDIR}/TrimmedReads/113*.fq.gz" # Don't forget the *
#need to merge the fastq files after trimming
#for

#Stop here 4/15/24 I processed the K. marxianus genome and i copied over all of the files
# Set the directory containing the sorted BAM files
#SORTED_BAM_DIR="/scratch/ry00555/McEachern/SortedBamFiles"

# Iterate over all BAM files in the directory
#  ml picard
# for bam_file in $SORTED_BAM_DIR/*.bam
# do
#     # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" .bam)
#
#     # Define the output file path
#     output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
#
#     # Run Picard to add or replace read groups
#     java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#     -I "$bam_file" \
#     -O "$output_file" \
#     -RGID 2 \
#     -RGLB lib1 \
#     -RGPL illumina \
#     -RGPU S34 \
#     -RGSM "${base_name%.*}"
# done
CountTSVsDIR="/scratch/ry00555/McEachern/CountTSVs"

#mkdir CountTSVs
#did this in command line
#   input_file="${SORTED_BAM_DIR}/*_output.bam"

#   ml SAMtools
#samtools index "$input_file"

#ml GATK
# for bam_file in ${SORTED_BAM_DIR}/*_output.bam; do
#     # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" _output.bam)
#     gatk CollectReadCounts \
#     -I "$bam_file" \
#     -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#     -L /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list \
#     --interval-merging-rule OVERLAPPING_ONLY \
#     -O "${CountTSVsDIR}/${base_name}.counts.tsv"
# done




# gatk CreateReadCountPanelOfNormals \
# -I ${CountTSVsDIR}/113-1-gDNA-CBS2359_merged.counts.tsv \
# -I ${CountTSVsDIR}/113-12-gDNA-7B520_merged.counts.tsv  \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# -O ${OUTDIR}/PanelofNormals/113_WT_Samples.pon.hdf5
#
# for count_files in $CountTSVsDIR/*.counts.tsv
# do
#
# #   # Get the base name of the counts file
# base_name=$(basename "$count_files" .counts.tsv)
# #  #   # Define the output file path
# gatk DenoiseReadCounts \
# -I "$count_files" \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# --count-panel-of-normals ${OUTDIR}/PanelofNormals/113_WT_Samples.pon.hdf5 \
# --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
# -denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
#
# done

ml R/3.6.2-foss-2019b
ml GATK/4.3.0.0-GCCcore-8.3.0-Java-11
# for copy_ratios in ${OUTDIR}/CopyRatios/*.standardizedCR.tsv
# do
#
# # # Get the base name of the counts file
#    base_name=$(basename "$copy_ratios" .standardizedCR.tsv)
# # #   # Define the output file path
# #
#  gatk PlotDenoisedCopyRatios \
#  --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
#  --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv  \
#  --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
#  --point-size-copy-ratio 1 \
#  --output-prefix ${base_name} \
#  --output ${OUTDIR}/PlotDenoisedCopyRatios
#  done
# #on local machine
#  scp -r ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/McEachern/PlotDenoisedCopyRatios /Users/rochelleyap/Desktop/McEchernBigWigs
for copy_ratios in ${OUTDIR}/CopyRatios/*.denoisedCR.tsv
do

base_name=$(basename "$copy_ratios" .denoisedCR.tsv)

gatk ModelSegments \
--denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv \
--output-prefix ${base_name} \
-O ${OUTDIR}/ModelSegments


gatk PlotModeledSegments \
--denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv \
--segments ${OUTDIR}/ModelSegments/${base_name}.modelFinal.seg \
--sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
       --point-size-copy-ratio 1 \
       --output-prefix ${base_name} \
       -O ${OUTDIR}/PlotModelSegments
done
### map the rest of Ailieen's samples from Run 113 copy all the files
