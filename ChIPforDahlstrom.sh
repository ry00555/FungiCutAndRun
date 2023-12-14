#!/bin/bash
#SBATCH --job-name=DahlstrohmChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../DahlstrohmChIP.%j.out
#SBATCH --error=../DahlstrohmChIP.%j.err

cd $SLURM_SUBMIT_DIR

#ml BWA/0.7.17-GCCcore-11.3.0
#ml SAMtools/0.1.20-GCC-11.3.0
#If incompatible modules use SAMtools/1.16.1-GCC-11.3.0 instead rec by Chelsea Taylor
ml deepTools/3.5.2-foss-2022a


ml Homer/4.11-foss-2022a
ml Perl/5.34.1-GCCcore-11.3.0
#ml Trim_Galore/0.6.7-GCCcore-11.2.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
#ml BEDTools/2.30.0-GCC-11.3.0
#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )


source ConfigforDahlstrom.txt

OUTDIR=/scratch/ry00555/DahlstromRun135/Output
#if [ ! -d $OUTDIR ]
#then
#mkdir -p $OUTDIR
#fi


# #process reads using trimGalore

 #trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
#
#FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
#mkdir "${OUTDIR}/SortedBamFiles"
#mkdir "${OUTDIR}/BigWigs"
#mkdir "${OUTDIR}/Peaks"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"
#
#Iterate over the files
#for f in $FILES
#do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-100]*_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	#read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	#bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	#bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

#
#bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
#samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

#Plot all reads
#bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#done



TAGDIR="$OUTDIR/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
PEAKSDIR="${OUTDIR}/Peaks"

#basic command
#makeTagDirectory ${TAGDIR}/129-43_ChIP_WT_input ${BAMDIR}/129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam
# List of sample names
#samples=(
#    "135-41_ChIP_WT_Anti6xHis_Rep1_S38"
    #"135-42_ChIP_delHcra_Anti6xHis_Rep1_S39"
#    "135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40"
#    "135-44_ChIP_WT_Input_Rep1_S41"
    #"135-45_ChIP_delHcra_Input_Rep1_S42"
#    "135-46_ChIP_hcra6xhis_Input_Rep1_S43"
#)

# Loop through each sample
#for sample in "${samples[@]}"; do
#    bam_file="${BAMDIR}/${sample}_L001_R1_001_val_1.fq.gz.bam"
#    tag_directory="${TAGDIR}/${sample}"
#    output_peak_file="${PEAKSDIR}/${sample}_peaks.txt"
#    input_bam="${BAMDIR}/${sample}_L001_R1_001_val_1.fq.gz.bam"

    # Run makeTagDirectory
#    makeTagDirectory "${tag_directory}" "${bam_file}"

    # Run findPeaks
  #  findPeaks "${tag_directory}" -style factor -region -size 150 -minDist 530 -o "${output_peak_file}" -i "${input_bam}"
#done
findPeaks "${TAGDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38" -style factor  -size 15 -o "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks.txt" -i "${TAGDIR}/135-44_ChIP_WT_Input_Rep1_S41"
findPeaks "${TAGDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40" -style factor -size 15  -o "${PEAKSDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40_peaks.txt" -i "${TAGDIR}/135-46_ChIP_hcra6xhis_Input_Rep1_S43"

annotatePeaks.pl "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks.txt" $GENOME -gtf /home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf > ${OUTDIR}/135-41_ChIP_WT.txt

annotatePeaks.pl "${PEAKSDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40_peaks.txt" $GENOME -gtf /home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf > ${OUTDIR}/135-43_hrcA6xhis.txt
