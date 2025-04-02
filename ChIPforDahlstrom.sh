#!/bin/bash
#SBATCH --job-name=DahlstrohmChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../DahlstrohmChIPRun147.%j.out
#SBATCH --error=../DahlstrohmChIPRun147.%j.err

cd $SLURM_SUBMIT_DIR

ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/0.1.20-GCC-11.3.0
#If incompatible modules use SAMtools/1.16.1-GCC-11.3.0 instead rec by Chelsea Taylor
ml deepTools/3.5.2-foss-2022a



#ml Trim_Galore/0.6.7-GCCcore-11.2.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
#ml BEDTools/2.30.0-GCC-11.3.0
#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )


source ConfigforDahlstrom.txt

OUTDIR="/scratch/ry00555/Run147/Dahlstrom"
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
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}
	#file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	#name=${file/%_S[1-190]*_L001_R1_001_val_1.fq.gz/}

# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	#read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 #	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
#	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"

#ml SAMtools
#ml BWA
#
#bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
#samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools
#ml deepTools
#Plot all reads
#bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#done

#mkdir $OUTDIR/MACSPeaks
PEAKDIR="${OUTDIR}/MACSPeaks"

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
#147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.bam
#147-105_ChIP_C3_6His_Input_Rep_1_Ped_48hrPDA.bam

#147-102_ChIP_C3_6His__IP_Rep_2_Ped_48hrPDA.bam
#147-106_ChIP_C3_6His__Input_Rep_2_Ped_48hrPDA.bam
 macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.bam -f BAMPE -n 147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-105_ChIP_C3_6His_Input_Rep_1_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

  macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-102_ChIP_C3_6His__IP_Rep_2_Ped_48hrPDA.bam -f BAMPE -n 147-102_ChIP_C3_6His__IP_Rep_2_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-106_ChIP_C3_6His__Input_Rep_2_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

  #  147-97_ChIP_C1_6His_Input_Rep_1_Ped_48hrPDA.bam
  # 147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA.bam

	# 147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA.bam
	 #147-98_ChIP_C1_6His_Input_Rep_2_Ped_48hrPDA.bam

	 macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA.bam -f BAMPE -n 147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-97_ChIP_C1_6His_Input_Rep_1_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

	 macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA.bam -f BAMPE -n 147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-98_ChIP_C1_6His_Input_Rep_2_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

#147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA.bam
#147-108_ChIP_C3_WT_Input_Rep_2_Ped_48hrPDA.bam

#147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA.bam
#147-107_ChIP_C3_WT_Input_Rep_1_Ped_48hrPDA.bam

#147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA.bam
#  147-99_ChIP_C1_WT_Input_Rep_1_Ped_48hrPDA.bam

#147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA.bam
#	147-100_ChIP_C1_WT_Input_Rep2_Ped_48hrPDA.bam
macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA.bam -f BAMPE -n 147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-99_ChIP_C1_WT_Input_Rep_1_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA.bam -f BAMPE -n 147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-100_ChIP_C1_WT_Input_Rep2_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA.bam -f BAMPE -n 147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-107_ChIP_C3_WT_Input_Rep_1_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50

macs3 callpeak -t ${OUTDIR}/SortedBamFiles/147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA.bam -f BAMPE -n 147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA -c ${OUTDIR}/SortedBamFiles/147-108_ChIP_C3_WT_Input_Rep_2_Ped_48hrPDA.bam --narrow -g 8298884 --broad-cutoff 0.1 --outdir ${OUTDIR}/MACSPeaks --min-length 50
#  for infile in $BAMDIR/*_Q30.bam
 #do
#    base=$(basename ${infile} _Q30.bam)
# #  Input=$BAMDIR/ ${infile} Input_Q30.bam
#  macs3 callpeak -t $infile -f BAMPE -n $base --narrow -g 8298884 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 50 --max-gap 500 #-c $Input
#  done


TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
PEAKSDIR="${OUTDIR}/HomerPeaks"
ml Homer/4.11-foss-2022a
ml Perl/5.34.1-GCCcore-11.3.0
#basic command

makeTagDirectory ${TAGDIR}/147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA ${BAMDIR}/147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-108_ChIP_C3_WT_Input_Rep_2_Ped_48hrPDA ${BAMDIR}/147-108_ChIP_C3_WT_Input_Rep_2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA ${BAMDIR}/147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-107_ChIP_C3_WT_Input_Rep_1_Ped_48hrPDA ${BAMDIR}/147-107_ChIP_C3_WT_Input_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA ${BAMDIR}/147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-99_ChIP_C1_WT_Input_Rep_1_Ped_48hrPDA ${BAMDIR}/147-99_ChIP_C1_WT_Input_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA ${BAMDIR}/147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-100_ChIP_C1_WT_Input_Rep2_Ped_48hrPDA ${BAMDIR}/147-100_ChIP_C1_WT_Input_Rep2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-97_ChIP_C1_6His_Input_Rep_1_Ped_48hrPDA ${BAMDIR}/147-97_ChIP_C1_6His_Input_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA ${BAMDIR}/147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA ${BAMDIR}/147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-98_ChIP_C1_6His_Input_Rep_2_Ped_48hrPDA ${BAMDIR}/147-98_ChIP_C1_6His_Input_Rep_2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA ${BAMDIR}/147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-105_ChIP_C3_6His_Input_Rep_1_Ped_48hrPDA ${BAMDIR}/147-105_ChIP_C3_6His_Input_Rep_1_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-102_ChIP_C3_6His_IP_Rep_2_Ped_48hrPDA ${BAMDIR}/147-102_ChIP_C3_6His__IP_Rep_2_Ped_48hrPDA.bam
makeTagDirectory ${TAGDIR}/147-106_ChIP_C3_6His_Input_Rep_2_Ped_48hrPDA ${BAMDIR}/147-106_ChIP_C3_6His__Input_Rep_2_Ped_48hrPDA.bam

findPeaks ${TAGDIR}/147-104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA -style factor -o ${PEAKSDIR}/104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA.txt -i ${TAGDIR}/147-108_ChIP_C3_WT_Input_Rep_2_Ped_48hrPDA
findPeaks ${TAGDIR}/147-103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA -style factor -o ${PEAKSDIR}/103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA.txt -i ${TAGDIR}/147-107_ChIP_C3_WT_Input_Rep_1_Ped_48hrPDA
findPeaks ${TAGDIR}/147-95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA -style factor -o ${PEAKSDIR}/95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA.txt -i ${TAGDIR}/147-99_ChIP_C1_WT_Input_Rep_1_Ped_48hrPDA
findPeaks ${TAGDIR}/147-96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA -style factor -o ${PEAKSDIR}/96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA.txt -i ${TAGDIR}/147-100_ChIP_C1_WT_Input_Rep2_Ped_48hrPDA
findPeaks ${TAGDIR}/147-93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA -style factor -o ${PEAKSDIR}/93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA.txt -i ${TAGDIR}/147-97_ChIP_C1_6His_Input_Rep_1_Ped_48hrPDA
findPeaks ${TAGDIR}/147-94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA -style factor -o ${PEAKSDIR}/94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA.txt -i ${TAGDIR}/147-98_ChIP_C1_6His_Input_Rep_2_Ped_48hrPDA
findPeaks ${TAGDIR}/147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA -style factor -o ${PEAKSDIR}/101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.txt -i ${TAGDIR}/147-105_ChIP_C3_6His_Input_Rep_1_Ped_48hrPDA

findPeaks ${TAGDIR}/147-101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA -style factor -o ${PEAKSDIR}/101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.txt -i ${TAGDIR}/147-106_ChIP_C3_6His_Input_Rep_2_Ped_48hrPDA

for peakfile in ${PEAKSDIR}/*.txt; do
    bedfile="${peakfile%.txt}.bed"
    pos2bed.pl "$peakfile" > "$bedfile"
done
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
#line by line commands that worked as the homer module has changed arguments.
#note that the WT antibody without input normalization used a peak size of 10 and tgat gave 1024 peaks. However 43_hrcA6xhis without input normalization peaks were called with a size of 100 not 10 and that gave 740. Ultimately we need the input normalized.
#findPeaks "135-41_ChIP_WT_Anti6xHis_Rep1_S38" -style factor -L 0 -LP 0.01 -P 0.01 -F 0 -fdr 0.01 -poisson 0.01 -size 100 -o "../Peaks/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks_wINPUT.txt" -i "135-44_ChIP_WT_Input_Rep1_S41"
#findPeaks "135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40" -style factor -L 0 -LP 0.01 -P 0.01 -F 0 -fdr 0.01 -poisson 0.01 -size 100 -o "../Peaks/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40_peaks_wINPUT.txt" -i "135-46_ChIP_hcra6xhis_Input_Rep1_S43"

#annotatePeaks.pl "135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks_wINPUT.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > ../135-41_ChIP_WT.txt
#!/bin/bash

# Define directories


# Run annotatePeaks.pl for each peak file
#annotatePeaks.pl "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks_wINPUT.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks_wINPUT_ann.txt"

annotatePeaks.pl "${PEAKSDIR}/104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA.txt" $GENOME -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/104_ChIP_C3_WT_IP_Rep_2_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/103_ChIP_C3_WT_IP_Rep_1_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/95_ChIP_C1_WT_IP_Rep1_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/96_ChIP_C1_WT_IP_Rep2_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/93_ChIP_C1_6His_IP_Rep1_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/94_ChIP_C1_6His_IP_Rep2_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/94_ChIP_C1_6His_IP_Rep_2_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/101_ChIP_C3_6His_IP_Rep_1_Ped_48hrPDA_ann.txt"
annotatePeaks.pl "${PEAKSDIR}/102_ChIP_C3_6His_IP_Rep_2_Ped_48hrPDA.txt" "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.fna" -gtf "/home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf" > "${PEAKSDIR}/102_ChIP_C3_6His_IP_Rep_2_Ped_48hrPDA_ann.txt"








#findPeaks "${TAGDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38" -style factor -o "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks.txt"
#-i "${TAGDIR}/135-44_ChIP_WT_Input_Rep1_S41"
#findPeaks "${TAGDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40" -style factor -o "${PEAKSDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40_peaks.txt"
#-i "${TAGDIR}/135-46_ChIP_hcra6xhis_Input_Rep1_S43"

#annotatePeaks.pl "${PEAKSDIR}/135-41_ChIP_WT_Anti6xHis_Rep1_S38_peaks_wINPUT.txt" $GENOME -gtf /home/ry00555/Research/Genomes/genomic.gtf > ${OUTDIR}/135-41_ChIP_WT.txt

#annotatePeaks.pl "${PEAKSDIR}/135-43_ChIP_hcra6xhis_Anti6xHis_Rep1_S40_peaks.txt" $GENOME -gtf /home/ry00555/Research/Genomes/GCA_019428685.1_ASM1942868v1_genomic.gtf > ${OUTDIR}/135-43_hrcA6xhis.txt
