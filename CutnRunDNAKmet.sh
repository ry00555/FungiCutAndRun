#!/bin/bash
#SBATCH --job-name=CUTandRun137
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapQual137.%j.out
#SBATCH --error=../MapQual137.%j.err
cd $SLURM_SUBMIT_DIR
OUTDIR="/scratch/ry00555/test"
#NCTRIMMED="/scratch/ry00555/OutputRun137/MapQual_CutNRun/TrimmedReads/"
FILES="/scratch/ry00555/test/trimmed/*_R1_001_val_1.fq.gz" #Don't forget the *

FASTQ="/scratch/ry00555/OutputRun137/CutandRun"

############## Take fastq files and align to Ncrassa genome to make sorted bam files
#test --- this ended up working
 #  ml STAR
 # for file in $FASTQ/*fastq\.gz;
 #  do
 #    if [[ $prefix ]]; then
 #          base=$(basename ${first} _R1_001.fastq.gz)
 #          sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/PE_trim_and_starNcrassa.sh -o $OUTDIR -n $base -m one $first $file
 #          prefix=
 #      else
 #         first=$file
 #          prefix=${file%%_*}
 #      fi
 #  done
 # ###Remove PCR duplicates
 # ml picard
 # module load SAMtools
 #  for infile in $OUTDIR/bams/*_q1.bam
 #  do
 #    base=$(basename ${infile} _q1.bam)
 #    java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/bams/"$base"_dupmetrics.txt -O $OUTDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
 #  done
#
ml BEDTools
# #Ncrassa nodups.bam conversion to bed files to call peaks and then pass through KMetSpikeIn Script
# for infile in $OUTDIR/bams/*nodups.bam
# do
# base=$(basename ${infile} .bam)
#    bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bed_files/$base.btb.bed
# done

#taking N crassa aligned beds into Kmet spike in script to make bedgraphs
for file in $OUTDIR/bams/*nodups.bam
   do
      base=$(basename "${file}" _nodups.bam)
  sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/kmet_spike.kd.sh $OUTDIR/KmetSpikeIn/bedgraphs $base $OUTDIR/trimmed/${base}*_R1_001_val_1.fq.gz \ $OUTDIR/trimmed/${base}*_R2_001_val_2.fq.gz $file bga $OUTDIR/genome/chrNameLength.txt
  done

  # sort bga files from Kmet spike in
     ml ucsc
      for infile in $OUTDIR/KmetSpikeIn/bedgraphs/*kmet.bga
        do
        base=$(basename ${infile} _kmet.bga)
            bedSort $infile $OUTDIR/KmetSpikeIn/bedgraphs/${base}.kmet_sort.bga
          done

  #       turning bedgraphs (normalized bga files) into bigwigs (bigwig files are for creating pictures)
  #     Kmet spike in of N crassa aligned
        ml ucsc
        for infile in $OUTDIR/KmetSpikeIn/bedgraphs/*.kmet_sort.bga
        do
         base=$(basename ${infile} .kmet_sort.bga)
        bedGraphToBigWig $infile $OUTDIR/ref/Ncrassa_ref/chrNameLength.txt $OUTDIR/KmetSpikeIn/BigWigs/${base}.KmetSpikeIn.bw
        done

############## Take fastq files and align to Ecoli genome to make sorted bam files
 ml STAR
  for file in $FASTQ/*fastq\.gz;
  do
    if [[ $prefix ]]; then
          base=$(basename ${first} _R1_001.fq.gz)
          sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/PE_trim_and_star_e_coli_RY.sh -o $OUTDIR/Ecoli_Aligned -n $base -m one $first $file
          prefix=
      else
          first=$file
          prefix=${file%%_*}
      fi
  done

  # ###Remove PCR duplicates
  ml picard
  module load SAMtools
  #
   for infile in $OUTDIR/Ecoli_Aligned/SortedBamFiles/*_ecoli_q1.bam
   do
     base=$(basename ${infile} _ecoli_q1.bam)
     java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/Ecoli_Aligned/SortedBamFiles/"$base"_dupmetrics.txt -O $OUTDIR/Ecoli_Aligned/SortedBamFiles/"$base"_ecoli_nodups.bam --REMOVE_DUPLICATES true
   done

   ml BEDTools
   # #Ecoli nodups.bam
#Feeding output of no duplicates to DNA Spike In normalization by taking the bam to convert to the bed and get a bedgraph
  for infile in $OUTDIR/Ecoli_Aligned/SortedBamFiles/*_ecoli_nodups.bam
 do
   base=$(basename ${infile} _ecoli_nodups.bam)
   bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/Ecoli_Aligned/bed_files/${base}_ecoli.btb.bed
   sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/DNAspike_in.kd.sh $OUTDIR/bed_files/${base}.btb.bed $OUTDIR/Ecoli_Aligned/bed_files/${base}_Ecoli.btb.bed 100000 bga "$OUTDIR/genome/chrNameLength.txt" 1 550 $OUTDIR/Ecoli_Aligned/bedgraphs/${base}Ecoli.norm.bga
    done

    # Taking bed files from no duplicates output into the DNA Spike In Normalization to make bed graphs
            ml BEDTools
            for f in $FILES
            do
               name=$(basename "${f}" _R1_001_val_1.fq.gz)
               sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/DNAspike_in.kd.sh $OUTDIR/bed_files/${base}.btb.bed $OUTDIR/Ecoli_Aligned/bed_files/${base}_Ecoli.btb.bed 100000 bga "$OUTDIR/genome/chrNameLength.txt" 1 550 $OUTDIR/bedgraphs/${name}.norm.bga
            done
    # Sorting DNA Spike In normalized bed graphs
              ml ucsc
               for infile in $OUTDIR/bedgraphs/*norm.bga
                do
               base=$(basename ${infile} .norm.bga)
                   bedSort $infile $OUTDIR/bedgraphs/${base}.norm_sort.bga
                done
    #Converting DNA Spike In Normalized bed graphs to big wigs
               ml ucsc

            for infile in $OUTDIR/bedgraphs/*.norm_sort.bga
              do
             base=$(basename ${infile} .norm_sort.bga)
             bedGraphToBigWig $infile $OUTDIR/genome/chrNameLength.txt $OUTDIR/BigWigs/${base}_DNASpikeinNorm.bw
            done

#sort the EColi aligned DNA spike in normalized bedgraph
    ml ucsc
     for infile in $OUTDIR/Ecoli_Aligned/bedgraphs/*Ecoli.norm.bga
      do
     base=$(basename ${infile} .norm.bga)
         bedSort $infile $OUTDIR/Ecoli_Aligned/bedgraphs/${base}EColi.norm_sort.bga
      done
      #convert the EColi aligned DNA spike in normalized bedgraph into a big wig

      ml ucsc
     for infile in $OUTDIR/Ecoli_Aligned/bedgraphs/*.norm_sort.bga
     do
     base=$(basename ${infile} .norm_sort.bga)
     bedGraphToBigWig $infile $OUTDIR/genome/chrNameLength.txt $OUTDIR/Ecoli_Aligned/BigWigs/${base}EColi_DNASpikeinNorm.bw
     done

#Here taking the EColi aligned bams with no duplicates to feed into Kmet Normalizaton
  # for file in $OUTDIR/Ecoli_Aligned/bams/*nodups.bam
  #    do
  #       base=$(basename "${file}" .bam)
  #   sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/kmet_spike.kd.sh $OUTDIR/Ecoli_Aligned/KMetSpikeInbedgraphs $base $OUTDIR/Ecoli_Aligned/TrimmedReads/${base}*R1_001_val_1.fq.gz \ $OUTDIR/Ecoli_Aligned/TrimmedReads/${base}*R2_001_val_2.fq.gz $file bga $OUTDIR/genome/chrNameLength.txt
  #   done
#Sort Ecoli aligned K met normalized bedgraphs
    # ml ucsc
    #  for infile in $OUTDIR/Ecoli_Aligned/KMetSpikeInbedgraphs/*norm.bga
    #   do
    #  base=$(basename ${infile} .norm.bga)
    #      bedSort $infile $OUTDIR/Ecoli_Aligned/KMetSpikeInbedgraphs/${base}.norm_sort.bga
    #   done
      #convert the EColi aligned Kmet spike in normalized bedgraph into a big wig

 #     ml ucsc
 #  for infile in $OUTDIR/Ecoli_Aligned/KMetSpikeInbedgraphs/*.norm_sort.bga
 #    do
 #   base=$(basename ${infile} .norm_sort.bga)
 #   bedGraphToBigWig $infile $OUTDIR/genome/chrNameLength.txt $OUTDIR/BigWigs/${base}EColi_KmetSpikeinNorm.bw
 # done
