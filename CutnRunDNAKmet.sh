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

#
# #aligning to ecoli genome
 ml STAR
  for file in $FASTQ/*fastq\.gz;
  do
    if [[ $prefix ]]; then
          base=$(basename ${first} _R1_001_val_1.fq.gz)
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
   for infile in $OUTDIR/bams/*_q1.bam
   do
     base=$(basename ${infile} _q1.bam)
     java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/bams/"$base"_dupmetrics.txt -O $OUTDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
   done
  #
  #
   for infile in $OUTDIR/Ecoli_Aligned/bams/*_ecoli_q1.bam
   do
     base=$(basename ${infile} _ecoli_q1.bam)
     java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/Ecoli_Aligned/bams/"$base"_dupmetrics.txt -O $OUTDIR/Ecoli_Aligned/bams/"$base"_ecoli_nodups.bam --REMOVE_DUPLICATES true
   done

   ml BEDTools
   # #Ncrassa nodups.bam
for infile in $OUTDIR/bams/*nodups.bam
 do
   base=$(basename ${infile} .bam)
      bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bed_files/$base.btb.bed
 done
   # #Ecoli nodups.bam
  for infile in $OUTDIR/Ecoli_Aligned/bams/*_ecoli_nodups.bam
 do
   base=$(basename ${infile} _ecoli_nodups.bam)
   bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/Ecoli_Aligned/bed_files/${base}_ecoli.btb.bed
    done

 ml BEDTools
for f in $FILES
 do
    name=$(basename "${f}" _R1_001_val_1.fq.gz)
    sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/DNAspike_in.kd.sh $OUTDIR/bed_files/${base}.btb.bed $OUTDIR/Ecoli_Aligned/bed_files/${base}_Ecoli.btb.bed 100000 bga "$OUTDIR/genome/chrNameLength.txt" 1 550 $OUTDIR/bedgraphs/${name}.norm.bga
 done

   ml ucsc
    for infile in $OUTDIR/bedgraphs/*norm.bga
     do
    base=$(basename ${infile} .norm.bga)
        bedSort $infile $OUTDIR/bedgraphs/${base}.norm_sort.bga
     done

    ml ucsc
 for infile in $OUTDIR/bedgraphs/*.norm_sort.bga
   do
  base=$(basename ${infile} .norm_sort.bga)
  bedGraphToBigWig $infile $OUTDIR/ref/Ncrassa_ref/chrNameLength.txt $OUTDIR/BigWigs/${base}_DNASpikeinNorm.bw
done

for file in $OUTDIR/bams/*nodups.bam
   do
      base=$(basename "${file}" .bam)
  sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/kmet_spike.kd.sh $OUTDIR/KmetSpikeIn/bedgraphs $base $OUTDIR/trimmed/${base}*R1_001_val_1.fq.gz \ $OUTDIR/trimmed/${base}*R2_001_val_2.fq.gz $file bga $OUTDIR/genome/chrNameLength.txt
  done

  for file in $OUTDIR/Ecoli_Aligned/bams/*nodups.bam
     do
        base=$(basename "${file}" .bam)
    sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/kmet_spike.kd.sh $OUTDIR/Ecoli_Aligned/bedgraphs $base $OUTDIR/trimmed/${base}*R1_001_val_1.fq.gz \ $OUTDIR/trimmed/${base}*R2_001_val_2.fq.gz $file bga $OUTDIR/genome/chrNameLength.txt
    done

    ml ucsc
     for infile in $OUTDIR/Ecoli_Aligned/bedgraphs/*norm.bga
      do
     base=$(basename ${infile} .norm.bga)
         bedSort $infile $OUTDIR/Ecoli_Aligned/bedgraphs/${base}.norm_sort.bga
      done

     ml ucsc
  for infile in $OUTDIR/Ecoli_Aligned/bedgraphs/*.norm_sort.bga
    do
   base=$(basename ${infile} .norm_sort.bga)
   bedGraphToBigWig $infile $OUTDIR/genome/chrNameLength.txt $OUTDIR/BigWigs/${base}EColi_DNASpikeinNorm.bw
 done
