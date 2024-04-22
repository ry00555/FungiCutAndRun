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
Working Directory = /home/ry00555/Bioinformatics/CrassaGenome
 scp -r ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/Bioinformatics/GATK/GATK109Bam /Users/ry00555/Desktop/BamFiles/Run109mus30mei3V2/GATKBamFiles

#Load these modules that are compatible with GATK version 4.3
ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
ml picard/2.27.4-Java-13.0.2
ml BWA/0.7.17-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load BEDOPS/2.4.39-foss-2019b
ml Bowtie2/2.4.1-GCC-8.3.0
ml R/3.6.2-foss-2019b

#This section will create the reference files needed for most GATK tool commands.
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  -R /scratch/ry00555/Genomes/CM002239.1_chr4.fasta \
  -O /scratch/ry00555/Genomes/CM002239.1_chr4.dict

java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
-I /scratch/ry00555/Genomes/crassachr4.bed \
-R /scratch/ry00555/Genomes/CM002239.1_chr4.fasta \
-SD /scratch/ry00555/Genomes/CM002239.1_chr4.dict \
-O /scratch/ry00555/Genomes/CM002239.1_chr4.interval_list


gatk PreprocessIntervals \
-R /scratch/ry00555/Genomes/CM002239.1_chr4.fasta \
-L /scratch/ry00555/Genomes/CM002239.1_chr4.interval_list \
--interval-merging-rule OVERLAPPING_ONLY \
--bin-length 10 \
--padding 0 \
-O /scratch/ry00555/Genomes/Crassachr4.preprocessed10_intervals.interval_list

gatk AnnotateIntervals \
 -R /scratch/ry00555/Genomes/CM002239.1_chr4.fasta \
 -L /scratch/ry00555/Genomes/Crassachr4.preprocessed10_intervals.interval_list \
 --interval-merging-rule OVERLAPPING_ONLY \
 -O /scratch/ry00555/Genomes/Crassachr4_preprocessed10_annotated_intervals.tsv

#109_58 is a wildtype sample. 109_59 is a ∆mus30, ∆mei3 sample. The bam files were aligned using MapCutandRun.sh using the reference fasta file, GCF_000182925.2.fasta, to keep consistent. The sorted bam files are in  /scratch/ry00555/OutputRun109/Run109Bam
#There is a way to create bam files from pair end reads fasta.gz files, but when mapped to IGV, they are empty.
#In MapCutandRun.sh, read group ID, sample name, library number are not annotated in the bam file. I ran GATK AddOrReplaceReadGroups after samtools index 109_nn_Genomic.bam files.
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I /scratch/ry00555/OutputRun109/V2/SortedBamFiles/109_59_Genomic.bam \
-O /scratch/ry00555/Bioinformatics/GATK/GATK109Bam/109_59_Genomic.bamoutput.bam \
-RGID 2 \
-RGLB lib1 \
-RGPL illumina \
-RGPU S34 \
-RGSM 109_59


gatk CollectReadCounts \
-I GATK109Bam/109_58_Genomic.bamoutput.bam \
-R /scratch/ry00555/Genomes/CM002239.1_chr4.fasta \
-L /scratch/ry00555/Genomes/Crassachr4.preprocessed100_intervals.interval_list  \
--interval-merging-rule OVERLAPPING_ONLY \
-O 109tsv/chr4_109_58preprocessed100cnv.counts.tsv

# only add arguments because this error occurs when I add multiple input files (2 WT samples) org.broadinstitute.hellbender.exceptions.GATKException: Could not create panel of normals.  It may be necessary to use stricter parameters for filtering.  For example, use a larger value of minimum-interval-median-percentile.
gatk CreateReadCountPanelOfNormals \
-I 109tsv/chr4_109_58preprocessed100cnv.counts.tsv \
--annotated-intervals /scratch/ry00555/Genomes/Crassachr4_preprocessed100_annotated_intervals.tsv \
-O 109PanelofNormals/chr4_109_WTpreprocessedcnv.pon.hdf5


gatk DenoiseReadCounts \
-I 109tsv/109_63preprocessed100cnv.counts.tsv \
--annotated-intervals /scratch/ry00555/Genomes/Crassachr4_preprocessed10_annotated_intervals.tsv \
--count-panel-of-normals 109PanelofNormals/109_58preprocessedcnv.pon.hdf5 \
--standardized-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.standardizedCR.tsv \
--denoised-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.denoisedCR.tsv

109_58v61preprocessed.denoisedCR.tsv	  109_60v58preprocessed.denoisedCR.tsv	    109_62v58preprocessed.denoisedCR.tsv      109_63v61preprocessed.denoisedCR.tsv
109_58v61preprocessed.standardizedCR.tsv  109_60v58preprocessed.standardizedCR.tsv  109_62v58preprocessed.standardizedCR.tsv  109_63v61preprocessed.standardizedCR.tsv
109_59v58preprocessed.denoisedCR.tsv	  109_60v61preprocessed.denoisedCR.tsv	    109_62v61preprocessed.denoisedCR.tsv      109_64v58preprocessed.denoisedCR.tsv
109_59v58preprocessed.standardizedCR.tsv  109_60v61preprocessed.standardizedCR.tsv  109_62v61preprocessed.standardizedCR.tsv  109_64v58preprocessed.standardizedCR.tsv
109_59v61preprocessed.denoisedCR.tsv	  109_61v58preprocessed.denoisedCR.tsv	    109_63v58preprocessed.denoisedCR.tsv      109_64v61preprocessed.denoisedCR.tsv
109_59v61preprocessed.standardizedCR.tsv  109_61v58preprocessed.standardizedCR.tsv  109_63v58preprocessed.standardizedCR.tsv  109_64v61preprocessed.standardizedCR.tsv

gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.standardizedCR.tsv \
--denoised-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.denoisedCR.tsv  \
--sequence-dictionary /scratch/ry00555/Genomes/CM002239.1_chr4.dict \
--point-size-copy-ratio 1 \
--output-prefix 63v58chr410point1 \
--output PlotDenoisedCopyRatios

scp -r ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/Bioinformatics/GATK/PlotModelSegments /Users/ry00555/Desktop/RochelleLabDesktopData/IGV/mus30xmei3/mus30Samples

gatk CollectAllelicCounts \
          -I GATK109Bam/109_58_Genomic.bamoutput.bam \
          -R /scratch/ry00555/Genomes/GCF_000182925.2.fasta \
          -L /scratch/ry00555/Genomes/Crassa.preprocessed100_intervals.interval_list \
          -O AllelicCounts/109_58.allelicCounts.tsv

gatk ModelSegments \
--denoised-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.denoisedCR.tsv \
--output-prefix 63v58chr410point1 \
-O ModelSegments

     gatk CallCopyRatioSegments \
          -I tumor.cr.seg \
          -O tumor.called.seg


 gatk PlotModeledSegments \
--denoised-copy-ratios CopyRatios/10chr4_109_63v58preprocessed.denoisedCR.tsv \
--segments ModelSegments/63v58chr410point1.modelFinal.seg \
         --sequence-dictionary /scratch/ry00555/Genomes/CM002239.1_chr4.dict\
         --point-size-copy-ratio 1 \
         --output-prefix 109_63v58chr410point1 \
         -O PlotModelSegments
