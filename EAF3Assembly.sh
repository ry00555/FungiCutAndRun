#!/bin/bash
#SBATCH --job-name=EAF3_assembly
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --output=EAF3_assembly.%j.out
#SBATCH --error=EAF3_assembly.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

ml SAMtools/1.23.1-GCC-13.3.0 StringTie/3.0.0-GCC-13.3.0 gffread/0.12.7-GCCcore-12.3.0 TransDecoder/5.7.1-GCC-12.3.0

GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
WTBAMs="/lustre2/scratch/ry00555/Run155/SortedBamFiles/"
mkdir -p $OUTDIR

WT1="$WTBAMs/155-93_RNA_WT__Rep1_S93_L002/155-93_RNA_WT__Rep1_S93_L002_Aligned.sortedByCoord.out.bam"
WT2="$WTBAMs/155-94_RNA_WT__Rep2_S94_L002/155-94_RNA_WT__Rep2_S94_L002_Aligned.sortedByCoord.out.bam"
WT3="$WTBAMs/155-95_RNA_WT__Rep3_S95_L002/155-95_RNA_WT__Rep3_S95_L002_Aligned.sortedByCoord.out.bam"

# Step 1: extract just the EAF3 locus from each BAM (with padding)
for bam in $WT1 $WT2 $WT3; do
    name=$(basename $bam .bam)
    samtools view -b $bam "CM002237.1:1934000-1937500" > $OUTDIR/${name}_EAF3.bam
    samtools index $OUTDIR/${name}_EAF3.bam
done

# Step 2: merge the 3 extracted BAMs
samtools merge -f -@ 8 \
    $OUTDIR/EAF3_WT_merged.bam \
    $OUTDIR/*_EAF3.bam
samtools index $OUTDIR/EAF3_WT_merged.bam

# Step 3: assemble transcripts
stringtie $OUTDIR/EAF3_WT_merged.bam \
    -G $GTF \
    -o $OUTDIR/EAF3_assembly.gtf \
    -p 8 -v

# Step 4: extract transcript sequences
gffread $OUTDIR/EAF3_assembly.gtf \
    -g $GENOME \
    -w $OUTDIR/EAF3_transcripts.fa

echo "Done — transcripts written to $OUTDIR/EAF3_transcripts.fa"

#figure out amino acid sequence using transdecoder identify long ORFs (minimum 50 aa to catch short isoform)
TransDecoder.LongOrfs -t  $OUTDIR/EAF3_transcripts.fa -m 50

# predict coding regions, keep all ORFs per transcript (not just best)
TransDecoder.Predict -t  $OUTDIR/EAF3_transcripts.fa

# copy final protein sequences to a clean output file
cp  $OUTDIR/EAF3_transcripts.fa.transdecoder.pep  $OUTDIR/EAF3_proteins.fa

echo "Done — protein sequences written to $OUTDIR/EAF3_proteins.fa"
echo ""
echo "Predicted proteins:" # this will say the protein length for each sequence made in the .out file and if any of them are incomplete proteins 
grep ">"  $OUTDIR/EAF3_proteins.fa
