#!/bin/bash
#SBATCH --job-name=EAF3_exonerate
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=00:30:00
#SBATCH --output=EAF3_exonerate.%j.out
#SBATCH --error=EAF3_exonerate.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

ml GCC/13.3.0 Exonerate/2.4.0-GCC-13.3.0 SAMtools/1.23.1-GCC-13.3.0

OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"



# Step 2: extract the genomic region spanning both NCU06788 and NCU06787
# Using wider window to make sure we capture everything
samtools faidx $GENOME "CM002237.1:1933959-1937714" > $OUTDIR/EAF3_fullgene.fa

# Step 3: run exonerate
exonerate \
    --model protein2genome \
    --query EAF3_NCU06788_NCU06787Concatenate.fa \
    --target $OUTDIR/EAF3_fullgene.fa \
    --showtargetgff yes \
    --showalignment yes \
    --showvulgar yes \
    --percent 50 \
    --score 100 \
    > $OUTDIR/EAF3_exonerate.out

echo "Done"

# Step 4: print the alignment
echo ""
echo "=== FULL OUTPUT ==="
cat $OUTDIR/EAF3_exonerate.out

# Step 5: print just the exon coordinates
echo ""
echo "=== EXON COORDINATES ==="
grep "exon\|cds" $OUTDIR/EAF3_exonerate.out | grep -v "^#"

# Step 6: extract the CDS sequence using the GFF output
echo ""
echo "=== EXTRACTING SPLICED CDS ==="
grep -v "^#" $OUTDIR/EAF3_exonerate.out | \
    grep $'^\t' | \
    awk '$3 == "cds"' > $OUTDIR/EAF3_exonerate_cds.gff

gffread $OUTDIR/EAF3_exonerate_cds.gff \
    -g $OUTDIR/EAF3_fullgene.fa \
    -x $OUTDIR/EAF3_spliced_cds.fa

echo ""
echo "=== SPLICED CDS SEQUENCE ==="
cat $OUTDIR/EAF3_spliced_cds.fa
