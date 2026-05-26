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

ml GCC/13.3.0 Exonerate/2.4.0-GCC-13.3.0 SAMtools/1.23.1-GCC-13.3.0 gffread/0.12.7-GCCcore-12.3.0

OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"

# Step 1: save the concatenated protein — written directly in script so path is guaranteed
cat > $OUTDIR/EAF3_full_protein.fa << 'EOF'
>EAF3_full
MAPSKTPQPPYSKDERVLCFHMEMLYEAKILDVQPTESGDGWSYKIHYKGWKSSWDDWVPQDRIRKLNDENKDLAQQLLAQYKQLQSGKAAKQPKKGGRPGGSDLSSARGSEERTAAGTTTQNNRNPRRARDFDLETVSGGVPFSMSRSMRSHGRDWDDDTPSRSSLISSFAACEARWSKLAEPADRLMRLPNKLTSTYVDRYGRTTFDDTYALGWTSTRPPPIGHNHPKIIRARNNGNWSEKLEAYYSELPPNYNIPLAGQACKKADEDWSSRSSTRVATPFRVDFGPFTILLEALAKIRKPREPTAGKKMWECRFTPKDVEKEKEDNFHNRPSIKLPLPDHVKALLVDDWENVTKNQQLVPIPHVHPVDEILKDYLEHERPNRVPESPQMDILEETVAGLREYFDRCLGRILLYRFERAQYHEQHLIWTAGTDEKHKSASDTYGAEHLARLLVSLPELVAQTNMDQQSVNRLREELIKFTNWFSRHTTKYFVSEYETPSQEYVDQARSV
EOF

# Verify the protein file was created
echo "Protein file:"
cat $OUTDIR/EAF3_full_protein.fa

# Step 2: extract genomic region
samtools faidx $GENOME "CM002237.1:1933959-1937714" > $OUTDIR/EAF3_fullgene.fa

# Verify genomic file was created
echo "Genomic file first line:"
head -1 $OUTDIR/EAF3_fullgene.fa

# Step 3: run exonerate
exonerate \
    --model protein2genome \
    --query $OUTDIR/EAF3_full_protein.fa \
    --target $OUTDIR/EAF3_fullgene.fa \
    --showtargetgff yes \
    --showalignment yes \
    --showvulgar yes \
    --percent 50 \
    --score 100 \
    > $OUTDIR/EAF3_exonerate.out

echo "Exonerate done"

# Step 4: print full output
cat $OUTDIR/EAF3_exonerate.out

# Step 5: extract exon coordinates
echo ""
echo "=== EXON COORDINATES ==="
grep -v "^#" $OUTDIR/EAF3_exonerate.out | grep "exon\|cds"

# Step 6: extract CDS GFF lines and get spliced sequence
grep -v "^#" $OUTDIR/EAF3_exonerate.out | \
    awk '$3 == "cds"' > $OUTDIR/EAF3_exonerate_cds.gff

echo ""
echo "CDS GFF lines:"
cat $OUTDIR/EAF3_exonerate_cds.gff

gffread $OUTDIR/EAF3_exonerate_cds.gff \
    -g $OUTDIR/EAF3_fullgene.fa \
    -x $OUTDIR/EAF3_spliced_cds.fa

echo ""
echo "=== SPLICED CDS ==="
cat $OUTDIR/EAF3_spliced_cds.fa
