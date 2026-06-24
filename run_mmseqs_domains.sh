#!/bin/bash
#SBATCH --job-name=seq_identity
#SBATCH --partition=inter_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=20:00:00
#SBATCH --output=seq_identity%j.out
#SBATCH --error=seq_identity%j.err

ml MMseqs2/18-8cc5c-gompi-2025a

DOMAIN_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/domain_sequences"

OUT="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/mmseqs_identity"

mkdir -p $OUT


#########################################
# DOMAIN SEQUENCE IDENTITY
#########################################

echo "Running domain MMseqs"


for fasta in $DOMAIN_DIR/*_fixed.fasta

do

domain=$(basename "$fasta" _fixed.fasta)

echo "Running $domain"


mmseqs createdb \
$fasta \
${OUT}/${domain}_db


mmseqs easy-search \
$fasta \
$fasta \
${OUT}/${domain}_identity.tsv \
tmp_${domain} \
--format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"


done



#########################################
# Combine domain identity files
#########################################


echo "Combining domain tables"


DOMAIN_OUT="${OUT}/all_domain_sequence_identity.tsv"


echo -e "query_id\ttarget_id\tdomain_pident\tdomain_alnlen\tdomain_qcov\tdomain_tcov\tdomain_evalue\tdomain_bits\tdomain" \
> $DOMAIN_OUT



for file in ${OUT}/*_identity.tsv

do

domain=$(basename "$file" _identity.tsv)


awk -v d="$domain" \
'BEGIN{OFS="\t"}
{
print $1,$2,$3,$4,$5,$6,$7,$8,d
}' "$file" >> $DOMAIN_OUT


done



#########################################
# Combine domain MMseqs2 results
#########################################

echo "Combining domain identity tables..."

OUTFILE="${OUT}/all_domain_sequence_identity.tsv"

echo -e "query_id\ttarget_id\tpident\talnlen\tqcov\ttcov\tevalue\tbits\tdomain" > $OUTFILE


for file in ${OUT}/*_identity.tsv
do

    domain=$(basename "$file" _identity.tsv)

    awk -v d="$domain" \
    'BEGIN{OFS="\t"}
    {print $1,$2,$3,$4,$5,$6,$7,$8,d}' \
    "$file" >> $OUTFILE

done


echo "Domain identity complete"



#########################################
# Extract full-length chromatin proteins
#########################################

#########################################
# Extract full length proteins by domain
#########################################

echo "Extracting full length proteins by domain..."


python3 <<'PY'


from Bio import SeqIO
import csv
import os


fasta="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/metadata/all_species.fasta"


domain_file="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/metadata/domain_occurrences_extended.tsv"


outdir="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/protein_sequences"


os.makedirs(outdir,exist_ok=True)



records={}


for r in SeqIO.parse(fasta,"fasta"):

    records[r.id]=r



domain_proteins={}



with open(domain_file) as f:

    for row in csv.DictReader(f,delimiter="\t"):


        for rid,rec in records.items():


            if f"|{row['accession']}|" in rid:


                new=rec[:]

                new.id=f"{row['gene']}_{row['accession']}"


                new.description=""


                domain=row["domain"]


                domain_proteins.setdefault(domain,{})[row["accession"]] = new


                break




for domain,proteins in domain_proteins.items():


    outfile=f"{outdir}/{domain}_proteins.fasta"


    SeqIO.write(
        proteins.values(),
        outfile,
        "fasta"
    )


    print(domain,len(proteins))


PY


echo "Full length domain-specific proteins complete"

#########################################
# Whole protein MMseqs2 by domain family
#########################################

PROTEIN_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/protein_sequences"


echo "Running full length protein MMseqs by domain"


mkdir -p ${OUT}/protein_results


for fasta in ${PROTEIN_DIR}/*.fasta

do

domain=$(basename "$fasta" .fasta)

echo "Running full protein identity for $domain"


mmseqs createdb \
$fasta \
${OUT}/protein_results/${domain}_db



mmseqs easy-search \
${OUT}/protein_results/${domain}_db \
${OUT}/protein_results/${domain}_db \
${OUT}/protein_results/${domain}_identity.tsv \
tmp_${domain}_protein \
--format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"


done


echo "Protein identity complete"

#########################################
# Combine protein identity tables
#########################################

echo "Combining protein identities"


PROTEIN_OUT="${OUT}/all_protein_sequence_identity.tsv"


echo -e "query_id\ttarget_id\tprotein_pident\tprotein_alnlen\tprotein_qcov\tprotein_tcov\tprotein_evalue\tprotein_bits\tdomain" > $PROTEIN_OUT



for file in ${OUT}/protein_results/*_identity.tsv

do

domain=$(basename "$file" _identity.tsv)


awk -v d="$domain" \
'BEGIN{OFS="\t"}
{
print $1,$2,$3,$4,$5,$6,$7,$8,d
}' "$file" >> $PROTEIN_OUT


done


echo $PROTEIN_OUT

#########################################
# Merge with FoldSeek master
#########################################


python3 <<'PY'


import pandas as pd

mmseqs_fullprotein_identity = pd.read_csv(
    "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/mmseqs_identity/all_protein_sequence_identity.tsv",
    sep="\t"
)


foldseek_master = foldseek_master.merge(
    mmseqs_fullprotein_identity,
    on=["query_id","target_id"],
    how="left"
)
