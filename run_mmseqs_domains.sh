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


BASE="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results"

FOLDSEEK_MASTER="${BASE}/annotated_hits_expanded.csv"

DOMAIN_FASTA="${BASE}/domain_sequences"

PROTEIN_FASTA="${BASE}/protein_sequences"

OUT="${BASE}/mmseqs_identity"

mkdir -p $OUT


##################################################
# 1. DOMAIN SEQUENCE IDENTITY
##################################################

echo "Extracting domain FASTAs..."



for fasta in ${DOMAIN_FASTA}/*.fasta

do

domain=$(basename "$fasta" .fasta)


echo "Running domain: $domain"



mmseqs createdb \
$fasta \
${OUT}/${domain}_domain_db



mmseqs easy-search \
${OUT}/${domain}_domain_db \
${OUT}/${domain}_domain_db \
${OUT}/${domain}_identity.tsv \
tmp_${domain} \
-e 1e-3 \
-c 0.5 \
--cov-mode 0 \
--format-output \
"query,target,pident,alnlen,qcov,tcov,evalue,bits"



done





#########################################
# Extract full-length proteins by domain
#########################################

echo "Extracting full-length proteins..."

python3 <<'PY'

from Bio import SeqIO
import csv
import os


fasta = "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/metadata/all_species.fasta"

domains = "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/metadata/domain_occurrences_extended.tsv"

outdir="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/protein_sequences"


os.makedirs(outdir,exist_ok=True)


# load proteins

records={}

for r in SeqIO.parse(fasta,"fasta"):

    records[r.id]=r



domain_proteins={}



with open(domains) as f:

    reader=csv.DictReader(f,delimiter="\t")

    for row in reader:


        accession=row["accession"]
        domain=row["domain"]
        gene=row["gene"]


        for rid,rec in records.items():

            if f"|{accession}|" in rid:


                new=rec[:]

                new.id=f"{gene}_{accession}"

                new.description=""


                domain_proteins.setdefault(
                    domain,
                    {}
                )[accession]=new


                break



# write one fasta per domain

for domain,proteins in domain_proteins.items():

    outfile=f"{outdir}/{domain}_proteins.fasta"


    SeqIO.write(
        proteins.values(),
        outfile,
        "fasta"
    )


    print(domain,len(proteins))


PY



echo "Full-length protein FASTAs created"



#########################################
# Whole protein identity by domain class
#########################################

echo "Running full protein MMseqs"


PROTEIN_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/protein_sequences"


mkdir -p ${OUT}/protein_results



for fasta in ${PROTEIN_DIR}/*_proteins.fasta

do


domain=$(basename "$fasta" _proteins.fasta)


echo "Running full protein identity: $domain"



# unique database per domain

mmseqs createdb \
$fasta \
${OUT}/protein_results/${domain}_protein_db



mmseqs easy-search \
${OUT}/protein_results/${domain}_protein_db \
${OUT}/protein_results/${domain}_protein_db \
${OUT}/protein_results/${domain}_protein_identity.tsv \
tmp_${domain}_protein \
-e 1e-5 \
-c 0.3 \
--cov-mode 0 \
--format-output \
"query,target,pident,alnlen,qcov,tcov,evalue,bits"



done



echo "Full protein MMseqs complete"



#########################################
# Combine full protein identity
#########################################


PROTEIN_OUT="${OUT}/all_protein_sequence_identity.tsv"


echo -e "query_id\ttarget_id\tprotein_pident\tprotein_alnlen\tprotein_qcov\tprotein_tcov\tprotein_evalue\tprotein_bits\tdomain" \
> $PROTEIN_OUT



for file in ${OUT}/protein_results/*_protein_identity.tsv

do


domain=$(basename "$file" _protein_identity.tsv)



awk -v d="$domain" \
'BEGIN{OFS="\t"}
{
print $1,$2,$3,$4,$5,$6,$7,$8,d
}' "$file" >> $PROTEIN_OUT



done



echo "Saved:"
echo $PROTEIN_OUT
#########################################
# Merge with FoldSeek master
#########################################


python3 <<'PY'


import pandas as pd


foldseek_master = foldseek_master.merge(
    mmseqs_fullprotein_identity,
    on=["query_id","target_id"],
    how="left"
)
