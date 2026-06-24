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

FULL_FASTA="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/metadata/all_species.fasta"

OUT="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/mmseqs_identity"

mkdir -p $OUT


#########################################
# DOMAIN SEQUENCE IDENTITY
#########################################

echo "Running domain MMseqs"


for fasta in $DOMAIN_DIR/*.fasta

do

domain=$(basename $(dirname $fasta))

echo "Running $domain"


mmseqs createdb \
$fasta \
${OUT}/${domain}_db


mmseqs easy-search \
${OUT}/${domain}_db \
${OUT}/${domain}_db \
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
# WHOLE PROTEIN MMSEQS
#########################################


echo "Running whole protein MMseqs"


mmseqs createdb \
$FULL_FASTA \
${OUT}/proteins_db



mmseqs easy-search \
${OUT}/proteins_db \
${OUT}/proteins_db \
${OUT}/protein_identity.tsv \
tmp_protein \
--format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"



echo "Whole protein complete"



#########################################
# Merge with FoldSeek master
#########################################


python3 <<'PY'


import pandas as pd



foldseek_file = "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/annotated_hits_expanded.csv"


domain_file = "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/mmseqs_identity/all_domain_sequence_identity.tsv"


protein_file = "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/mmseqs_identity/protein_identity.tsv"



print("Reading FoldSeek master")


foldseek_master = pd.read_csv(
    foldseek_file
)



#########################################
# Domain identity
#########################################


print("Reading domain identity")


domain_identity = pd.read_csv(
    domain_file,
    sep="\t"
)



foldseek_master = foldseek_master.merge(

    domain_identity,

    on=[
        "query_id",
        "target_id"
    ],

    how="left"

)



#########################################
# Protein identity
#########################################


print("Reading protein identity")


protein_identity = pd.read_csv(

    protein_file,

    sep="\t",

    header=None,

    names=[

        "query_accession",

        "target_accession",

        "protein_pident",

        "protein_alnlen",

        "protein_qcov",

        "protein_tcov",

        "protein_evalue",

        "protein_bits"

    ]

)



foldseek_master = foldseek_master.merge(

    protein_identity,

    on=[

        "query_accession",

        "target_accession"

    ],

    how="left"

)



#########################################
# Save
#########################################


outfile="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/annotated_hits_expanded_with_sequence_identity.csv"



foldseek_master.to_csv(

    outfile,

    index=False

)



print("DONE")
print(outfile)

print(foldseek_master.head())


PY
