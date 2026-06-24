#!/bin/bash
#SBATCH --job-name=chromatin_domains
#SBATCH --partition=inter_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=chromatin_domains_%j.out
#SBATCH --error=chromatin_domains_%j.err


set -euo pipefail


##################################################
# MODULES
##################################################

module load HMMER/3.4-gompi-2024a
module load Foldseek/10-941cd33-GPU
module load MAFFT/7.505-GCC-11.3.0-with-extensions


##################################################
# PATHS
##################################################

BASE=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek

OUT=${BASE}/chromatin_domain_results

cd ${BASE}

#mkdir -p \
#${OUT}/metadata \
#${OUT}/domains \
#${OUT}/domain_sequences \
#${OUT}/domain_hmms \
#${OUT}/domains_extracted \
#${OUT}/foldseek \
#${OUT}/tmp



##################################################
# PFAM
##################################################

#PFAM=${OUT}/Pfam-A.hmm

#if [ ! -f ${PFAM}.h3i ]
#then

#    wget -q \
#        -O ${PFAM}.gz \
#        https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

#    gunzip ${PFAM}.gz

#    hmmpress ${PFAM}

#fi


##################################################
# MAKE CHROMATIN-ONLY PFAM DATABASE
##################################################

CHROM_PFAM=${OUT}/metadata/chromatin_Pfam.hmm

DOMAIN_LIST=${OUT}/metadata/chromatin_pfams.txt

# if [ ! -f "${CHROM_PFAM}.h3i" ]
# then
#
# echo "Creating chromatin Pfam database"
#
# cat > "${DOMAIN_LIST}" << EOF
# BAH
# CBAH
# PHD
# PHD_2
# PHD_3
# PHD_4
# PHD_5
# PHD_like
# PHD_JHD1
# PHD_ash2p_like
# zf-PHD-like
# PWWP
# MUM1-like_PWWP
# Chromo
# Chromo_2
# Chromo_shadow
# Bromodomain
# MRG
# PHF12_MRG_bd
# MBT
# Tudor_2
# Tudor_3
# Tudor_4
# Tudor_5
# 53-BP1_Tudor
# Crb2_Tudor
# SGF29_Tudor
# SMN_Tudor
# JmjC
# JmjC_2
# JmjN
# SET
# N-SET
# Pre-SET
# preSET_CXC
# SET_assoc
# SET7_N
# DOT1
# WD40
# EOF
#
#
# echo "Indexing Pfam"
#
# hmmfetch --index "${PFAM}"
#
#
# echo "Extracting chromatin Pfam models"
#
# hmmfetch \
#     -f \
#     "${PFAM}" \
#     "${DOMAIN_LIST}" \
#     > "${CHROM_PFAM}"
#
#
# echo "Pressing chromatin Pfam database"
#
# hmmpress "${CHROM_PFAM}"
#
#
# echo "Chromatin Pfam database created"
#
# fi
##################################################
# COMBINE PROTEOMES
##################################################

ALL_FASTA=${OUT}/metadata/all_species.fasta

#> ${ALL_FASTA}


#for f in ${BASE}/*_proteome.fasta

#do

#species=$(basename ${f} _proteome.fasta)

#echo "Adding ${species}"

#sed "s/^>/>${species}|/" ${f} >> ${ALL_FASTA}

#done



##################################################
# PFAM SCAN
##################################################

#hmmscan \
#--cpu ${SLURM_CPUS_PER_TASK} \
#--domtblout ${OUT}/metadata/chromatin_hits.domtbl \
#-E 1e-3 \
#--domE 1e-3 \
#${CHROM_PFAM} \
#${ALL_FASTA}



##################################################
# PARSE PFAM DOMAINS
##################################################

python3 <<'PY'

import csv

DOMAINS={

"BAH":"BAH",
"CBAH":"BAH",

"PHD":"PHD",
"PHD_2":"PHD",
"PHD_3":"PHD",
"PHD_4":"PHD",
"PHD_5":"PHD",
"PHD_like":"PHD",
"PHD_JHD1":"PHD",
"PHD_ash2p_like":"PHD",
"zf-PHD-like":"PHD",

"PWWP":"PWWP",
"MUM1-like_PWWP":"PWWP",

"Chromo":"Chromo",
"Chromo_2":"Chromo",
"Chromo_shadow":"ChromoShadow",

"Bromodomain":"Bromodomain",

"MRG":"MRG",
"PHF12_MRG_bd":"MRG",

"MBT":"MBT",

"JmjC":"JmjC",
"JmjC_2":"JmjC",
"JmjN":"JmjN",

"SET":"SET",
"N-SET":"SET",
"Pre-SET":"SET",
"preSET_CXC":"SET",
"SET_assoc":"SET",
"SET7_N":"SET",

"DOT1":"DOT1",
"53-BP1_Tudor":"Tudor",
"Crb2_Tudor":"Tudor",
"SGF29_Tudor":"Tudor",
"SMN_Tudor":"Tudor",
"Tudor_2":"Tudor",
"Tudor_3":"Tudor",
"Tudor_4":"Tudor",
"Tudor_5":"Tudor",
"WD40":"WD40"
}



inp="chromatin_domain_results/metadata/chromatin_hits.domtbl"

out="chromatin_domain_results/metadata/domain_occurrences.tsv"



rows=[]


for line in open(inp):

    if line.startswith("#"):
        continue


    x=line.split()


    if len(x)<23:
        continue


    pfam=x[0]


    if pfam not in DOMAINS:
        continue


    protein=x[3]


    parts=protein.split("|")


    if len(parts) < 4:
        continue

    species=parts[0]
    accession=parts[2]
    gene=parts[3]

    rows.append([

        species,
        accession,
        gene,
        DOMAINS[pfam],
        x[17],
        x[18],
        x[6]

    ])



with open(out,"w") as f:

    w=csv.writer(f,delimiter="\t")

    w.writerow([
    "species",
    "accession",
    "gene",
    "domain",
    "start",
    "end",
    "evalue"
    ])

    w.writerows(rows)



print("PFAM domains:",len(rows))


PY



##################################################
# EXTRACT DOMAIN SEQUENCES
##################################################

python3 <<'PY'


from Bio import SeqIO
import csv


fasta="chromatin_domain_results/metadata/all_species.fasta"

dom="chromatin_domain_results/metadata/domain_occurrences.tsv"


records={}

for r in SeqIO.parse(fasta,"fasta"):

    records[r.id]=r


out={}


with open(dom) as f:

    for row in csv.DictReader(f,delimiter="\t"):


        key=f"{row['species']}|sp|{row['accession']}|{row['gene']}"


        if key not in records:
            continue


        seq=records[key].seq


        start=int(row["start"])
        end=int(row["end"])


        domain=row["domain"]


        new=records[key][start-1:end]


        new.id=f"{row['gene']}_{domain}"


        out.setdefault(domain,[]).append(new)



for domain,seqs in out.items():

    SeqIO.write(
    seqs,
    f"chromatin_domain_results/domain_sequences/{domain}.fasta",
    "fasta"
    )



PY



##################################################
# BUILD DOMAIN HMMs
##################################################

for f in ${OUT}/domain_sequences/*.fasta

do

domain=$(basename $f .fasta)

nseq=$(grep -c "^>" $f)

if [ $nseq -lt 2 ]
then
    echo "Skipping ${domain} (only ${nseq} sequence)"
    continue
fi

mafft \
--auto \
$f \
> ${OUT}/domain_sequences/${domain}.aln.fasta


hmmbuild \
${OUT}/domain_hmms/${domain}.hmm \
${OUT}/domain_sequences/${domain}.aln.fasta


done




##################################################
# SEARCH FUNGAL DISTANT RELATIVES
##################################################

FUNGAL_FASTA=${OUT}/metadata/fungal_proteomes.fasta


> ${FUNGAL_FASTA}


for f in \
${BASE}/ncr_proteome.fasta \
${BASE}/fgr_proteome.fasta \
${BASE}/mgr_proteome.fasta \
${BASE}/cne_proteome.fasta \
${BASE}/zt_proteome.fasta


do

species=$(basename ${f} _proteome.fasta)

sed "s/^>/>${species}|/" ${f} >> ${FUNGAL_FASTA}


done



for hmm in ${OUT}/domain_hmms/*.hmm

do


domain=$(basename ${hmm} .hmm)


hmmsearch \
--cpu ${SLURM_CPUS_PER_TASK} \
--domtblout ${OUT}/metadata/${domain}_fungal_hits.domtbl \
-E 1e-3 \
${hmm} \
${FUNGAL_FASTA}



done




##################################################
# MERGE ALL DOMAIN HITS
##################################################

python3 <<'PY'

import csv
from pathlib import Path


base=Path(
"chromatin_domain_results"
)


out=base/"metadata/domain_occurrences_extended.tsv"


rows=[]


with open(base/"metadata/domain_occurrences.tsv") as f:

    rows=list(csv.reader(f,delimiter="\t"))



for domtbl in Path(base/"metadata").glob("*_fungal_hits.domtbl"):


    domain=domtbl.name.replace("_fungal_hits.domtbl","")


    for line in open(domtbl):


        if line.startswith("#"):
            continue


        x=line.split()


        if len(x)<23:
            continue


        prot=x[0]

        parts=prot.split("|")


        rows.append([

        parts[0],
        parts[2],
        parts[3],
        domain,
        x[17],
        x[18],
        x[6]

        ])



with open(out,"w") as f:

    w=csv.writer(f,delimiter="\t")

    w.writerow([
    "species",
    "accession",
    "gene",
    "domain",
    "start",
    "end",
    "evalue"
    ])

    w.writerows(rows)



print("Extended hits:",len(rows))

PY




##################################################
# DOMAIN ARCHITECTURE
##################################################

python3 <<'PY'


from collections import defaultdict


d=defaultdict(list)


for line in open(
"chromatin_domain_results/metadata/domain_occurrences_extended.tsv"
):


    if line.startswith("species"):
        continue


    s,a,g,dom,*_=line.strip().split("\t")


    d[(s,g)].append(dom)



with open(
"chromatin_domain_results/metadata/domain_architecture.tsv",
"w"
) as out:


    out.write("species\tgene\tarchitecture\n")


    for k,v in d.items():

        out.write(
        f"{k[0]}\t{k[1]}\t{';'.join(sorted(set(v)))}\n"
        )


PY




##################################################
# EXTRACT STRUCTURAL DOMAINS
##################################################

python3 <<'PY'


import csv
from pathlib import Path
import re


ROOT=Path(
"/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek"
)


CIFS=list(ROOT.glob("cif_*"))

print("Indexing CIF files...")

cif_index={}

for c in CIFS:
    for f in c.rglob("*.cif"):
        cif_index[f.name]=f

print("Indexed",len(cif_index),"CIFs")

out=ROOT/"chromatin_domain_results/domains_extracted"


out.mkdir(exist_ok=True)



for row in csv.DictReader(
open(ROOT/"chromatin_domain_results/metadata/domain_occurrences_extended.tsv"),
delimiter="\t"
):


    accession=row["accession"]

    if not accession:
        continue


        cif=None

        for name,path in cif_index.items():

            if accession in name:
                cif=path
                break

    if cif is None:
        continue



    atoms=[]


    for line in open(cif):

        if not line.startswith("ATOM"):
            continue


        x=line.split()


        if len(x)<12:
            continue


        if x[3]!="CA":
            continue


        try:

            res=int(x[8])

            if int(row["start"]) <= res <= int(row["end"]):

                atoms.append(line)


        except:
            pass



    if len(atoms)<10:
        continue



    domain=row["domain"]


    d=out/domain

    d.mkdir(exist_ok=True)



    name=f"{row['gene']}_{domain}_{row['start']}-{row['end']}.pdb"


    with open(d/name,"w") as f:

        f.writelines(atoms)

        f.write("END\n")


print("done")


PY




##################################################
# FOLDSEEK
##################################################

for d in ${OUT}/domains_extracted/*

do


domain=$(basename ${d})


n=$(find ${d} -name "*.pdb" | wc -l)

if [ ${n} -lt 2 ]
then
continue
fi



foldseek easy-search \
${d} \
${d} \
${OUT}/foldseek/${domain}_allvall.tsv \
${OUT}/tmp/${domain} \
--format-output \
"query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,rmsd" \
--threads ${SLURM_CPUS_PER_TASK}


done



echo "COMPLETE"
echo ${OUT}
