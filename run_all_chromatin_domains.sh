#!/bin/bash
#SBATCH --job-name=chromatin_domains
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
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



##################################################
# PATHS
##################################################

BASE=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek

OUT=${BASE}/chromatin_domain_results

mkdir -p ${OUT}/{metadata,domain_lists,domains,foldseek,tmp}





##################################################
# PFAM DATABASE
##################################################

PFAM=${OUT}/Pfam-A.hmm


if [ ! -f ${PFAM}.h3i ]

then

wget -q \
-O ${PFAM}.gz \
https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip ${PFAM}.gz

hmmpress ${PFAM}

fi





##################################################
# COMBINE ALL PROTEOMES
##################################################

echo "Combining proteomes"


ALL_FASTA=${OUT}/metadata/all_species.fasta


> ${ALL_FASTA}


for f in ${BASE}/*_proteome.fasta

do

species=$(basename ${f} _proteome.fasta)

echo "Adding ${species}"


sed "s/^>/>${species}|/" ${f} >> ${ALL_FASTA}


done





##################################################
# HMMER
##################################################

echo "Running hmmscan"



hmmscan \
--cpu ${SLURM_CPUS_PER_TASK} \
--domtblout ${OUT}/metadata/pfam_hits.domtbl \
-E 1e-3 \
--domE 1e-3 \
${PFAM} \
${ALL_FASTA}





##################################################
# PARSE HMMER DOMAINS
##################################################

echo "Parsing chromatin domains"



python3 <<'PY'


import csv



DOMAINS={


"BAH":"BAH",

"Chromo":"Chromo",

"Chromo_shadow":"ChromoShadow",

"MRG":"MRG",

"PWWP":"PWWP",

"PHD":"PHD",

"Tudor":"Tudor",

"MBT":"MBT",

"Bromodomain":"Bromodomain",


"ADD":"ADD",

"zf-CXXC":"CXXC",

"CW":"CW",



"SET":"SET",

"DOT1":"DOT1",


"HAT":"HAT_MYST",

"GNAT":"HAT_GNAT",



"JmjC":"JmjC",

"JmjN":"JmjN",



"WD40":"WD40"


}



inp="chromatin_domain_results/metadata/pfam_hits.domtbl"


out="chromatin_domain_results/metadata/domain_occurrences.tsv"



hits=[]



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



    species=parts[0]



    accession=parts[2]



    gene=parts[3]



    hits.append([

    species,

    accession,

    gene,

    DOMAINS[pfam],

    x[19],

    x[20],

    x[11]

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



    w.writerows(hits)



print("Domain hits:",len(hits))

print("Proteins:",len(set(x[2] for x in hits)))



PY





##################################################
# CREATE DOMAIN PROTEIN LISTS
##################################################

echo "Creating protein lists"



for d in $(cut -f4 ${OUT}/metadata/domain_occurrences.tsv | tail -n +2 | sort -u)

do


grep -w ${d} \
${OUT}/metadata/domain_occurrences.tsv \
| cut -f3 \
| sort -u \
> ${OUT}/domain_lists/${d}.txt



done





##################################################
# DOMAIN ARCHITECTURE
##################################################

python3 <<'PY'


from collections import defaultdict



d=defaultdict(list)



for line in open(
"chromatin_domain_results/metadata/domain_occurrences.tsv"
):


    if line.startswith("species"):
        continue



    s,a,g,dom,*_=line.strip().split("\t")


    d[(s,g)].append(dom)



with open(
"chromatin_domain_results/metadata/domain_architectures.tsv",
"w"
) as out:


    out.write(
    "species\tgene\tarchitecture\n"
    )


    for (s,g),v in d.items():


        out.write(

        f"{s}\t{g}\t{';'.join(sorted(set(v)))}\n"

        )


PY





##################################################
# FIND CIFS
##################################################

echo "Linking structures"



while IFS=$'\t' read species accession gene domain start end evalue

do


if [ "$gene" = "gene" ]; then
continue
fi



mkdir -p ${OUT}/domains/${domain}/${species}



find ${BASE}/cif_* \
-name "*${accession}*.cif" \
-exec ln -sf {} \
${OUT}/domains/${domain}/${species}/ \;


done < ${OUT}/metadata/domain_occurrences.tsv




##################################################
# EXTRACT DOMAIN STRUCTURES FROM HMMER BOUNDARIES
##################################################

echo "Extracting domain-only structures"


DOMAIN_PDB=${OUT}/domains_extracted

mkdir -p ${DOMAIN_PDB}



python3 <<'PY'


import csv
from pathlib import Path
import re


BASE=Path(
"chromatin_domain_results"
)

MANIFEST_DIR=BASE/"metadata"
DOMAIN_TABLE=MANIFEST_DIR/"domain_occurrences.tsv"

CIF_ROOT=BASE/"domains"

OUT=BASE/"domains_extracted"



def read_cif_ca(path):

    atoms=[]

    lines=open(path).read().splitlines()

    atom_started=False
    cols={}

    for i,line in enumerate(lines):

        if line.startswith("_atom_site."):

            atom_started=True
            key=line.strip().replace("_atom_site.","")
            cols[key]=len(cols)


        elif atom_started and line.startswith("ATOM"):

            x=line.split()

            try:

                atom_name=x[cols["label_atom_id"]]

                if atom_name != "CA":
                    continue


                residue=int(
                    x[cols["label_seq_id"]]
                )


                atoms.append(
                    (
                    residue,
                    float(x[cols["Cartn_x"]]),
                    float(x[cols["Cartn_y"]]),
                    float(x[cols["Cartn_z"]]),
                    x[cols["label_comp_id"]]
                    )
                )


            except:
                continue


    return atoms





def write_pdb(atoms,out):


    lines=[]


    for i,a in enumerate(atoms,1):

        res,x,y,z,name=a


        lines.append(
        f"ATOM  {i:5d}  CA  {name:3s} A{res:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )


    lines.append("END")


    out.write_text(
        "\n".join(lines)
    )





written=0


with open(DOMAIN_TABLE) as f:


    for row in csv.DictReader(f,delimiter="\t"):


        gene=row["gene"]
        domain=row["domain"]

        start=int(row["start"])
        end=int(row["end"])



        # find CIF
        matches=list(
            CIF_ROOT.rglob("*.cif")
        )

        matches=[
            m for m in matches
            if row["gene"] in m.name
            or row["accession"] in m.name
        ]


        if not matches:
            continue



        cif=matches[0]



        atoms=read_cif_ca(cif)



        domain_atoms=[
            a for a in atoms
            if start-5 <= a[0] <= end+5
        ]



        if len(domain_atoms)<10:
            continue



        outdir=OUT/domain
        outdir.mkdir(
            parents=True,
            exist_ok=True
        )


        safe=re.sub(
            "[^A-Za-z0-9_-]",
            "_",
            gene
        )


        outfile=outdir / (
        f"{safe}_{domain}_{start}-{end}.pdb"
        )


        write_pdb(
            domain_atoms,
            outfile
        )


        written+=1



print(
"Domain structures written:",
written
)

if written == 0:
    raise RuntimeError(
    "No domain PDBs extracted. Check CIF parser or accession matching."
    )
PY




##################################################
# FOLDSEEK DOMAIN-ONLY COMPARISONS
##################################################

echo "Running FoldSeek on extracted domains"



for d in ${DOMAIN_PDB}/*

do


domain=$(basename ${d})


n=$(find ${d} -name "*.pdb" | wc -l)


if [ ${n} -lt 2 ]

then

echo "Skip ${domain}: ${n}"

continue

fi



echo "FoldSeek ${domain}"



foldseek easy-search \
${d} \
${d} \
${OUT}/foldseek/${domain}_domain_allvall.tsv \
${OUT}/tmp/${domain}_domain \
--format-output \
"query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,rmsd" \
--alignment-type 1 \
--threads ${SLURM_CPUS_PER_TASK}



done



echo "DOMAIN LEVEL FOLDSEEK COMPLETE"

echo ${OUT}
