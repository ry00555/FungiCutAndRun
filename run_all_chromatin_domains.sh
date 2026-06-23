#!/bin/bash
#SBATCH --job-name=chromatin_domains_all
#SBATCH --partition=inter_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=chromatin_domains_%j.out
#SBATCH --error=chromatin_domains_%j.err


set -euo pipefail


########################################################
# MODULES
########################################################

module load HMMER/3.4-gompi-2024a
module load Foldseek/10-941cd33-GPU



########################################################
# PATHS
########################################################

BASE=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek

OUT=${BASE}/chromatin_domain_results

PFAM=${OUT}/Pfam-A.hmm


mkdir -p \
${OUT}/metadata \
${OUT}/domain_lists \
${OUT}/domains \
${OUT}/foldseek \
${OUT}/tmp



########################################################
# SPECIES METADATA
########################################################


declare -A SPECIES


# Fungi

SPECIES[ncr]="Neurospora_crassa,Fungi"
SPECIES[cne]="Cryptococcus_neoformans,Fungi"
SPECIES[fgr]="Fusarium_graminearum,Fungi"
SPECIES[mgr]="Magnaporthe_oryzae,Fungi"
SPECIES[zt]="Zymoseptoria_tritici,Fungi"

SPECIES[yeast]="Saccharomyces_cerevisiae,Fungi"
SPECIES[schpo]="Schizosaccharomyces_pombe,Fungi"



# Metazoa

SPECIES[human]="Homo_sapiens,Metazoa"
SPECIES[mouse]="Mus_musculus,Metazoa"
SPECIES[zebrafish]="Danio_rerio,Metazoa"

SPECIES[drome]="Drosophila_melanogaster,Metazoa"
SPECIES[caeel]="Caenorhabditis_elegans,Metazoa"



# Plantae

SPECIES[arath]="Arabidopsis_thaliana,Plantae"




########################################################
# PFAM DOMAINS
########################################################


declare -A DOMAINS


DOMAINS[BAH]="PF01426"

DOMAINS[Chromo]="PF00385"

DOMAINS[ChromoShadow]="PF01393"

DOMAINS[MRG]="PF05521"

DOMAINS[PWWP]="PF00855"

DOMAINS[PHD]="PF00628"

DOMAINS[Tudor]="PF00567"

DOMAINS[MBT]="PF02820"

DOMAINS[Bromodomain]="PF00439"


DOMAINS[ADD]="PF17980"

DOMAINS[CXXC]="PF02008"

DOMAINS[CW]="PF07496"



DOMAINS[SET]="PF00856"

DOMAINS[DOT1]="PF08161"


DOMAINS[HAT_MYST]="PF01853"

DOMAINS[HAT_GNAT]="PF00583"



DOMAINS[JmjC]="PF02373"

DOMAINS[JmjN]="PF02375"



DOMAINS[WD40]="PF00400"





########################################################
# DOWNLOAD PFAM
########################################################


echo "Preparing Pfam database"



if [ ! -f ${PFAM}.h3i ]

then

wget -q \
-O ${PFAM}.gz \
https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz


gunzip ${PFAM}.gz


hmmpress ${PFAM}


else

echo "Pfam already prepared"

fi





########################################################
# COMBINE ALL PROTEOMES
########################################################


echo "Combining all FASTA files"


ALL_FASTA=${OUT}/metadata/all_species.fasta


> ${ALL_FASTA}



for f in ${BASE}/*_proteome.fasta

do


species=$(basename ${f} _proteome.fasta)


echo "Adding ${species}"


sed "s/^>/>${species}|/" ${f} >> ${ALL_FASTA}



done





########################################################
# RUN HMMER
########################################################


echo "Running HMMER"


hmmscan \
--cpu ${SLURM_CPUS_PER_TASK} \
--domtblout ${OUT}/metadata/pfam_hits.domtbl \
-E 1e-3 \
--domE 1e-3 \
${PFAM} \
${ALL_FASTA}





########################################################
# PARSE DOMAIN HITS
########################################################


echo "Parsing domains"



python3 <<'EOF'


import csv



PFAM={


"PF01426":"BAH",

"PF00385":"Chromo",

"PF01393":"ChromoShadow",

"PF05521":"MRG",

"PF00855":"PWWP",

"PF00628":"PHD",

"PF00567":"Tudor",

"PF02820":"MBT",

"PF00439":"Bromodomain",


"PF17980":"ADD",

"PF02008":"CXXC",

"PF07496":"CW",


"PF00856":"SET",

"PF08161":"DOT1",


"PF01853":"HAT_MYST",

"PF00583":"HAT_GNAT",


"PF02373":"JmjC",

"PF02375":"JmjN",


"PF00400":"WD40"


}



inp="chromatin_domain_results/metadata/pfam_hits.domtbl"


out="chromatin_domain_results/metadata/domain_occurrences.tsv"



rows=[]



with open(inp) as f:


    for line in f:


        if line.startswith("#"):
            continue


        x=line.split()


        if len(x)<23:
            continue



        pfam=x[1].split(".")[0]


        if pfam not in PFAM:
            continue



        protein=x[3]


        species=protein.split("|")[0]

        gene=protein.split("|")[1]



        rows.append([

        species,

        gene,

        PFAM[pfam],

        x[19],

        x[20],

        x[6]

        ])




with open(out,"w") as f:


    w=csv.writer(f,delimiter="\t")


    w.writerow([

    "species",

    "protein",

    "domain",

    "start",

    "end",

    "evalue"

    ])



    w.writerows(rows)



print("domain hits:",len(rows))


EOF





########################################################
# DOMAIN PROTEIN LISTS
########################################################


echo "Creating domain protein lists"



cut -f3 ${OUT}/metadata/domain_occurrences.tsv \
| tail -n +2 \
| sort -u \
| while read domain


do


grep -w ${domain} \
${OUT}/metadata/domain_occurrences.tsv \
| cut -f2 \
| sort -u \
> ${OUT}/domain_lists/${domain}.txt



done





########################################################
# DOMAIN ARCHITECTURES
########################################################


python3 <<'EOF'


from collections import defaultdict


d={}



with open(
"chromatin_domain_results/metadata/domain_occurrences.tsv"
) as f:


    next(f)


    for line in f:


        s,p,dom,*_=line.strip().split("\t")


        d.setdefault(
        (s,p),
        []
        ).append(dom)



with open(
"chromatin_domain_results/metadata/domain_architectures.tsv",
"w"
) as out:


    out.write(
    "species\tprotein\tdomains\n"
    )


    for (s,p),v in d.items():


        out.write(

        f"{s}\t{p}\t{';'.join(sorted(set(v)))}\n"

        )


EOF





########################################################
# LINK CIF STRUCTURES BY DOMAIN
########################################################


echo "Finding CIF structures"



while IFS=$'\t' read species protein domain start end evalue

do


if [ "$protein" == "protein" ]; then
continue
fi



mkdir -p ${OUT}/domains/${domain}/${species}



find ${BASE}/cif_* \
-name "*${protein}*.cif" \
-exec ln -sf {} \
${OUT}/domains/${domain}/${species}/ \;



done < ${OUT}/metadata/domain_occurrences.tsv





########################################################
# FOLDSEEK PER DOMAIN CLASS
########################################################


echo "Running FoldSeek"



for d in ${OUT}/domains/*


do


domain=$(basename $d)


n=$(find ${d} -name "*.cif" | wc -l)



if [ ${n} -lt 2 ]

then

echo "Skipping ${domain}: ${n} structures"

continue

fi




echo "FoldSeek ${domain}"



foldseek easy-search \
${d} \
${d} \
${OUT}/foldseek/${domain}_allvall.tsv \
${OUT}/tmp/${domain} \
--format-output \
"query,target,evalue,bits,alntmscore,rmsd" \
--threads ${SLURM_CPUS_PER_TASK}



done





echo "====================================="
echo " COMPLETE"
echo " Results:"
echo ${OUT}
echo "====================================="
