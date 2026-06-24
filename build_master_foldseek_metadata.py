#!/usr/bin/env python3

import pandas as pd
import glob
import os
import re


#########################################
# PATHS
#########################################

from pathlib import Path


# FoldSeek data location
FOLDSEEK = Path(
    "/scratch/ry00555/RNASeqPaper2026/"
    "Proteome/StructuralSimilarity/FoldSeek"
)


FASTA = FOLDSEEK / "all_species.fasta"


DOMAIN_META = (
    FOLDSEEK /
    "chromatin_domain_results" /
    "metadata" /
    "domain_occurrences_extended.tsv"
)


PDB_DIR = (
    FOLDSEEK /
    "chromatin_domain_results" /
    "domains_extracted"
)


OUTPUT = (
    FOLDSEEK /
    "annotated_hits_expanded.csv"
)

#########################################
# 1. Parse all_species.fasta
#########################################

print("Reading FASTA metadata...")


fasta_rows=[]


with open(FASTA) as f:

    for line in f:

        if line.startswith(">"):

            h=line.strip()[1:]


            accession=re.search(
                r'\|([A-Z0-9]+)\|',
                h
            )

            entry=re.search(
                r'\|([A-Z0-9]+_[A-Z]+)',
                h
            )

            gene=re.search(
                r'GN=([^\s]+)',
                h
            )

            organism=re.search(
                r'OS=(.*?) OX=',
                h
            )

            description=h.split(" ",1)[1]


            fasta_rows.append({

                "accession":
                    accession.group(1)
                    if accession else None,

                "uniprot_entry":
                    entry.group(1)
                    if entry else None,

                "gene":
                    gene.group(1)
                    if gene else None,

                "organism":
                    organism.group(1)
                    if organism else None,

                "description":
                    description

            })


fasta_meta=pd.DataFrame(fasta_rows)



#########################################
# 2. Parse domain extracted structures
#########################################

print("Reading domain PDB IDs...")


domain_rows=[]


for pdb in glob.glob(
    PDB_DIR+"/**/*.pdb",
    recursive=True
):

    filename=os.path.basename(pdb)

    fid=filename.replace(".pdb","")


    # Example:
    # V5IPW0_NEUCR_V5IPW0_WD40_311-343

    parts=fid.split("_")


    accession=parts[0]

    domain=parts[-2]


    coords=parts[-1]


    domain_rows.append({

        "foldseek_id":fid,

        "accession":accession,

        "domain":domain,

        "coords":coords

    })



domain_map=pd.DataFrame(domain_rows)



#########################################
# 3. Merge FASTA metadata
#########################################

print("Merging FASTA annotations...")


domain_map=domain_map.merge(
    fasta_meta,
    on="accession",
    how="left"
)



#########################################
# 4. Add domain information
#########################################

print("Adding domain metadata...")


domains=pd.read_csv(
    DOMAIN_META,
    sep="\t"
)


domains=domains.rename(
columns={
"accession":"accession"
}
)


domain_map=domain_map.merge(
    domains[
        [
        "accession",
        "domain",
        "start",
        "end",
        "evalue"
        ]
    ],
    on=[
        "accession",
        "domain"
    ],
    how="left"
)



#########################################
# 5. Read ALL FoldSeek results
#########################################

print("Reading FoldSeek results...")


glob.glob(
    str(FOLDSEEK / "*allvall.tsv")
)



results=[]


for f in allvall_files:


    domain=os.path.basename(f)\
        .replace("_allvall.tsv","")


    x=pd.read_csv(
        f,
        sep="\t",
        header=None
    )


    x.columns=[

        "query_id",
        "target_id",

        "alnlen",

        "qstart",
        "qend",

        "tstart",
        "tend",

        "evalue",

        "bits",

        "qtmscore",

        "ttmscore",

        "lddt",

        "rmsd",

        "fident"

    ]


    x["query_domain"]=domain

    x["target_domain"]=domain


    results.append(x)



foldseek=pd.concat(results)



#########################################
# 6. Annotate query proteins
#########################################

print("Annotating queries...")


foldseek=foldseek.merge(

    domain_map,

    left_on="query_id",

    right_on="foldseek_id",

    how="left"

)


foldseek=foldseek.rename(

columns={

"gene":"query_gene",

"organism":"query_organism",

"accession":"query_accession",

"description":"query_description",

"uniprot_entry":"query_entry"

}

)



#########################################
# 7. Annotate targets
#########################################

print("Annotating targets...")


foldseek=foldseek.merge(

    domain_map,

    left_on="target_id",

    right_on="foldseek_id",

    how="left",

    suffixes=("","_target")

)



foldseek=foldseek.rename(

columns={

"gene":"target_gene",

"organism":"target_organism",

"accession":"target_accession",

"description":"target_description",

"uniprot_entry":"target_entry"

}

)



#########################################
# 8. Add biology fields
#########################################

print("Adding comparison labels...")


foldseek["mean_tmscore"]=(

    foldseek.qtmscore +

    foldseek.ttmscore

)/2



foldseek["same_domain"]=(

    foldseek.query_domain ==

    foldseek.target_domain

)



foldseek["same_modification"]=False



foldseek["cross_kingdom"]=(

    foldseek.query_organism !=

    foldseek.target_organism

)



# placeholder until your modification table

foldseek["query_modification"]="unknown"

foldseek["target_modification"]="unknown"



#########################################
# 9. Output
#########################################

final=foldseek[

[

"query_gene",

"target_gene",

"query_organism",

"target_organism",

"query_modification",

"target_modification",

"query_domain",

"target_domain",

"mean_tmscore",

"qtmscore",

"ttmscore",

"lddt",

"rmsd",

"fident",

"evalue",

"same_modification",

"same_domain",

"cross_kingdom",

"query_accession",

"target_accession",

"query_entry",

"target_entry",

"query_description",

"target_description"

]

]



final.to_csv(
    OUTPUT,
    index=False
)


print("\nDONE")
print(final.head())
print(
"Rows:",
len(final)
)
