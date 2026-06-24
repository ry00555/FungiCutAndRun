#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd
from pathlib import Path


#########################################
# Paths
#########################################

BASE = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek")

FASTA = BASE / "chromatin_domain_results/metadata/all_species.fasta"

PDB_DIR = BASE / "chromatin_domain_results/domains_extracted"

FOLDSEEK_DIR = BASE / "chromatin_domain_results/foldseek"

OUT = BASE / "annotated_hits_expanded.csv"


#########################################
# 1. Parse FASTA metadata
#########################################

print("Reading FASTA metadata...")

fasta_rows=[]

with open(FASTA) as f:

    for line in f:

        if line.startswith(">"):

            header=line.strip()

            accession = re.search(r"\|sp\|([^|]+)",header)

            if accession:
                accession=accession.group(1)
            else:
                continue


            gene = re.search(r"GN=([^\s]+)",header)

            organism = re.search(r"OS=(.*?) OX=",header)

            description = header.split(" ",1)[1].split(" OS=")[0]


            fasta_rows.append({

                "accession":accession,

                "gene":
                    gene.group(1) if gene else accession,

                "description":
                    description,

                "organism":
                    organism.group(1) if organism else "unknown"

            })


fasta_meta=pd.DataFrame(fasta_rows)

print("FASTA proteins:",len(fasta_meta))


#########################################
# Parse domain extracted structures
#########################################

print("Reading domain PDBs...")

domain_rows=[]

for pdb in glob.glob(
    str(PDB_DIR / "**" / "*.pdb"),
    recursive=True
):

    filename=os.path.basename(pdb).replace(".pdb","")

    parts=filename.split("_")

    # Example:
    # JADE1_MOUSE_Q6ZPI0_PHD_207-251

    gene = parts[0]

    species = parts[1]

    accession = parts[2]

    domain = parts[3]

    coords = parts[4]


    domain_rows.append({

        "foldseek_id": filename,

        "gene_from_pdb": gene,

        "species_code": species,

        "accession": accession,

        "domain": domain,

        "coords": coords

    })


domain_map=pd.DataFrame(domain_rows)

print("Domain structures:", len(domain_map))


#########################################
# Merge FASTA metadata
#########################################

domain_map = domain_map.merge(
    fasta_meta,
    on="accession",
    how="left"
)

print(domain_map.head())

#########################################
# 3. Read FoldSeek all-vs-all results
#########################################

print("Reading FoldSeek TSVs...")


tsvs=glob.glob(
    str(FOLDSEEK_DIR / "*_allvall.tsv")
)


all_hits=[]


for tsv in tsvs:


    domain=os.path.basename(tsv).replace("_allvall.tsv","")


    df=pd.read_csv(
        tsv,
        sep="\t",
        header=None
    )


    df.columns=[

        "query_id",
        "target_id",
        "mean_tmscore",
        "aln_len",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits",
        "qtmscore",
        "ttmscore",
        "lddt",
        "rmsd"

    ]


    df["domain"]=domain


    all_hits.append(df)



foldseek=pd.concat(all_hits,ignore_index=True)


print("FoldSeek comparisons:",len(foldseek))


#########################################
# 4. Annotate query and target proteins
#########################################

print("Annotating query/target IDs...")


hits = foldseek.merge(
    domain_map.add_prefix("query_"),
    left_on="query_id",
    right_on="query_foldseek_id",
    how="left"
)


hits = hits.merge(
    domain_map.add_prefix("target_"),
    left_on="target_id",
    right_on="target_foldseek_id",
    how="left"
)



#########################################
# Rename final fields
#########################################

hits_final = hits.rename(columns={

    "query_gene_from_pdb":"query_gene",
    "target_gene_from_pdb":"target_gene",

    "query_organism":"query_organism",
    "target_organism":"target_organism",

    "query_domain":"query_domain",
    "target_domain":"target_domain"

})


#########################################
# Add classifications
#########################################

hits_final["same_domain"] = (
    hits_final["query_domain"] ==
    hits_final["target_domain"]
)


hits_final["cross_species"] = (
    hits_final["query_organism"] !=
    hits_final["target_organism"]
)



#########################################
# Save
#########################################

hits_final.to_csv(
    "/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/annotated_hits_expanded.csv",
    index=False
)


print(hits_final.head())
