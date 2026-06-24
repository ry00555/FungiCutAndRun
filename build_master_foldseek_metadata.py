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
# 2. Parse extracted domain PDB files
#########################################

print("Reading domain PDBs...")


pdb_rows=[]


pdb_files = glob.glob(
    str(PDB_DIR / "**" / "*.pdb"),
    recursive=True
)


print("Found PDB files:",len(pdb_files))


for pdb in pdb_files:


    filename=os.path.basename(pdb)

    fid=filename.replace(".pdb","")


    # Example:
    # V5IPW0_NEUCR_V5IPW0_WD40_311-343

    m=re.match(
        r"(.+?)_([A-Za-z0-9_]+)_(\d+-\d+)$",
        fid
    )


    if not m:
        continue


    prefix=m.group(1)

    domain=m.group(2)

    coords=m.group(3)


    accession=prefix.split("_")[0]


    start,end=coords.split("-")


    pdb_rows.append({

        "foldseek_id":fid,

        "accession":accession,

        "domain":domain,

        "start":int(start),

        "end":int(end)

    })


domain_map=pd.DataFrame(pdb_rows)


print("Domain structures:",len(domain_map))


#########################################
# Merge FASTA annotations
#########################################

print("Merging FASTA annotations...")


domain_map=domain_map.merge(
    fasta_meta,
    on="accession",
    how="left"
)



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


query_map=domain_map.rename(
    columns={
        "foldseek_id":"query_id",
        "gene":"query_gene",
        "organism":"query_organism",
        "domain":"query_domain"
    }
)


target_map=domain_map.rename(
    columns={
        "foldseek_id":"target_id",
        "gene":"target_gene",
        "organism":"target_organism",
        "domain":"target_domain"
    }
)



hits=foldseek.merge(
    query_map[
        [
        "query_id",
        "query_gene",
        "query_organism",
        "query_domain"
        ]
    ],
    on="query_id",
    how="left"
)



hits=hits.merge(
    target_map[
        [
        "target_id",
        "target_gene",
        "target_organism",
        "target_domain"
        ]
    ],
    on="target_id",
    how="left"
)



#########################################
# 5. Add classifications
#########################################


hits["same_domain"] = (
    hits.query_domain ==
    hits.target_domain
)


hits["cross_species"] = (
    hits.query_organism !=
    hits.target_organism
)



#########################################
# 6. Save
#########################################


hits.to_csv(
    OUT,
    index=False
)


print("DONE")
print("Saved:",OUT)

print(hits.head())
