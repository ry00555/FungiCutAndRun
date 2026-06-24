#!/usr/bin/env python3
"""
build_master_foldseek_metadata.py  — FIXED VERSION
----------------------------------------------------
Reads the outputs of run_all_chromatin_domains.sh and builds a fully
annotated FoldSeek hit table.

Fixes from original:
  1. FoldSeek TSVs read with header=None (they have no header row)
  2. Column count verified against actual FoldSeek format-output string
  3. FASTA parser handles both sp| (Swiss-Prot) and tr| (TrEMBL) entries
     so fungal proteins get metadata
  4. PDB filename parsing handles both naming conventions:
       GENE_SPECIES_ACC_DOMAIN_COORDS.pdb  (5 parts)
       GENE_ACC_DOMAIN_START-END.pdb       (4 parts)
  5. Species/organism derived from FASTA when missing from PDB filename
  6. Multi-domain proteins flagged correctly
  7. Cross-kingdom flag added
"""

import os, glob, re
import pandas as pd
from pathlib import Path

BASE = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek")
OUT  = BASE / "chromatin_domain_results"

FASTA      = OUT / "metadata/all_species.fasta"
PDB_DIR    = OUT / "domains_extracted"
FOLDSEEK_DIR = OUT / "foldseek"
OUTFILE    = BASE / "annotated_hits_expanded.csv"

# ── Kingdom mapping ───────────────────────────────────────────────────
KINGDOM = {
    "Homo sapiens":               "Metazoa",
    "Mus musculus":               "Metazoa",
    "Danio rerio":                "Metazoa",
    "Drosophila melanogaster":    "Metazoa",
    "Caenorhabditis elegans":     "Metazoa",
    "Arabidopsis thaliana":       "Plantae",
    "Saccharomyces cerevisiae":   "Fungi",
    "Schizosaccharomyces pombe":  "Fungi",
    "Neurospora crassa":          "Fungi",
    "Fusarium graminearum":       "Fungi",
    "Magnaporthe oryzae":         "Fungi",
    "Zymoseptoria tritici":       "Fungi",
    "Cryptococcus neoformans":    "Fungi",
}

# ── Step 1: Parse FASTA metadata ─────────────────────────────────────
# Handles both:
#   >sp|ACC|NAME Description OS=Organism OX=... GN=Gene PE=... SV=...
#   >tr|ACC|NAME Description OS=Organism OX=... GN=Gene PE=... SV=...
#   >species|sp|ACC|NAME ... (species-prefixed format from the pipeline)
print("Reading FASTA metadata...")
fasta_rows = []
with open(FASTA) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line.strip()

        # Extract accession — handles both |sp| and |tr|
        acc_m = re.search(r"\|(?:sp|tr)\|([^|]+)", header)
        if not acc_m:
            # Try plain accession at start
            acc_m = re.search(r">(?:\w+\|)?(\w+)\|", header)
        if not acc_m:
            continue
        accession = acc_m.group(1)

        gene_m = re.search(r"\bGN=([^\s]+)", header)
        org_m  = re.search(r"\bOS=(.*?)\s+OX=", header)
        desc_m = re.search(r"[A-Z0-9]+\s+(.*?)\s+OS=", header)

        organism = org_m.group(1).strip() if org_m else "unknown"

        fasta_rows.append({
            "accession":   accession,
            "gene":        gene_m.group(1) if gene_m else accession,
            "description": desc_m.group(1).strip() if desc_m else "",
            "organism":    organism,
            "kingdom":     KINGDOM.get(organism, "unknown"),
        })

fasta_meta = pd.DataFrame(fasta_rows).drop_duplicates("accession")
print(f"  FASTA proteins: {len(fasta_meta)}")
print(f"  Organisms: {sorted(fasta_meta['organism'].unique())[:5]}...")

# ── Step 2: Parse domain PDB filenames ───────────────────────────────
# Two possible naming conventions from the extraction script:
#   GENE_SPECIES_ACC_DOMAIN_START-END.pdb  (5 underscore-parts)
#   GENE_ACC_DOMAIN_START-END.pdb          (4 underscore-parts)
print("\nReading domain PDB files...")
domain_rows = []
for pdb in glob.glob(str(PDB_DIR / "**" / "*.pdb"), recursive=True):
    filename = os.path.basename(pdb).replace(".pdb", "")
    parts    = filename.split("_")

    if len(parts) >= 5:
        # GENE_SPECIES_ACC_DOMAIN_COORDS
        gene      = parts[0]
        species   = parts[1]
        accession = parts[2]
        domain    = parts[3]
        coords    = parts[4]
    elif len(parts) == 4:
        # GENE_ACC_DOMAIN_COORDS
        gene      = parts[0]
        species   = ""
        accession = parts[1]
        domain    = parts[2]
        coords    = parts[3]
    else:
        continue

    domain_rows.append({
        "foldseek_id":   filename,
        "gene_from_pdb": gene,
        "species_code":  species,
        "accession":     accession,
        "domain":        domain,
        "coords":        coords,
        "pdb_path":      pdb,
    })

domain_map = pd.DataFrame(domain_rows)
print(f"  Domain structures: {len(domain_map)}")

# Merge FASTA metadata onto domain map
domain_map = domain_map.merge(fasta_meta, on="accession", how="left")

# Fill gene from FASTA where PDB name was just accession
mask = domain_map["gene_from_pdb"] == domain_map["accession"]
domain_map.loc[mask, "gene_from_pdb"] = domain_map.loc[mask, "gene"]

# Flag multi-domain proteins
domain_counts = domain_map.groupby("accession")["domain"].nunique().reset_index()
domain_counts.columns = ["accession", "n_domains"]
all_domains  = domain_map.groupby("accession")["domain"].apply(
    lambda x: ";".join(sorted(set(x)))).reset_index()
all_domains.columns = ["accession", "all_domains"]
domain_map = domain_map.merge(domain_counts, on="accession", how="left")
domain_map = domain_map.merge(all_domains,   on="accession", how="left")
domain_map["is_multidomain"] = domain_map["n_domains"] > 1

print(f"  Multi-domain proteins: {domain_map[domain_map['is_multidomain']]['accession'].nunique()}")
print(domain_map[["accession","gene_from_pdb","domain","organism","kingdom"]].head(3))

# ── Step 3: Read FoldSeek TSVs ────────────────────────────────────────
# format-output used in run_all_chromatin_domains.sh:
#   "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,rmsd"
# That is 14 columns — read with header=None
FOLDSEEK_COLS = [
    "query_id","target_id","fident","alnlen",
    "qstart","qend","tstart","tend",
    "evalue","bits","alntmscore","qtmscore","ttmscore","rmsd"
]

print("\nReading FoldSeek TSVs...")
tsvs = glob.glob(str(FOLDSEEK_DIR / "*_allvall.tsv"))
print(f"  Found {len(tsvs)} TSV files")

all_hits = []
for tsv in sorted(tsvs):
    domain = os.path.basename(tsv).replace("_allvall.tsv","")
    try:
        df = pd.read_csv(tsv, sep="\t", header=None, names=FOLDSEEK_COLS)
        df["domain"] = domain
        # Compute mean TM-score
        df["mean_tmscore"] = ((df["qtmscore"].astype(float) +
                               df["ttmscore"].astype(float)) / 2).round(4)
        all_hits.append(df)
        print(f"  {domain}: {len(df):,} hits")
    except Exception as e:
        print(f"  FAILED {tsv}: {e}")

if not all_hits:
    print("ERROR: No FoldSeek TSVs found or all failed to load")
    print(f"  Checked: {FOLDSEEK_DIR}")
    raise SystemExit(1)

foldseek = pd.concat(all_hits, ignore_index=True)
print(f"\n  Total FoldSeek comparisons: {len(foldseek):,}")

# ── Step 4: Annotate query and target ────────────────────────────────
print("\nAnnotating query and target proteins...")

# Build lookup: foldseek_id -> metadata row
lookup = domain_map.set_index("foldseek_id")

def get_meta(foldseek_id: str, prefix: str) -> dict:
    """Look up metadata for a FoldSeek query or target ID."""
    row = lookup.loc[foldseek_id] if foldseek_id in lookup.index else None
    if row is None:
        return {
            f"{prefix}gene":        foldseek_id,
            f"{prefix}accession":   "",
            f"{prefix}domain":      "",
            f"{prefix}organism":    "unknown",
            f"{prefix}kingdom":     "unknown",
            f"{prefix}is_multidomain": False,
            f"{prefix}all_domains": "",
            f"{prefix}coords":      "",
        }
    # Handle case where lookup returns multiple rows (same foldseek_id)
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    return {
        f"{prefix}gene":         row.get("gene_from_pdb", ""),
        f"{prefix}accession":    row.get("accession", ""),
        f"{prefix}domain":       row.get("domain", ""),
        f"{prefix}organism":     row.get("organism", "unknown"),
        f"{prefix}kingdom":      row.get("kingdom", "unknown"),
        f"{prefix}is_multidomain": row.get("is_multidomain", False),
        f"{prefix}all_domains":  row.get("all_domains", ""),
        f"{prefix}coords":       row.get("coords", ""),
    }

# Apply metadata — vectorised where possible, fallback to row-by-row
hits_meta = foldseek.copy()

# Merge query metadata
q_meta = domain_map[["foldseek_id","gene_from_pdb","accession","domain",
                       "organism","kingdom","is_multidomain","all_domains",
                       "coords"]].copy()
q_meta.columns = ["query_id","query_gene","query_accession","query_domain",
                   "query_organism","query_kingdom","query_is_multidomain",
                   "query_all_domains","query_coords"]
q_meta = q_meta.drop_duplicates("query_id")

t_meta = domain_map[["foldseek_id","gene_from_pdb","accession","domain",
                       "organism","kingdom","is_multidomain","all_domains",
                       "coords"]].copy()
t_meta.columns = ["target_id","target_gene","target_accession","target_domain",
                   "target_organism","target_kingdom","target_is_multidomain",
                   "target_all_domains","target_coords"]
t_meta = t_meta.drop_duplicates("target_id")

hits_meta = hits_meta.merge(q_meta, on="query_id",  how="left")
hits_meta = hits_meta.merge(t_meta, on="target_id", how="left")

# ── Step 5: Add classification columns ────────────────────────────────
hits_meta["same_domain"]   = hits_meta["query_domain"] == hits_meta["target_domain"]
hits_meta["cross_species"] = hits_meta["query_organism"] != hits_meta["target_organism"]
hits_meta["cross_kingdom"] = hits_meta["query_kingdom"] != hits_meta["target_kingdom"]
hits_meta["is_self_hit"]   = hits_meta["query_id"] == hits_meta["target_id"]

# Confidence bin
def tm_bin(tm):
    try:
        v = float(tm)
        if v > 0.6: return "high(>0.6)"
        if v > 0.5: return "medium(0.5-0.6)"
        return "low(0.4-0.5)"
    except: return "unknown"
hits_meta["tm_confidence"] = hits_meta["mean_tmscore"].apply(tm_bin)

# ── Step 6: Save ──────────────────────────────────────────────────────
hits_meta.to_csv(str(OUTFILE), index=False)
print(f"\nSaved: {OUTFILE}  ({len(hits_meta):,} rows)")

# ── Summary ───────────────────────────────────────────────────────────
print("\n" + "="*65)
print("SUMMARY")
print("="*65)
print(f"Total hits          : {len(hits_meta):,}")
nself = hits_meta["is_self_hit"].sum()
print(f"Self hits           : {nself:,}")
hits_ns = hits_meta[~hits_meta["is_self_hit"]]
print(f"Non-self hits       : {len(hits_ns):,}")
print(f"TM > 0.5            : {(hits_ns['mean_tmscore']>0.5).sum():,}")
print(f"TM > 0.6            : {(hits_ns['mean_tmscore']>0.6).sum():,}")
print(f"Cross-kingdom TM>0.5: {((hits_ns['cross_kingdom']) & (hits_ns['mean_tmscore']>0.5)).sum():,}")
print(f"Cross-domain TM>0.5 : {((~hits_ns['same_domain']) & (hits_ns['mean_tmscore']>0.5)).sum():,}")

print("\nHits per domain:")
for dom, n in hits_ns.groupby("domain").size().sort_values(ascending=False).items():
    print(f"  {dom:<20} {n:>8,}")

print("\nTop cross-kingdom structural hits (TM>0.5, non-self):")
xk = hits_ns[(hits_ns["cross_kingdom"]) & (hits_ns["mean_tmscore"]>0.5)]
xk_top = xk.sort_values("mean_tmscore", ascending=False).head(20)
for _, r in xk_top.iterrows():
    print(f"  {r.get('query_gene','?'):<12} ({r.get('query_domain','?'):<10}) vs "
          f"{r.get('target_gene','?'):<12} ({r.get('target_domain','?'):<10}) "
          f"TM={r['mean_tmscore']:.3f}  "
          f"{r.get('query_organism','?').split()[0]} vs {r.get('target_organism','?').split()[0]}")

print("\nDone.")
