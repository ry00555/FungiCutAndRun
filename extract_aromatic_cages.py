#!/usr/bin/env python3
"""
extract_aromatic_cages.py
--------------------------
Extracts aromatic cage residues from chromatin reader domain structures.
These residues form the methyl-lysine binding pocket and their conservation
reflects functional conservation of mark readout.

Aromatic cage definitions (literature-based):
  Chromo/EAF3 : 3-residue cage (Y/W/F triad), positions relative to domain
  BAH/EPR-1   : aromatic cage in binding cleft
  WD40/EED    : Y/W cage (EED: Y148, Y306, W364, Y365 in human)
  PHD         : aromatic cage for H3K4me3 (W/Y/F triad)
  Tudor       : aromatic cage (4 residues)
  PWWP        : aromatic cage (W residue + 2-3 supporting aromatics)
  MRG/EAF3   : Chromo-like aromatic cage

For each domain instance, scans for aromatic residues (F, Y, W) within
the domain boundary and clusters them by proximity to identify cage.

Output: aromatic_cage_residues.csv
  gene, organism, domain, pfam_start, pfam_end,
  cage_residue_pos, cage_residue_aa, cage_x, cage_y, cage_z,
  n_aromatic_in_domain, cage_density
"""
import csv, re
from pathlib import Path
from itertools import combinations
import numpy as np

SCRATCH  = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
BOUNDS   = SCRATCH / "pfam/domain_boundaries.csv"
MANIFEST = SCRATCH / "pdb_lists/structure_manifest.tsv"
OUT_CSV  = SCRATCH / "results/aromatic_cage_residues.csv"
OUT_SUMMARY = SCRATCH / "results/aromatic_cage_summary.csv"

# Aromatic amino acids that form the cage
AROMATICS = {"PHE", "TYR", "TRP", "F", "Y", "W"}

# Domain types to analyze
TARGET_DOMAINS = {"Chromo", "BAH", "WD40", "PHD", "Tudor", "PWWP", "MRG"}

# Load gene -> structure path
gene_to_struct = {}
with open(MANIFEST) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        src  = row.get("structure_file","").strip()
        gene = row.get("protein_gene_name","").strip()
        org  = row.get("organism","").strip()
        mod  = row.get("modification_read","").strip()
        if src and src != "MISSING" and gene:
            p = Path(src)
            if p.exists():
                gene_to_struct[gene] = (p, org, mod)

def parse_all_atoms(path):
    """Parse all heavy atoms from PDB or CIF.
    Returns list of (resnum, resname, atomname, x, y, z, chain)"""
    atoms = []
    if path.suffix == ".cif":
        col_map, cols = {}, []
        for line in path.read_text().splitlines():
            if line.startswith("_atom_site."):
                col = line.strip().lstrip("_atom_site.")
                col_map[col] = len(cols); cols.append(col)
            elif cols and line.startswith("ATOM"):
                parts = line.split()
                if len(parts) < len(cols): continue
                try:
                    resname  = parts[col_map.get("label_comp_id", 5)]
                    atomname = parts[col_map.get("label_atom_id", 3)].strip('"')
                    resnum   = int(parts[col_map.get("label_seq_id", 8)])
                    chain    = parts[col_map.get("label_asym_id", 6)]
                    x = float(parts[col_map.get("Cartn_x", 10)])
                    y = float(parts[col_map.get("Cartn_y", 11)])
                    z = float(parts[col_map.get("Cartn_z", 12)])
                    atoms.append((resnum, resname, atomname, x, y, z, chain))
                except (ValueError, IndexError, KeyError): continue
    else:
        for line in path.read_text().splitlines():
            if not line.startswith("ATOM"): continue
            try:
                atomname = line[12:16].strip()
                resname  = line[17:20].strip()
                chain    = line[21]
                resnum   = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((resnum, resname, atomname, x, y, z, chain))
            except ValueError: continue
    return atoms

def get_aromatic_residues_in_domain(atoms, env_start, env_end, padding=5):
    """Get all aromatic residues within domain boundaries."""
    seen = {}
    for resnum, resname, atomname, x, y, z, chain in atoms:
        if chain != "A": continue
        if not (env_start - padding <= resnum <= env_end + padding): continue
        if resname not in {"PHE","TYR","TRP"}: continue
        # Use CA or ring atom for position
        if atomname in {"CA","CZ","CE1","NE1","CG"}:
            if resnum not in seen:
                seen[resnum] = (resname, x, y, z)
    return seen  # {resnum: (resname, x, y, z)}

def find_aromatic_cluster(aromatic_residues, max_dist=8.0, min_cage=2):
    """
    Find the densest cluster of aromatic residues — the aromatic cage.
    Returns list of (resnum, resname, x, y, z) for cage residues.
    """
    if len(aromatic_residues) < min_cage:
        # Return in same format as main path: (resnum, resname, x, y, z)
        return [(k, v[0], v[1], v[2], v[3]) for k, v in aromatic_residues.items()]

    residues = list(aromatic_residues.items())
    coords   = np.array([(v[1], v[2], v[3]) for _, v in residues])

    # Find residue with most neighbors within max_dist
    best_center = 0
    best_count  = 0
    for i, (_, v) in enumerate(residues):
        c = np.array([v[1], v[2], v[3]])
        dists = np.linalg.norm(coords - c, axis=1)
        count = np.sum(dists <= max_dist) - 1  # exclude self
        if count > best_count:
            best_count  = count
            best_center = i

    # Return all aromatics within max_dist of best center
    center_coord = coords[best_center]
    dists = np.linalg.norm(coords - center_coord, axis=1)
    cage = [(residues[i][0], residues[i][1][0],
             residues[i][1][1], residues[i][1][2], residues[i][1][3])
            for i in np.where(dists <= max_dist)[0]]
    return cage

# Process all domain boundaries
print(f"Loading domain boundaries...")
with open(BOUNDS) as f:
    bounds_rows = list(csv.DictReader(f))

print(f"  {len(bounds_rows)} domain instances")
target_rows = [r for r in bounds_rows if r["domain"] in TARGET_DOMAINS]
print(f"  {len(target_rows)} in target domains: {TARGET_DOMAINS}")

cage_records   = []
summary_records = []
errors = skipped = written = 0

for row in target_rows:
    gene      = row["gene"]
    domain    = row["domain"]
    env_start = int(row["env_start"])
    env_end   = int(row["env_end"])

    if gene not in gene_to_struct:
        skipped += 1
        continue

    struct_path, organism, modification = gene_to_struct[gene]

    try:
        atoms = parse_all_atoms(struct_path)
    except Exception as e:
        errors += 1
        continue

    aromatic_residues = get_aromatic_residues_in_domain(atoms, env_start, env_end)
    n_aromatic = len(aromatic_residues)
    domain_len = env_end - env_start + 1
    cage_density = n_aromatic / domain_len if domain_len > 0 else 0

    cage = find_aromatic_cluster(aromatic_residues)

    # Per-residue records
    for resnum, resname, x, y, z in cage:
        cage_records.append({
            "gene":               gene,
            "organism":           organism,
            "modification":       modification,
            "domain":             domain,
            "domain_start":       env_start,
            "domain_end":         env_end,
            "domain_length":      domain_len,
            "cage_resnum":        resnum,
            "cage_resname":       resname,
            "cage_x":             round(x, 3),
            "cage_y":             round(y, 3),
            "cage_z":             round(z, 3),
            "n_aromatic_domain":  n_aromatic,
            "cage_size":          len(cage),
            "cage_density":       round(cage_density, 4),
        })

    # Per-domain summary
    cage_comp = "".join(sorted(r for _, r, *_ in cage))
    summary_records.append({
        "gene":              gene,
        "organism":          organism,
        "modification":      modification,
        "domain":            domain,
        "domain_start":      env_start,
        "domain_end":        env_end,
        "domain_length":     domain_len,
        "n_aromatic_domain": n_aromatic,
        "cage_size":         len(cage),
        "cage_density":      round(cage_density, 4),
        "cage_composition":  cage_comp,  # e.g. "FWYY"
        "n_phe":             sum(1 for _, r, *_ in cage if r=="PHE"),
        "n_tyr":             sum(1 for _, r, *_ in cage if r=="TYR"),
        "n_trp":             sum(1 for _, r, *_ in cage if r=="TRP"),
    })
    written += 1

print(f"\nProcessed: {written} domain instances")
print(f"Skipped (no structure): {skipped}")
print(f"Errors: {errors}")

# Write outputs
with open(OUT_CSV, "w", newline="") as f:
    if cage_records:
        w = csv.DictWriter(f, fieldnames=list(cage_records[0].keys()))
        w.writeheader(); w.writerows(cage_records)
print(f"Cage residues -> {OUT_CSV}  ({len(cage_records)} rows)")

with open(OUT_SUMMARY, "w", newline="") as f:
    if summary_records:
        w = csv.DictWriter(f, fieldnames=list(summary_records[0].keys()))
        w.writeheader(); w.writerows(summary_records)
print(f"Cage summary  -> {OUT_SUMMARY}  ({len(summary_records)} rows)")
