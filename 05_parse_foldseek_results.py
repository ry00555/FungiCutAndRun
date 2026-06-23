#!/usr/bin/env python3
"""
05_parse_foldseek_results.py
-----------------------------
Parse FoldSeek all-vs-all output:
  1. Build a filename-stem -> gene_name reverse map from structure_manifest.tsv
     (e.g. "CBX7_Homo_AF3" -> "CBX7")
  2. Filter hits by TM-score / e-value thresholds
  3. Join with master annotated CSV metadata
  4. Build N×N TM-score similarity matrix
  5. Export annotated_hits.csv and tmscore_matrix.csv

Outputs (written to SCRATCH_DIR)
---------------------------------
  results/annotated_hits.csv
  results/tmscore_matrix.csv
  logs/05_parse.log

Usage:
    python 05_parse_foldseek_results.py          # uses all defaults
"""

import csv, sys, re, argparse, logging
from pathlib import Path
from collections import defaultdict

import numpy as np

SCRATCH = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")

FOLDSEEK_COLS = [
    "query","target","fident","alnlen","mismatch","gapopen",
    "qstart","qend","tstart","tend","evalue","bits",
    "lddt","alntmscore","qtmscore","ttmscore","rmsd","prob",
]


def setup_logging(logdir: Path) -> logging.Logger:
    logdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(logdir / "05_parse.log"),
        ],
    )
    return logging.getLogger(__name__)


def build_stem_to_gene(manifest_path: str, log) -> dict[str, str]:
    """
    Read structure_manifest.tsv and build a mapping from PDB filename stem
    to protein_gene_name. Also adds entries for chain/model variants that
    FoldSeek extracts from multi-model RCSB PDB files.

    Examples:
        PWWP1_Homo_AF.cif       -> "PWWP1_Homo_AF" -> "PWWP1"
        CBX7_Homo_RCSB_1x4p.pdb -> various chain variants -> "CBX7"
    """
    stem_to_gene: dict[str, str] = {}
    if not Path(manifest_path).exists():
        log.warning(f"Manifest not found: {manifest_path} — gene names will fall back to stem")
        return stem_to_gene
    with open(manifest_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row.get("protein_gene_name", "").strip()
            src  = row.get("structure_file", "").strip()
            if not gene or not src or src == "MISSING":
                continue
            stem = Path(src).stem
            stem_to_gene[stem] = gene
            # Also map chain/model variants FoldSeek extracts from multi-model PDBs:
            # e.g. CBX7_Homo_RCSB_1x4p_MODEL_1_A, CBX7_Homo_RCSB_1x4p_A, etc.
            for suffix_pat in [
                r"_MODEL_\d+_[A-Za-z0-9]+$",  # _MODEL_1_A
                r"_[A-Z]$",                      # _A _B _C (chain letter)
                r"_\d+$",                        # _1 _2 (chain number)
            ]:
                # pre-register a few common variants
                for n in range(1, 25):
                    for c in list("ABCDEFGHIJKLMNOPQRSTUVWXYZ") + [str(i) for i in range(1,10)]:
                        stem_to_gene[f"{stem}_MODEL_{n}_{c}"] = gene
                        stem_to_gene[f"{stem}_{c}"] = gene
    # Also scan pdbs dir to catch RCSB files not listed in manifest
    # e.g. 53BP1_Homo_RCSB_8s0e.pdb -> gene=53BP1
    import re as _re
    pdbs_dir = Path(manifest_path).parent / "pdbs"
    if pdbs_dir.exists():
        for pdb_file in pdbs_dir.iterdir():
            file_stem = pdb_file.stem
            if file_stem not in stem_to_gene:
                # Extract gene name by stripping species/source suffix
                m = _re.match(r"^(.+?)_(Homo|Mus_|Dros|Neur|Arab|Caen|Dani|Schi|Cryp|Magn|Zymo|Fusa)", file_stem)
                if m:
                    # Check if the gene part matches a known gene
                    gene_part = m.group(1)
                    # Find the gene in stem_map values
                    for known_stem, gene in list(stem_to_gene.items()):
                        if gene == gene_part:
                            stem_to_gene[file_stem] = gene
                            break
    log.info(f"Loaded {len(stem_to_gene)} stem->gene mappings from manifest+pdbs")
    return stem_to_gene


def pdb_to_gene(pdb_name: str, stem_map: dict[str, str]) -> str:
    """
    Convert a FoldSeek query/target field to a gene name.
    RCSB multi-chain PDBs have chain/model suffixes (_A, _2, _MODEL_1_A)
    that need to be stripped before manifest lookup.
    Falls back to the raw stem if not found in the manifest map.
    """
    stem = Path(pdb_name).stem
    if stem in stem_map:
        return stem_map[stem]
    # Strip chain/model suffixes progressively until we find a match
    # Handles: _MODEL_1_A, _A, _2, _MODEL_20_B etc.
    candidate = re.sub(r"(_MODEL_\d+)?_[A-Za-z0-9]+$", "", stem)
    while candidate and candidate != stem:
        if candidate in stem_map:
            return stem_map[candidate]
        stem = candidate
        candidate = re.sub(r"(_MODEL_\d+)?_[A-Za-z0-9]+$", "", stem)
    return pdb_name.split("/")[-1].split(".")[0]  # last resort: raw stem


def load_metadata(csv_path: str) -> dict:
    meta = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            meta[row["protein_gene_name"].strip()] = row
    return meta


def parse_hits(hits_path: str, tmscore_thresh: float,
               evalue_thresh: float, log) -> list[dict]:
    hits = []
    with open(hits_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < len(FOLDSEEK_COLS):
                continue
            hit = dict(zip(FOLDSEEK_COLS, row))
            try:
                qtm     = float(hit["qtmscore"])
                ttm     = float(hit["ttmscore"])
                ev      = float(hit["evalue"])
                mean_tm = (qtm + ttm) / 2
            except (ValueError, KeyError):
                continue
            if mean_tm < tmscore_thresh:
                continue
            hit["mean_tmscore"] = round(mean_tm, 4)
            hits.append(hit)
    log.info(f"Loaded {len(hits)} hits (TM≥{tmscore_thresh}, e≤{evalue_thresh})")
    return hits


def load_clusters(cluster_path: str) -> dict:
    mapping = {}
    with open(cluster_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                mapping[row[1]] = row[0]
    return mapping


def build_matrix(hits: list[dict], proteins: list[str],
                 stem_map: dict[str, str]) -> np.ndarray:
    idx = {p: i for i, p in enumerate(proteins)}
    n   = len(proteins)
    mat = np.zeros((n, n), dtype=np.float32)
    np.fill_diagonal(mat, 1.0)
    for hit in hits:
        q = pdb_to_gene(hit["query"],  stem_map)
        t = pdb_to_gene(hit["target"], stem_map)
        if q in idx and t in idx:
            s = hit["mean_tmscore"]
            mat[idx[q], idx[t]] = s
            mat[idx[t], idx[q]] = s
    return mat


def annotate_hits(hits: list[dict], meta: dict,
                  clusters: dict, stem_map: dict) -> list[dict]:
    out = []
    for hit in hits:
        q  = pdb_to_gene(hit["query"],  stem_map)
        t  = pdb_to_gene(hit["target"], stem_map)
        if q == t:   # skip same-gene hits (different models/chains)
            continue
        qm = meta.get(q, {})
        tm = meta.get(t, {})
        out.append({
            "query_gene":          q,
            "target_gene":         t,
            "query_organism":      qm.get("organism", "unknown"),
            "target_organism":     tm.get("organism", "unknown"),
            "query_modification":  qm.get("modification_read", "unknown"),
            "target_modification": tm.get("modification_read", "unknown"),
            "query_domain":        qm.get("reader_domain", "unknown"),
            "target_domain":       tm.get("reader_domain", "unknown"),
            "mean_tmscore":        hit["mean_tmscore"],
            "qtmscore":            hit.get("qtmscore", ""),
            "ttmscore":            hit.get("ttmscore", ""),
            "lddt":                hit.get("lddt", ""),
            "rmsd":                hit.get("rmsd", ""),
            "fident":              hit.get("fident", ""),
            "evalue":              hit.get("evalue", ""),
            "cluster_rep":         clusters.get(hit["query"], ""),
            "same_modification":   (qm.get("modification_read","") ==
                                    tm.get("modification_read","")),
            "same_domain":         (qm.get("reader_domain","") ==
                                    tm.get("reader_domain","")),
            "cross_kingdom":       (qm.get("organism","") !=
                                    tm.get("organism","")),
        })
    return out


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hits",
        default=str(SCRATCH / "results/allvall_results.tsv"))
    p.add_argument("--clusters",
        default=str(SCRATCH / "results/clusters_cluster.tsv"))
    p.add_argument("--meta",
        default=str(SCRATCH / "master_reader_annotated_v4.csv"))
    p.add_argument("--manifest",
        default=str(SCRATCH / "pdb_lists/structure_manifest.tsv"),
        help="structure_manifest.tsv — used to map filename stems to gene names")
    p.add_argument("--out",
        default=str(SCRATCH / "results/annotated_hits.csv"))
    p.add_argument("--matrix",
        default=str(SCRATCH / "results/tmscore_matrix.csv"))
    p.add_argument("--logdir",
        default=str(SCRATCH / "logs"))
    p.add_argument("--tmscore", type=float, default=0.4)
    p.add_argument("--evalue",  type=float, default=0.001)
    return p.parse_args()


def main():
    args = parse_args()
    log  = setup_logging(Path(args.logdir))

    log.info(f"Scratch dir : {SCRATCH}")
    log.info(f"Hits        : {args.hits}")
    log.info(f"Manifest    : {args.manifest}")
    log.info(f"Meta        : {args.meta}")

    # Build the filename-stem -> gene reverse map FIRST
    stem_map = build_stem_to_gene(args.manifest, log)

    meta     = load_metadata(args.meta)
    clusters = load_clusters(args.clusters)
    proteins = list(meta.keys())
    hits     = parse_hits(args.hits, args.tmscore, args.evalue, log)
    annotated = annotate_hits(hits, meta, clusters, stem_map)

    # Warn if we have unresolved gene names (fallback stems)
    unknown_q = {r["query_gene"]  for r in annotated if r["query_organism"]  == "unknown"}
    unknown_t = {r["target_gene"] for r in annotated if r["target_organism"] == "unknown"}
    unresolved = unknown_q | unknown_t
    if unresolved:
        log.warning(
            f"{len(unresolved)} gene names not found in metadata "
            f"(metadata will show 'unknown'). First 5: "
            f"{sorted(unresolved)[:5]}"
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if annotated:
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(annotated[0].keys()))
            w.writeheader()
            w.writerows(annotated)
        log.info(f"Annotated hits : {out_path}  ({len(annotated)} rows)")
    else:
        log.warning("No hits above thresholds — check --tmscore / --evalue")

    mat = build_matrix(hits, proteins, stem_map)
    with open(args.matrix, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + proteins)
        for i, prot in enumerate(proteins):
            w.writerow([prot] + [f"{v:.4f}" for v in mat[i]])
    log.info(f"TM-score matrix: {args.matrix}  ({len(proteins)}×{len(proteins)})")

    cluster_map = defaultdict(list)
    for member, rep in clusters.items():
        cluster_map[rep].append(member)
    log.info("\nTop structural clusters (≥3 members):")
    for rep, members in sorted(cluster_map.items(), key=lambda x: -len(x[1])):
        if len(members) >= 3:
            rg = pdb_to_gene(rep, stem_map)
            rm = meta.get(rg, {})
            log.info(
                f"  {rg} ({rm.get('organism','?')}, "
                f"{rm.get('reader_domain','?')}): {len(members)} members"
            )

    log.info("Done.")


if __name__ == "__main__":
    main()
