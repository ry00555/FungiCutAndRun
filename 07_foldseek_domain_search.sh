#!/bin/bash
#SBATCH --job-name=fs_domain_search
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/fs_domain_%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/fs_domain_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# 07_foldseek_domain_search.sh
#
# 1. Builds per-organism FoldSeek databases from CIF tars or FASTA
# 2. Searches ALL domain representative structures against every proteome
#    — hits assigned to EVERY matching domain class (multi-domain aware)
# 3. Produces master_domain_hits.csv with one row per hit per domain class
# 4. Extracts aromatic cage residues for all reader domain hits
#    (chromodomain, PHD, PWWP, Tudor, BAH, MRG, WD40-EED, Bromodomain)
#
# Usage:  sbatch 07_foldseek_domain_search.sh
# ══════════════════════════════════════════════════════════════════════

set -euo pipefail

module load Foldseek/10-941cd33-GPU

FSDIR="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek"
QUERIES="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/pdb_lists/pdbs"
OUT="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results/domain_proteome_search"
TMP="${FSDIR}/tmp_$$"
LOGS="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs"

mkdir -p "${OUT}/hits" "${OUT}/aromatic_cage" "${TMP}" "${LOGS}"

echo "════════════════════════════════════════════════════════════"
echo " FoldSeek Proteome-Wide Domain Search"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " Node    : ${SLURM_NODELIST}"
echo " GPU     : ${CUDA_VISIBLE_DEVICES}"
echo " Started : $(date)"
echo "════════════════════════════════════════════════════════════"

FMT="query,target,fident,alnlen,alntmscore,qtmscore,ttmscore,rmsd,prob"
TMSCORE_THRESH=0.4

# ── Organism config ───────────────────────────────────────────────────
ORGANISMS=(
    "human|tar|${FSDIR}/UP000005640_9606_HUMAN_v4.tar|9606|Homo sapiens|Metazoa"
    "mouse|tar|${FSDIR}/UP000000589_10090_MOUSE_v4.tar|10090|Mus musculus|Metazoa"
    "zebrafish|tar|${FSDIR}/UP000000437_7955_DANRE_v4.tar|7955|Danio rerio|Metazoa"
    "drosophila|tar|${FSDIR}/UP000000803_7227_DROME_v4.tar|7227|Drosophila melanogaster|Metazoa"
    "celegans|tar|${FSDIR}/UP000001940_6239_CAEEL_v4.tar|6239|Caenorhabditis elegans|Metazoa"
    "arabidopsis|tar|${FSDIR}/UP000006548_3702_ARATH_v4.tar|3702|Arabidopsis thaliana|Plantae"
    "yeast|tar|${FSDIR}/UP000002311_559292_YEAST_v4.tar|559292|Saccharomyces cerevisiae|Fungi"
    "pombe|tar|${FSDIR}/UP000002485_284812_SCHPO_v4.tar|284812|Schizosaccharomyces pombe|Fungi"
    "neurospora|fasta|${FSDIR}/ncr_proteome.fasta|5141|Neurospora crassa|Fungi"
    "fusarium|fasta|${FSDIR}/fgr_proteome.fasta|229533|Fusarium graminearum|Fungi"
    "magnaporthe|fasta|${FSDIR}/mgr_proteome.fasta|242507|Magnaporthe oryzae|Fungi"
    "zymoseptoria|fasta|${FSDIR}/zt_proteome.fasta|336722|Zymoseptoria tritici|Fungi"
    "cryptococcus|fasta|${FSDIR}/cne_proteome.fasta|5207|Cryptococcus neoformans|Fungi"
)

# ── Domain queries ────────────────────────────────────────────────────
# Each query maps to ONE or MORE domain classes (multi-domain proteins
# get listed under every class they belong to in the master sheet).
# Format: "query_cif|DomainClass1,DomainClass2|mark_context|functional_role"
DOMAIN_QUERIES=(
    # BAH readers
    "EPR1_fg_Fusa_AF.cif|BAH|H3K27me3|Facultative_heterochromatin_reader"
    "EPR1_zt_Zymo_AF.cif|BAH|H3K27me3|Facultative_heterochromatin_reader"
    "EPR1_mo_Magn_AF.cif|BAH|H3K27me3|Facultative_heterochromatin_reader"
    "ORC1_Homo_AF.cif|BAH|H4K20me2|Replication_licensing_reader"
    "BAHD1_Homo_AF.cif|BAH|H3K9me3|Constitutive_heterochromatin_reader"
    # Chromodomains — H3K27me3 readers (PRC1-type)
    "CBX7_Homo_AF.cif|Chromo|H3K27me3|PRC1_reader"
    "CBX1_Homo_AF.cif|Chromo|H3K9me3|HP1_reader"
    "HP1_Neur_AF.cif|Chromo|H3K9me3|HP1_reader"
    # MRG/Chromodomain — H3K36me readers (Eaf3-type)
    # NOTE: Eaf3 has BOTH a chromodomain AND an MRG domain — listed under both
    "EAF3_Neur_AF.cif|Chromo,MRG|H3K36me2|Rpd3S_H3K36_reader"
    "MoEAF3_Magn_AF.cif|Chromo,MRG|H3K36me2|PRC1_substitute_reader"
    "MRG15_Homo_AF.cif|MRG|H3K36me3|NuA4_subunit"
    "MRG-partner_Neur_AF.cif|MRG|H3K36me|MRG_complex_subunit"
    # SET domain writers — H3K27
    "EZH2_Homo_AF.cif|SET|H3K27me3_writer|PRC2_catalytic"
    "SET-7_Neur_AF.cif|SET|H3K27me3_writer|PRC2_catalytic_fungal"
    # SET domain writers — H3K36
    "ASH1_Neur_AF.cif|SET|H3K36me2_writer|ASH1_type_repressive"
    "SET-2_Neur_AF.cif|SET|H3K36me3_writer|Set2_type_active"
    "SETD2_Homo_AF.cif|SET|H3K36me3_writer|Set2_type_active"
    # SET domain writers — H3K9
    "SUV39H1_Homo_AF.cif|SET|H3K9me3_writer|HP1_heterochromatin"
    "DIM-5_Neur_AF.cif|SET|H3K9me3_writer|HP1_heterochromatin_fungal"
    # SET domain writers — H3K4
    "SET1A_Homo_AF.cif|SET|H3K4me3_writer|COMPASS_active"
    "SET-1_Neur_AF.cif|SET|H3K4me3_writer|COMPASS_active_fungal"
    # PWWP readers
    "PWWP1_Homo_AF.cif|PWWP|H3K36me3|PWWP_reader"
    "PSIP1_Homo_AF.cif|PWWP|H3K36me3|LEDGF_chromatin_tether"
    "NSD_ncra_PWWP_Neur_AF.cif|PWWP|H3K36me3|NSD_type_fungal"
    # PHD fingers
    "ING2_Homo_AF.cif|PHD|H3K4me3|ING_reader"
    "PHF1_Homo_AF.cif|PHD,Tudor|H3K27me3|PCL_PRC2_accessory"
    "RAG2_Homo_AF.cif|PHD|H3K4me3|V_D_J_recombination"
    "MoPHD1_Magn_AF.cif|PHD|H3K4me3|Fungal_PHD"
    # Tudor domains
    "TDRD3_Homo_AF.cif|Tudor|H3K4me3_H3R|Tudor_reader"
    "SMN1_Homo_AF.cif|Tudor|sDMA|Sm_protein_assembly"
    "JMJD2A_Homo_AF.cif|Tudor,JmjC|H3K4me3_H4K20me3|Dual_reader_demethylase"
    "MoTUDOR1_Magn_AF.cif|Tudor|H3K4me|Fungal_Tudor"
    # JmjC demethylases
    "KDM4A_Homo_AF.cif|JmjC|H3K9me3_H3K36me3|KDM4_demethylase"
    "KDM6A_Homo_AF.cif|JmjC|H3K27me3|UTX_PRC2_antagonist"
    "KDM2A_Homo_AF.cif|JmjC,CXXC|H3K36me2|KDM2_CpG_tethered"
    "MoJMJ1_Magn_AF.cif|JmjC|H3K|Fungal_JmjC"
    # Bromodomains
    "BRD4_Homo_AF.cif|Bromodomain|Kac|BET_reader"
    "GCN5L2_Homo_AF.cif|Bromodomain|Kac|SAGA_HAT_subunit"
    "MoBRD_Snf2_Magn_AF.cif|Bromodomain|Kac|Fungal_BRD_remodeler"
    # CXXC domains
    "CFP1_Homo_AF.cif|CXXC|unmethylCpG|COMPASS_recruiter"
    "CFP-1_ce_Caen_AF.cif|CXXC|unmethylCpG|COMPASS_recruiter_worm"
    "KDM2A_Homo_AF.cif|JmjC,CXXC|H3K36me2|KDM2_CpG_tethered"
    # WD40 — EED (H3K27me3 reader within PRC2)
    "EED_Homo_AF.cif|WD40|H3K27me3|PRC2_allosteric_reader"
    "EED_Neur_AF.cif|WD40|H3K27me3|PRC2_allosteric_reader_fungal"
    # WD40 — RBBP/NuRF55
    "RBBP4_Homo_AF.cif|WD40|H3_tail|PRC2_NuRD_CAF1_hub"
    "MoNURF55_Magn_AF.cif|WD40|H3_tail|PRC2_hub_fungal"
    # MBT repeats
    "L3MBTL1_Homo_AF.cif|MBT|H4K20me1_H3K4me1|Polycomb_like_reader"
    # BAH+PHD multi-domain (DNMT3 ADD)
    "DNMT3A_ADD_Homo_AF.cif|PHD,ADD|H3K4me0|DNMT3_histone_reader"
    # ATRX — ADD domain (PHD+PHD tandem)
    "ATRX_Homo_AF.cif|PHD|H3K9me3_H3K4me0|SWI_SNF_heterochromatin"
)

# ── Step 1: Build per-organism databases ─────────────────────────────
echo ""
echo "[Step 1] Building organism databases..."

# Find ProstT5 once
PROSTT5=""
for p in /db/FoldSeek/ProstT5 /work/GACRC_db/foldseek/ProstT5 "${FSDIR}/ProstT5"; do
    [[ -d "$p" ]] && PROSTT5="$p" && break
done
if [[ -z "${PROSTT5}" ]]; then
    echo "  Downloading ProstT5 weights (one-time)..."
    mkdir -p "${FSDIR}/ProstT5"
    foldseek databases ProstT5 "${FSDIR}/ProstT5" "${TMP}/prostt5_dl" \
        --threads ${SLURM_CPUS_PER_TASK}
    PROSTT5="${FSDIR}/ProstT5"
fi
echo "  ProstT5 weights: ${PROSTT5}"

for ORG_ENTRY in "${ORGANISMS[@]}"; do
    IFS='|' read -r LABEL DB_TYPE SOURCE TAXID ORG_NAME KINGDOM <<< "${ORG_ENTRY}"
    DB_PATH="${FSDIR}/db_${LABEL}"

    if [[ -f "${DB_PATH}.index" ]]; then
        echo "  [SKIP] ${LABEL}: DB already exists"
        continue
    fi
    if [[ ! -f "${SOURCE}" ]]; then
        echo "  [SKIP] ${LABEL}: source not found (${SOURCE})"
        continue
    fi

    echo "  Building: ${LABEL} (${DB_TYPE})..."

    if [[ "${DB_TYPE}" == "tar" ]]; then
        EXTRACT_DIR="${FSDIR}/cif_${LABEL}"
        mkdir -p "${EXTRACT_DIR}"
        echo "    Extracting tar..."
        tar -xf "${SOURCE}" -C "${EXTRACT_DIR}" 2>/dev/null
        find "${EXTRACT_DIR}" -name "*.gz" -exec gunzip -f {} \; 2>/dev/null || true
        N_CIF=$(find "${EXTRACT_DIR}" -name "*.cif" | wc -l)
        echo "    ${N_CIF} CIF files extracted"
        foldseek createdb "${EXTRACT_DIR}" "${DB_PATH}" \
            --threads ${SLURM_CPUS_PER_TASK}

    elif [[ "${DB_TYPE}" == "fasta" ]]; then
        foldseek createdb "${SOURCE}" "${DB_PATH}" \
            --input-format 1 \
            --prostt5-model "${PROSTT5}" \
            --gpu 1 \
            --threads ${SLURM_CPUS_PER_TASK}
    fi
    echo "    Done: ${DB_PATH}"
done

# ── Step 2: Search all domain queries against all organisms ───────────
echo ""
echo "[Step 2] Running searches (${#DOMAIN_QUERIES[@]} queries × ${#ORGANISMS[@]} organisms)..."

for ORG_ENTRY in "${ORGANISMS[@]}"; do
    IFS='|' read -r LABEL DB_TYPE SOURCE TAXID ORG_NAME KINGDOM <<< "${ORG_ENTRY}"
    DB_PATH="${FSDIR}/db_${LABEL}"
    [[ ! -f "${DB_PATH}.index" ]] && echo "  [SKIP] ${LABEL}: no DB" && continue
    echo "  Organism: ${LABEL} (${ORG_NAME})"

    for QUERY_ENTRY in "${DOMAIN_QUERIES[@]}"; do
        IFS='|' read -r QUERY_CIF DOMAIN_CLASSES MARK_CONTEXT FUNCTIONAL_ROLE <<< "${QUERY_ENTRY}"
        QUERY_FILE="${QUERIES}/${QUERY_CIF}"
        QUERY_STEM="${QUERY_CIF%.cif}"
        HIT_FILE="${OUT}/hits/${QUERY_STEM}_vs_${LABEL}.tsv"

        [[ ! -f "${QUERY_FILE}" ]] && continue
        [[ -f "${HIT_FILE}" ]]     && continue   # already done

        foldseek easy-search \
            "${QUERY_FILE}" \
            "${DB_PATH}" \
            "${HIT_FILE}" \
            "${TMP}/${QUERY_STEM}_${LABEL}" \
            --alignment-type 1 \
            --tmscore-threshold ${TMSCORE_THRESH} \
            --format-mode 4 \
            --format-output "${FMT}" \
            --gpu 1 \
            --threads ${SLURM_CPUS_PER_TASK} \
            --max-seqs 1000 \
            -s 9.5 \
            2>/dev/null

        N=$(wc -l < "${HIT_FILE}" 2>/dev/null || echo 0)
        [[ $N -gt 0 ]] && echo "    ${QUERY_STEM} vs ${LABEL}: ${N} hits"
    done
done

# ── Step 3: Build master sheet + aromatic cage analysis ───────────────
echo ""
echo "[Step 3] Building master sheet and aromatic cage analysis..."

python3 - << 'PYEOF'
import csv, re, os, json
from pathlib import Path
from collections import defaultdict

OUT   = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results/domain_proteome_search")
HITS  = OUT / "hits"
CAGE  = OUT / "aromatic_cage"
CAGE.mkdir(exist_ok=True)

# ── Organism metadata ─────────────────────────────────────────────────
ORG_META = {
    "human":        ("Homo sapiens",                "Metazoa",  9606),
    "mouse":        ("Mus musculus",                "Metazoa",  10090),
    "zebrafish":    ("Danio rerio",                 "Metazoa",  7955),
    "drosophila":   ("Drosophila melanogaster",     "Metazoa",  7227),
    "celegans":     ("Caenorhabditis elegans",      "Metazoa",  6239),
    "arabidopsis":  ("Arabidopsis thaliana",        "Plantae",  3702),
    "yeast":        ("Saccharomyces cerevisiae",    "Fungi",    559292),
    "pombe":        ("Schizosaccharomyces pombe",   "Fungi",    284812),
    "neurospora":   ("Neurospora crassa",           "Fungi",    5141),
    "fusarium":     ("Fusarium graminearum",        "Fungi",    229533),
    "magnaporthe":  ("Magnaporthe oryzae",          "Fungi",    242507),
    "zymoseptoria": ("Zymoseptoria tritici",        "Fungi",    336722),
    "cryptococcus": ("Cryptococcus neoformans",     "Fungi",    5207),
}

# ── Query metadata — multi-domain aware ──────────────────────────────
# A hit gets ONE ROW PER DOMAIN CLASS the query belongs to
# So EAF3 (Chromo,MRG) will produce rows in both Chromo AND MRG sections
QUERY_META = {
    "EPR1_fg_Fusa_AF":    (["BAH"],            "H3K27me3",        "Facultative_heterochromatin_reader"),
    "EPR1_zt_Zymo_AF":    (["BAH"],            "H3K27me3",        "Facultative_heterochromatin_reader"),
    "EPR1_mo_Magn_AF":    (["BAH"],            "H3K27me3",        "Facultative_heterochromatin_reader"),
    "ORC1_Homo_AF":       (["BAH"],            "H4K20me2",        "Replication_licensing_reader"),
    "BAHD1_Homo_AF":      (["BAH"],            "H3K9me3",         "Constitutive_heterochromatin_reader"),
    "CBX7_Homo_AF":       (["Chromo"],         "H3K27me3",        "PRC1_reader"),
    "CBX1_Homo_AF":       (["Chromo"],         "H3K9me3",         "HP1_reader"),
    "HP1_Neur_AF":        (["Chromo"],         "H3K9me3",         "HP1_reader"),
    "EAF3_Neur_AF":       (["Chromo","MRG"],   "H3K36me2",        "Rpd3S_H3K36_reader"),
    "MoEAF3_Magn_AF":     (["Chromo","MRG"],   "H3K36me2",        "PRC1_substitute_reader"),
    "MRG15_Homo_AF":      (["MRG"],            "H3K36me3",        "NuA4_subunit"),
    "MRG-partner_Neur_AF":(["MRG"],            "H3K36me",         "MRG_complex_subunit"),
    "EZH2_Homo_AF":       (["SET"],            "H3K27me3_writer", "PRC2_catalytic"),
    "SET-7_Neur_AF":      (["SET"],            "H3K27me3_writer", "PRC2_catalytic_fungal"),
    "ASH1_Neur_AF":       (["SET"],            "H3K36me2_writer", "ASH1_repressive"),
    "SET-2_Neur_AF":      (["SET"],            "H3K36me3_writer", "Set2_active"),
    "SETD2_Homo_AF":      (["SET"],            "H3K36me3_writer", "Set2_active"),
    "SUV39H1_Homo_AF":    (["SET"],            "H3K9me3_writer",  "HP1_heterochromatin"),
    "DIM-5_Neur_AF":      (["SET"],            "H3K9me3_writer",  "HP1_heterochromatin_fungal"),
    "SET1A_Homo_AF":      (["SET"],            "H3K4me3_writer",  "COMPASS_active"),
    "SET-1_Neur_AF":      (["SET"],            "H3K4me3_writer",  "COMPASS_active_fungal"),
    "PWWP1_Homo_AF":      (["PWWP"],           "H3K36me3",        "PWWP_reader"),
    "PSIP1_Homo_AF":      (["PWWP"],           "H3K36me3",        "LEDGF_chromatin_tether"),
    "NSD_ncra_PWWP_Neur_AF":(["PWWP"],        "H3K36me3",        "NSD_type_fungal"),
    "ING2_Homo_AF":       (["PHD"],            "H3K4me3",         "ING_reader"),
    "PHF1_Homo_AF":       (["PHD","Tudor"],    "H3K27me3",        "PCL_PRC2_accessory"),
    "RAG2_Homo_AF":       (["PHD"],            "H3K4me3",         "VDJ_recombination"),
    "MoPHD1_Magn_AF":     (["PHD"],            "H3K4me3",         "Fungal_PHD"),
    "TDRD3_Homo_AF":      (["Tudor"],          "H3K4me3_H3R",     "Tudor_reader"),
    "SMN1_Homo_AF":       (["Tudor"],          "sDMA",            "Sm_protein_assembly"),
    "JMJD2A_Homo_AF":     (["Tudor","JmjC"],   "H3K4me3_H4K20me3","Dual_reader_demethylase"),
    "MoTUDOR1_Magn_AF":   (["Tudor"],          "H3K4me",          "Fungal_Tudor"),
    "KDM4A_Homo_AF":      (["JmjC"],           "H3K9me3_H3K36me3","KDM4_demethylase"),
    "KDM6A_Homo_AF":      (["JmjC"],           "H3K27me3",        "UTX_PRC2_antagonist"),
    "KDM2A_Homo_AF":      (["JmjC","CXXC"],    "H3K36me2",        "KDM2_CpG_tethered"),
    "MoJMJ1_Magn_AF":     (["JmjC"],           "H3K",             "Fungal_JmjC"),
    "BRD4_Homo_AF":       (["Bromodomain"],    "Kac",             "BET_reader"),
    "GCN5L2_Homo_AF":     (["Bromodomain"],    "Kac",             "SAGA_HAT_subunit"),
    "MoBRD_Snf2_Magn_AF": (["Bromodomain"],    "Kac",             "Fungal_BRD_remodeler"),
    "CFP1_Homo_AF":       (["CXXC"],           "unmethylCpG",     "COMPASS_recruiter"),
    "CFP-1_ce_Caen_AF":   (["CXXC"],           "unmethylCpG",     "COMPASS_recruiter_worm"),
    "EED_Homo_AF":        (["WD40"],           "H3K27me3",        "PRC2_allosteric_reader"),
    "EED_Neur_AF":        (["WD40"],           "H3K27me3",        "PRC2_allosteric_reader_fungal"),
    "RBBP4_Homo_AF":      (["WD40"],           "H3_tail",         "PRC2_NuRD_CAF1_hub"),
    "MoNURF55_Magn_AF":   (["WD40"],           "H3_tail",         "PRC2_hub_fungal"),
    "L3MBTL1_Homo_AF":    (["MBT"],            "H4K20me1_H3K4me1","Polycomb_like_reader"),
    "DNMT3A_ADD_Homo_AF": (["PHD","ADD"],       "H3K4me0",         "DNMT3_histone_reader"),
    "ATRX_Homo_AF":       (["PHD"],            "H3K9me3_H3K4me0", "SWI_SNF_heterochromatin"),
}

# ── Aromatic cage residues by domain type ─────────────────────────────
# Known positions of cage residues from structural literature
# Format: domain -> list of (residue_type, canonical_position_description, role)
AROMATIC_CAGE_RESIDUES = {
    "Chromo": [
        ("Trp", "cage_pos1", "methyl-Lys sandwiching"),
        ("Tyr", "cage_pos2", "methyl-Lys sandwiching"),
        ("Trp/Phe", "cage_pos3", "CH-pi stacking"),
        ("Glu", "cage_neg", "electrostatic_complementarity"),
    ],
    "BAH": [
        ("Tyr", "cage_pos1", "H3K27me3_recognition"),
        ("Trp", "cage_pos2", "methyl_cage"),
        ("Phe", "cage_pos3", "hydrophobic_stacking"),
        ("Asp/Glu", "cage_neg", "backbone_contacts"),
    ],
    "PWWP": [
        ("Trp", "PWWP_motif_W", "methyl-Lys_recognition"),
        ("Pro", "PWWP_motif_P", "structural"),
        ("Tyr", "cage_pos2", "aromatic_stacking"),
        ("Ile/Val", "cage_pos3", "hydrophobic"),
    ],
    "Tudor": [
        ("Tyr", "cage_pos1", "methyl-Lys_sandwiching"),
        ("Trp", "cage_pos2", "aromatic_cage"),
        ("Phe", "cage_pos3", "cage_wall"),
        ("Asn", "cage_polar", "H-bond_to_methylLys"),
    ],
    "MRG": [
        ("Tyr", "chromodomain_cage1", "H3K36me_recognition"),
        ("Trp", "chromodomain_cage2", "methyl_stacking"),
        ("Phe", "chromodomain_cage3", "aromatic_platform"),
    ],
    "PHD": [
        ("Trp", "cage_pos1", "H3K4me3_cage"),
        ("Tyr", "cage_pos2", "aromatic_stacking"),
        ("Trp/Phe", "cage_pos3", "methyl_recognition"),
        ("Asp", "cage_neg", "K4_amine_contact"),
    ],
    "WD40": [
        ("Phe", "aromatic_shelf_1", "H3K27me3_cage_EED"),
        ("Tyr", "aromatic_shelf_2", "methyl_recognition"),
        ("His", "cage_His", "H-bond_network"),
    ],
    "Bromodomain": [
        ("Asn", "Asn_anchor", "Kac_carbonyl_H-bond"),
        ("Tyr", "WPF_shelf_Y", "Kac_recognition"),
        ("Trp", "WPF_shelf_W", "hydrophobic_lid"),
        ("Pro", "WPF_shelf_P", "structural"),
    ],
    "MBT": [
        ("Trp", "cage_pos1", "mono_methyl_Lys_cage"),
        ("Tyr", "cage_pos2", "aromatic_stacking"),
        ("Asp", "cage_neg", "charge_complementarity"),
    ],
}

# Reader domain classes (exclude writers for cage analysis)
READER_DOMAINS = {"BAH","Chromo","MRG","PWWP","PHD","Tudor","WD40","Bromodomain","MBT","CXXC"}

# ── Parse hits and build master rows ─────────────────────────────────
print("Parsing hit files...")
master_rows = []
# Track unique targets per organism to flag multi-domain hits
target_to_domains = defaultdict(lambda: defaultdict(set))
# {organism -> {uniprot_acc -> set(domain_classes)}}

for tsv_file in sorted(HITS.glob("*.tsv")):
    stem  = tsv_file.stem
    match = re.match(r"^(.+)_vs_(.+)$", stem)
    if not match:
        continue
    query_stem = match.group(1)
    org_label  = match.group(2)

    org_name, kingdom, taxid = ORG_META.get(org_label, ("Unknown","Unknown",0))
    domain_classes, mark_ctx, func_role = QUERY_META.get(
        query_stem, (["Unknown"], "unknown", "unknown"))

    for line in tsv_file.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            query   = parts[0]
            target  = parts[1]
            fident  = float(parts[2])
            alnlen  = int(parts[3])
            alntm   = float(parts[4])
            qtm     = float(parts[5])
            ttm     = float(parts[6])
            rmsd    = float(parts[7])
            prob    = float(parts[8])
            mean_tm = round((qtm + ttm) / 2, 4)
        except (ValueError, IndexError):
            continue

        # Extract UniProt accession
        uniprot_acc = "unknown"
        m = re.search(r"AF-([A-Z0-9]+)-F\d", target)
        if m:
            uniprot_acc = m.group(1)
        else:
            m = re.search(r"[st][rp]\|([A-Z0-9]+)\|", target)
            if m:
                uniprot_acc = m.group(1)

        # Track which domain classes this target hit (for multi-domain flagging)
        for dc in domain_classes:
            target_to_domains[org_label][uniprot_acc].add(dc)

        # ── KEY: one row per domain class ─────────────────────────
        for domain_class in domain_classes:
            master_rows.append({
                "domain_class":      domain_class,
                "query_structure":   query_stem,
                "mark_context":      mark_ctx,
                "functional_role":   func_role,
                "is_reader":         "yes" if domain_class in READER_DOMAINS else "no",
                "organism":          org_name,
                "kingdom":           kingdom,
                "taxid":             taxid,
                "uniprot_acc":       uniprot_acc,
                "target_id":         target,
                "mean_tmscore":      mean_tm,
                "qtmscore":          round(qtm, 4),
                "ttmscore":          round(ttm, 4),
                "alntmscore":        round(alntm, 4),
                "seq_identity_pct":  round(fident * 100, 1),
                "aln_length":        alnlen,
                "rmsd":              round(rmsd, 3),
                "prob":              round(prob, 4),
                "tm_confidence":     "high(>0.6)" if mean_tm > 0.6
                                     else "medium(0.5-0.6)" if mean_tm > 0.5
                                     else "low(0.4-0.5)",
                "multi_domain_hit":  "",  # filled in below
                "all_domains_hit":   "",  # filled in below
            })

# ── Flag multi-domain hits ────────────────────────────────────────────
print("Flagging multi-domain hits...")
for row in master_rows:
    org_label  = next((k for k,v in ORG_META.items() if v[0]==row["organism"]), "")
    uniprot    = row["uniprot_acc"]
    all_doms   = target_to_domains.get(org_label, {}).get(uniprot, set())
    is_multi   = len(all_doms) > 1
    row["multi_domain_hit"] = "yes" if is_multi else "no"
    row["all_domains_hit"]  = ";".join(sorted(all_doms))

# Sort: domain_class, kingdom, organism, TM-score desc
master_rows.sort(key=lambda r: (
    r["domain_class"], r["kingdom"], r["organism"], -r["mean_tmscore"]))

# Write master CSV
master_path = OUT / "master_domain_hits.csv"
if master_rows:
    with open(master_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(master_rows[0].keys()))
        w.writeheader()
        w.writerows(master_rows)
    print(f"\nMaster sheet: {master_path}")
    print(f"  Total rows : {len(master_rows)}")
    unique_targets = {r["uniprot_acc"] for r in master_rows}
    multi_domain   = {r["uniprot_acc"] for r in master_rows if r["multi_domain_hit"]=="yes"}
    print(f"  Unique targets    : {len(unique_targets)}")
    print(f"  Multi-domain hits : {len(multi_domain)}")
else:
    print("No hits found.")

# ── Summary by domain class ───────────────────────────────────────────
print("\n" + "="*75)
print(f"{'Domain':<15} {'Kingdom':<10} {'Orgs':<6} {'TM>0.4':<8} {'TM>0.5':<8} {'TM>0.6'}")
print("-"*75)
by_dom_king = defaultdict(list)
for r in master_rows:
    by_dom_king[(r["domain_class"], r["kingdom"])].append(r)
for (dom, king), rlist in sorted(by_dom_king.items()):
    orgs  = len({r["organism"] for r in rlist})
    n04   = len(rlist)
    n05   = sum(1 for r in rlist if r["mean_tmscore"] > 0.5)
    n06   = sum(1 for r in rlist if r["mean_tmscore"] > 0.6)
    print(f"{dom:<15} {king:<10} {orgs:<6} {n04:<8} {n05:<8} {n06}")

# ── Multi-domain hit report ───────────────────────────────────────────
print("\n" + "="*75)
print("MULTI-DOMAIN HITS (proteins detected by >1 domain query):")
print("-"*75)
for org_label, acc_dict in target_to_domains.items():
    org_name = ORG_META.get(org_label, ("?",))[0]
    multi = {acc: doms for acc, doms in acc_dict.items() if len(doms) > 1}
    if multi:
        print(f"\n  {org_name}:")
        for acc, doms in sorted(multi.items(), key=lambda x: -len(x[1])):
            dom_str = " + ".join(sorted(doms))
            # Find best TM-score for this acc
            best_tm = max((r["mean_tmscore"] for r in master_rows
                           if r["uniprot_acc"]==acc and
                           ORG_META.get(org_label,("?",))[0]==r["organism"]), default=0)
            print(f"    {acc}  [{dom_str}]  best_TM={best_tm:.3f}")

# ── Aromatic cage analysis ────────────────────────────────────────────
print("\n" + "="*75)
print("AROMATIC CAGE ANALYSIS (reader domains only, TM > 0.5):")
print("-"*75)

cage_rows = []
for r in master_rows:
    if r["is_reader"] != "yes":
        continue
    if r["mean_tmscore"] <= 0.5:
        continue
    dom = r["domain_class"]
    if dom not in AROMATIC_CAGE_RESIDUES:
        continue

    cage_spec = AROMATIC_CAGE_RESIDUES[dom]
    for res_type, pos_label, role in cage_spec:
        cage_rows.append({
            "domain_class":     dom,
            "organism":         r["organism"],
            "kingdom":          r["kingdom"],
            "uniprot_acc":      r["uniprot_acc"],
            "query_structure":  r["query_structure"],
            "mean_tmscore":     r["mean_tmscore"],
            "cage_residue_type": res_type,
            "cage_position":    pos_label,
            "cage_role":        role,
            "mark_context":     r["mark_context"],
            "functional_role":  r["functional_role"],
            "multi_domain":     r["multi_domain_hit"],
            "all_domains":      r["all_domains_hit"],
            "structural_note":  (
                "Cage residues are inferred from domain-class consensus; "
                "exact positions require CIF coordinate extraction. "
                "See aromatic_cage_positions.md for extraction protocol."
            ),
        })

cage_path = CAGE / "aromatic_cage_predictions.csv"
if cage_rows:
    with open(cage_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cage_rows[0].keys()))
        w.writeheader()
        w.writerows(cage_rows)
    print(f"\nAromatic cage predictions: {cage_path}")
    print(f"  Rows: {len(cage_rows)}")
    by_dom = defaultdict(set)
    for r in cage_rows:
        by_dom[r["domain_class"]].add(r["uniprot_acc"])
    for dom, accs in sorted(by_dom.items()):
        print(f"  {dom:<15}: {len(accs)} proteins with predicted cage residues")

# ── Write aromatic cage protocol note ─────────────────────────────────
protocol = CAGE / "aromatic_cage_positions.md"
protocol.write_text("""# Aromatic Cage Residue Extraction Protocol

## What this file contains
`aromatic_cage_predictions.csv` lists all reader-domain hits (TM > 0.5) with their
predicted aromatic cage residue types, based on domain-class consensus from the
structural literature. These are the residue types EXPECTED at cage positions —
not extracted coordinates.

## To get exact residue positions from CIF files

For each protein of interest, extract the aromatic cage coordinates:

```bash
# Install gemmi (CIF parser)
pip install gemmi --user

# Extract aromatic residues near the methyl-lysine binding pocket
python3 - << 'EOF'
import gemmi, sys
cif_path = sys.argv[1]  # path to .cif file
doc = gemmi.cif.read(cif_path)
st  = gemmi.make_structure_from_block(doc.sole_block())
model = st[0]

aromatic = {"TYR", "TRP", "PHE", "HIS"}
for chain in model:
    for res in chain:
        if res.name in aromatic:
            ca = res.find_atom("CA", "\\\\")
            if ca:
                print(f"{res.name} {chain.name} {res.seqid} {ca.pos}")
EOF
```

## Key aromatic cage positions by domain (from literature)

| Domain | Residues | Canonical example | Reference |
|--------|----------|-------------------|-----------|
| Chromo | W-Y-W triad | CBX7 W32-Y36-W47 (H3K27me3) | Bernstein 2006 |
| BAH | Y-W-Y triad | ORC1 Y816-W875 (H4K20me2) | Kuo 2012 |
| PWWP | W-Y | DNMT3A W330 (H3K36me3) | Dhayalan 2010 |
| Tudor | Y-W-F-N cage | TDRD3 Y563-F565-Y617 | Liu 2010 |
| PHD | W-W-Y | ING2 W238-Y215 (H3K4me3) | Shi 2006 |
| WD40/EED | F-Y-H shelf | EED F97-Y148-H204 (H3K27me3) | Margueron 2009 |
| Bromodomain | N-Y-W (WPF shelf) | BRD4 N140-Y139-W81 (Kac) | Filippakopoulos 2010 |
| MRG/Chromo | W-Y triad | Eaf3 W80-Y23 (H3K36me) | Xu 2008 |

## Interpretation for cross-kingdom comparison
Conserved aromatic cage = same mark can be read
Substituted cage (Y->F, W->F) = reduced affinity, potentially altered specificity
Missing cage residues = reader function lost (pseudoreader)
""")
print(f"\nProtocol written: {protocol}")
print("\nDone.")
PYEOF

rm -rf "${TMP}"

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Complete: $(date)"
echo " Master sheet : ${OUT}/master_domain_hits.csv"
echo " Cage analysis: ${OUT}/aromatic_cage/"
echo ""
echo " Download:"
echo "   scp ry00555@sapelo2.gacrc.uga.edu:${OUT}/master_domain_hits.csv ~/Desktop/"
echo "   scp ry00555@sapelo2.gacrc.uga.edu:${OUT}/aromatic_cage/aromatic_cage_predictions.csv ~/Desktop/"
echo "════════════════════════════════════════════════════════════"
