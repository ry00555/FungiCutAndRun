#!/usr/bin/env python3
"""
07a_uniprot_domain_pull.py
---------------------------
Query UniProt for ALL proteins containing chromatin reader/writer domains
across every organism of interest. This gives us the ground-truth list of
domain-containing proteins BEFORE any structure comparison.

For each organism × domain combination, queries:
  https://rest.uniprot.org/uniprotkb/search?query=...

Outputs:
  domain_proteins_master.csv     — all domain-containing proteins, all orgs
  domain_proteins_by_org/        — per-organism CSVs
  accessions_for_structure.txt   — flat list of UniProt accs to fetch CIFs for

Usage:
  python3 07a_uniprot_domain_pull.py
  (run on login node — API calls only, ~30 min)
"""

import csv, json, re, time, sys
from pathlib import Path
from collections import defaultdict
import requests

SCRATCH = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
OUT     = SCRATCH / "results/proteome_domain_pull"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "domain_proteins_by_org").mkdir(exist_ok=True)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# ── Organisms: taxid -> label ─────────────────────────────────────────
ORGANISMS = {
    9606:   ("Homo sapiens",               "Metazoa"),
    10090:  ("Mus musculus",               "Metazoa"),
    7955:   ("Danio rerio",                "Metazoa"),
    7227:   ("Drosophila melanogaster",    "Metazoa"),
    6239:   ("Caenorhabditis elegans",     "Metazoa"),
    3702:   ("Arabidopsis thaliana",       "Plantae"),
    559292: ("Saccharomyces cerevisiae",   "Fungi"),
    284812: ("Schizosaccharomyces pombe",  "Fungi"),
    5141:   ("Neurospora crassa",          "Fungi"),
    229533: ("Fusarium graminearum",       "Fungi"),
    242507: ("Magnaporthe oryzae",         "Fungi"),
    336722: ("Zymoseptoria tritici",       "Fungi"),
    5207:   ("Cryptococcus neoformans",    "Fungi"),
}

# ── Domain definitions ────────────────────────────────────────────────
# Each entry: domain_name -> (pfam_acc, interpro_acc, description)
# We search by Pfam ID which is most reliable
DOMAINS = {
    # Readers
    "BAH":          ("PF01428", "IPR001025", "Bromo-Adjacent Homology domain"),
    "Chromo":       ("PF00385", "IPR023780", "Chromodomain"),
    "Chromo_shadow":("PF01393", "IPR008251", "Chromodomain shadow"),
    "MRG":          ("PF05521", "IPR008574", "MRG domain"),
    "PWWP":         ("PF00855", "IPR000313", "PWWP domain"),
    "PHD":          ("PF00628", "IPR019787", "PHD finger"),
    "Tudor":        ("PF00567", "IPR001206", "Tudor domain"),
    "Bromodomain":  ("PF00439", "IPR001487", "Bromodomain"),
    "WD40":         ("PF00400", "IPR001680", "WD40 repeat"),
    "MBT":          ("PF02820", "IPR004092", "MBT repeat"),
    "CXXC":         ("PF02008", "IPR002857", "CXXC zinc finger"),
    "ADD":          ("PF06827", "IPR010536", "ADD domain"),
    # Writers — SET domain methyltransferases (all share SET Pfam)
    "SET":          ("PF00856", "IPR001214", "SET domain"),
    # Demethylases
    "JmjC":         ("PF02373", "IPR003347", "JmjC domain"),
    "JmjN":         ("PF02375", "IPR003349", "JmjN domain"),
    # LSD1-type demethylase
    "amine_oxidase":("PF01593", "IPR002937", "Amine oxidase (LSD1-type)"),
}

# ── Fields to retrieve from UniProt ──────────────────────────────────
FIELDS = [
    "accession",
    "gene_names",
    "protein_name",
    "organism_name",
    "organism_id",
    "length",
    "reviewed",       # Swiss-Prot vs TrEMBL
    "xref_pfam",      # Pfam domain list
    "xref_interpro",  # InterPro entries
    "go_p",           # GO biological process
    "cc_function",    # function annotation
    "ft_domain",      # domain feature annotations (positions)
]


def query_uniprot(taxid: int, pfam_acc: str, session: requests.Session,
                  page_size: int = 500) -> list[dict]:
    """
    Retrieve all UniProt entries for a given taxid + Pfam domain.
    Handles pagination automatically.
    """
    params = {
        "query":  f"taxonomy_id:{taxid} AND database:pfam AND xref_pfam:{pfam_acc}",
        "format": "json",
        "fields": ",".join(FIELDS),
        "size":   page_size,
    }

    all_results = []
    url = UNIPROT_SEARCH

    while url:
        try:
            r = session.get(url, params=params if url == UNIPROT_SEARCH else None,
                           timeout=30)
            r.raise_for_status()
            data = r.json()
            results = data.get("results", [])
            all_results.extend(results)

            # Check for next page via Link header
            link = r.headers.get("Link", "")
            next_url = None
            if 'rel="next"' in link:
                m = re.search(r'<([^>]+)>;\s*rel="next"', link)
                if m:
                    next_url = m.group(1)
            url    = next_url
            params = None   # only used for first request
            time.sleep(0.2)

        except Exception as e:
            print(f"    WARNING: API error for taxid={taxid} pfam={pfam_acc}: {e}")
            break

    return all_results


def parse_entry(entry: dict, domain_name: str, pfam_acc: str,
                org_name: str, kingdom: str, taxid: int) -> dict:
    """Parse a UniProt JSON entry into a flat row."""

    acc = entry.get("primaryAccession", "")

    # Gene name
    genes = entry.get("genes", [])
    gene_name = ""
    if genes:
        gn = genes[0].get("geneName", {})
        gene_name = gn.get("value", "")
        if not gene_name and genes[0].get("synonyms"):
            gene_name = genes[0]["synonyms"][0].get("value", "")

    # Protein name
    prot_desc = entry.get("proteinDescription", {})
    rec_name  = prot_desc.get("recommendedName", {})
    prot_name = rec_name.get("fullName", {}).get("value", "")
    if not prot_name:
        sub_names = prot_desc.get("submittedNames", [])
        if sub_names:
            prot_name = sub_names[0].get("fullName", {}).get("value", "")

    # Sequence length
    seq_len = entry.get("sequence", {}).get("length", 0)

    # Reviewed status
    reviewed = "yes" if entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)" else "no"

    # ALL Pfam domains in this protein (not just the query one)
    all_pfam = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "Pfam":
            pid  = ref.get("id", "")
            pname = next((p["value"] for p in ref.get("properties", [])
                         if p["key"] == "EntryName"), pid)
            all_pfam.append(f"{pid}({pname})")

    # All InterPro
    all_ipr = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "InterPro":
            iid   = ref.get("id", "")
            iname = next((p["value"] for p in ref.get("properties", [])
                         if p["key"] == "EntryName"), "")
            all_ipr.append(f"{iid}({iname})")

    # Domain feature positions (from sequence annotation)
    domain_positions = []
    for feature in entry.get("features", []):
        if feature.get("type") in ("Domain", "Zinc finger", "Repeat"):
            desc  = feature.get("description", "")
            start = feature.get("location", {}).get("start", {}).get("value", "")
            end   = feature.get("location", {}).get("end", {}).get("value", "")
            if start and end:
                domain_positions.append(f"{desc}:{start}-{end}")

    # GO biological process
    go_terms = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            props = {p["key"]: p["value"] for p in ref.get("properties", [])}
            term  = props.get("GoTerm", "")
            if term.startswith("P:"):
                go_terms.append(term[2:])

    # Function comment
    func = ""
    for cc in entry.get("comments", []):
        if cc.get("commentType") == "FUNCTION":
            texts = cc.get("texts", [])
            if texts:
                func = texts[0].get("value", "")[:200]
                break

    # Detect chromatin-relevant co-occurring domains
    chromatin_pfam_ids = {d.split("(")[0] for d in all_pfam}
    known_chromatin_pfam = {v[0] for v in DOMAINS.values()}
    cooccurring = [
        d.split("(")[1].rstrip(")") for d in all_pfam
        if d.split("(")[0] in known_chromatin_pfam
        and d.split("(")[0] != pfam_acc
    ]

    return {
        "uniprot_acc":          acc,
        "gene_name":            gene_name,
        "protein_name":         prot_name,
        "organism":             org_name,
        "kingdom":              kingdom,
        "taxid":                taxid,
        "query_domain":         domain_name,
        "query_pfam":           pfam_acc,
        "seq_length":           seq_len,
        "reviewed":             reviewed,
        "alphafold_id":         f"AF-{acc}-F1",
        "all_pfam_domains":     ";".join(all_pfam),
        "all_interpro":         ";".join(all_ipr[:10]),
        "domain_positions":     ";".join(domain_positions),
        "cooccurring_chromatin_domains": ";".join(cooccurring),
        "n_chromatin_domains":  len(cooccurring) + 1,
        "is_multidomain":       "yes" if cooccurring else "no",
        "go_biological_process":";".join(go_terms[:5]),
        "function_annotation":  func,
    }


def main():
    session = requests.Session()
    session.headers.update({
        "User-Agent": "ChromatinDomainPuller/1.0 (ry00555@uga.edu)",
        "Accept": "application/json",
    })

    all_rows  = []
    all_accs  = set()
    stats     = defaultdict(lambda: defaultdict(int))

    total_queries = len(ORGANISMS) * len(DOMAINS)
    query_n = 0

    for taxid, (org_name, kingdom) in ORGANISMS.items():
        org_rows = []
        print(f"\n{'='*65}")
        print(f"Organism: {org_name} (taxid {taxid}, {kingdom})")
        print(f"{'='*65}")

        for domain_name, (pfam_acc, ipr_acc, description) in DOMAINS.items():
            query_n += 1
            print(f"  [{query_n}/{total_queries}] {domain_name} ({pfam_acc})...",
                  end=" ", flush=True)

            results = query_uniprot(taxid, pfam_acc, session)

            if not results:
                print(f"0 hits")
                continue

            rows = [parse_entry(e, domain_name, pfam_acc, org_name, kingdom, taxid)
                    for e in results]

            print(f"{len(rows)} proteins "
                  f"({sum(1 for r in rows if r['is_multidomain']=='yes')} multidomain)")

            org_rows.extend(rows)
            all_rows.extend(rows)
            all_accs.update(r["uniprot_acc"] for r in rows)
            stats[org_name][domain_name] = len(rows)

            time.sleep(0.3)

        # Write per-organism file
        org_safe = org_name.replace(" ", "_")
        org_path = OUT / "domain_proteins_by_org" / f"{org_safe}_domain_proteins.csv"
        if org_rows:
            with open(org_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(org_rows[0].keys()))
                w.writeheader()
                # deduplicate within organism by (acc, domain)
                seen = set()
                for r in org_rows:
                    key = (r["uniprot_acc"], r["query_domain"])
                    if key not in seen:
                        seen.add(key)
                        w.writerow(r)
            print(f"  -> {org_path}")

    # Write global master sheet (deduplicated by acc+domain)
    master_path = OUT / "domain_proteins_master.csv"
    seen_global = set()
    deduped = []
    for r in all_rows:
        key = (r["uniprot_acc"], r["query_domain"])
        if key not in seen_global:
            seen_global.add(key)
            deduped.append(r)

    deduped.sort(key=lambda r: (r["domain_class"] if "domain_class" in r
                                else r["query_domain"],
                                r["kingdom"], r["organism"],
                                r["gene_name"]))

    if deduped:
        with open(master_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(deduped[0].keys()))
            w.writeheader()
            w.writerows(deduped)

    # Write flat accession list for CIF fetching
    acc_path = OUT / "accessions_for_structure.txt"
    acc_path.write_text("\n".join(sorted(all_accs)) + "\n")

    # Print summary
    print("\n" + "="*65)
    print("SUMMARY")
    print("="*65)
    print(f"Total unique proteins   : {len(all_accs)}")
    print(f"Total rows (acc×domain) : {len(deduped)}")
    print(f"Master sheet            : {master_path}")
    print(f"Accession list          : {acc_path}")

    print(f"\n{'Organism':<35} {'Domain':<15} {'Count'}")
    print("-"*60)
    for org in sorted(stats.keys()):
        for dom, n in sorted(stats[org].items()):
            if n > 0:
                print(f"  {org:<33} {dom:<15} {n}")

    multidomain = [r for r in deduped if r["is_multidomain"] == "yes"]
    print(f"\nMulti-domain proteins   : {len(multidomain)}")
    print("\nTop multi-domain hits:")
    by_acc = defaultdict(list)
    for r in multidomain:
        by_acc[(r["uniprot_acc"], r["organism"])].append(r["query_domain"])
    top = sorted(by_acc.items(), key=lambda x: -len(x[1]))[:20]
    for (acc, org), doms in top:
        print(f"  {acc:<12} {org:<35} {'+'.join(sorted(set(doms)))}")

    print("\nDone. Next: run 07b_fetch_domain_cifs.sh to get structures,")
    print("then 07c_foldseek_allvall.sh for the comparison.")


if __name__ == "__main__":
    main()
