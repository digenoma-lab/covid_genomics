#!/usr/bin/env python3

import os
import re
import json
import glob
import argparse
from pathlib import Path
import pandas as pd

# =========================
# SARS-CoV-2 reference model
# MN908947.3 / Wuhan-Hu-1
# Coordinates are 1-based inclusive.
# =========================

GENES = {
    "1ab": {"start": 266, "end": 21555},
    "s":   {"start": 21563, "end": 25384},
    "3a":  {"start": 25393, "end": 26220},
    "e":   {"start": 26245, "end": 26472},
    "m":   {"start": 26523, "end": 27191},
    "6":   {"start": 27202, "end": 27387},
    "7a":  {"start": 27394, "end": 27759},
    "7b":  {"start": 27756, "end": 27887},
    "8":   {"start": 27894, "end": 28259},
    "n":   {"start": 28274, "end": 29533},
    "10":  {"start": 29558, "end": 29674},
}

GENE_ALIASES = {
    # canonical
    "1ab": "1ab",
    "s": "s",
    "e": "e",
    "m": "m",
    "n": "n",
    "3a": "3a",
    "6": "6",
    "7a": "7a",
    "7b": "7b",
    "8": "8",
    "10": "10",

    # common aliases in constellations
    "spike": "s",
    "orf3a": "3a",
    "orf6": "6",
    "orf7a": "7a",
    "orf7b": "7b",
    "orf8": "8",
    "orf10": "10",
    "envelope": "e",
    "membrane": "m",
    "nucleocapsid": "n",

    # polyprotein label
    "orf1ab": "1ab",
}

# Standard genetic code
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# Regex patterns
NUC_SUB_RE = re.compile(r"^nuc:([ACGT])(\d+)([ACGT])$", re.IGNORECASE)
NUC_DEL_RE = re.compile(r"^nuc:([ACGT]+)(\d+)-$", re.IGNORECASE)
NUC_INS_RE = re.compile(r"^nuc:(\d+)\+([ACGT]+)$", re.IGNORECASE)

AA_SUB_RE = re.compile(r"^([A-Za-z0-9]+):([A-Z\*])(\d+)([A-Z\*])$", re.IGNORECASE)
AA_DEL_RE = re.compile(r"^([A-Za-z0-9]+):([A-Z]+)(\d+)-$", re.IGNORECASE)

# explicit nucleotide deletion syntax in many constellations
DEL_LEN_RE = re.compile(r"^del:(\d+):(\d+)$", re.IGNORECASE)

def normalize_gene_name(gene: str) -> str:
    g = gene.strip().lower()
    return GENE_ALIASES.get(g, g)

def load_fasta_single(path: str) -> str:
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)

def infer_lineage_name(obj: dict, filename: str) -> str:
    variant = obj.get("variant", {})
    pango = variant.get("Pango_lineages", [])
    if isinstance(pango, list) and pango:
        return pango[0]

    mrca = variant.get("mrca_lineage")
    if mrca:
        return mrca

    stem = Path(filename).stem
    return stem[1:] if stem.startswith("c") else stem

def aa_index_to_genome_positions(gene: str, aa_index: int):
    """
    aa_index is 1-based.
    Returns the 3 genomic positions for that codon.
    """
    g = GENES[gene]
    codon_start = g["start"] + (aa_index - 1) * 3
    return [codon_start, codon_start + 1, codon_start + 2]

def get_ref_codon(ref: str, gene: str, aa_index: int):
    pos = aa_index_to_genome_positions(gene, aa_index)
    codon = "".join(ref[p - 1] for p in pos)
    return codon, pos

def best_alt_codons(ref_codon: str, alt_aa: str):
    """
    Return all codons encoding alt_aa with minimal Hamming distance to ref_codon.
    """
    candidates = AA_TO_CODONS.get(alt_aa, [])
    if not candidates:
        return []

    dists = [(sum(a != b for a, b in zip(ref_codon, c)), c) for c in candidates]
    min_dist = min(d for d, _ in dists)
    return [c for d, c in dists if d == min_dist]

def codon_to_nt_events(ref_positions, ref_codon: str, alt_codon: str):
    events = []
    for p, r, a in zip(ref_positions, ref_codon, alt_codon):
        if r != a:
            events.append({"pos": p, "ref": r, "alt": a})
    return events

def deletion_event_from_ref(ref: str, pos: int, length: int):
    """
    Return a deletion event using 1-based genomic coordinates.
    """
    ref_nt = ref[pos - 1: pos - 1 + length]
    if len(ref_nt) != length:
        return None
    return {"pos": pos, "ref": ref_nt, "alt": "-"}

def parse_site_to_nt_events(site: str, ref: str):
    """
    Returns list of nt-event dicts:
      {"pos": int, "ref": str, "alt": str, "raw_site": str, "event_group": str, "kind": str}
    """

    site = site.strip()

    # -------------------------
    # Direct nucleotide SNV
    # e.g. nuc:C913T
    # -------------------------
    m = NUC_SUB_RE.match(site)
    if m:
        refb, pos, altb = m.groups()
        return [{
            "pos": int(pos),
            "ref": refb.upper(),
            "alt": altb.upper(),
            "raw_site": site,
            "event_group": site,
            "kind": "nuc_snv"
        }]

    # -------------------------
    # Direct nucleotide deletion
    # e.g. nuc:AT21991-
    # -------------------------
    m = NUC_DEL_RE.match(site)
    if m:
        refb, pos = m.groups()
        return [{
            "pos": int(pos),
            "ref": refb.upper(),
            "alt": "-",
            "raw_site": site,
            "event_group": site,
            "kind": "nuc_del"
        }]

    # -------------------------
    # Direct nucleotide insertion
    # e.g. nuc:22205+GAGCCAGAA
    # -------------------------
    m = NUC_INS_RE.match(site)
    if m:
        pos, altb = m.groups()
        return [{
            "pos": int(pos),
            "ref": "-",
            "alt": f"+{altb.upper()}",
            "raw_site": site,
            "event_group": site,
            "kind": "nuc_ins"
        }]

    # -------------------------
    # Explicit deletion by start:length
    # e.g. del:28271:1
    # -------------------------
    m = DEL_LEN_RE.match(site)
    if m:
        pos, length = m.groups()
        pos = int(pos)
        length = int(length)
        ev = deletion_event_from_ref(ref, pos, length)
        if ev is None:
            return []
        ev.update({
            "raw_site": site,
            "event_group": site,
            "kind": "del_len"
        })
        return [ev]

    # -------------------------
    # Amino-acid substitution
    # e.g. spike:N501Y, ORF8:E106*, s:A67V
    # -------------------------
    m = AA_SUB_RE.match(site)
    if m:
        gene, ref_aa, aa_idx, alt_aa = m.groups()
        gene = normalize_gene_name(gene)
        aa_idx = int(aa_idx)
        ref_aa = ref_aa.upper()
        alt_aa = alt_aa.upper()

        if gene not in GENES:
            return []

        ref_codon, positions = get_ref_codon(ref, gene, aa_idx)
        translated = CODON_TABLE.get(ref_codon, "?")

        # Do not hard fail, but skip clearly inconsistent/unsupported codons
        if translated == "?":
            return []

        alt_codons = best_alt_codons(ref_codon, alt_aa)
        if not alt_codons:
            return []

        out = []
        for alt_codon in alt_codons:
            events = codon_to_nt_events(positions, ref_codon, alt_codon)
            for ev in events:
                ev.update({
                    "raw_site": site,
                    "event_group": site,
                    "kind": "aa_sub",
                })
            out.extend(events)
        return out

    # -------------------------
    # Amino-acid deletion
    # e.g. s:HV69-, orf1ab:SGF3675-
    # -------------------------
    m = AA_DEL_RE.match(site)
    if m:
        gene, aa_string, aa_idx = m.groups()
        gene = normalize_gene_name(gene)
        aa_idx = int(aa_idx)
        aa_string = aa_string.upper()

        if gene not in GENES:
            return []

        n_aa = len(aa_string)
        start_nt = aa_index_to_genome_positions(gene, aa_idx)[0]
        del_len = n_aa * 3
        ref_nt = ref[start_nt - 1: start_nt - 1 + del_len]

        if len(ref_nt) != del_len:
            return []

        return [{
            "pos": start_nt,
            "ref": ref_nt,
            "alt": "-",
            "raw_site": site,
            "event_group": site,
            "kind": "aa_del"
        }]

    # unsupported syntax
    return []

def parse_constellation_json(path: str, ref: str):
    with open(path) as fh:
        obj = json.load(fh)

    lineage = infer_lineage_name(obj, path)
    sites = obj.get("sites", [])

    rows = []
    skipped = []

    for site in sites:
        events = parse_site_to_nt_events(site, ref)
        if not events:
            skipped.append(site)
            continue

        for ev in events:
            rows.append({
                "lineage": lineage,
                "pos": ev["pos"],
                "ref": ev["ref"],
                "alt": ev["alt"],
                "is_defining": 1,
                "raw_site": ev["raw_site"],
                "event_group": ev["event_group"],
                "kind": ev["kind"],
                "source_file": os.path.basename(path),
            })

    return rows, skipped

def main():
    ap = argparse.ArgumentParser(
        description="Build lineage_defs.tsv from SARS-CoV-2 constellation JSON files"
    )
    ap.add_argument("--input-dir", required=True, help="Folder with constellation JSON files")
    ap.add_argument("--reference-fasta", required=True, help="Reference FASTA for MN908947.3")
    ap.add_argument("--output", default="lineage_defs.tsv", help="Minimal output TSV")
    ap.add_argument("--full-output", default="lineage_defs.full.tsv", help="Full annotated TSV")
    ap.add_argument("--skipped-output", default="skipped_sites.tsv", help="Skipped/unparsed sites")
    ap.add_argument(
        "--allowed-lineages",
        default=None,
        help="Optional text file with one lineage per line to retain"
    )
    args = ap.parse_args()

    ref = load_fasta_single(args.reference_fasta)

    allowed = None
    if args.allowed_lineages:
        with open(args.allowed_lineages) as fh:
            allowed = {x.strip() for x in fh if x.strip()}

    all_rows = []
    skipped_rows = []

    json_files = sorted(glob.glob(os.path.join(args.input_dir, "*.json")))
    if not json_files:
        raise SystemExit(f"No JSON files found in {args.input_dir}")

    for path in json_files:
        rows, skipped = parse_constellation_json(path, ref)

        if allowed is not None:
            rows = [r for r in rows if r["lineage"] in allowed]
            skipped = skipped  # keep skipped report unfiltered

        all_rows.extend(rows)

        for s in skipped:
            skipped_rows.append({
                "source_file": os.path.basename(path),
                "raw_site": s
            })

    if not all_rows:
        raise SystemExit("No lineage definitions recovered.")

    df = pd.DataFrame(all_rows).sort_values(
        ["lineage", "pos", "ref", "alt", "raw_site"]
    )

    # full output with provenance
    df.to_csv(args.full_output, sep="\t", index=False)

    # minimal output for downstream Bal-style script
    minimal = (
        df[["lineage", "pos", "ref", "alt", "is_defining"]]
        .drop_duplicates()
        .sort_values(["lineage", "pos", "ref", "alt"])
    )
    minimal.to_csv(args.output, sep="\t", index=False)

    print(f"[OK] Full output: {args.full_output} ({len(df)} rows)")
    print(f"[OK] Minimal output: {args.output} ({len(minimal)} rows)")

    if skipped_rows:
        skipped_df = pd.DataFrame(skipped_rows)
        skipped_df.to_csv(args.skipped_output, sep="\t", index=False)
        print(f"[INFO] Skipped sites: {args.skipped_output} ({len(skipped_df)} rows)")
    else:
        print("[INFO] No skipped sites")

if __name__ == "__main__":
    main()
