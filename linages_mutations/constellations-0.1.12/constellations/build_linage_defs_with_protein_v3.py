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

    "orf1ab": "1ab",
    "orf1a": "orf1a",
    "orf1b": "orf1b",
}

# ORF1a / ORF1b coordinates in amino-acid space converted to genomic coordinates.
# ORF1a starts at nt 266.
# ORF1b translation starts after the programmed frameshift and is conventionally
# indexed separately from amino acid 1.
ORF1A_START = 266
ORF1B_START = 13468  # genomic coordinate used for amino-acid numbering convenience

# NSP coordinates in amino-acid space within ORF1ab polyprotein
# Using standard SARS-CoV-2 nsp boundaries:
NSP_RANGES = {
    "nsp1":  (1, 180),
    "nsp2":  (181, 818),
    "nsp3":  (819, 2763),
    "nsp4":  (2764, 3263),
    "nsp5":  (3264, 3569),
    "nsp6":  (3570, 3859),
    "nsp7":  (3860, 3942),
    "nsp8":  (3943, 4140),
    "nsp9":  (4141, 4253),
    "nsp10": (4254, 4392),
    "nsp12": (4393, 5324),
    "nsp13": (5325, 5925),
    "nsp14": (5926, 6452),
    "nsp15": (6453, 6798),
    "nsp16": (6799, 7096),
}

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

NUC_SUB_RE = re.compile(r"^nuc:([ACGT])(\d+)([ACGT])$", re.IGNORECASE)
NUC_DEL_RE = re.compile(r"^nuc:([ACGT]+)(\d+)-$", re.IGNORECASE)
NUC_INS_RE = re.compile(r"^nuc:(\d+)\+([ACGT]+)$", re.IGNORECASE)

AA_SUB_RE = re.compile(r"^([A-Za-z0-9]+):([A-Z\*])(\d+)([A-Z\*])$", re.IGNORECASE)
AA_DEL_RE = re.compile(r"^([A-Za-z0-9]+):([A-Z]+)(\d+)-$", re.IGNORECASE)
DEL_LEN_RE = re.compile(r"^del:(\d+):(\d+)$", re.IGNORECASE)

# special compound events
RG203KR_RE = re.compile(r"^n:RG203KR$", re.IGNORECASE)

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

def orf1a_aa_to_genome_positions(aa_index: int):
    codon_start = ORF1A_START + (aa_index - 1) * 3
    return [codon_start, codon_start + 1, codon_start + 2]

def orf1b_aa_to_genome_positions(aa_index: int):
    codon_start = ORF1B_START + (aa_index - 1) * 3
    return [codon_start, codon_start + 1, codon_start + 2]

def nsp_aa_to_orf1ab_aa(nsp: str, aa_index: int):
    start_aa, end_aa = NSP_RANGES[nsp]
    polyprotein_aa = start_aa + aa_index - 1
    if polyprotein_aa > end_aa:
        return None
    return polyprotein_aa

def aa_index_to_genome_positions(gene: str, aa_index: int):
    if gene in GENES:
        g = GENES[gene]
        codon_start = g["start"] + (aa_index - 1) * 3
        return [codon_start, codon_start + 1, codon_start + 2]

    if gene == "orf1a":
        return orf1a_aa_to_genome_positions(aa_index)

    if gene == "orf1b":
        return orf1b_aa_to_genome_positions(aa_index)

    if gene.startswith("nsp"):
        poly_aa = nsp_aa_to_orf1ab_aa(gene, aa_index)
        if poly_aa is None:
            return None
        codon_start = ORF1A_START + (poly_aa - 1) * 3
        return [codon_start, codon_start + 1, codon_start + 2]

    return None

def get_ref_codon(ref: str, gene: str, aa_index: int):
    pos = aa_index_to_genome_positions(gene, aa_index)
    if pos is None:
        return None, None
    codon = "".join(ref[p - 1] for p in pos)
    return codon, pos

def best_alt_codons(ref_codon: str, alt_aa: str):
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
    ref_nt = ref[pos - 1: pos - 1 + length]
    if len(ref_nt) != length:
        return None
    return {"pos": pos, "ref": ref_nt, "alt": "-"}

def special_compound_event(site: str, ref: str):
    if RG203KR_RE.match(site):
        # Standard SARS-CoV-2 N RG203KR event corresponds to:
        # 28881 G>A, 28882 G>A, 28883 G>C
        return [
            {"pos": 28881, "ref": "G", "alt": "A", "raw_site": site, "event_group": site, "kind": "compound"},
            {"pos": 28882, "ref": "G", "alt": "A", "raw_site": site, "event_group": site, "kind": "compound"},
            {"pos": 28883, "ref": "G", "alt": "C", "raw_site": site, "event_group": site, "kind": "compound"},
        ]
    return None

def parse_site_to_nt_events(site: str, ref: str):
    site = site.strip()

    special = special_compound_event(site, ref)
    if special is not None:
        return special

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

    m = AA_SUB_RE.match(site)
    if m:
        gene, ref_aa, aa_idx, alt_aa = m.groups()
        gene = normalize_gene_name(gene)
        aa_idx = int(aa_idx)
        ref_aa = ref_aa.upper()
        alt_aa = alt_aa.upper()

        ref_codon, positions = get_ref_codon(ref, gene, aa_idx)
        if ref_codon is None or positions is None:
            return []

        translated = CODON_TABLE.get(ref_codon, "?")
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

    m = AA_DEL_RE.match(site)
    if m:
        gene, aa_string, aa_idx = m.groups()
        gene = normalize_gene_name(gene)
        aa_idx = int(aa_idx)
        aa_string = aa_string.upper()

        positions = aa_index_to_genome_positions(gene, aa_idx)
        if positions is None:
            return []

        start_nt = positions[0]
        del_len = len(aa_string) * 3
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
    ap = argparse.ArgumentParser(description="Build lineage_defs.tsv from constellation JSON files")
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--output", default="lineage_defs.tsv")
    ap.add_argument("--full-output", default="lineage_defs.full.tsv")
    ap.add_argument("--skipped-output", default="skipped_sites.tsv")
    ap.add_argument("--allowed-lineages", default=None)
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

    df.to_csv(args.full_output, sep="\t", index=False)

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
