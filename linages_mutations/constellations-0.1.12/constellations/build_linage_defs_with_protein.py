#!/usr/bin/env python3

import os
import re
import json
import glob
import argparse
from pathlib import Path
from itertools import product
import pandas as pd

# =========================
# SARS-CoV-2 reference model
# NC_045512.2 / Wuhan-Hu-1
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

# Parsers
NUC_SUB_RE = re.compile(r"^nuc:([ACGT])(\d+)([ACGT])$")
AA_SUB_RE  = re.compile(r"^([A-Za-z0-9]+):([A-Z\*])(\d+)([A-Z\*])$")
AA_DEL_RE  = re.compile(r"^([A-Za-z0-9]+):([A-Z]+)(\d+)-$")

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def load_fasta_single(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)

def infer_lineage_name(obj, filename):
    variant = obj.get("variant", {})
    pango = variant.get("Pango_lineages", [])
    if isinstance(pango, list) and pango:
        return pango[0]
    mrca = variant.get("mrca_lineage")
    if mrca:
        return mrca
    stem = Path(filename).stem
    return stem[1:] if stem.startswith("c") else stem

def cds_seq(ref, gene):
    g = GENES[gene]
    return ref[g["start"]-1:g["end"]]

def aa_index_to_genome_positions(gene, aa_index):
    """
    aa_index is 1-based.
    Returns 3 genomic positions (1-based) for that codon.
    """
    g = GENES[gene]
    codon_start = g["start"] + (aa_index - 1) * 3
    return [codon_start, codon_start + 1, codon_start + 2]

def get_ref_codon(ref, gene, aa_index):
    pos = aa_index_to_genome_positions(gene, aa_index)
    codon = "".join(ref[p-1] for p in pos)
    return codon, pos

def best_alt_codons(ref_codon, alt_aa):
    """
    Return all minimal-Hamming-distance codons encoding alt_aa.
    """
    candidates = AA_TO_CODONS.get(alt_aa, [])
    if not candidates:
        return []
    dists = [(sum(a != b for a, b in zip(ref_codon, c)), c) for c in candidates]
    mind = min(d for d, _ in dists)
    return [c for d, c in dists if d == mind]

def codon_to_nt_events(ref_positions, ref_codon, alt_codon):
    events = []
    for p, r, a in zip(ref_positions, ref_codon, alt_codon):
        if r != a:
            events.append({"pos": p, "ref": r, "alt": a})
    return events

def parse_site_to_nt_events(site, ref):
    """
    Returns list of nt-event dicts:
      {"pos": int, "ref": str, "alt": str, "raw_site": str, "event_group": str}
    """
    # nuc direct
    m = NUC_SUB_RE.match(site)
    if m:
        refb, pos, altb = m.groups()
        return [{
            "pos": int(pos), "ref": refb, "alt": altb,
            "raw_site": site, "event_group": site, "kind": "nuc"
        }]

    # aa substitution like s:N501Y, 8:Q27*
    m = AA_SUB_RE.match(site)
    if m:
        gene, ref_aa, aa_idx, alt_aa = m.groups()
        gene = gene.lower()
        aa_idx = int(aa_idx)
        if gene not in GENES:
            return []

        ref_codon, positions = get_ref_codon(ref, gene, aa_idx)
        translated = CODON_TABLE.get(ref_codon, "?")
        # Do not hard fail; constellations can occasionally refer to mature peptides in ORF1ab.
        alt_codons = best_alt_codons(ref_codon, alt_aa)
        out = []
        for alt_codon in alt_codons:
            events = codon_to_nt_events(positions, ref_codon, alt_codon)
            for ev in events:
                ev.update({"raw_site": site, "event_group": site, "kind": "aa_sub"})
            out.extend(events)
        return out

    # aa deletion like s:HV69-, s:Y144-, 1ab:SGF3675-
    m = AA_DEL_RE.match(site)
    if m:
        gene, aa_string, aa_idx = m.groups()
        gene = gene.lower()
        aa_idx = int(aa_idx)
        if gene not in GENES:
            return []

        n_aa = len(aa_string)
        start_nt = aa_index_to_genome_positions(gene, aa_idx)[0]
        del_len = n_aa * 3
        ref_nt = ref[start_nt-1:start_nt-1+del_len]
        return [{
            "pos": start_nt,
            "ref": ref_nt,
            "alt": "-",
            "raw_site": site,
            "event_group": site,
            "kind": "aa_del"
        }]

    return []

def parse_constellation_json(path, ref):
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True, help="Folder with constellation JSON files")
    ap.add_argument("--reference-fasta", required=True, help="NC_045512.2 / Wuhan-Hu-1 FASTA")
    ap.add_argument("--output", default="lineage_defs.tsv")
    ap.add_argument("--full-output", default="lineage_defs.full.tsv")
    ap.add_argument("--skipped-output", default="skipped_sites.tsv")
    ap.add_argument("--allowed-lineages", default=None,
                    help="Optional text file, one lineage per line")
    args = ap.parse_args()

    ref = load_fasta_single(args.reference_fasta)

    allowed = None
    if args.allowed_lineages:
        with open(args.allowed_lineages) as fh:
            allowed = {x.strip() for x in fh if x.strip()}

    all_rows = []
    skipped_rows = []

    for path in sorted(glob.glob(os.path.join(args.input_dir, "*.json"))):
        rows, skipped = parse_constellation_json(path, ref)
        if allowed is not None:
            rows = [r for r in rows if r["lineage"] in allowed]
        all_rows.extend(rows)
        for s in skipped:
            skipped_rows.append({"source_file": os.path.basename(path), "raw_site": s})

    df = pd.DataFrame(all_rows)
    if df.empty:
        raise SystemExit("No lineage definitions recovered.")

    # full table
    df = df.sort_values(["lineage", "pos", "ref", "alt", "raw_site"])
    df.to_csv(args.full_output, sep="\t", index=False)

    # minimal table for your Bal-style script
    minimal = (
        df[["lineage", "pos", "ref", "alt", "is_defining"]]
        .drop_duplicates()
        .sort_values(["lineage", "pos", "ref", "alt"])
    )
    minimal.to_csv(args.output, sep="\t", index=False)

    print(f"[OK] full table: {args.full_output} ({len(df)} rows)")
    print(f"[OK] minimal table: {args.output} ({len(minimal)} rows)")

    if skipped_rows:
        skipped = pd.DataFrame(skipped_rows)
        skipped.to_csv(args.skipped_output, sep="\t", index=False)
        print(f"[INFO] skipped sites: {args.skipped_output} ({len(skipped)} rows)")

if __name__ == "__main__":
    main()
