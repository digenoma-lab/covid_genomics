#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd


# =========================
# Helpers
# =========================

def normalize_alt(alt: str) -> str:
    """
    Keep the same representation used in your mutation tables and lineage_defs:
      SNV: A / C / G / T
      deletion: -
      insertion: +SEQ
    """
    return str(alt).strip().upper()


def normalize_ref(ref: str) -> str:
    return str(ref).strip().upper()


def mutation_key(pos, ref, alt):
    return f"{int(pos)}:{normalize_ref(ref)}>{normalize_alt(alt)}"


def load_lineage_defs(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"lineage", "pos", "ref", "alt"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in lineage_defs: {missing}")

    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].map(normalize_ref)
    df["alt"] = df["alt"].map(normalize_alt)
    df["mut_key"] = df.apply(lambda r: mutation_key(r["pos"], r["ref"], r["alt"]), axis=1)

    if "is_defining" not in df.columns:
        df["is_defining"] = 1

    df = df.drop_duplicates(subset=["lineage", "mut_key"]).copy()
    return df


def load_sample_variants(
    path: str,
    min_total_dp: int = 100,
    min_alt_dp: int = 20,
    pass_only: bool = True,
    snv_only: bool = False,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    required = {"POS", "REF", "ALT", "ALT_FREQ", "ALT_DP", "TOTAL_DP", "PASS"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {missing}")

    df["POS"] = df["POS"].astype(int)
    df["REF"] = df["REF"].map(normalize_ref)
    df["ALT"] = df["ALT"].map(normalize_alt)
    df["PASS"] = df["PASS"].astype(str).str.upper()
    df["ALT_FREQ"] = df["ALT_FREQ"].astype(float)
    df["ALT_DP"] = df["ALT_DP"].astype(float)
    df["TOTAL_DP"] = df["TOTAL_DP"].astype(float)

    if pass_only:
        df = df[df["PASS"] == "TRUE"]

    df = df[df["TOTAL_DP"] >= min_total_dp]
    df = df[df["ALT_DP"] >= min_alt_dp]

    if snv_only:
        df = df[
            (df["REF"].str.len() == 1) &
            (df["ALT"].str.len() == 1) &
            (~df["ALT"].str.startswith("+")) &
            (df["ALT"] != "-")
        ]

    df["mut_key"] = df.apply(lambda r: mutation_key(r["POS"], r["REF"], r["ALT"]), axis=1)

    # keep highest ALT_FREQ if duplicates exist
    df = df.sort_values(["mut_key", "ALT_FREQ", "ALT_DP"], ascending=[True, False, False])
    df = df.drop_duplicates(subset=["mut_key"]).copy()

    return df


def lineage_variant_sets(lineage_defs: pd.DataFrame):
    """
    Returns:
      lineage_to_muts: lineage -> set(mut_key)
      mut_to_lineages: mut_key -> set(lineages)
    """
    lineage_to_muts = {}
    for lineage, sub in lineage_defs.groupby("lineage"):
        lineage_to_muts[lineage] = set(sub["mut_key"])

    mut_to_lineages = defaultdict(set)
    for lineage, muts in lineage_to_muts.items():
        for m in muts:
            mut_to_lineages[m].add(lineage)

    return lineage_to_muts, mut_to_lineages


def build_sample_maps(sample_df: pd.DataFrame):
    af_map = dict(zip(sample_df["mut_key"], sample_df["ALT_FREQ"]))
    alt_dp_map = dict(zip(sample_df["mut_key"], sample_df["ALT_DP"]))
    total_dp_map = dict(zip(sample_df["mut_key"], sample_df["TOTAL_DP"]))
    return af_map, alt_dp_map, total_dp_map


# =========================
# Major lineage
# =========================

def score_major_lineages(
    sample_df: pd.DataFrame,
    lineage_to_muts: dict,
    major_af_min: float = 0.6,
):
    """
    Major lineage score = overlap between high-AF sample mutations and lineage-defining mutations.
    """
    high_af = set(sample_df.loc[sample_df["ALT_FREQ"] >= major_af_min, "mut_key"])
    rows = []

    for lineage, muts in lineage_to_muts.items():
        overlap = high_af & muts
        rows.append({
            "lineage": lineage,
            "n_defining_sites": len(muts),
            "n_major_supported": len(overlap),
            "frac_major_supported": len(overlap) / len(muts) if len(muts) else 0.0,
            "major_supported_mutations": sorted(overlap),
        })

    out = pd.DataFrame(rows).sort_values(
        ["n_major_supported", "frac_major_supported", "n_defining_sites", "lineage"],
        ascending=[False, False, False, True]
    ).reset_index(drop=True)
    return out


# =========================
# Minor lineage
# =========================

def score_minor_lineages(
    sample_df: pd.DataFrame,
    lineage_to_muts: dict,
    major_lineage: str,
    minor_af_min: float = 0.05,
    minor_af_max: float = 0.5,
    min_minor_support: int = 3,
    min_minor_median_af: float = 0.08,
):
    """
    Minor lineage score uses mutations specific to candidate lineage relative to the major lineage.
    """
    af_map, alt_dp_map, total_dp_map = build_sample_maps(sample_df)

    major_set = lineage_to_muts[major_lineage]

    minor_df = sample_df[
        (sample_df["ALT_FREQ"] >= minor_af_min) &
        (sample_df["ALT_FREQ"] <= minor_af_max)
    ].copy()
    observed_minor = set(minor_df["mut_key"])

    rows = []

    for lineage, muts in lineage_to_muts.items():
        if lineage == major_lineage:
            continue

        candidate_specific = muts - major_set
        supported = sorted(candidate_specific & observed_minor)

        if supported:
            afs = [af_map[m] for m in supported]
            alt_dps = [alt_dp_map[m] for m in supported]
            total_dps = [total_dp_map[m] for m in supported]
            median_af = float(pd.Series(afs).median())
            mean_af = float(pd.Series(afs).mean())
        else:
            afs = []
            alt_dps = []
            total_dps = []
            median_af = 0.0
            mean_af = 0.0

        passes = (len(supported) >= min_minor_support and median_af >= min_minor_median_af)

        rows.append({
            "candidate_minor_lineage": lineage,
            "n_specific_sites": len(candidate_specific),
            "n_minor_supported": len(supported),
            "minor_median_af": median_af,
            "minor_mean_af": mean_af,
            "minor_supported_mutations": supported,
            "minor_supported_alt_dp": alt_dps,
            "minor_supported_total_dp": total_dps,
            "passes_minor_rule": passes,
        })

    out = pd.DataFrame(rows).sort_values(
        ["passes_minor_rule", "n_minor_supported", "minor_median_af", "minor_mean_af", "candidate_minor_lineage"],
        ascending=[False, False, False, False, True]
    ).reset_index(drop=True)

    return out


# =========================
# Optional: confidence tier
# =========================

def confidence_tier(n_minor_supported: int, minor_median_af: float) -> str:
    if n_minor_supported >= 6 and minor_median_af >= 0.12:
        return "high"
    if n_minor_supported >= 4 and minor_median_af >= 0.08:
        return "medium"
    if n_minor_supported >= 3 and minor_median_af >= 0.05:
        return "low"
    return "none"


# =========================
# Single sample analysis
# =========================

def analyze_sample(
    sample_path: str,
    lineage_defs: pd.DataFrame,
    min_total_dp: int = 100,
    min_alt_dp: int = 20,
    pass_only: bool = True,
    snv_only: bool = False,
    major_af_min: float = 0.6,
    minor_af_min: float = 0.05,
    minor_af_max: float = 0.5,
    min_minor_support: int = 3,
    min_minor_median_af: float = 0.08,
):
    sample_name = Path(sample_path).stem

    sample_df = load_sample_variants(
        sample_path,
        min_total_dp=min_total_dp,
        min_alt_dp=min_alt_dp,
        pass_only=pass_only,
        snv_only=snv_only,
    )

    if sample_df.empty:
        return {
            "sample": sample_name,
            "n_filtered_variants": 0,
            "major_lineage": None,
            "major_n_supported": 0,
            "major_frac_supported": 0.0,
            "co_infection_flag": False,
            "minor_lineage": None,
            "minor_n_supported": 0,
            "minor_median_af": 0.0,
            "confidence_tier": "none",
        }, pd.DataFrame(), pd.DataFrame(), sample_df

    lineage_to_muts, _ = lineage_variant_sets(lineage_defs)

    major_scores = score_major_lineages(
        sample_df=sample_df,
        lineage_to_muts=lineage_to_muts,
        major_af_min=major_af_min,
    )

    best_major = major_scores.iloc[0].to_dict()
    major_lineage = best_major["lineage"]

    minor_scores = score_minor_lineages(
        sample_df=sample_df,
        lineage_to_muts=lineage_to_muts,
        major_lineage=major_lineage,
        minor_af_min=minor_af_min,
        minor_af_max=minor_af_max,
        min_minor_support=min_minor_support,
        min_minor_median_af=min_minor_median_af,
    )

    if len(minor_scores):
        best_minor = minor_scores.iloc[0].to_dict()
        co_flag = bool(best_minor["passes_minor_rule"])
        minor_lineage = best_minor["candidate_minor_lineage"]
        minor_n = int(best_minor["n_minor_supported"])
        minor_med_af = float(best_minor["minor_median_af"])
        tier = confidence_tier(minor_n, minor_med_af)
        minor_supported_mutations = best_minor["minor_supported_mutations"]
    else:
        co_flag = False
        minor_lineage = None
        minor_n = 0
        minor_med_af = 0.0
        tier = "none"
        minor_supported_mutations = []

    summary = {
        "sample": sample_name,
        "n_filtered_variants": len(sample_df),
        "major_lineage": major_lineage,
        "major_n_supported": int(best_major["n_major_supported"]),
        "major_frac_supported": float(best_major["frac_major_supported"]),
        "co_infection_flag": co_flag,
        "minor_lineage": minor_lineage,
        "minor_n_supported": minor_n,
        "minor_median_af": minor_med_af,
        "confidence_tier": tier,
        "minor_supported_mutations": "|".join(minor_supported_mutations),
    }

    return summary, major_scores, minor_scores, sample_df


# =========================
# Batch mode
# =========================

def main():
    ap = argparse.ArgumentParser(description="Call putative SARS-CoV-2 co-infections from sample variant tables")
    ap.add_argument("--variants-dir", required=True, help="Folder with per-sample TSV mutation tables")
    ap.add_argument("--lineage-defs", required=True, help="lineage_defs.tsv")
    ap.add_argument("--glob", default="*.tsv", help="Glob for sample files inside variants-dir")
    ap.add_argument("--outdir", default="coinfection_calls", help="Output folder")

    ap.add_argument("--min-total-dp", type=int, default=100)
    ap.add_argument("--min-alt-dp", type=int, default=20)
    ap.add_argument("--no-pass-filter", action="store_true")
    ap.add_argument("--snv-only", action="store_true")

    ap.add_argument("--major-af-min", type=float, default=0.6)
    ap.add_argument("--minor-af-min", type=float, default=0.05)
    ap.add_argument("--minor-af-max", type=float, default=0.5)
    ap.add_argument("--min-minor-support", type=int, default=3)
    ap.add_argument("--min-minor-median-af", type=float, default=0.08)

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lineage_defs = load_lineage_defs(args.lineage_defs)

    sample_paths = sorted(glob.glob(os.path.join(args.variants_dir, args.glob)))
    if not sample_paths:
        raise SystemExit(f"No sample files found in {args.variants_dir} matching {args.glob}")

    summaries = []

    for sample_path in sample_paths:
        sample_name = Path(sample_path).stem
        print(f"[INFO] analyzing {sample_name}")

        try:
            summary, major_scores, minor_scores, sample_df = analyze_sample(
                sample_path=sample_path,
                lineage_defs=lineage_defs,
                min_total_dp=args.min_total_dp,
                min_alt_dp=args.min_alt_dp,
                pass_only=not args.no_pass_filter,
                snv_only=args.snv_only,
                major_af_min=args.major_af_min,
                minor_af_min=args.minor_af_min,
                minor_af_max=args.minor_af_max,
                min_minor_support=args.min_minor_support,
                min_minor_median_af=args.min_minor_median_af,
            )
            summaries.append(summary)

            # save per-sample details
            sample_prefix = os.path.join(args.outdir, sample_name)
            major_scores.to_csv(f"{sample_prefix}.major_scores.tsv", sep="\t", index=False)
            minor_scores.to_csv(f"{sample_prefix}.minor_scores.tsv", sep="\t", index=False)
            sample_df.to_csv(f"{sample_prefix}.filtered_variants.tsv", sep="\t", index=False)

        except Exception as e:
            summaries.append({
                "sample": sample_name,
                "n_filtered_variants": 0,
                "major_lineage": None,
                "major_n_supported": 0,
                "major_frac_supported": 0.0,
                "co_infection_flag": False,
                "minor_lineage": None,
                "minor_n_supported": 0,
                "minor_median_af": 0.0,
                "confidence_tier": "error",
                "minor_supported_mutations": str(e),
            })

    cohort = pd.DataFrame(summaries).sort_values(
        ["co_infection_flag", "confidence_tier", "minor_n_supported", "minor_median_af", "sample"],
        ascending=[False, True, False, False, True]
    )

    cohort.to_csv(os.path.join(args.outdir, "cohort_coinfection_summary.tsv"), sep="\t", index=False)

    print(f"[OK] wrote {os.path.join(args.outdir, 'cohort_coinfection_summary.tsv')}")


if __name__ == "__main__":
    main()
