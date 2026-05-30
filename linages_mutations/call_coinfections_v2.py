#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd


def normalize_alt(alt: str) -> str:
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


def build_lineage_maps(lineage_defs: pd.DataFrame):
    lineage_to_muts = {}
    for lineage, sub in lineage_defs.groupby("lineage"):
        lineage_to_muts[lineage] = set(sub["mut_key"])

    mut_to_lineages = defaultdict(set)
    for lineage, muts in lineage_to_muts.items():
        for m in muts:
            mut_to_lineages[m].add(lineage)

    return lineage_to_muts, mut_to_lineages


def filter_lineage_defs_by_sharing(lineage_defs: pd.DataFrame, max_lineages_per_site: int):
    """
    Remove mutations shared by too many lineages.
    """
    _, mut_to_lineages = build_lineage_maps(lineage_defs)
    keep_mut = {m for m, ls in mut_to_lineages.items() if len(ls) <= max_lineages_per_site}
    return lineage_defs[lineage_defs["mut_key"].isin(keep_mut)].copy()


def load_sample_variants(
    path: str,
    lineage_defs: pd.DataFrame,
    min_total_dp: int = 100,
    min_alt_dp: int = 50,
    pass_only: bool = True,
    major_af_min: float = 0.7,
    keep_only_lineage_positions: bool = True,
    top_alt_per_pos: bool = True,
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

    if keep_only_lineage_positions:
        lineage_positions = set(lineage_defs["pos"].astype(int))
        df = df[df["POS"].isin(lineage_positions)].copy()

    if top_alt_per_pos and not df.empty:
        # keep highest AF ALT per genomic position
        df = df.sort_values(["POS", "ALT_FREQ", "ALT_DP"], ascending=[True, False, False])
        df = df.drop_duplicates(subset=["POS"], keep="first").copy()

    # separate major-supporting and minor-range calls later, but keep all filtered lineage-position variants
    return df


def build_sample_maps(sample_df: pd.DataFrame):
    af_map = dict(zip(sample_df["mut_key"], sample_df["ALT_FREQ"]))
    alt_dp_map = dict(zip(sample_df["mut_key"], sample_df["ALT_DP"]))
    total_dp_map = dict(zip(sample_df["mut_key"], sample_df["TOTAL_DP"]))
    pos_to_mut = dict(zip(sample_df["POS"], sample_df["mut_key"]))
    return af_map, alt_dp_map, total_dp_map, pos_to_mut


def score_major_lineages(sample_df: pd.DataFrame, lineage_to_muts: dict, major_af_min: float = 0.7):
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


def score_minor_lineages(
    sample_df: pd.DataFrame,
    lineage_to_muts: dict,
    major_lineage: str,
    minor_af_min: float = 0.10,
    minor_af_max: float = 0.40,
    min_minor_support: int = 5,
    min_minor_median_af: float = 0.12,
):
    af_map, alt_dp_map, total_dp_map, _ = build_sample_maps(sample_df)

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


def confidence_tier(n_minor_supported: int, minor_median_af: float) -> str:
    if n_minor_supported >= 8 and minor_median_af >= 0.15:
        return "high"
    if n_minor_supported >= 6 and minor_median_af >= 0.12:
        return "medium"
    if n_minor_supported >= 5 and minor_median_af >= 0.10:
        return "low"
    return "none"


def analyze_sample(
    sample_path: str,
    lineage_defs: pd.DataFrame,
    lineage_to_muts: dict,
    min_total_dp: int = 100,
    min_alt_dp: int = 50,
    pass_only: bool = True,
    snv_only: bool = False,
    keep_only_lineage_positions: bool = True,
    top_alt_per_pos: bool = True,
    major_af_min: float = 0.7,
    major_min_supported: int = 8,
    major_min_frac: float = 0.4,
    minor_af_min: float = 0.10,
    minor_af_max: float = 0.40,
    min_minor_support: int = 5,
    min_minor_median_af: float = 0.12,
):
    sample_name = Path(sample_path).stem

    sample_df = load_sample_variants(
        path=sample_path,
        lineage_defs=lineage_defs,
        min_total_dp=min_total_dp,
        min_alt_dp=min_alt_dp,
        pass_only=pass_only,
        major_af_min=major_af_min,
        keep_only_lineage_positions=keep_only_lineage_positions,
        top_alt_per_pos=top_alt_per_pos,
        snv_only=snv_only,
    )

    if sample_df.empty:
        return {
            "sample": sample_name,
            "n_filtered_variants": 0,
            "major_lineage": None,
            "major_n_supported": 0,
            "major_frac_supported": 0.0,
            "major_status": "no_data",
            "co_infection_flag": False,
            "minor_lineage": None,
            "minor_n_supported": 0,
            "minor_median_af": 0.0,
            "confidence_tier": "none",
            "minor_supported_mutations": "",
        }, pd.DataFrame(), pd.DataFrame(), sample_df

    major_scores = score_major_lineages(
        sample_df=sample_df,
        lineage_to_muts=lineage_to_muts,
        major_af_min=major_af_min,
    )

    best_major = major_scores.iloc[0].to_dict()
    major_lineage = best_major["lineage"]
    major_n_supported = int(best_major["n_major_supported"])
    major_frac_supported = float(best_major["frac_major_supported"])

    if major_n_supported < major_min_supported or major_frac_supported < major_min_frac:
        summary = {
            "sample": sample_name,
            "n_filtered_variants": len(sample_df),
            "major_lineage": None,
            "major_n_supported": major_n_supported,
            "major_frac_supported": major_frac_supported,
            "major_status": "unassigned_low_support",
            "co_infection_flag": False,
            "minor_lineage": None,
            "minor_n_supported": 0,
            "minor_median_af": 0.0,
            "confidence_tier": "none",
            "minor_supported_mutations": "",
        }
        return summary, major_scores, pd.DataFrame(), sample_df

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
        "major_n_supported": major_n_supported,
        "major_frac_supported": major_frac_supported,
        "major_status": "assigned",
        "co_infection_flag": co_flag,
        "minor_lineage": minor_lineage if co_flag else None,
        "minor_n_supported": minor_n if co_flag else 0,
        "minor_median_af": minor_med_af if co_flag else 0.0,
        "confidence_tier": tier if co_flag else "none",
        "minor_supported_mutations": "|".join(minor_supported_mutations) if co_flag else "",
    }

    return summary, major_scores, minor_scores, sample_df


def main():
    ap = argparse.ArgumentParser(description="Strict caller for putative SARS-CoV-2 co-infections")
    ap.add_argument("--variants-dir", required=True, help="Folder with per-sample TSV mutation tables")
    ap.add_argument("--lineage-defs", required=True, help="lineage_defs.tsv")
    ap.add_argument("--glob", default="*.tsv", help="Glob for sample files")
    ap.add_argument("--outdir", default="coinfection_calls_strict", help="Output folder")

    ap.add_argument("--min-total-dp", type=int, default=100)
    ap.add_argument("--min-alt-dp", type=int, default=50)
    ap.add_argument("--no-pass-filter", action="store_true")
    ap.add_argument("--snv-only", action="store_true")

    ap.add_argument("--major-af-min", type=float, default=0.7)
    ap.add_argument("--major-min-supported", type=int, default=8)
    ap.add_argument("--major-min-frac", type=float, default=0.4)

    ap.add_argument("--minor-af-min", type=float, default=0.10)
    ap.add_argument("--minor-af-max", type=float, default=0.40)
    ap.add_argument("--min-minor-support", type=int, default=5)
    ap.add_argument("--min-minor-median-af", type=float, default=0.12)

    ap.add_argument("--max-lineages-per-site", type=int, default=3,
                    help="Drop lineage-defining mutations shared by more than this many lineages")
    ap.add_argument("--disable-lineage-position-filter", action="store_true",
                    help="Do not restrict sample variants to lineage-defining positions")
    ap.add_argument("--disable-top-alt-per-pos", action="store_true",
                    help="Do not keep only top ALT per position")

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lineage_defs = load_lineage_defs(args.lineage_defs)
    lineage_defs = filter_lineage_defs_by_sharing(
        lineage_defs=lineage_defs,
        max_lineages_per_site=args.max_lineages_per_site,
    )
    lineage_to_muts, _ = build_lineage_maps(lineage_defs)

    sample_paths = sorted(glob.glob(os.path.join(args.variants_dir, args.glob)))
    if not sample_paths:
        raise SystemExit(f"No sample files found in {args.variants_dir} matching {args.glob}")

    summaries = []

    lineage_defs.to_csv(os.path.join(args.outdir, "lineage_defs.filtered.tsv"), sep="\t", index=False)

    for sample_path in sample_paths:
        sample_name = Path(sample_path).stem
        print(f"[INFO] analyzing {sample_name}")

        try:
            summary, major_scores, minor_scores, sample_df = analyze_sample(
                sample_path=sample_path,
                lineage_defs=lineage_defs,
                lineage_to_muts=lineage_to_muts,
                min_total_dp=args.min_total_dp,
                min_alt_dp=args.min_alt_dp,
                pass_only=not args.no_pass_filter,
                snv_only=args.snv_only,
                keep_only_lineage_positions=not args.disable_lineage_position_filter,
                top_alt_per_pos=not args.disable_top_alt_per_pos,
                major_af_min=args.major_af_min,
                major_min_supported=args.major_min_supported,
                major_min_frac=args.major_min_frac,
                minor_af_min=args.minor_af_min,
                minor_af_max=args.minor_af_max,
                min_minor_support=args.min_minor_support,
                min_minor_median_af=args.min_minor_median_af,
            )
            summaries.append(summary)

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
                "major_status": "error",
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
