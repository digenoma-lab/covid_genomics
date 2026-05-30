#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd


# =========================
# Basic normalization
# =========================

def normalize_ref(x: str) -> str:
    return str(x).strip().upper()

def normalize_alt(x: str) -> str:
    return str(x).strip().upper()

def mut_key(pos, ref, alt) -> str:
    return f"{int(pos)}:{normalize_ref(ref)}>{normalize_alt(alt)}"

def site_key(pos, ref) -> str:
    return f"{int(pos)}:{normalize_ref(ref)}"


# =========================
# Optional lineage grouping
# =========================

def lineage_group(lineage: str) -> str:
    lineage = str(lineage)
    if lineage == "B.1.617.2" or lineage.startswith("AY."):
        return "Delta"
    if lineage == "B.1.1.529" or lineage.startswith("BA."):
        return "Omicron"
    if lineage == "P.1":
        return "Gamma"
    if lineage == "C.37":
        return "Lambda"
    if lineage == "B.1.1.7":
        return "Alpha"
    if lineage == "B.1.351":
        return "Beta"
    if lineage in {"B.1.427", "B.1.429"}:
        return "Epsilon"
    if lineage == "B.1.525":
        return "Eta"
    if lineage == "B.1.526":
        return "Iota"
    if lineage == "B.1.621":
        return "Mu"
    if lineage == "P.2":
        return "Zeta"
    if lineage == "P.3":
        return "Theta"
    if lineage == "B.1.617.1":
        return "Kappa"
    return lineage


# =========================
# Load lineage definitions
# =========================

def load_lineage_defs(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"lineage", "pos", "ref", "alt"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in lineage_defs: {missing}")

    df["lineage"] = df["lineage"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].map(normalize_ref)
    df["alt"] = df["alt"].map(normalize_alt)

    if "is_defining" not in df.columns:
        df["is_defining"] = 1

    df["mut_key"] = df.apply(lambda r: mut_key(r["pos"], r["ref"], r["alt"]), axis=1)
    df["site_key"] = df.apply(lambda r: site_key(r["pos"], r["ref"]), axis=1)
    df["group"] = df["lineage"].map(lineage_group)

    return df.drop_duplicates(subset=["lineage", "mut_key"]).copy()


def filter_lineage_defs_by_site_sharing(lineage_defs: pd.DataFrame, max_lineages_per_site: int) -> pd.DataFrame:
    """
    Remove sites shared by too many exact lineages.
    Sharing is counted at the site level, not mut_key level.
    """
    site_to_lineages = (
        lineage_defs.groupby("site_key")["lineage"]
        .apply(lambda x: set(x))
        .to_dict()
    )
    keep_sites = {s for s, ls in site_to_lineages.items() if len(ls) <= max_lineages_per_site}
    return lineage_defs[lineage_defs["site_key"].isin(keep_sites)].copy()


def build_lineage_site_maps(lineage_defs: pd.DataFrame, label_col: str = "lineage"):
    """
    Returns:
      label_to_sites: label -> site_key -> set(mut_key allowed for that site)
      label_to_site_set: label -> set(site_keys)
      mut_to_labels: mut_key -> set(labels)
      site_to_labels: site_key -> set(labels)
    """
    label_to_sites = defaultdict(lambda: defaultdict(set))
    label_to_site_set = defaultdict(set)
    mut_to_labels = defaultdict(set)
    site_to_labels = defaultdict(set)

    for _, r in lineage_defs.iterrows():
        label = r[label_col]
        s = r["site_key"]
        m = r["mut_key"]
        label_to_sites[label][s].add(m)
        label_to_site_set[label].add(s)
        mut_to_labels[m].add(label)
        site_to_labels[s].add(label)

    return label_to_sites, label_to_site_set, mut_to_labels, site_to_labels


# =========================
# Load sample variants
# =========================

def load_sample_variants(
    path: str,
    lineage_defs: pd.DataFrame,
    min_total_dp: int = 100,
    min_alt_dp: int = 50,
    pass_only: bool = True,
    snv_only: bool = False,
    keep_only_lineage_sites: bool = True,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    required = {"POS", "REF", "ALT", "ALT_FREQ", "ALT_DP", "TOTAL_DP", "PASS"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {missing}")

    df["POS"] = df["POS"].astype(int)
    df["REF"] = df["REF"].map(normalize_ref)
    df["ALT"] = df["ALT"].map(normalize_alt)
    df["ALT_FREQ"] = df["ALT_FREQ"].astype(float)
    df["ALT_DP"] = df["ALT_DP"].astype(float)
    df["TOTAL_DP"] = df["TOTAL_DP"].astype(float)
    df["PASS"] = df["PASS"].astype(str).str.upper()

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

    df["mut_key"] = df.apply(lambda r: mut_key(r["POS"], r["REF"], r["ALT"]), axis=1)
    df["site_key"] = df.apply(lambda r: site_key(r["POS"], r["REF"]), axis=1)

    if keep_only_lineage_sites:
        lineage_sites = set(lineage_defs["site_key"])
        df = df[df["site_key"].isin(lineage_sites)].copy()

    # keep duplicate ALTs if distinct; only collapse exact duplicates
    df = df.sort_values(["mut_key", "ALT_FREQ", "ALT_DP"], ascending=[True, False, False])
    df = df.drop_duplicates(subset=["mut_key"], keep="first").copy()

    return df


# =========================
# Sample site support maps
# =========================

def build_sample_support_maps(sample_df: pd.DataFrame):
    """
    sample_site_to_muts: site_key -> set(mut_key observed in sample at that site)
    mut_to_af / mut_to_altdp / mut_to_totaldp
    """
    sample_site_to_muts = defaultdict(set)
    mut_to_af = {}
    mut_to_altdp = {}
    mut_to_totaldp = {}

    for _, r in sample_df.iterrows():
        m = r["mut_key"]
        s = r["site_key"]
        sample_site_to_muts[s].add(m)
        mut_to_af[m] = float(r["ALT_FREQ"])
        mut_to_altdp[m] = float(r["ALT_DP"])
        mut_to_totaldp[m] = float(r["TOTAL_DP"])

    return sample_site_to_muts, mut_to_af, mut_to_altdp, mut_to_totaldp


# =========================
# Site-level support
# =========================

def site_supported(sample_site_to_muts, allowed_muts_for_site, af_map, af_min=None, af_max=None):
    """
    Returns matching observed mut_keys at the site after AF bounds.
    """
    observed = sample_site_to_muts.get(next(iter({k for k in []}), None), set())  # unused safeguard
    matches = []
    # allowed_muts_for_site is a set(mut_key)
    # need the sample site muts for the corresponding site from caller
    # so caller must pass actual observed set
    raise RuntimeError("site_supported should not be called directly")


def matching_muts_at_site(observed_site_muts, allowed_site_muts, af_map, af_min=None, af_max=None):
    matches = []
    for m in observed_site_muts:
        if m not in allowed_site_muts:
            continue
        af = af_map[m]
        if af_min is not None and af < af_min:
            continue
        if af_max is not None and af > af_max:
            continue
        matches.append(m)
    return matches


# =========================
# Major lineage scoring by site
# =========================

def score_major_labels(
    sample_df: pd.DataFrame,
    label_to_sites: dict,
    label_to_site_set: dict,
    major_af_min: float = 0.70,
):
    sample_site_to_muts, af_map, altdp_map, totaldp_map = build_sample_support_maps(sample_df)

    rows = []
    for label, site_map in label_to_sites.items():
        supported_sites = []
        supported_mutations = []

        for s, allowed_muts in site_map.items():
            obs = sample_site_to_muts.get(s, set())
            matches = matching_muts_at_site(
                observed_site_muts=obs,
                allowed_site_muts=allowed_muts,
                af_map=af_map,
                af_min=major_af_min,
                af_max=None,
            )
            if matches:
                supported_sites.append(s)
                # keep highest-AF matching mut for reporting
                best = max(matches, key=lambda m: af_map[m])
                supported_mutations.append(best)

        n_def_sites = len(label_to_site_set[label])
        n_supported = len(supported_sites)
        frac = n_supported / n_def_sites if n_def_sites else 0.0

        rows.append({
            "label": label,
            "n_defining_sites": n_def_sites,
            "n_major_supported_sites": n_supported,
            "frac_major_supported_sites": frac,
            "major_supported_sites": sorted(supported_sites),
            "major_supported_mutations": sorted(supported_mutations),
        })

    out = pd.DataFrame(rows).sort_values(
        ["n_major_supported_sites", "frac_major_supported_sites", "n_defining_sites", "label"],
        ascending=[False, False, False, True]
    ).reset_index(drop=True)

    return out


# =========================
# Minor lineage scoring by site
# =========================

def score_minor_labels(
    sample_df: pd.DataFrame,
    label_to_sites: dict,
    label_to_site_set: dict,
    major_label: str,
    minor_af_min: float = 0.10,
    minor_af_max: float = 0.40,
    min_minor_support_sites: int = 5,
    min_minor_median_af: float = 0.12,
):
    sample_site_to_muts, af_map, altdp_map, totaldp_map = build_sample_support_maps(sample_df)

    major_sites = label_to_site_set[major_label]
    rows = []

    for label, site_map in label_to_sites.items():
        if label == major_label:
            continue

        candidate_specific_sites = set(site_map.keys()) - set(major_sites)

        supported_sites = []
        supported_mutations = []
        supported_af = []
        supported_altdp = []
        supported_totaldp = []

        for s in candidate_specific_sites:
            allowed_muts = site_map[s]
            obs = sample_site_to_muts.get(s, set())
            matches = matching_muts_at_site(
                observed_site_muts=obs,
                allowed_site_muts=allowed_muts,
                af_map=af_map,
                af_min=minor_af_min,
                af_max=minor_af_max,
            )
            if matches:
                best = max(matches, key=lambda m: af_map[m])
                supported_sites.append(s)
                supported_mutations.append(best)
                supported_af.append(af_map[best])
                supported_altdp.append(altdp_map[best])
                supported_totaldp.append(totaldp_map[best])

        n_specific = len(candidate_specific_sites)
        n_supported = len(supported_sites)
        median_af = float(pd.Series(supported_af).median()) if supported_af else 0.0
        mean_af = float(pd.Series(supported_af).mean()) if supported_af else 0.0

        passes = (n_supported >= min_minor_support_sites and median_af >= min_minor_median_af)

        rows.append({
            "label": label,
            "n_specific_sites": n_specific,
            "n_minor_supported_sites": n_supported,
            "minor_median_af": median_af,
            "minor_mean_af": mean_af,
            "minor_supported_sites": sorted(supported_sites),
            "minor_supported_mutations": sorted(supported_mutations),
            "minor_supported_alt_dp": supported_altdp,
            "minor_supported_total_dp": supported_totaldp,
            "passes_minor_rule": passes,
        })

    out = pd.DataFrame(rows).sort_values(
        ["passes_minor_rule", "n_minor_supported_sites", "minor_median_af", "minor_mean_af", "label"],
        ascending=[False, False, False, False, True]
    ).reset_index(drop=True)

    return out


# =========================
# Confidence tier
# =========================

def confidence_tier(n_sites: int, median_af: float) -> str:
    if n_sites >= 8 and median_af >= 0.15:
        return "high"
    if n_sites >= 6 and median_af >= 0.12:
        return "medium"
    if n_sites >= 5 and median_af >= 0.10:
        return "low"
    return "none"


# =========================
# Single-sample analysis
# =========================

def analyze_sample(
    sample_path: str,
    lineage_defs: pd.DataFrame,
    label_col: str = "lineage",
    min_total_dp: int = 100,
    min_alt_dp: int = 50,
    pass_only: bool = True,
    snv_only: bool = False,
    keep_only_lineage_sites: bool = True,
    major_af_min: float = 0.70,
    major_min_supported_sites: int = 8,
    major_min_frac: float = 0.40,
    minor_af_min: float = 0.10,
    minor_af_max: float = 0.40,
    min_minor_support_sites: int = 5,
    min_minor_median_af: float = 0.12,
):
    sample_name = Path(sample_path).stem

    sample_df = load_sample_variants(
        path=sample_path,
        lineage_defs=lineage_defs,
        min_total_dp=min_total_dp,
        min_alt_dp=min_alt_dp,
        pass_only=pass_only,
        snv_only=snv_only,
        keep_only_lineage_sites=keep_only_lineage_sites,
    )

    if sample_df.empty:
        return {
            "sample": sample_name,
            "n_filtered_variants": 0,
            "major_label": None,
            "major_n_supported_sites": 0,
            "major_frac_supported_sites": 0.0,
            "major_status": "no_data",
            "co_infection_flag": False,
            "minor_label": None,
            "minor_n_supported_sites": 0,
            "minor_median_af": 0.0,
            "confidence_tier": "none",
            "minor_supported_mutations": "",
        }, pd.DataFrame(), pd.DataFrame(), sample_df

    defs = lineage_defs.copy()
    defs["label"] = defs[label_col].astype(str)

    label_to_sites, label_to_site_set, _, _ = build_lineage_site_maps(defs, label_col="label")

    major_scores = score_major_labels(
        sample_df=sample_df,
        label_to_sites=label_to_sites,
        label_to_site_set=label_to_site_set,
        major_af_min=major_af_min,
    )

    best_major = major_scores.iloc[0].to_dict()
    major_label = best_major["label"]
    major_n_supported = int(best_major["n_major_supported_sites"])
    major_frac_supported = float(best_major["frac_major_supported_sites"])

    if major_n_supported < major_min_supported_sites or major_frac_supported < major_min_frac:
        summary = {
            "sample": sample_name,
            "n_filtered_variants": len(sample_df),
            "major_label": None,
            "major_n_supported_sites": major_n_supported,
            "major_frac_supported_sites": major_frac_supported,
            "major_status": "unassigned_low_support",
            "co_infection_flag": False,
            "minor_label": None,
            "minor_n_supported_sites": 0,
            "minor_median_af": 0.0,
            "confidence_tier": "none",
            "minor_supported_mutations": "",
        }
        return summary, major_scores, pd.DataFrame(), sample_df

    minor_scores = score_minor_labels(
        sample_df=sample_df,
        label_to_sites=label_to_sites,
        label_to_site_set=label_to_site_set,
        major_label=major_label,
        minor_af_min=minor_af_min,
        minor_af_max=minor_af_max,
        min_minor_support_sites=min_minor_support_sites,
        min_minor_median_af=min_minor_median_af,
    )

    if len(minor_scores):
        best_minor = minor_scores.iloc[0].to_dict()
        co_flag = bool(best_minor["passes_minor_rule"])
        minor_label = best_minor["label"]
        minor_n = int(best_minor["n_minor_supported_sites"])
        minor_med_af = float(best_minor["minor_median_af"])
        tier = confidence_tier(minor_n, minor_med_af)
        minor_supported_mutations = best_minor["minor_supported_mutations"]
    else:
        co_flag = False
        minor_label = None
        minor_n = 0
        minor_med_af = 0.0
        tier = "none"
        minor_supported_mutations = []

    summary = {
        "sample": sample_name,
        "n_filtered_variants": len(sample_df),
        "major_label": major_label,
        "major_n_supported_sites": major_n_supported,
        "major_frac_supported_sites": major_frac_supported,
        "major_status": "assigned",
        "co_infection_flag": co_flag,
        "minor_label": minor_label if co_flag else None,
        "minor_n_supported_sites": minor_n if co_flag else 0,
        "minor_median_af": minor_med_af if co_flag else 0.0,
        "confidence_tier": tier if co_flag else "none",
        "minor_supported_mutations": "|".join(minor_supported_mutations) if co_flag else "",
    }

    return summary, major_scores, minor_scores, sample_df


# =========================
# Main
# =========================

def main():
    ap = argparse.ArgumentParser(
        description="Site-based SARS-CoV-2 co-infection caller using lineage-defining mutations"
    )
    ap.add_argument("--variants-dir", required=True, help="Folder with per-sample TSV mutation tables")
    ap.add_argument("--lineage-defs", required=True, help="lineage_defs.tsv")
    ap.add_argument("--glob", default="*.tsv", help="Glob for sample files")
    ap.add_argument("--outdir", default="coinfection_calls_site_based", help="Output folder")

    ap.add_argument("--label-mode", choices=["lineage", "group"], default="lineage",
                    help="Score exact lineages or broad variant groups")

    ap.add_argument("--min-total-dp", type=int, default=100)
    ap.add_argument("--min-alt-dp", type=int, default=50)
    ap.add_argument("--no-pass-filter", action="store_true")
    ap.add_argument("--snv-only", action="store_true")
    ap.add_argument("--disable-lineage-site-filter", action="store_true")

    ap.add_argument("--major-af-min", type=float, default=0.70)
    ap.add_argument("--major-min-supported-sites", type=int, default=8)
    ap.add_argument("--major-min-frac", type=float, default=0.40)

    ap.add_argument("--minor-af-min", type=float, default=0.10)
    ap.add_argument("--minor-af-max", type=float, default=0.40)
    ap.add_argument("--min-minor-support-sites", type=int, default=5)
    ap.add_argument("--min-minor-median-af", type=float, default=0.12)

    ap.add_argument("--max-lineages-per-site", type=int, default=5,
                    help="Drop sites shared by more than this many exact lineages")

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lineage_defs = load_lineage_defs(args.lineage_defs)
    lineage_defs = filter_lineage_defs_by_site_sharing(
        lineage_defs=lineage_defs,
        max_lineages_per_site=args.max_lineages_per_site
    )
    lineage_defs.to_csv(os.path.join(args.outdir, "lineage_defs.filtered.tsv"), sep="\t", index=False)

    sample_paths = sorted(glob.glob(os.path.join(args.variants_dir, args.glob)))
    if not sample_paths:
        raise SystemExit(f"No files found in {args.variants_dir} matching {args.glob}")

    summaries = []

    for sample_path in sample_paths:
        sample_name = Path(sample_path).stem
        print(f"[INFO] analyzing {sample_name}")

        try:
            summary, major_scores, minor_scores, sample_df = analyze_sample(
                sample_path=sample_path,
                lineage_defs=lineage_defs,
                label_col=("group" if args.label_mode == "group" else "lineage"),
                min_total_dp=args.min_total_dp,
                min_alt_dp=args.min_alt_dp,
                pass_only=not args.no_pass_filter,
                snv_only=args.snv_only,
                keep_only_lineage_sites=not args.disable_lineage_site_filter,
                major_af_min=args.major_af_min,
                major_min_supported_sites=args.major_min_supported_sites,
                major_min_frac=args.major_min_frac,
                minor_af_min=args.minor_af_min,
                minor_af_max=args.minor_af_max,
                min_minor_support_sites=args.min_minor_support_sites,
                min_minor_median_af=args.min_minor_median_af,
            )
            summaries.append(summary)

            prefix = os.path.join(args.outdir, sample_name)
            major_scores.to_csv(f"{prefix}.major_scores.tsv", sep="\t", index=False)
            minor_scores.to_csv(f"{prefix}.minor_scores.tsv", sep="\t", index=False)
            sample_df.to_csv(f"{prefix}.filtered_variants.tsv", sep="\t", index=False)

        except Exception as e:
            summaries.append({
                "sample": sample_name,
                "n_filtered_variants": 0,
                "major_label": None,
                "major_n_supported_sites": 0,
                "major_frac_supported_sites": 0.0,
                "major_status": "error",
                "co_infection_flag": False,
                "minor_label": None,
                "minor_n_supported_sites": 0,
                "minor_median_af": 0.0,
                "confidence_tier": "error",
                "minor_supported_mutations": str(e),
            })

    cohort = pd.DataFrame(summaries).sort_values(
        ["co_infection_flag", "confidence_tier", "minor_n_supported_sites", "minor_median_af", "sample"],
        ascending=[False, True, False, False, True]
    )
    cohort.to_csv(os.path.join(args.outdir, "cohort_coinfection_summary.tsv"), sep="\t", index=False)
    print(f"[OK] wrote {os.path.join(args.outdir, 'cohort_coinfection_summary.tsv')}")


if __name__ == "__main__":
    main()
