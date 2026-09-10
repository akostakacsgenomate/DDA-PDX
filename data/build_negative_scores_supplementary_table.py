# -*- coding: utf-8 -*-
"""Build Supplementary Table: Negative_scores_

Recovers the manuscript claim that every validation molecular profile had
at least one PCT MTA with a negative DDA score (mean 7.3 negatives/profile).

This statistic is computed from the full DDA compound x profile score matrix
restricted to the PCT MTA panel, NOT from the 1,151 outcome-linked S3 rows.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent
S3_XLSX = DATA / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the negative-DDA-score supplementary table from the full "
            "compound–profile score matrix (not redistributed with this repository)."
        )
    )
    parser.add_argument(
        "--compounds-csv",
        type=Path,
        required=True,
        help="CSV of the full DDA compound × profile score matrix (LEVEL, COMPOUND, REPORT_ID).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for the generated Excel workbook (default: data/).",
    )
    return parser.parse_args()


def profile_table(df: pd.DataFrame, profile_ids: list[str] | None = None) -> pd.DataFrame:
    work = df.copy()
    if profile_ids is not None:
        work = work[work["REPORT_ID"].isin(profile_ids)]
    rows: list[dict] = []
    for pid, group in work.groupby("REPORT_ID", sort=True):
        n_neg = int((group["LEVEL"] < 0).sum())
        n_pos = int((group["LEVEL"] > 0).sum())
        n_zero = int((group["LEVEL"] == 0).sum())
        rows.append(
            {
                "PCT_ID": pid,
                "n_PCT_MTA_compounds_scored": int(len(group)),
                "n_negative_DDA_scores": n_neg,
                "n_positive_DDA_scores": n_pos,
                "n_zero_DDA_scores": n_zero,
                "has_at_least_one_negative": n_neg >= 1,
                "most_negative_DDA_score": float(group["LEVEL"].min()),
                "most_positive_DDA_score": float(group["LEVEL"].max()),
            }
        )
    return pd.DataFrame(rows).sort_values("PCT_ID").reset_index(drop=True)


def main() -> Path:
    args = parse_args()
    compounds_path = args.compounds_csv.expanduser().resolve()
    if not compounds_path.is_file():
        raise FileNotFoundError(
            f"Compound–profile score matrix not found: {compounds_path}. "
            "This file is not included in the public reproducibility package."
        )
    out_dir = args.output_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y_%m_%d")
    out_xlsx = out_dir / f"Negative_scores_Supplementary_Table_{stamp}.xlsx"

    comp = pd.read_csv(compounds_path)
    comp["LEVEL"] = pd.to_numeric(comp["LEVEL"], errors="coerce")
    comp["COMPOUND_U"] = comp["COMPOUND"].astype(str).str.strip().str.upper()
    comp["REPORT_ID"] = comp["REPORT_ID"].astype(str).str.strip()

    s3 = pd.read_excel(S3_XLSX, sheet_name="PCT_DRUG_RESPONSE")
    s3["PCT_ID"] = s3["PCT_ID"].astype(str).str.strip()
    s3["COMPOUND_U"] = s3["COMPOUND"].astype(str).str.strip().str.upper()
    s3["DDA_score"] = pd.to_numeric(s3["DDA score"], errors="coerce")

    is_chemo = s3["On_label"].astype(str).str.lower().eq("chemo")
    mta_compounds = sorted(s3.loc[~is_chemo, "COMPOUND_U"].dropna().unique())
    val_profiles = sorted(s3["PCT_ID"].unique())

    mta_scores = comp[comp["COMPOUND_U"].isin(mta_compounds)].copy()
    mta_scores = mta_scores[mta_scores["LEVEL"].notna()].copy()

    per_178 = profile_table(mta_scores, val_profiles)
    per_192 = profile_table(mta_scores, None)

    mta_val = mta_scores[mta_scores["REPORT_ID"].isin(val_profiles)].copy()
    comp_freq = (
        mta_val.assign(is_neg=mta_val["LEVEL"] < 0)
        .groupby("COMPOUND_U", as_index=False)
        .agg(
            n_profiles_scored=("REPORT_ID", "nunique"),
            n_profiles_negative=("is_neg", "sum"),
            mean_DDA_score=("LEVEL", "mean"),
            median_DDA_score=("LEVEL", "median"),
            min_DDA_score=("LEVEL", "min"),
            max_DDA_score=("LEVEL", "max"),
        )
        .sort_values(["n_profiles_negative", "COMPOUND_U"], ascending=[False, True])
        .reset_index(drop=True)
        .rename(columns={"COMPOUND_U": "COMPOUND"})
    )
    comp_freq["pct_profiles_negative"] = (
        100.0 * comp_freq["n_profiles_negative"] / comp_freq["n_profiles_scored"]
    ).round(1)

    neg_long = mta_val.loc[mta_val["LEVEL"] < 0, ["REPORT_ID", "COMPOUND", "LEVEL"]].copy()
    neg_long = neg_long.rename(columns={"REPORT_ID": "PCT_ID", "LEVEL": "DDA_score"})
    neg_long["COMPOUND"] = neg_long["COMPOUND"].astype(str)
    neg_long = neg_long.sort_values(["PCT_ID", "DDA_score", "COMPOUND"]).reset_index(drop=True)

    s3_scored = s3[s3["DDA_score"].notna()].copy()
    s3_per = (
        s3_scored.groupby("PCT_ID")
        .agg(
            n_outcome_linked_scored_treatments=("DDA_score", "size"),
            n_negative_in_S3=("DDA_score", lambda s: int((s < 0).sum())),
        )
        .reset_index()
    )
    s3_per = pd.DataFrame({"PCT_ID": val_profiles}).merge(s3_per, on="PCT_ID", how="left")
    s3_per["n_outcome_linked_scored_treatments"] = (
        s3_per["n_outcome_linked_scored_treatments"].fillna(0).astype(int)
    )
    s3_per["n_negative_in_S3"] = s3_per["n_negative_in_S3"].fillna(0).astype(int)
    s3_per["has_at_least_one_negative_in_S3"] = s3_per["n_negative_in_S3"] >= 1

    mean_178 = float(per_178["n_negative_DDA_scores"].mean())
    mean_192 = float(per_192["n_negative_DDA_scores"].mean())
    mean_s3 = float(s3_per["n_negative_in_S3"].mean())

    summary = pd.DataFrame(
        [
            {
                "cohort": "Validation molecular profiles (n=178)",
                "definition": (
                    "PCT MTA compounds (n=23 non-chemo agents in PCT panel) "
                    "scored by DDA for each profile"
                ),
                "n_profiles": int(len(per_178)),
                "n_profiles_with_ge1_negative": int(per_178["has_at_least_one_negative"].sum()),
                "pct_profiles_with_ge1_negative": round(
                    100.0 * per_178["has_at_least_one_negative"].mean(), 1
                ),
                "mean_n_negatives_per_profile": round(mean_178, 1),
                "mean_n_negatives_per_profile_exact": mean_178,
                "median_n_negatives_per_profile": float(
                    per_178["n_negative_DDA_scores"].median()
                ),
                "min_n_negatives_per_profile": int(per_178["n_negative_DDA_scores"].min()),
                "max_n_negatives_per_profile": int(per_178["n_negative_DDA_scores"].max()),
                "total_negative_assignments": int(per_178["n_negative_DDA_scores"].sum()),
                "mean_PCT_MTAs_scored_per_profile": round(
                    float(per_178["n_PCT_MTA_compounds_scored"].mean()), 2
                ),
            },
            {
                "cohort": "All DDA-scored molecular profiles in compound export (n=192)",
                "definition": (
                    "Same PCT MTA compound set; includes profiles beyond the "
                    "178 validation tumors"
                ),
                "n_profiles": int(len(per_192)),
                "n_profiles_with_ge1_negative": int(per_192["has_at_least_one_negative"].sum()),
                "pct_profiles_with_ge1_negative": round(
                    100.0 * per_192["has_at_least_one_negative"].mean(), 1
                ),
                "mean_n_negatives_per_profile": round(mean_192, 1),
                "mean_n_negatives_per_profile_exact": mean_192,
                "median_n_negatives_per_profile": float(
                    per_192["n_negative_DDA_scores"].median()
                ),
                "min_n_negatives_per_profile": int(per_192["n_negative_DDA_scores"].min()),
                "max_n_negatives_per_profile": int(per_192["n_negative_DDA_scores"].max()),
                "total_negative_assignments": int(per_192["n_negative_DDA_scores"].sum()),
                "mean_PCT_MTAs_scored_per_profile": round(
                    float(per_192["n_PCT_MTA_compounds_scored"].mean()), 2
                ),
            },
            {
                "cohort": "S3 outcome-linked validation treatments only (contrast)",
                "definition": (
                    "Only the 1,151 scored PCT treatments with experimental outcomes "
                    "in Supplementary PCT_DRUG_RESPONSE"
                ),
                "n_profiles": int(len(s3_per)),
                "n_profiles_with_ge1_negative": int(
                    s3_per["has_at_least_one_negative_in_S3"].sum()
                ),
                "pct_profiles_with_ge1_negative": round(
                    100.0 * s3_per["has_at_least_one_negative_in_S3"].mean(), 1
                ),
                "mean_n_negatives_per_profile": round(mean_s3, 1),
                "mean_n_negatives_per_profile_exact": mean_s3,
                "median_n_negatives_per_profile": float(s3_per["n_negative_in_S3"].median()),
                "min_n_negatives_per_profile": int(s3_per["n_negative_in_S3"].min()),
                "max_n_negatives_per_profile": int(s3_per["n_negative_in_S3"].max()),
                "total_negative_assignments": int(s3_per["n_negative_in_S3"].sum()),
                "mean_PCT_MTAs_scored_per_profile": round(
                    float(s3_per["n_outcome_linked_scored_treatments"].mean()), 2
                ),
            },
        ]
    )

    readme = pd.DataFrame(
        {
            "item": [
                "Table title",
                "Manuscript claim addressed",
                "Why not recoverable from S3 alone",
                "Definition of negative score",
                "PCT MTA compound set",
                "Primary cohort for the 7.3 claim",
                "Data source (full scores)",
                "Data source (validation IDs / compound list)",
            ],
            "content": [
                "Negative_scores_ - per-profile counts of PCT MTA drugs with negative DDA scores",
                (
                    "Every molecular profile had >=1 PCT drug with a negative DDA score; "
                    "mean = 7.3 negatives per profile"
                ),
                (
                    "S3 lists only outcome-linked scored treatments (n=1,151; ~6.5/tumor). "
                    "The claim counts DDA scores for the full PCT MTA panel on each "
                    "molecular profile, including agents not tested experimentally in that tumor."
                ),
                (
                    "DDA score (LEVEL) < 0, indicating net resistance / negative evidence "
                    "outweighs sensitivity evidence for that compound-profile pair"
                ),
                ", ".join(mta_compounds) + f" (n={len(mta_compounds)} non-chemo PCT compounds)",
                (
                    f"178 validation PCT molecular profiles; mean negatives/profile = "
                    f"{mean_178:.1f} (exact {mean_178:.6f})"
                ),
                str(compounds_path),
                str(S3_XLSX) + " - sheet PCT_DRUG_RESPONSE",
            ],
        }
    )

    methods = pd.DataFrame(
        {
            "section": [
                "Definition",
                "Denominator",
                "Statistic",
                "Relation to Extended Data Fig. 1a",
                "Relation to S3",
            ],
            "text": [
                (
                    "A PCT molecularly targeted agent (MTA) was counted as negative-scoring "
                    "for a molecular profile if its DDA score (LEVEL) was strictly less than zero."
                ),
                (
                    "For each of the 178 validation molecular profiles, DDA scores were taken "
                    "for all non-chemotherapy PCT compounds in the study drug panel (n=23 MTAs), "
                    "regardless of whether that compound was experimentally tested in that tumor. "
                    "Chemotherapy (On_label = chemo) was excluded."
                ),
                (
                    "We report (i) the proportion of profiles with >=1 negative-scoring PCT MTA "
                    "and (ii) the mean number of negative-scoring PCT MTAs per profile. In the "
                    "validation cohort these values are 100% and 7.3, respectively."
                ),
                (
                    "Extended Data Fig. 1a displays case-level DDA scores for PCT drugs; this "
                    "table quantifies the negative-score subset underlying the "
                    "resistance-mechanism statement in Results."
                ),
                (
                    "Supplementary PCT_DRUG_RESPONSE (S3) contains only treatments with "
                    "experimental outcomes (1,151 scored monotherapies). Negative-score "
                    "prevalence cannot be recovered from S3 alone because most PCT MTA x "
                    "profile scores lack matched PCT outcome rows."
                ),
            ],
        }
    )

    blurb = pd.DataFrame(
        {
            "use": ["Results (clarifying sentence)", "Methods (definition)"],
            "text": [
                (
                    "Considering DDA scores for the full PCT MTA panel on each molecular "
                    "profile (not restricted to outcome-linked treatments), every validation "
                    "profile had at least one negative-scoring PCT MTA (mean 7.3 negatives "
                    "per profile; Supplementary Table Negative_scores_)."
                ),
                (
                    "Negative-scoring PCT drugs were defined as non-chemotherapy PCT compounds "
                    "with DDA score < 0 for a given molecular profile. Counts were computed "
                    "from the complete DDA compound-profile score matrix for the 23 PCT MTAs "
                    "across the 178 validation profiles (Supplementary Table Negative_scores_)."
                ),
            ],
        }
    )

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)
        per_178.to_excel(writer, sheet_name="Per_profile_validation_178", index=False)
        neg_long.to_excel(writer, sheet_name="Negative_assignments_long", index=False)
        comp_freq.to_excel(writer, sheet_name="Compound_frequency", index=False)
        s3_per.to_excel(writer, sheet_name="Contrast_S3_outcome_linked", index=False)
        methods.to_excel(writer, sheet_name="Methods_note", index=False)
        blurb.to_excel(writer, sheet_name="Suggested_manuscript_text", index=False)
        per_192.to_excel(writer, sheet_name="Per_profile_all_192", index=False)

    print(f"Wrote {out_xlsx}")
    print(
        summary[
            [
                "cohort",
                "n_profiles",
                "n_profiles_with_ge1_negative",
                "mean_n_negatives_per_profile",
            ]
        ].to_string(index=False)
    )
    return out_xlsx


if __name__ == "__main__":
    main()
