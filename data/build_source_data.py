"""Build Nature Cancer statistics Source Data Excel workbooks.

One workbook per figure, saved under data/ (same folder as the input tables).
Each workbook has one sheet per panel letter (a, b, c, ...) with the
underlying source values for that panel (not derived summary statistics).

Usage
-----
python data/build_source_data.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent
REPO = DATA.parent
SRC = REPO / "src"
OUT = DATA

PRIMARY_XLSX = DATA / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
PRIMARY_SHEET = "PCT_DRUG_RESPONSE"
VARIANTS_SHEET = "VARIANTS"

BOX_DIR = SRC / "variant_counts"
MOA_DIR = SRC / "moa_distribution"

BOX_CLASSIFICATION_LABELS = {
    "DRIVER": "Driver",
    "VUS": "VUS",
    "NON_DRIVER": "Non-Driver",
    "All": "All",
}

ON_LABEL_COL_ALIASES = {"on_label", "onlabel", "on label", "on_label_utility"}


def _finalize_source_df(df: pd.DataFrame | None) -> pd.DataFrame:
    """Fill empty On_label as N/A, rename display columns, round numerics."""
    if df is None:
        return pd.DataFrame({"note": ["No rows"]})
    out = df.copy()
    if out.empty:
        return pd.DataFrame({"note": ["No rows"]})

    out = _standardize_source_colnames(out)

    for col in list(out.columns):
        key = str(col).strip().lower().replace("-", "_")
        if key in ON_LABEL_COL_ALIASES or key.endswith("on_label"):
            s = out[col]
            empty = s.isna() | s.astype(str).str.strip().isin(
                {"", "nan", "None", "<NA>", "NaN"}
            )
            out[col] = s.astype(object)
            out.loc[empty, col] = "N/A"

    skip_round = {
        "n",
        "n_models",
        "n_targets",
        "count",
        "rank_bin",
        "rank_pooled",
        "rank_int",
        "rank_within_model",
        "n_pairs",
        "level",
        "score",
        "n_Rank1_no_onlabel",
        "n_SC",
    }
    for col in out.columns:
        if not pd.api.types.is_numeric_dtype(out[col]):
            continue
        name = str(col).lower()
        if name in skip_round or name.startswith("n_") or pd.api.types.is_integer_dtype(
            out[col]
        ):
            continue
        out[col] = out[col].round(4)
    return out


def _standardize_source_colnames(df: pd.DataFrame) -> pd.DataFrame:
    """Rename stratum→Cohort and TTD/TimeToDouble→PFS for source tables."""
    rename: dict[str, str] = {}
    for col in df.columns:
        name = str(col)
        new = name
        if new.lower() == "stratum":
            new = "Cohort"
        elif new == "TimeToDouble":
            new = "PFS"
        else:
            # rank_TTD → rank_PFS, diff_TTD → diff_PFS, etc.
            new = new.replace("TTD", "PFS").replace("ttd", "PFS")
        if new != name and new not in rename.values() and new not in df.columns:
            rename[name] = new
        elif new != name and name not in rename:
            # Collision: keep original if target already exists
            if new not in df.columns:
                rename[name] = new
    return df.rename(columns=rename) if rename else df


def _center_align_worksheet(ws) -> None:
    from openpyxl.styles import Alignment

    center = Alignment(horizontal="center", vertical="center", wrap_text=False)
    for row in ws.iter_rows(min_row=1, max_row=ws.max_row, max_col=ws.max_column):
        for cell in row:
            cell.alignment = center


def _write_workbook(path: Path, sheets: dict[str, pd.DataFrame]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with pd.ExcelWriter(tmp, engine="openpyxl") as writer:
        for name, df in sheets.items():
            safe = name[:31]
            out = _finalize_source_df(df)
            out.to_excel(writer, sheet_name=safe, index=False)
            _center_align_worksheet(writer.sheets[safe])
    try:
        tmp.replace(path)
        return path
    except PermissionError:
        alt = path.with_name(path.stem + "_updated" + path.suffix)
        if alt.exists():
            alt.unlink()
        tmp.replace(alt)
        print(
            f"WARNING: locked file {path.name}; wrote {alt.name} instead "
            f"(close Excel and rename to replace)."
        )
        return alt


def _read_any(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    if path.suffix.lower() in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        frames = []
        for sheet in xl.sheet_names:
            part = pd.read_excel(path, sheet_name=sheet)
            if len(xl.sheet_names) > 1:
                part.insert(0, "source_sheet", sheet)
            frames.append(part)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return pd.read_csv(path)


def _latest_matching(directory: Path, pattern: str) -> Path | None:
    if not directory.exists():
        return None
    matches = sorted(
        directory.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True
    )
    return matches[0] if matches else None


def _latest_dir(directory: Path, pattern: str) -> Path | None:
    if not directory.exists():
        return None
    matches = sorted(
        [p for p in directory.glob(pattern) if p.is_dir()],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return matches[0] if matches else None


def _resolve_dm_dir() -> Path | None:
    dm = SRC / "data_modeling"
    return dm if dm.is_dir() else None


def _resolve_dcr_dir() -> Path | None:
    dcr = SRC / "dcr_orr_db"
    return dcr if dcr.is_dir() else None


def _resolve_km_ranking_dir() -> Path | None:
    km = SRC / "kaplan_meier" / "ranking"
    return km if km.is_dir() else None


def _resolve_hr_ranking_dir() -> Path | None:
    ranking = SRC / "hazard_ratios"
    return ranking if ranking.is_dir() else None


def _resolve_cu_dir() -> Path | None:
    cu = SRC / "clinical_utility_ranking"
    return cu if cu.is_dir() else None


def _resolve_loo_dir() -> Path | None:
    loo = SRC / "data_modeling" / "leave_one_out"
    return loo if loo.is_dir() else None


def _load_primary() -> pd.DataFrame:
    df = pd.read_excel(PRIMARY_XLSX, sheet_name=PRIMARY_SHEET)
    rename = {
        "PCT_ID": "Model",
        "DDA score": "LEVEL",
        "TimeToDouble (days)": "TimeToDouble",
        "TARGET": "Treatment target",
        "BestResponse (%)": "BestResponse",
        "BestAvgResponse (%)": "BestAvgResponse",
        "Response BestResponse": "Response_Best_Response",
        "Response BestAvgResponse": "Response_BestAvgResponse",
    }
    keep = {k: v for k, v in rename.items() if k in df.columns and v not in df.columns}
    return df.rename(columns=keep)


def _load_sc_chemo_rows(primary: pd.DataFrame) -> pd.DataFrame:
    """Standard-of-care (SC) rows: On_label == chemo from the primary workbook."""
    if "On_label" not in primary.columns:
        return pd.DataFrame()
    chemo = primary[
        primary["On_label"].astype(str).str.strip().str.lower() == "chemo"
    ].copy()
    if chemo.empty:
        return chemo
    cols = [
        c
        for c in ["Model", "COMPOUND", "LEVEL", "TimeToDouble", "On_label"]
        if c in chemo.columns
    ]
    out = chemo[cols].copy()
    out["Rank"] = "SC"
    out["GROUP"] = "SC"
    return out.reset_index(drop=True)


def _load_ranked_modeling_data() -> pd.DataFrame:
    dm = _resolve_dm_dir()
    if dm is None:
        return pd.DataFrame()
    path = dm / "data_all_ranks_normal_ranks.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    rename = {
        "REPORT_ID": "Model",
        "SCORE": "LEVEL",
        "TIME_TO_DOUBLE": "TimeToDouble",
        "ON_LABEL": "On_label",
        "RESPONSE_BESTAVGRESPONSE": "Response_BestAvgResponse",
        "RESPONSE_BEST_RESPONSE": "Response_Best_Response",
        "rank_pooled": "Rank",
    }
    keep = {k: v for k, v in rename.items() if k in df.columns}
    out = df.rename(columns=keep)
    if "Rank" in out.columns:
        out["Rank"] = out["Rank"].map(
            lambda r: "7-10" if int(r) >= 7 else str(int(r))
        )
    return out


def _load_cu_model_table() -> pd.DataFrame:
    cu = _resolve_cu_dir()
    if cu is None:
        return pd.DataFrame()
    path = _latest_matching(cu, "clinical_utility_model_table.xlsx")
    if path is None:
        path = _latest_matching(cu, "clinical_utility_model_table_*.xlsx")
    return _read_any(path)


def _fig2_site_counts(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work = work[
        ~work["Treatment target"].astype(str).str.contains(
            r"\bchemotherapy\b", case=False, na=False
        )
    ]
    work = work[
        ~work["Treatment target"].astype(str).str.contains(
            r"\bTubulin\b", case=False, na=False
        )
    ]
    work["Primary Tumor Site"] = work["Primary Tumor Site"].replace(
        {"UNKNOWN": "unknown"}
    )
    grouped = (
        work.groupby("Primary Tumor Site", dropna=False)["Model"]
        .nunique()
        .reset_index()
        .rename(columns={"Model": "n_models"})
    )
    total = grouped["n_models"].sum()
    grouped["percentage"] = (grouped["n_models"] / total) * 100.0
    return grouped.sort_values("Primary Tumor Site").reset_index(drop=True)


def _fig2_score_points(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        c
        for c in [
            "Model",
            "COMPOUND",
            "Treatment target",
            "LEVEL",
            "Primary Tumor Site",
            "BestResponse",
            "BestAvgResponse",
            "Response_Best_Response",
            "Response_BestAvgResponse",
            "TimeToDouble",
            "On_label",
            "Approved",
        ]
        if c in df.columns
    ]
    out = df[cols].copy()
    return out.sort_values(
        ["LEVEL", "Model", "COMPOUND"], ascending=[False, True, True]
    ).reset_index(drop=True)


def _resolve_variants_path() -> Path | None:
    if not PRIMARY_XLSX.exists():
        return None
    try:
        if VARIANTS_SHEET in pd.ExcelFile(PRIMARY_XLSX).sheet_names:
            return PRIMARY_XLSX
    except Exception:
        return None
    return None


def _fig2_boxplot_data() -> pd.DataFrame:
    """Per-model alteration counts by classification (boxplot source values)."""
    variants_path = _resolve_variants_path()
    if variants_path is None:
        return pd.DataFrame()

    variants = pd.read_excel(variants_path, sheet_name=VARIANTS_SHEET)
    data = variants[variants["CLASSIFICATION"] != "GENOMIC_MARKER"].copy()
    data["CLASSIFICATION"] = data["CLASSIFICATION"].replace({"VUS_IN_DRIVER": "VUS"})

    wes_data = (
        data.groupby(["PCT_ID", "CLASSIFICATION"]).size().reset_index(name="COUNT")
    )
    all_data = data.groupby(["PCT_ID"]).size().reset_index(name="COUNT")
    all_data["CLASSIFICATION"] = "All"
    box_plot_data = pd.concat([wes_data, all_data], ignore_index=True)
    box_plot_data["CLASSIFICATION"] = box_plot_data["CLASSIFICATION"].replace(
        BOX_CLASSIFICATION_LABELS
    )
    box_plot_data = box_plot_data.rename(columns={"PCT_ID": "Model"})
    return box_plot_data[["Model", "CLASSIFICATION", "COUNT"]].sort_values(
        ["Model", "CLASSIFICATION"]
    ).reset_index(drop=True)


def _fig2_moa_counts() -> pd.DataFrame:
    counts_path = _latest_matching(MOA_DIR, "Fig2c_MoA_piechart_counts_*.xlsx")
    return _read_any(counts_path)


def _assign_rank_groups(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_rank"] = (
        out.groupby("Model", group_keys=False)["LEVEL"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )
    out["Rank"] = out["_rank"].apply(lambda r: str(r) if r <= 6 else "7-10")
    return out.drop(columns=["_rank"])


def _fig3_panel_abc(ranked: pd.DataFrame, response_col: str | None) -> pd.DataFrame:
    """Patient/treatment-level rows underlying Fig3a–c stacked bars."""
    cols = ["Model", "COMPOUND", "LEVEL", "Rank", "TimeToDouble"]
    if response_col and response_col in ranked.columns:
        cols.append(response_col)
    keep = [c for c in cols if c in ranked.columns]
    out = ranked[keep].copy()
    if response_col is None and "TimeToDouble" in out.columns:
        out["DurableBenefit"] = out["TimeToDouble"].apply(
            lambda x: "DB"
            if pd.notna(x) and float(x) > 42
            else ("NonDB" if pd.notna(x) else pd.NA)
        )
    return out.sort_values(["Rank", "Model", "COMPOUND"]).reset_index(drop=True)


def _fig3_waterfall_source(primary: pd.DataFrame, response_cols: list[str]) -> pd.DataFrame:
    ranked = _assign_rank_groups(primary.dropna(subset=["LEVEL"]).copy())
    cols = ["Model", "COMPOUND", "LEVEL", "Rank", *response_cols]
    keep = [c for c in cols if c in ranked.columns]
    return (
        ranked[keep]
        .sort_values(["Rank", "LEVEL"], ascending=[True, False])
        .reset_index(drop=True)
    )


def _km_survival_source() -> pd.DataFrame:
    """Underlying survival times by rank group."""
    ranked = _load_ranked_modeling_data()
    if ranked.empty or "TimeToDouble" not in ranked.columns:
        return pd.DataFrame()
    cols = [
        c
        for c in ["Model", "COMPOUND", "LEVEL", "Rank", "TimeToDouble"]
        if c in ranked.columns
    ]
    return ranked[cols].sort_values(["Rank", "Model"]).reset_index(drop=True)


def _fig3_g_pairwise_hr(hr_dir: Path | None) -> pd.DataFrame:
    """Pairwise HR vs rank 1 only; numeric P and FDR, no significance labels."""
    if hr_dir is None:
        return pd.DataFrame()
    df = _read_any(hr_dir / "pairwise_hr.xlsx")
    if df.empty:
        return df
    comp_col = next(
        (c for c in df.columns if str(c).lower().replace(" ", "") == "groupcomparison"),
        None,
    )
    if comp_col is not None:
        df = df[df[comp_col].astype(str).str.match(r"^1\s+vs\.?\s", case=False, na=False)]
    drop_cols = [
        c
        for c in df.columns
        if str(c).lower() in {
            "fdr_significance_changed",
            "raw_significance",
            "fdr_values",
        }
    ]
    return df.drop(columns=drop_cols, errors="ignore").reset_index(drop=True)


def _fig3_h_scam_input(dm: Path | None) -> pd.DataFrame:
    """SCAM/GAM raw inputs (Fig3h); drop rank_pooled."""
    if dm is None:
        return pd.DataFrame()
    df = _read_any(dm / "data_all_ranks_normal_ranks.csv")
    if df.empty:
        return df
    return df.drop(columns=["rank_pooled"], errors="ignore")


def _fig3_i_effect_sizes(dm: Path | None) -> pd.DataFrame:
    """LMM effect sizes (Fig3i): calculated p_value and p_adj_fdr only."""
    if dm is None:
        return pd.DataFrame()
    df = _read_any(dm / "effect_sizes_all_ranks_vs_top_rank.csv")
    if df.empty:
        return df
    drop_cols = [
        c
        for c in df.columns
        if str(c).lower()
        in {"p_label", "p_label_plot", "p_adj_bonf", "p_adj_bonferroni"}
        or "bonf" in str(c).lower()
    ]
    return df.drop(columns=drop_cols, errors="ignore")


def _fig4_rank_plus_sc(
    ranked: pd.DataFrame, sc: pd.DataFrame, ranks: set[str]
) -> pd.DataFrame:
    """KM source for Fig4a/b: selected rank rows + SC chemo rows (no On_label)."""
    frames: list[pd.DataFrame] = []
    if not ranked.empty:
        part = ranked[ranked["Rank"].astype(str).isin(ranks)].copy()
        part["GROUP"] = part["Rank"].astype(str)
        cols = [
            c
            for c in ["Model", "COMPOUND", "LEVEL", "Rank", "TimeToDouble", "GROUP"]
            if c in part.columns
        ]
        frames.append(part[cols])
    if not sc.empty:
        cols = [
            c
            for c in ["Model", "COMPOUND", "LEVEL", "Rank", "TimeToDouble", "GROUP"]
            if c in sc.columns
        ]
        frames.append(sc[cols])
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["GROUP", "Model", "COMPOUND"], kind="stable").reset_index(
        drop=True
    )


def _import_clinical_utility():
    import sys

    src_str = str(SRC)
    if src_str not in sys.path:
        sys.path.insert(0, src_str)
    import clinical_utility_ranking as cu

    return cu


def _fig4_c_onlabel_status() -> pd.DataFrame:
    """Fig4c stacked-bar source: rank × On-label / Off-label / Experimental status."""
    try:
        cu = _import_clinical_utility()
    except Exception:
        return pd.DataFrame()
    ranked = cu._prepare_ranked_onlabel_patient_table(PRIMARY_XLSX, PRIMARY_SHEET)
    if ranked.empty or "ON_LABEL_UTILITY" not in ranked.columns:
        return pd.DataFrame()
    work = ranked[ranked["ON_LABEL_UTILITY"].isin(cu.ON_LABEL_UTILITY_ORDER)].copy()
    work["On_label"] = work["ON_LABEL_UTILITY"].map(cu.SHORT_SUMMARY_RANK_UTILITY_LABELS)
    rank_col = "RANK" if "RANK" in work.columns else "Rank"
    work["Rank"] = work[rank_col].astype(str)
    cols = [
        c
        for c in ["Model", "COMPOUND", "LEVEL", "Rank", "On_label"]
        if c in work.columns
    ]
    return work[cols].sort_values(
        ["Rank", "Model", "COMPOUND"], kind="stable"
    ).reset_index(drop=True)


def _fig4_d_no_onlabel_vs_sc(
    ranked: pd.DataFrame, model_table: pd.DataFrame, sc: pd.DataFrame
) -> pd.DataFrame:
    """Raw KM data for Fig4d: no-on-label rank-1 treatments + SC chemo rows."""
    frames: list[pd.DataFrame] = []

    if not model_table.empty and "scenario" in model_table.columns:
        s1 = model_table[
            (model_table["scenario"].astype(str) == "scenario1_no_onlabel")
            & model_table["rank1_TTD"].notna()
        ].copy()
        if not s1.empty:
            left = pd.DataFrame(
                {
                    "Model": s1["Model"].astype(str).values,
                    "COMPOUND": s1["rank1_COMPOUND"].values
                    if "rank1_COMPOUND" in s1.columns
                    else pd.NA,
                    "LEVEL": pd.NA,
                    "Rank": "1",
                    "TimeToDouble": pd.to_numeric(s1["rank1_TTD"], errors="coerce"),
                    "GROUP": "Rank1_no_onlabel",
                }
            )
            # Prefer full treatment row LEVEL from ranked table when available.
            if not ranked.empty:
                r1 = ranked[ranked["Rank"].astype(str) == "1"][
                    [c for c in ["Model", "COMPOUND", "LEVEL", "TimeToDouble"] if c in ranked.columns]
                ].copy()
                r1 = r1[r1["Model"].astype(str).isin(set(s1["Model"].astype(str)))]
                if not r1.empty:
                    left = r1.copy()
                    left["Rank"] = "1"
                    left["GROUP"] = "Rank1_no_onlabel"
            frames.append(left)

    if not sc.empty:
        # Prefer matched SC for scenario-1 models when chemo_TTD is available.
        if not model_table.empty and "chemo_TTD" in model_table.columns:
            s1_models = set(
                model_table.loc[
                    model_table["scenario"].astype(str) == "scenario1_no_onlabel",
                    "Model",
                ].astype(str)
            )
            matched = model_table[
                model_table["Model"].astype(str).isin(s1_models)
                & model_table["chemo_TTD"].notna()
            ].copy()
            if not matched.empty:
                right = pd.DataFrame(
                    {
                        "Model": matched["Model"].astype(str).values,
                        "COMPOUND": pd.NA,
                        "LEVEL": pd.NA,
                        "Rank": "SC",
                        "TimeToDouble": pd.to_numeric(
                            matched["chemo_TTD"], errors="coerce"
                        ),
                        "GROUP": "SC",
                    }
                )
                # Attach chemo compound names from SC table when possible.
                if "COMPOUND" in sc.columns:
                    chemo_map = (
                        sc.groupby("Model", as_index=True)["COMPOUND"]
                        .agg(lambda s: s.iloc[0] if len(s) else pd.NA)
                        .to_dict()
                    )
                    right["COMPOUND"] = right["Model"].map(chemo_map)
                frames.append(right)
            else:
                frames.append(
                    sc[
                        [
                            c
                            for c in [
                                "Model",
                                "COMPOUND",
                                "LEVEL",
                                "Rank",
                                "TimeToDouble",
                                "GROUP",
                            ]
                            if c in sc.columns
                        ]
                    ]
                )
        else:
            frames.append(
                sc[
                    [
                        c
                        for c in [
                            "Model",
                            "COMPOUND",
                            "LEVEL",
                            "Rank",
                            "TimeToDouble",
                            "GROUP",
                        ]
                        if c in sc.columns
                    ]
                ]
            )

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out = out.dropna(subset=["TimeToDouble"])
    return out.sort_values(["GROUP", "Model"], kind="stable").reset_index(drop=True)


def _fig4_e_paired() -> pd.DataFrame:
    """Fig4e: paired rank vs SC PFS rows; drop flags, ratios, and on_label_negative."""
    paired_long = _latest_matching(SRC / "paired_tests", "*_paired_data_long.csv")
    df = _read_any(paired_long)
    if df.empty:
        return df
    drop_cols = [
        c
        for c in df.columns
        if str(c).lower()
        in {"improved", "worse", "tie", "diff_pfs", "ratio_pfs", "pct_change_pfs"}
        or str(c).lower().replace("_", "") in {"diffpfs", "ratiopfs", "pctchangepfs"}
    ]
    df = df.drop(columns=drop_cols, errors="ignore")
    skip = {"on_label_negative", "on-label-negative"}
    for col in ("SUBGROUP", "Cohort"):
        if col in df.columns:
            df = df[~df[col].astype(str).str.strip().str.lower().isin(skip)]
    return df.reset_index(drop=True)


def _q1_combined_plot_labels(stratum: str, labels: dict[str, str]) -> tuple[str, str, str]:
    """Return (x-axis comparison, rank-1 arm, on-label arm) matching Fig4f ticks."""
    full = labels.get(stratum, stratum)
    if "\n" in full:
        left, right = full.split("\n", 1)
        rank1_arm = left.replace(" vs", "").strip()
        onlabel_arm = right.strip()
    elif " vs " in full:
        rank1_arm, onlabel_arm = full.split(" vs ", 1)
        rank1_arm = rank1_arm.strip()
        onlabel_arm = onlabel_arm.strip()
    else:
        rank1_arm, onlabel_arm = "Rank 1 (off/exp)", full
    return onlabel_arm, rank1_arm, onlabel_arm


def _fig4_f_paired(model_table: pd.DataFrame) -> pd.DataFrame:
    """Fig4f: Rank 1 listed once, then each compared on-label pool with its own label."""
    if model_table.empty:
        return pd.DataFrame()
    try:
        cu = _import_clinical_utility()
    except Exception:
        return pd.DataFrame()

    ranked = cu._prepare_ranked_onlabel_patient_table(PRIMARY_XLSX, PRIMARY_SHEET)
    rank1_arm = "Rank 1 (off/exp)"
    rank1_parts: list[pd.DataFrame] = []
    pool_parts: list[pd.DataFrame] = []
    pool_order: list[str] = []
    for stratum in cu.Q1_COMBINED_POOL_ORDER:
        paired = cu._q1_combined_stratum_paired(
            stratum, model_table, ranked, data_column=cu.DEFAULT_DATA_COLUMN
        )
        if paired is None or paired.empty:
            continue
        _, _, onlabel_arm = _q1_combined_plot_labels(
            stratum, cu.Q1_COMBINED_SHORT_LABELS
        )
        rank1_parts.append(
            pd.DataFrame(
                {
                    "Model": paired["Model"],
                    "PFS": pd.to_numeric(paired["rank_TTD"], errors="coerce"),
                    "arm": rank1_arm,
                }
            )
        )
        pool_parts.append(
            pd.DataFrame(
                {
                    "Model": paired["Model"],
                    "PFS": pd.to_numeric(paired["chemo_TTD"], errors="coerce"),
                    "arm": onlabel_arm,
                }
            )
        )
        pool_order.append(onlabel_arm)
    if not rank1_parts:
        return pd.DataFrame()
    rank1 = (
        pd.concat(rank1_parts, ignore_index=True)
        .dropna(subset=["PFS"])
        .drop_duplicates(subset=["Model"], keep="first")
        .sort_values("Model", kind="stable")
        .reset_index(drop=True)
    )
    pools = pd.concat(pool_parts, ignore_index=True).dropna(subset=["PFS"])
    arm_order = {name: i for i, name in enumerate(pool_order)}
    pools["_ord"] = pools["arm"].map(arm_order)
    pools = (
        pools.sort_values(["_ord", "Model"], kind="stable")
        .drop(columns=["_ord"])
        .reset_index(drop=True)
    )
    return pd.concat([rank1, pools], ignore_index=True)


def _load_non_assigned_rows() -> pd.DataFrame:
    """Non-assigned KM rows from data/non_assigned_merged.xlsx."""
    df = _read_any(DATA / "non_assigned_merged.xlsx")
    if df.empty:
        return df
    if "COMPOUND" not in df.columns and "Treatment" in df.columns:
        df = df.rename(columns={"Treatment": "COMPOUND"})
    if "TimeToDouble" not in df.columns:
        return pd.DataFrame()
    out = pd.DataFrame(
        {
            "Model": df["Model"].astype(str) if "Model" in df.columns else pd.NA,
            "COMPOUND": df["COMPOUND"] if "COMPOUND" in df.columns else pd.NA,
            "LEVEL": df["LEVEL"] if "LEVEL" in df.columns else pd.NA,
            "Rank": "Non_assigned",
            "TimeToDouble": pd.to_numeric(df["TimeToDouble"], errors="coerce"),
            "GROUP": "Non_assigned",
        }
    )
    return out.dropna(subset=["TimeToDouble"]).reset_index(drop=True)


def _fig4_g_nonassigned_vs_710(ranked: pd.DataFrame) -> pd.DataFrame:
    """Fig4g KM source: rank 7-10 plus Non-assigned rows."""
    frames: list[pd.DataFrame] = []
    if not ranked.empty:
        bottom = ranked[ranked["Rank"].astype(str) == "7-10"].copy()
        bottom["GROUP"] = "7-10"
        cols = [
            c
            for c in ["Model", "COMPOUND", "LEVEL", "Rank", "TimeToDouble", "GROUP"]
            if c in bottom.columns
        ]
        frames.append(bottom[cols])
    nona = _load_non_assigned_rows()
    if not nona.empty:
        frames.append(nona)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(
        ["GROUP", "Model", "COMPOUND"], kind="stable"
    ).reset_index(drop=True)


def build_fig2(primary: pd.DataFrame) -> Path:
    path = OUT / "Source_Data_Fig2.xlsx"
    points = _fig2_score_points(primary)
    wf_e = points[
        [
            c
            for c in [
                "Model",
                "COMPOUND",
                "LEVEL",
                "BestResponse",
                "Response_Best_Response",
                "Primary Tumor Site",
            ]
            if c in points.columns
        ]
    ].copy()
    wf_f = points[
        [
            c
            for c in [
                "Model",
                "COMPOUND",
                "LEVEL",
                "BestAvgResponse",
                "Response_BestAvgResponse",
                "Primary Tumor Site",
            ]
            if c in points.columns
        ]
    ].copy()
    sheets = {
        "a": _fig2_site_counts(primary),
        "b": _fig2_boxplot_data(),
        "c": _fig2_moa_counts(),
        "d": points,
        "e": wf_e,
        "f": wf_f,
    }
    return _write_workbook(path, sheets)


def build_fig3(primary: pd.DataFrame) -> Path:
    path = OUT / "Source_Data_Fig3.xlsx"
    ranked = _load_ranked_modeling_data()
    if ranked.empty:
        ranked = _assign_rank_groups(primary.dropna(subset=["LEVEL"]).copy())

    dm = _resolve_dm_dir()
    hr_dir = _resolve_hr_ranking_dir()

    panel_a = _fig3_panel_abc(ranked, "Response_Best_Response")
    panel_b = _fig3_panel_abc(ranked, "Response_BestAvgResponse")
    panel_c = _fig3_panel_abc(ranked, None)
    panel_c = panel_c[
        [
            c
            for c in [
                "Model",
                "COMPOUND",
                "LEVEL",
                "Rank",
                "TimeToDouble",
                "DurableBenefit",
            ]
            if c in panel_c.columns
        ]
    ]

    # Fig3h = SCAM/GAM raw inputs; Fig3i = LMM effect sizes (FDR only).
    sheets: dict[str, pd.DataFrame] = {
        "a": panel_a,
        "b": panel_b,
        "c": panel_c,
        "d": _fig3_waterfall_source(
            primary, ["BestResponse", "Response_Best_Response"]
        ),
        "e": _fig3_waterfall_source(
            primary, ["BestAvgResponse", "Response_BestAvgResponse"]
        ),
        "f": _km_survival_source(),
        "g": _fig3_g_pairwise_hr(hr_dir),
        "h": _fig3_h_scam_input(dm),
        "i": _fig3_i_effect_sizes(dm),
    }
    return _write_workbook(path, sheets)


def build_fig4(primary: pd.DataFrame) -> Path:
    path = OUT / "Source_Data_Fig4.xlsx"
    ranked = _load_ranked_modeling_data()
    if ranked.empty:
        ranked = _assign_rank_groups(primary.dropna(subset=["LEVEL"]).copy())
    sc = _load_sc_chemo_rows(primary)
    model_table = _load_cu_model_table()

    sheets: dict[str, pd.DataFrame] = {
        "a": _fig4_rank_plus_sc(ranked, sc, {"1"}),
        "b": _fig4_rank_plus_sc(ranked, sc, {"7-10"}),
        "c": _fig4_c_onlabel_status(),
        "d": _fig4_d_no_onlabel_vs_sc(ranked, model_table, sc),
        "e": _fig4_e_paired(),
        "f": _fig4_f_paired(model_table),
        "g": _fig4_g_nonassigned_vs_710(ranked),
    }

    return _write_workbook(path, sheets)


def build_ed_fig1(primary: pd.DataFrame) -> Path:
    path = OUT / "Source_Data_ED_Fig1.xlsx"
    points = _fig2_score_points(primary)
    ranked = points.dropna(subset=["LEVEL"]).copy()
    ranked["rank_within_model"] = (
        ranked.groupby("Model")["LEVEL"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    sheets = {
        "a": points,
        "b": ranked.sort_values(["Model", "rank_within_model"]).reset_index(drop=True),
    }
    return _write_workbook(path, sheets)


def build_ed_fig2() -> Path:
    path = OUT / "Source_Data_ED_Fig2.xlsx"
    loo = _resolve_loo_dir()
    dm = _resolve_dm_dir()
    sheets = {
        "a": _read_any(loo / "SCAM_GAM_Final_table.xlsx") if loo else pd.DataFrame(),
        "b": _read_any(dm / "data_all_ranks_normal_ranks.csv") if dm else pd.DataFrame(),
    }
    if sheets["a"].empty and loo is not None:
        all_dir = loo / "ALL"
        sheets["a"] = (
            _read_any(_latest_matching(all_dir, "*.csv"))
            if all_dir.exists()
            else pd.DataFrame()
        )
    return _write_workbook(path, sheets)


def write_inventory(written: list[Path]) -> Path:
    inv = OUT / "source_data_build_inventory.md"
    lines = [
        "# Nature Cancer statistics source data - inventory",
        "",
        f"Output folder: `{OUT.relative_to(REPO)}`",
        "",
        "Requirement: Excel format, one file per relevant figure.",
        "Each workbook uses panel-letter sheet names (`a`, `b`, `c`, ...).",
        "Sheets contain underlying source values for each panel (not derived summaries).",
        "Empty On_label values are written as N/A; numeric values are rounded to 4 decimals.",
        "",
        "## Packaged workbooks",
        "",
    ]
    for p in written:
        lines.append(f"- `{p.name}`")
    lines.extend(
        [
            "",
            "Regenerate workbooks with:",
            "",
            "```bash",
            "python data/build_source_data.py",
            "```",
            "",
        ]
    )
    inv.write_text("\n".join(lines), encoding="utf-8")
    return inv


def main() -> None:
    primary = _load_primary()
    written = [
        build_fig2(primary),
        build_fig3(primary),
        build_fig4(primary),
        build_ed_fig1(primary),
        build_ed_fig2(),
    ]
    inv = write_inventory(written)
    for p in written + [inv]:
        print(f"Wrote {p} ({p.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
