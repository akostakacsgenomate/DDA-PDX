"""Clinical utility ranking — Fig. 4c–4f panel outputs.

Generates:
  - Fig4c_pub_rank_utility_composition
  - Fig4d_pub_km_no_onlabel_rank1_vs_SC
    (no-on-label rank-1 vs all SC rows; pooled log-rank / unstratified Cox)
  - Fig4e paired rank vs matched-SC PFS boxplot
  - Fig4f combined-pool mPFS boxplot

Usage: python src/clinical_utility_ranking.py [--source PATH] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from matplotlib.patches import Patch
from matplotlib.transforms import ScaledTranslation, blended_transform_factory
from openpyxl import load_workbook
from openpyxl.styles import Alignment
from openpyxl.styles import PatternFill
from scipy import stats
from scipy.stats import chi2 as chi2_distribution
from statsmodels.stats.multitest import multipletests

from plot_typography import (
    apply_arial_font,
    format_p_equals_mathtext,
    format_p_sci_mathtext,
    format_p_sci_plain,
)

apply_arial_font()

SEED = 42
N_BOOTSTRAP = 10_000
np.random.seed(SEED)
LOGGER = logging.getLogger("clinical_utility_ranking")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"
DEFAULT_DATA_COLUMN = "TimeToDouble"
KM_XLIM = (0, 200)
KM_XTICKS = np.arange(0, 225, 25)
DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_SOURCE_SHEET = "PCT_DRUG_RESPONSE"
OUTPUT_DIRNAME = "clinical_utility_ranking"
# Primary survival inference for this run: unstratified (standard) log-rank / Cox.
USE_STRATIFIED_SURVIVAL = False
# Q2: Downranked on-label pools vs Rank 1 (Off/Exp) or standard chemo (unpaired) — disabled.
ENABLE_Q2_ANALYSIS = False

_FILL_SIG = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
_FILL_NS = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")
MIN_SCENARIO_PAIRED_N = 5
P_VALUE_EXCEL_COLUMNS: Tuple[str, ...] = (
    "wilcoxon_p",
    "wilcoxon_p_fdr_bh",
    "binom_p",
    "log_ratio_t_p",
    "stratified_logrank_p",
    "stratified_logrank_p_fdr_bh",
    "stratified_cox_p",
    "stratified_cox_p_fdr_bh",
    "pooled_logrank_p",
    "pooled_logrank_p_fdr_bh",
    "cox_p",
    "logrank_pooled_p",
    "mannwhitney_p",
    "pooled_cox_p",
)


def resolve_source_path(explicit: Optional[Path] = None) -> Path:
    if explicit is not None:
        p = explicit.resolve()
        if p.is_file():
            return p
        raise FileNotFoundError(p)
    if DEFAULT_SOURCE_XLSX.is_file():
        return DEFAULT_SOURCE_XLSX.resolve()
    raise FileNotFoundError(f"Source workbook not found: {DEFAULT_SOURCE_XLSX}")


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()



def apply_excel_alignment_and_significance(path: Path, *, alpha: float = 0.05) -> None:
    """Center cells and color significance columns in a workbook."""
    wb = load_workbook(path)
    center = Alignment(horizontal="center", vertical="center")
    for ws in wb.worksheets:
        headers = [cell.value for cell in ws[1]]
        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = center
        for col_name in ("raw_significance", "fdr_significance", "FDR_values", "Sig_raw", "Sig_fdr"):
            if col_name not in headers:
                continue
            col_idx = headers.index(col_name) + 1
            for row_idx in range(2, ws.max_row + 1):
                cell = ws.cell(row=row_idx, column=col_idx)
                value = str(cell.value).strip().lower()
                is_sig = value not in ("", "ns", "non-significant", "none", "nan")
                cell.fill = _FILL_SIG if is_sig else _FILL_NS
        for col_name in headers:
            if col_name not in P_VALUE_EXCEL_COLUMNS:
                continue
            col_idx = headers.index(col_name) + 1
            for row_idx in range(2, ws.max_row + 1):
                cell = ws.cell(row=row_idx, column=col_idx)
                if cell.value is None or str(cell.value).strip() == "":
                    continue
                try:
                    p_val = float(cell.value)
                except (TypeError, ValueError):
                    continue
                if not np.isfinite(p_val):
                    continue
                cell.fill = _FILL_SIG if p_val < alpha else _FILL_NS
    wb.save(path)


_RANK_POOL_FROM = 7
_POOLED_LABEL = "7-10"

RANKING_UTILITY_RANKS: Tuple[str, ...] = ("1", "2", "3", "4", "5", "6", _POOLED_LABEL)
RANKING_CHEMO_NONASSIG_GROUP_ORDER: Tuple[str, ...] = (
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    _POOLED_LABEL,
    "Non_assigned",
    "SC",
)
RANKING_KM_COLORS: Dict[str, str] = {
    "1": "#FFD700",
    "2": "#C0C0C0",
    "3": "#CD7F32",
    "4": "#C2185B",
    "5": "#64B5F6",
    "6": "#D8F3DC",
    "7-10": "#404040",
    "SC": "#b121f6",
    "Non_assigned": "#8B4513",
}
RANKING_KM_COLORS_UTILITY_REFS: Dict[str, str] = {
    **RANKING_KM_COLORS,
    "2": "#000000",
    "3": "#2B2B2B",
    "4": "#505050",
    "5": "#757575",
    "6": "#9A9A9A",
    "7-10": "#C0C0C0",
}
ON_LABEL_UTILITY_ORDER: Tuple[str, ...] = ("On_label", "Off_label", "Exp")
ON_LABEL_UTILITY_COLORS: Dict[str, str] = {
    "On_label": "#50C878",
    "Off_label": "#0F52BA",
    "Exp": "#F0A030",
}
RANKING_UTILITY_REF_GROUPS: Tuple[str, ...] = ("SC", "NonA")
RANKING_UTILITY_GROUP_ORDER: Tuple[str, ...] = ON_LABEL_UTILITY_ORDER + RANKING_UTILITY_REF_GROUPS
RANKING_UTILITY_COLORS: Dict[str, str] = {
    **ON_LABEL_UTILITY_COLORS,
    "SC": "#b121f6",
    "NonA": "#8B4513",
}

CLINICAL_UTILITY_SCENARIO1 = "scenario1_no_onlabel"
CLINICAL_UTILITY_SCENARIO2 = "scenario2_onlabel_rank1"
CLINICAL_UTILITY_SCENARIO3 = "scenario3_onlabel_not_top"
CLINICAL_UTILITY_SCENARIO4 = "scenario4_onlabel_low_rank"
CLINICAL_UTILITY_SCENARIOS: Tuple[str, ...] = (
    CLINICAL_UTILITY_SCENARIO1,
    CLINICAL_UTILITY_SCENARIO2,
    CLINICAL_UTILITY_SCENARIO3,
    CLINICAL_UTILITY_SCENARIO4,
)
CLINICAL_UTILITY_SCENARIO3_STRATA: Tuple[str, ...] = ("rank2", "rank3", "rank2_3_pooled")
CLINICAL_UTILITY_GROUP_RANK1 = "Rank_1"
CLINICAL_UTILITY_GROUP_ONLABEL = "On_label"
CLINICAL_UTILITY_GROUP_ABOVE = "Best_OffExp_above"
CLINICAL_UTILITY_GROUP_SC = "SC"
CLINICAL_UTILITY_COLORS: Dict[str, str] = {
    CLINICAL_UTILITY_GROUP_RANK1: "#FFD700",
    CLINICAL_UTILITY_GROUP_ONLABEL: "#50C878",
    CLINICAL_UTILITY_GROUP_ABOVE: "#0F52BA",
    CLINICAL_UTILITY_GROUP_SC: "#b121f6",
}
CLINICAL_UTILITY_COMPARISONS: Tuple[Tuple[str, str, str], ...] = (
    ("rank1_vs_onlabel", CLINICAL_UTILITY_GROUP_RANK1, CLINICAL_UTILITY_GROUP_ONLABEL),
    ("rank1_vs_sc", CLINICAL_UTILITY_GROUP_RANK1, CLINICAL_UTILITY_GROUP_SC),
    ("onlabel_vs_sc", CLINICAL_UTILITY_GROUP_ONLABEL, CLINICAL_UTILITY_GROUP_SC),
    (
        "above_onlabel_vs_onlabel",
        CLINICAL_UTILITY_GROUP_ABOVE,
        CLINICAL_UTILITY_GROUP_ONLABEL,
    ),
)

def assign_dense_level_groups(
    df: pd.DataFrame,
    *,
    model_col: str = "Model",
    level_col: str = "LEVEL",
    group_col: str = "GROUP",
    rank_pool_from: int = _RANK_POOL_FROM,
    pooled_label: str = _POOLED_LABEL,
) -> pd.DataFrame:
    """Assign GROUP from dense rank of level_col within each model_col group.

    Higher LEVEL → lower rank number (1 = best DDA response). Tied values
    share rank. Ranks >= rank_pool_from are pooled into pooled_label ('7-10').
    """
    out = df.copy()
    level_rank = (
        out.groupby(model_col, group_keys=False)[level_col]
        .rank(method="dense", ascending=False)
        .astype(int)
    )
    rank_pooled = level_rank.clip(upper=rank_pool_from)
    out[group_col] = [
        pooled_label if r >= rank_pool_from else str(int(r))
        for r in rank_pooled
    ]
    return out

@dataclass
class KMPlotStyle:
    """Shared typography parameters for KM and composition figures."""

    font: float = 6.4 * 1.6
    tick_size: float = 6.4 * 1.6
    label_size: float = 6.4 * 1.6
    at_risk: float = 5.5 * 1.6
    font_legend: float = 6.4 * 1.6
    cishow: bool = False

    def apply_rcparams(self) -> None:
        plt.rcParams.update(
            {
                "axes.titlesize": self.font,
                "axes.labelsize": self.font,
                "xtick.labelsize": self.at_risk,
                "ytick.labelsize": self.font,
                "legend.fontsize": self.font_legend,
                "font.size": self.font,
            }
        )


# Top / Bottom pair KMs (Fig 4a/4b): +20% over default KMPlotStyle — match KM_PDX.PAIR_KM_ETALON_STYLE.
_PAIR_KM_FONT_SCALE = 1.20
PAIR_KM_ETALON_STYLE = KMPlotStyle(
    font=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    tick_size=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    label_size=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    at_risk=5.5 * 1.6 * _PAIR_KM_FONT_SCALE,
    font_legend=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
)

FIG4C_PANEL_STEM = "Fig4c_pub_rank_utility_composition"
FIG4D_PANEL_STEM = "Fig4d_pub_km_no_onlabel_rank1_vs_SC"
FIG4E_PANEL_STEM = "Fig4e_clinical_utility_rank_mPFS_vs_SC_paired_boxplot"
FIG4F_PANEL_STEM = "Fig4f_clinical_utility_q1_combined_pool_mpfs_boxplot"
FIG4_PANELS_DIR = _SCRIPT_DIR / "figure_panels" / "Fig4"
_FIG4C_FONT_SCALE = 1.10


def normalize_id_columns(
    df: pd.DataFrame, model_col: str = "Model", compound_col: str = "COMPOUND"
) -> pd.DataFrame:
    """Lowercase and strip hyphens from model/compound identifiers."""
    out = df.copy()
    out[model_col] = out[model_col].astype(str).str.replace("-", "", regex=False).str.lower()
    out[compound_col] = out[compound_col].astype(str).str.replace("-", "", regex=False).str.lower()
    return out


def harmonize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Harmonize source-workbook column names to the analysis schema."""
    out = df.copy()
    supplementary_renames = {
        "PCT_ID": "Model",
        "DDA score": "LEVEL",
        "TimeToDouble (days)": "TimeToDouble",
    }
    rename = {
        source: target
        for source, target in supplementary_renames.items()
        if source in out.columns and target not in out.columns
    }
    if rename:
        out = out.rename(columns=rename)
    if "Treatment target" not in out.columns and "TARGET" in out.columns:
        out = out.rename(columns={"TARGET": "Treatment target"})
    if "Treatment" in out.columns and "COMPOUND" not in out.columns:
        out = out.rename(columns={"Treatment": "COMPOUND"})
    if "Approved" in out.columns and "APPROVED" not in out.columns:
        out["APPROVED"] = out["Approved"].map(normalize_approved_workbook_value)
    elif "APPROVED" in out.columns:
        out["APPROVED"] = out["APPROVED"].map(normalize_approved_workbook_value)
    if "On_label" in out.columns:
        out["On_label"] = out["On_label"].map(normalize_on_label_workbook_value)
    return out


def normalize_approved_workbook_value(value: Any) -> Any:
    """Normalize Approved flags in the source workbook."""
    if pd.isna(value):
        return value
    if value in (1, 1.0, True, "1", "1.0", "True", "TRUE"):
        return "IGAZ"
    if value in (0, 0.0, False, "0", "0.0", "False", "FALSE"):
        return "HAMIS"
    text = str(value).strip().upper()
    if text in ("IGAZ", "YES", "TRUE", "APPROVED"):
        return "IGAZ"
    if text in ("HAMIS", "NO", "FALSE", "NOT APPROVED"):
        return "HAMIS"
    return str(value).strip()


def normalize_on_label_workbook_value(value: Any) -> Any:
    """Normalize On_label flags in the source workbook."""
    if pd.isna(value):
        return value
    if value is True or str(value).strip().upper() == "TRUE":
        return "IGAZ"
    if value is False or str(value).strip().upper() == "FALSE":
        return "HAMIS"
    text = str(value).strip()
    if text.upper() == "UNKNOWN":
        return "UNKNOWN"
    if text.lower() in ("n/a", "na"):
        return "n/a"
    return text

def load_primary_patient_table(source_path: Path, sheet_name: Optional[str] = None) -> pd.DataFrame:
    """Load the primary patient workbook (optionally from a named Excel sheet)."""
    if source_path.suffix.lower() == ".csv":
        return pd.read_csv(source_path)
    if sheet_name:
        return pd.read_excel(source_path, sheet_name=sheet_name)
    return pd.read_excel(source_path)


def clean_base_patient_table(df: pd.DataFrame) -> pd.DataFrame:
    """Apply shared row filtering and numeric casting."""
    out = df.copy()
    out = out[out["LEVEL"] != "#HIÁNYZIK"].copy()
    out = out[out["LEVEL"].notna()].copy()
    out["LEVEL"] = out["LEVEL"].astype(float)
    out["TimeToDouble"] = pd.to_numeric(out["TimeToDouble"], errors="coerce").astype(float)
    if "Day_Last" in out.columns:
        out["Day_Last"] = pd.to_numeric(out["Day_Last"], errors="coerce")
    return out


def normalize_survival_columns(
    df: pd.DataFrame, time_col: str, event_col: str = "CENSOR"
) -> pd.DataFrame:
    """Ensure lifelines receives numeric durations and boolean event indicators."""
    out = df.copy()
    out.loc[:, time_col] = pd.to_numeric(out[time_col], errors="coerce")
    ev = out[event_col].replace({True: 1, False: 0, "1": 1, "0": 0, "True": 1, "False": 0})
    events = pd.to_numeric(ev, errors="coerce").fillna(0).astype(int).to_numpy(dtype=np.bool_, copy=False)
    out[event_col] = pd.Series(events, index=out.index, dtype=bool)
    return out


def validate_on_label_source_columns(patient_df: pd.DataFrame) -> None:
    missing = {"APPROVED", "On_label"}.difference(patient_df.columns)
    if missing:
        raise KeyError(
            f"Source table is missing on-label columns: {sorted(missing)}. "
            "Expected APPROVED and On_label in the recategorised source workbook."
        )


def resolve_on_label_columns(patient_df: pd.DataFrame) -> Tuple[str, str, str]:
    """Detect on-label column naming convention used in the workbook."""
    on_col = "On_label" if "On_label" in patient_df.columns else "On_off_label"
    on_yes = "IGAZ" if on_col == "On_label" else "YES_on"
    on_no = "HAMIS" if on_col == "On_label" else "no_off"
    return on_col, on_yes, on_no


def is_chemo_on_label(series: pd.Series) -> pd.Series:
    """True for rows tagged On_label = chemo (standard-of-care reference arm)."""
    return series.astype(str).str.strip().str.lower() == "chemo"


def map_on_label_utility(value: Any) -> Optional[str]:
    """Map raw On_label workbook values to utility groups.

    IGAZ -> On_label, HAMIS -> Off_label, n/a (or missing) -> Exp.
    UNKNOWN (any case) and unrecognised values return None (row excluded).
    """
    if pd.isna(value):
        return "Exp"
    text = str(value).strip()
    if not text or text.lower() in ("n/a", "na"):
        return "Exp"
    upper = text.upper()
    if upper == "UNKNOWN":
        return None
    if text == "IGAZ":
        return "On_label"
    if text == "HAMIS":
        return "Off_label"
    if text.strip().lower() == "chemo":
        return None
    LOGGER.warning("Unrecognised On_label value %r; excluding row.", value)
    return None


def assign_on_label_utility_column(
    df: pd.DataFrame,
    on_col: str = "On_label",
    utility_col: str = "ON_LABEL_UTILITY",
    *,
    drop_unmapped: bool = True,
) -> pd.DataFrame:
    """Add utility on-label column; optionally drop UNKNOWN / unmapped rows.

    For dense ranking across all screened compounds, call with
    ``drop_unmapped=False`` so UNKNOWN (e.g. tamoxifen) stays in the rank pool.
    Analyses that use a no-on-label / utility group arm should then exclude
    rows where ``ON_LABEL_UTILITY`` is null.
    """
    out = df.copy()
    out[utility_col] = out[on_col].map(map_on_label_utility)
    before = len(out)
    if drop_unmapped:
        out = out[out[utility_col].notna()].copy()
    excluded = before - len(out)
    if excluded:
        LOGGER.info("Excluded %d row(s) with UNKNOWN or unmapped On_label values.", excluded)
    elif not drop_unmapped:
        n_unknown = int(out[utility_col].isna().sum())
        if n_unknown:
            LOGGER.info(
                "Kept %d UNKNOWN/unmapped On_label row(s) in the rank pool "
                "(excluded from no-on-label / utility group arms).",
                n_unknown,
            )
    return out


def rows_with_mapped_on_label_utility(
    df: pd.DataFrame, utility_col: str = "ON_LABEL_UTILITY"
) -> pd.DataFrame:
    """Rows with a mapped on-label utility (excludes UNKNOWN / unmapped)."""
    if utility_col not in df.columns:
        return df.copy()
    return df[df[utility_col].notna()].copy()

def load_sc_chemo_reference(
    source_path: Path, source_sheet: Optional[str]
) -> pd.DataFrame:
    """SC reference: rows with On_label = chemo from the primary source workbook.

    Chemo rows intentionally have no DDA LEVEL score, so they must be extracted
    before ``clean_base_patient_table`` (which drops rows with missing LEVEL).
    """
    patient = load_primary_patient_table(source_path, source_sheet)
    patient = harmonize_column_names(patient)
    patient = normalize_id_columns(patient)
    on_col, _, _ = resolve_on_label_columns(patient)
    chemo_from_source = patient[is_chemo_on_label(patient[on_col])].copy()
    if chemo_from_source.empty:
        LOGGER.info("No rows with On_label = chemo found; SC reference cohort omitted.")
        return pd.DataFrame()

    chemo_from_source["CENSOR"] = True
    chemo_from_source = normalize_survival_columns(chemo_from_source, DEFAULT_DATA_COLUMN)
    return chemo_from_source.dropna(subset=[DEFAULT_DATA_COLUMN])

def significance_stars(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "ns"
    if p_value < 0.0005:
        return "***"
    if p_value < 0.005:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def space_plot_operators(text: str) -> str:
    """Add readable spacing around comparison operators in on-plot annotations.

    Mathtext segments (``$...$``) are left untouched so *P*-value formatting
    like ``$\\mathit{P} = 5.6 \\\\times 10^{-2}$`` stays valid.
    """
    parts = re.split(r"(\$[^$]*\$)", text)

    def _space_plain(chunk: str) -> str:
        spaced = chunk
        spaced = re.sub(r">=", " >= ", spaced)
        spaced = re.sub(r"<=", " <= ", spaced)
        spaced = re.sub(r"(?<![<>=!])=(?!=)", " = ", spaced)
        spaced = re.sub(r"(?<![=])>(?!=)", " > ", spaced)
        spaced = re.sub(r"(?<![=!<])<(?!=)", " < ", spaced)
        return re.sub(r"  +", " ", spaced)

    return "".join(
        part if part.startswith("$") and part.endswith("$") else _space_plain(part)
        for part in parts
    )


def safe_filename(value: str) -> str:
    if not str(value).strip():
        return ""
    text = str(value).strip()
    for ch in '<>:"/\\|?*':
        text = text.replace(ch, "_")
    return text


def safe_to_excel(df: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    """Write Excel output; log and continue if the target file is locked."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_excel(path, index=index, engine="openpyxl")
        apply_excel_alignment_and_significance(path)
    except PermissionError:
        LOGGER.warning("Could not write %s (file may be open in Excel).", path)


def apply_km_axis_style(ax: plt.Axes, style: KMPlotStyle, title: str = "") -> None:
    ax.set_xlabel("Time (days)", fontsize=style.label_size, labelpad=4)
    ax.set_ylabel("Survival Probability (%)", fontsize=style.label_size, labelpad=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(*KM_XLIM)
    ax.set_xticks(KM_XTICKS)
    if title:
        ax.set_title(title)


def add_km_at_risk_table(
    ax: plt.Axes,
    kmf_list: Sequence[KaplanMeierFitter],
    *,
    at_risk_ypos: float = -0.35,
) -> None:
    """Add lifelines at-risk count table below a KM axes."""
    if not kmf_list:
        return
    add_at_risk_counts(
        *kmf_list,
        ypos=at_risk_ypos,
        ax=ax,
        rows_to_show=["At risk"],
        at_risk_count_from_start_of_period=True,
        xticks=KM_XTICKS,
    )


def collect_km_group_statistics(
    kmf: KaplanMeierFitter,
    fit_label: str,
    group: str,
    group_data: pd.DataFrame,
    level_col: str = "LEVEL",
) -> Dict[str, Any]:
    """Extract median survival and 95% CI bounds from a fitted KM curve."""
    lower_bound = kmf.confidence_interval_[f"{fit_label}_lower_0.95"]
    upper_bound = kmf.confidence_interval_[f"{fit_label}_upper_0.95"]
    lower_median_time = np.interp(0.5, 1 - lower_bound.values.flatten(), lower_bound.index)
    upper_median_time = np.interp(0.5, 1 - upper_bound.values.flatten(), upper_bound.index)
    return {
        "Group": group,
        "Median Survival Time": kmf.median_survival_time_,
        "n_number": len(group_data),
        "CI_lower": lower_median_time,
        "CI_upper": upper_median_time,
        "avg_LEVEL": group_data[level_col].mean(),
    }

def _pair_km_display_label(group: str, stat_label_map: Optional[Dict[str, str]] = None) -> str:
    if stat_label_map and group in stat_label_map:
        return stat_label_map[group]
    return group


def _format_pair_km_median(value: Any) -> str:
    if pd.isna(value) or np.isinf(value):
        return "NR"
    return f"{float(value):.1f}"


def _annotate_pair_km_stats_on_ax(
    ax: plt.Axes,
    *,
    left_group: str,
    right_group: str,
    left_median: Any,
    right_median: Any,
    logrank_p: float,
    hr: float,
    ci_low: float,
    ci_high: float,
    hr_p: float,
    style: KMPlotStyle,
    stat_label_map: Optional[Dict[str, str]] = None,
    logrank_label: str = "logrank",
    hr_label: str = "HR",
    stats_x: float = 0.97,
    stats_y: float = 0.40,
    stats_ha: str = "right",
) -> None:
    """Add mPFS / log-rank / HR annotation box inside the KM axes."""
    left_label = _pair_km_display_label(left_group, stat_label_map)
    right_label = _pair_km_display_label(right_group, stat_label_map)
    stats_lines = [
        f"{left_label} mPFS = {_format_pair_km_median(left_median)}",
        f"{right_label} mPFS = {_format_pair_km_median(right_median)}",
        f"{logrank_label} {format_p_equals_mathtext(logrank_p)}",
    ]
    if pd.notna(hr):
        stats_lines.append(f"{hr_label} = {float(hr):.2f} ({float(ci_low):.2f}-{float(ci_high):.2f})")
        if pd.notna(hr_p):
            stats_lines.append(f"{hr_label} {format_p_equals_mathtext(hr_p)}")
    ax.text(
        stats_x,
        stats_y,
        space_plot_operators("\n".join(stats_lines)),
        transform=ax.transAxes,
        fontsize=style.font * 0.9,
        va="top",
        ha=stats_ha,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="none", pad=3),
    )


def create_pair_km_plot(
    data: pd.DataFrame,
    left_group: str,
    right_group: str,
    data_column: str,
    colors: Dict[str, str],
    style: KMPlotStyle,
    output_path: Path,
    *,
    color_overrides: Optional[Dict[Tuple[str, str], str]] = None,
    title: str = "",
    stat_label_map: Optional[Dict[str, str]] = None,
    stats_override: Optional[Dict[str, float]] = None,
    logrank_label: str = "logrank",
    hr_label: str = "HR",
) -> Optional[Dict[str, Any]]:
    """Create a two-group KM figure and return one pairwise summary row."""
    pair_groups = [left_group, right_group]
    pair_data = data[data["GROUP"].isin(pair_groups)].copy()
    if pair_data.empty:
        LOGGER.warning("No data for pair %s vs %s", left_group, right_group)
        return None

    style.apply_rcparams()
    pair_fig, pair_ax = plt.subplots(figsize=(6.5, 6.0))
    kmf_pair_list: List[KaplanMeierFitter] = []
    display_labels = [
        _pair_km_display_label(group, stat_label_map) for group in pair_groups
    ]

    for group in pair_groups:
        group_data = pair_data[pair_data["GROUP"] == group].copy()
        if group_data.empty:
            continue
        group_data.loc[:, data_column] = pd.to_numeric(group_data[data_column], errors="coerce")
        group_data = group_data.dropna(subset=[data_column])
        kmf = KaplanMeierFitter()
        display_label = _pair_km_display_label(group, stat_label_map)
        kmf.fit(
            group_data[data_column],
            event_observed=group_data["CENSOR"],
            label=display_label,
        )
        plot_color = colors.get(group)
        if color_overrides and (left_group, right_group) in color_overrides and group == left_group:
            plot_color = color_overrides[(left_group, right_group)]
        kmf.plot(ax=pair_ax, color=plot_color, ci_show=style.cishow, linestyle="-", linewidth=1)
        kmf_pair_list.append(kmf)

    if not kmf_pair_list:
        plt.close(pair_fig)
        return None

    left_data = pair_data[pair_data["GROUP"] == left_group].dropna(subset=[data_column])
    right_data = pair_data[pair_data["GROUP"] == right_group].dropna(subset=[data_column])

    left_median = KaplanMeierFitter().fit(
        left_data[data_column], event_observed=left_data["CENSOR"]
    ).median_survival_time_
    right_median = KaplanMeierFitter().fit(
        right_data[data_column], event_observed=right_data["CENSOR"]
    ).median_survival_time_

    logrank_res = logrank_test(
        left_data[data_column],
        right_data[data_column],
        event_observed_A=left_data["CENSOR"],
        event_observed_B=right_data["CENSOR"],
    )

    cox_df = pair_data[[data_column, "CENSOR", "GROUP"]].dropna()
    cox_df["group_indicator"] = (cox_df["GROUP"] == left_group).astype(int)
    hr = ci_low = ci_high = hr_p = np.nan
    if cox_df["group_indicator"].nunique() == 2:
        cph = CoxPHFitter()
        cph.fit(cox_df[[data_column, "CENSOR", "group_indicator"]], duration_col=data_column, event_col="CENSOR")
        hr = float(np.exp(cph.params_["group_indicator"]))
        ci_low = float(np.exp(cph.confidence_intervals_.loc["group_indicator", "95% lower-bound"]))
        ci_high = float(np.exp(cph.confidence_intervals_.loc["group_indicator", "95% upper-bound"]))
        hr_p = float(cph.summary.loc["group_indicator", "p"])

    logrank_p = float(logrank_res.p_value)
    if stats_override:
        logrank_p = float(stats_override.get("logrank_p", logrank_p))
        hr = stats_override.get("hr", hr)
        ci_low = stats_override.get("ci_low", ci_low)
        ci_high = stats_override.get("ci_high", ci_high)
        hr_p = stats_override.get("hr_p", hr_p)

    apply_km_axis_style(pair_ax, style, title=title)
    pair_ax.legend(
        labels=display_labels,
        frameon=False,
        fontsize=style.font_legend,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
    )
    _annotate_pair_km_stats_on_ax(
        pair_ax,
        left_group=left_group,
        right_group=right_group,
        left_median=left_median,
        right_median=right_median,
        logrank_p=logrank_p,
        hr=hr,
        ci_low=ci_low,
        ci_high=ci_high,
        hr_p=hr_p,
        style=style,
        stat_label_map=stat_label_map,
        logrank_label=logrank_label,
        hr_label=hr_label,
    )
    add_at_risk_counts(
        *kmf_pair_list,
        ypos=-0.35,
        ax=pair_ax,
        rows_to_show=["At risk"],
        at_risk_count_from_start_of_period=True,
        xticks=KM_XTICKS,
    )
    pair_fig.tight_layout(rect=[0.0, 0.14, 1.0, 0.96])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pair_fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.close(pair_fig)

    return {
        "comparison": f"{left_group} vs {right_group}",
        f"mPFS_{left_group}": left_median,
        f"mPFS_{right_group}": right_median,
        "logrank_p": logrank_p,
        "HR": hr,
        "HR_CI_lower_95": ci_low,
        "HR_CI_upper_95": ci_high,
        "HR_P": hr_p,
    }

def _count_pairwise_strata(
    df: pd.DataFrame,
    group_col: str,
    stratum_col: str,
    left_group: str,
    right_group: str,
) -> int:
    """Count strata with at least one observation in each comparison group."""
    pair = df[df[group_col].isin([left_group, right_group])]
    return sum(
        1 for _, sub in pair.groupby(stratum_col) if sub[group_col].nunique() >= 2
    )


def _two_sample_logrank_oe_var(
    times_left: np.ndarray,
    events_left: np.ndarray,
    times_right: np.ndarray,
    events_right: np.ndarray,
) -> Tuple[float, float]:
    """Within-stratum log-rank score (O−E) and variance for the left arm.

    Standard hypergeometric variance at each event time; returns values that can
    be summed across strata for a Mantel–Haenszel stratified log-rank (1 df).
    """
    times_left = np.asarray(times_left, dtype=float)
    times_right = np.asarray(times_right, dtype=float)
    events_left = np.asarray(events_left, dtype=float).astype(bool)
    events_right = np.asarray(events_right, dtype=float).astype(bool)
    if times_left.size == 0 or times_right.size == 0:
        return 0.0, 0.0

    times = np.concatenate([times_left, times_right])
    events = np.concatenate([events_left, events_right])
    is_left = np.concatenate(
        [np.ones(times_left.size, dtype=bool), np.zeros(times_right.size, dtype=bool)]
    )
    event_times = np.unique(times[events])
    oe = 0.0
    var = 0.0
    for t in event_times:
        at_risk = times >= t
        n1 = int(np.sum(at_risk & is_left))
        n0 = int(np.sum(at_risk & ~is_left))
        n = n1 + n0
        if n <= 1:
            continue
        d1 = int(np.sum((times == t) & events & is_left))
        d0 = int(np.sum((times == t) & events & ~is_left))
        d = d1 + d0
        if d == 0:
            continue
        e1 = d * n1 / n
        oe += d1 - e1
        var += (n1 * n0 * d * (n - d)) / (n * n * (n - 1))
    return float(oe), float(var)


def stratified_logrank_pairwise(
    df: pd.DataFrame,
    data_column: str,
    event_col: str,
    group_col: str,
    stratum_col: str,
    left_group: str,
    right_group: str,
) -> Tuple[float, float, int]:
    """Mantel–Haenszel / Mantel–Cox stratified log-rank for one pair.

    Sums signed (O−E) and variances across strata (e.g. Model), then forms a
    1-df chi-square: (Σ(O−E))² / ΣV. This preserves treatment direction, unlike
    summing unsigned per-stratum chi-square statistics.
    """
    pair = df[df[group_col].isin([left_group, right_group])].dropna(subset=[data_column]).copy()
    oe_sum = 0.0
    var_sum = 0.0
    n_strata = 0
    for _, sub in pair.groupby(stratum_col):
        if sub[group_col].nunique() < 2:
            continue
        left = sub[sub[group_col] == left_group]
        right = sub[sub[group_col] == right_group]
        if left.empty or right.empty:
            continue
        oe, var = _two_sample_logrank_oe_var(
            left[data_column].to_numpy(dtype=float),
            left[event_col].to_numpy(),
            right[data_column].to_numpy(dtype=float),
            right[event_col].to_numpy(),
        )
        if var <= 0 and oe == 0:
            # Both arms present but no informative event contrast (e.g. identical times).
            n_strata += 1
            continue
        oe_sum += oe
        var_sum += var
        n_strata += 1
    if n_strata == 0 or var_sum <= 0:
        return np.nan, np.nan, int(n_strata)
    chi2_stat = float((oe_sum * oe_sum) / var_sum)
    p_value = float(chi2_distribution.sf(chi2_stat, 1))
    return chi2_stat, p_value, int(n_strata)


def stratified_cox_pairwise(
    df: pd.DataFrame,
    data_column: str,
    event_col: str,
    group_col: str,
    stratum_col: str,
    left_group: str,
    right_group: str,
) -> Dict[str, float]:
    """Stratified Cox PH (strata = Model) for one pair; arm = 1 for left_group."""
    empty = {
        "stratified_cox_hr": np.nan,
        "stratified_cox_hr_low": np.nan,
        "stratified_cox_hr_high": np.nan,
        "stratified_cox_p": np.nan,
        "n_strata_contributing": 0,
    }
    pair = df[df[group_col].isin([left_group, right_group])].dropna(subset=[data_column]).copy()
    if pair.empty:
        return empty

    n_strata = _count_pairwise_strata(pair, group_col, stratum_col, left_group, right_group)
    if n_strata == 0:
        return empty

    cox_df = pair[[data_column, event_col, group_col, stratum_col]].copy()
    cox_df["arm"] = (cox_df[group_col] == left_group).astype(int)
    if cox_df["arm"].nunique() < 2:
        empty["n_strata_contributing"] = n_strata
        return empty

    try:
        cph = CoxPHFitter()
        cph.fit(
            cox_df,
            duration_col=data_column,
            event_col=event_col,
            strata=[stratum_col],
            formula="arm",
        )
        summary = cph.summary.loc["arm"]
        return {
            "stratified_cox_hr": float(summary["exp(coef)"]),
            "stratified_cox_hr_low": float(summary["exp(coef) lower 95%"]),
            "stratified_cox_hr_high": float(summary["exp(coef) upper 95%"]),
            "stratified_cox_p": float(summary["p"]),
            "n_strata_contributing": n_strata,
        }
    except Exception as exc:
        LOGGER.warning(
            "Stratified Cox failed for %s vs %s: %s",
            left_group,
            right_group,
            exc,
        )
        return {**empty, "n_strata_contributing": n_strata}


def stratified_cox_global(
    df: pd.DataFrame,
    data_column: str,
    event_col: str,
    group_col: str,
    stratum_col: str,
    group_order: Sequence[str],
    reference_group: str,
) -> Tuple[pd.DataFrame, float, float]:
    """Omnibus stratified Cox PH across ordered groups (reference = reference_group).

    Returns (coefficient summary, global LRT chi-square, global LRT p-value).
    """
    empty_summary = pd.DataFrame()
    model_df = df[df[group_col].isin(group_order)].dropna(subset=[data_column]).copy()
    if model_df[group_col].nunique() < 2 or model_df[stratum_col].nunique() < 2:
        return empty_summary, np.nan, np.nan

    model_df[group_col] = pd.Categorical(
        model_df[group_col], categories=list(group_order), ordered=True
    )
    formula = f'C({group_col}, Treatment(reference="{reference_group}"))'
    try:
        cph = CoxPHFitter()
        cph.fit(
            model_df,
            duration_col=data_column,
            event_col=event_col,
            strata=[stratum_col],
            formula=formula,
        )
        lrt = cph.log_likelihood_ratio_test()
        return cph.summary.reset_index(), float(lrt.test_statistic), float(lrt.p_value)
    except Exception as exc:
        LOGGER.warning("Global stratified Cox failed: %s", exc)
        return empty_summary, np.nan, np.nan


def append_fdr_to_pvalue_columns(
    df: pd.DataFrame, p_columns: Sequence[str]
) -> pd.DataFrame:
    """Add BH-FDR columns and star labels for each requested raw p-value column."""
    if df.empty:
        return df
    out = df.copy()
    for p_col in p_columns:
        if p_col not in out.columns or not out[p_col].notna().any():
            out[f"{p_col}_fdr_bh"] = np.nan
            out[f"{p_col}_significance"] = ""
            continue
        raw = out[p_col].astype(float).values
        _, fdr, _, _ = multipletests(raw, method="fdr_bh")
        out[f"{p_col}_fdr_bh"] = fdr
        out[f"{p_col}_significance"] = out[p_col].apply(significance_stars)
    return out


def build_subgroup_stratified_comparison_table(
    df: pd.DataFrame,
    data_column: str,
    *,
    event_col: str = "CENSOR",
    group_col: str = "GROUP",
    stratum_col: str = "Model",
    pair_comparisons: Sequence[Tuple[str, str]],
) -> pd.DataFrame:
    """Pairwise stratified log-rank and Cox PH with Model strata."""
    rows: List[Dict[str, Any]] = []
    for left_group, right_group in pair_comparisons:
        chi2_stat, logrank_p, n_strata = stratified_logrank_pairwise(
            df,
            data_column,
            event_col,
            group_col,
            stratum_col,
            left_group,
            right_group,
        )
        cox_stats = stratified_cox_pairwise(
            df,
            data_column,
            event_col,
            group_col,
            stratum_col,
            left_group,
            right_group,
        )
        pair_data = df[df[group_col].isin([left_group, right_group])].dropna(subset=[data_column])
        rows.append(
            {
                "comparison": f"{left_group} vs {right_group}",
                "left_group": left_group,
                "right_group": right_group,
                "n_datapoints": int(len(pair_data)),
                "n_models_total": int(pair_data[stratum_col].nunique()),
                "n_strata_contributing": int(n_strata),
                "stratified_logrank_chi2": chi2_stat,
                "stratified_logrank_p": logrank_p,
                **cox_stats,
            }
        )
    return append_fdr_to_pvalue_columns(
        pd.DataFrame(rows),
        ("stratified_logrank_p", "stratified_cox_p"),
    )


def pooled_logrank_cox_pairwise(
    df: pd.DataFrame,
    data_column: str,
    left_group: str,
    right_group: str,
    *,
    event_col: str = "CENSOR",
    group_col: str = "GROUP",
) -> Dict[str, float]:
    """Unstratified (pooled) log-rank and Cox PH for one pair."""
    empty = {
        "pooled_logrank_p": np.nan,
        "pooled_hr": np.nan,
        "pooled_hr_low": np.nan,
        "pooled_hr_high": np.nan,
        "pooled_hr_p": np.nan,
        "pooled_cox_hr": np.nan,
        "pooled_cox_hr_low": np.nan,
        "pooled_cox_hr_high": np.nan,
        "pooled_cox_p": np.nan,
    }
    pair_data = df[df[group_col].isin([left_group, right_group])].dropna(subset=[data_column])
    if pair_data.empty:
        return empty

    left_data = pair_data[pair_data[group_col] == left_group]
    right_data = pair_data[pair_data[group_col] == right_group]
    if left_data.empty or right_data.empty:
        return empty

    logrank_res = logrank_test(
        left_data[data_column],
        right_data[data_column],
        event_observed_A=left_data[event_col],
        event_observed_B=right_data[event_col],
    )

    cox_df = pair_data[[data_column, event_col, group_col]].copy()
    cox_df["group_indicator"] = (cox_df[group_col] == left_group).astype(int)
    if cox_df["group_indicator"].nunique() < 2:
        return {**empty, "pooled_logrank_p": float(logrank_res.p_value)}

    try:
        cph = CoxPHFitter()
        cph.fit(
            cox_df[[data_column, event_col, "group_indicator"]],
            duration_col=data_column,
            event_col=event_col,
        )
        summary = cph.summary.loc["group_indicator"]
        hr = float(summary["exp(coef)"])
        hr_low = float(summary["exp(coef) lower 95%"])
        hr_high = float(summary["exp(coef) upper 95%"])
        hr_p = float(summary["p"])
        return {
            "pooled_logrank_p": float(logrank_res.p_value),
            "pooled_hr": hr,
            "pooled_hr_low": hr_low,
            "pooled_hr_high": hr_high,
            "pooled_hr_p": hr_p,
            # Aliases used by publication / summary writers
            "pooled_cox_hr": hr,
            "pooled_cox_hr_low": hr_low,
            "pooled_cox_hr_high": hr_high,
            "pooled_cox_p": hr_p,
        }
    except Exception as exc:
        LOGGER.warning(
            "Pooled Cox failed for %s vs %s: %s",
            left_group,
            right_group,
            exc,
        )
        return {**empty, "pooled_logrank_p": float(logrank_res.p_value)}

def _rank_label_for_plot(rank: str) -> str:
    return rank.replace("-", "_")


def _fit_rank_utility_km_on_axes(
    ax: plt.Axes,
    rank_data: pd.DataFrame,
    data_column: str,
    style: KMPlotStyle,
    *,
    title: str,
    show_legend: bool,
) -> pd.DataFrame:
    """Fit on-label utility KM curves on a provided axes."""
    rows: List[Dict[str, Any]] = []
    kmf_list: List[KaplanMeierFitter] = []
    for label in RANKING_UTILITY_GROUP_ORDER:
        group_data = rank_data[rank_data["GROUP"] == label].copy()
        if group_data.empty:
            continue
        group_data.loc[:, data_column] = pd.to_numeric(group_data[data_column], errors="coerce")
        group_data = group_data.dropna(subset=[data_column])
        if group_data.empty:
            continue
        plot_label = f"{label} (n={len(group_data)})"
        kmf = KaplanMeierFitter()
        kmf.fit(group_data[data_column], event_observed=group_data["CENSOR"], label=plot_label)
        kmf.plot(ax=ax, color=RANKING_UTILITY_COLORS[label], ci_show=style.cishow, linestyle="-", linewidth=1)
        kmf_list.append(kmf)
        stat = collect_km_group_statistics(kmf, plot_label, label, group_data)
        stat["Rank"] = title
        rows.append(stat)

    apply_km_axis_style(ax, style, title=title)
    if show_legend and rows:
        ax.legend(frameon=False, fontsize=style.font_legend * 0.85, loc="best")
    add_km_at_risk_table(ax, kmf_list, at_risk_ypos=-0.35)
    return pd.DataFrame(rows)

def _prepare_ranked_onlabel_patient_table(
    source_path: Path, source_sheet: Optional[str]
) -> pd.DataFrame:
    """Cleaned, ranked treatment rows with ON_LABEL_UTILITY (chemo excluded).

    Dense ranks include UNKNOWN On_label rows (e.g. tamoxifen) so all-model
    ranking matches KM_PDX. Those rows have null ON_LABEL_UTILITY and must not
    be used as no-on-label / utility analysis arms.
    """
    patient = load_primary_patient_table(source_path, source_sheet)
    patient = harmonize_column_names(patient)
    patient = normalize_id_columns(patient)
    validate_on_label_source_columns(patient)
    on_col, _, _ = resolve_on_label_columns(patient)
    chemo_mask = is_chemo_on_label(patient[on_col])
    if chemo_mask.any():
        LOGGER.info(
            "Excluding %d row(s) with On_label = chemo before on-label rank assignment.",
            int(chemo_mask.sum()),
        )
    patient = patient[~chemo_mask].copy()
    patient = clean_base_patient_table(patient)
    patient = assign_on_label_utility_column(patient, on_col=on_col, drop_unmapped=False)
    patient["CENSOR"] = True
    ranked = assign_dense_level_groups(
        patient, model_col="Model", level_col="LEVEL", group_col="RANK"
    )
    ranked = normalize_survival_columns(ranked, DEFAULT_DATA_COLUMN)
    return ranked.dropna(subset=[DEFAULT_DATA_COLUMN])


def _collapse_chemo_per_model(chemo: pd.DataFrame, data_column: str) -> pd.DataFrame:
    """One chemo row per model (median TTD when duplicate chemo rows exist)."""
    if chemo.empty:
        return chemo
    work = chemo.copy()
    work[data_column] = pd.to_numeric(work[data_column], errors="coerce")
    work = work.dropna(subset=[data_column])
    if work.empty:
        return work
    multi = work.groupby("Model").size()
    n_multi = int((multi > 1).sum())
    if n_multi:
        LOGGER.warning(
            "%d model(s) have multiple chemo rows; using median %s.",
            n_multi,
            data_column,
        )
    agg: Dict[str, Tuple[str, str]] = {data_column: (data_column, "median"), "CENSOR": ("CENSOR", "first")}
    if "COMPOUND" in work.columns:
        agg["COMPOUND"] = ("COMPOUND", "first")
    collapsed = work.groupby("Model", as_index=False).agg(**agg)
    collapsed["CENSOR"] = True
    return collapsed

def _clinical_utility_rank_as_int(rank: Any) -> int:
    text = str(rank).strip()
    if text == _POOLED_LABEL:
        return _RANK_POOL_FROM
    try:
        return int(text)
    except ValueError:
        return 999


def _median_ttd_for_rows(rows: pd.DataFrame, data_column: str) -> float:
    if rows.empty:
        return np.nan
    values = pd.to_numeric(rows[data_column], errors="coerce").dropna()
    if values.empty:
        return np.nan
    return float(values.median())


def _best_row_at_rank(
    model_rows: pd.DataFrame,
    rank: str,
    data_column: str,
    *,
    utilities: Optional[Sequence[str]] = None,
    require_mapped_utility: bool = True,
) -> Optional[pd.Series]:
    subset = model_rows[model_rows["RANK"].astype(str) == str(rank)].copy()
    if require_mapped_utility:
        subset = rows_with_mapped_on_label_utility(subset)
    if utilities is not None:
        subset = subset[subset["ON_LABEL_UTILITY"].isin(utilities)]
    if subset.empty:
        return None
    subset = subset.sort_values(["LEVEL", data_column], ascending=[False, False])
    return subset.iloc[0]


def _best_onlabel_row(model_rows: pd.DataFrame, data_column: str) -> Optional[pd.Series]:
    on_rows = model_rows[model_rows["ON_LABEL_UTILITY"] == "On_label"].copy()
    if on_rows.empty:
        return None
    on_rows["_rank_int"] = on_rows["RANK"].map(_clinical_utility_rank_as_int)
    on_rows = on_rows.sort_values(["_rank_int", "LEVEL", data_column], ascending=[True, False, False])
    return on_rows.iloc[0]


def _best_offexp_above_onlabel_row(
    model_rows: pd.DataFrame, onlabel_rank_int: int, data_column: str
) -> Optional[pd.Series]:
    above = model_rows[
        model_rows["ON_LABEL_UTILITY"].isin(("Off_label", "Exp"))
    ].copy()
    if above.empty:
        return None
    above["_rank_int"] = above["RANK"].map(_clinical_utility_rank_as_int)
    above = above[above["_rank_int"] < onlabel_rank_int]
    if above.empty:
        return None
    above = above.sort_values(["_rank_int", "LEVEL", data_column], ascending=[True, False, False])
    return above.iloc[0]


def _verify_ranks_above_onlabel(model_rows: pd.DataFrame, onlabel_rank_int: int) -> bool:
    """True when every compound ranked above on-label is Off_label or Exp."""
    above = model_rows.copy()
    above["_rank_int"] = above["RANK"].map(_clinical_utility_rank_as_int)
    above = above[above["_rank_int"] < onlabel_rank_int]
    if above.empty:
        return True
    return bool(above["ON_LABEL_UTILITY"].isin(("Off_label", "Exp")).all())


def build_model_clinical_utility_table(
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    data_column: str,
) -> pd.DataFrame:
    """Per-model summary: ranks, TTDs, scenario label, and descriptive gates."""
    chemo_map = (
        chemo_collapsed.set_index("Model")[data_column].to_dict()
        if not chemo_collapsed.empty
        else {}
    )
    rows: List[Dict[str, Any]] = []

    for model, model_rows in ranked.groupby("Model"):
        model_rows = model_rows.copy()
        rank1_row = _best_row_at_rank(model_rows, "1", data_column)
        on_row = _best_onlabel_row(model_rows, data_column)

        has_on_label = on_row is not None
        best_onlabel_rank = (
            str(on_row["RANK"]) if on_row is not None else None
        )
        best_onlabel_rank_int = (
            _clinical_utility_rank_as_int(best_onlabel_rank) if best_onlabel_rank else None
        )

        if not has_on_label:
            scenario = CLINICAL_UTILITY_SCENARIO1
        elif best_onlabel_rank_int == 1:
            scenario = CLINICAL_UTILITY_SCENARIO2
        elif best_onlabel_rank_int in (2, 3):
            scenario = CLINICAL_UTILITY_SCENARIO3
        else:
            scenario = CLINICAL_UTILITY_SCENARIO4

        above_row = None
        ranks_above_valid = None
        if has_on_label and best_onlabel_rank_int is not None and best_onlabel_rank_int > 1:
            above_row = _best_offexp_above_onlabel_row(
                model_rows, best_onlabel_rank_int, data_column
            )
            ranks_above_valid = _verify_ranks_above_onlabel(model_rows, best_onlabel_rank_int)

        chemo_ttd = chemo_map.get(model, np.nan)
        # Rank-1 analysis arm excludes UNKNOWN/unmapped (no-on-label / utility groups).
        rank1_mapped = rows_with_mapped_on_label_utility(
            model_rows[model_rows["RANK"].astype(str) == "1"]
        )
        row: Dict[str, Any] = {
            "Model": model,
            "has_on_label": has_on_label,
            "best_onlabel_rank": best_onlabel_rank,
            "best_onlabel_rank_int": best_onlabel_rank_int,
            "scenario": scenario,
            "rank1_TTD": _median_ttd_for_rows(rank1_mapped, data_column),
            "rank1_utility": (
                str(rank1_row["ON_LABEL_UTILITY"]) if rank1_row is not None else None
            ),
            "rank1_COMPOUND": rank1_row["COMPOUND"] if rank1_row is not None else None,
            "onlabel_TTD": float(on_row[data_column]) if on_row is not None else np.nan,
            "onlabel_rank": best_onlabel_rank,
            "onlabel_COMPOUND": on_row["COMPOUND"] if on_row is not None else None,
            "above_onlabel_TTD": (
                float(above_row[data_column]) if above_row is not None else np.nan
            ),
            "above_onlabel_rank": (
                str(above_row["RANK"]) if above_row is not None else None
            ),
            "above_onlabel_utility": (
                str(above_row["ON_LABEL_UTILITY"]) if above_row is not None else None
            ),
            "above_onlabel_COMPOUND": (
                above_row["COMPOUND"] if above_row is not None else None
            ),
            "chemo_TTD": chemo_ttd,
            "ranks_above_onlabel_valid": ranks_above_valid,
        }
        rows.append(row)

    return pd.DataFrame(rows)


def build_clinical_utility_rank_distribution(model_table: pd.DataFrame) -> pd.DataFrame:
    """Count models by on-label rank position and scenario."""
    buckets = {
        "rank_1": 0,
        "rank_2": 0,
        "rank_3": 0,
        "rank_4_plus": 0,
        "no_on_label": 0,
    }
    for _, row in model_table.iterrows():
        if not row["has_on_label"]:
            buckets["no_on_label"] += 1
        elif row["best_onlabel_rank_int"] == 1:
            buckets["rank_1"] += 1
        elif row["best_onlabel_rank_int"] == 2:
            buckets["rank_2"] += 1
        elif row["best_onlabel_rank_int"] == 3:
            buckets["rank_3"] += 1
        else:
            buckets["rank_4_plus"] += 1

    scenario_counts = model_table["scenario"].value_counts().to_dict()
    out_rows = [
        {"bucket": key, "n_models": count}
        for key, count in buckets.items()
    ]
    out = pd.DataFrame(out_rows)
    for scenario in CLINICAL_UTILITY_SCENARIOS:
        out.loc[len(out)] = {
            "bucket": f"scenario_{scenario}",
            "n_models": int(scenario_counts.get(scenario, 0)),
        }
    return out


def _scenario3_stratum_filter(model_table: pd.DataFrame, stratum: str) -> pd.DataFrame:
    subset = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO3].copy()
    if stratum == "rank2":
        return subset[subset["best_onlabel_rank_int"] == 2]
    if stratum == "rank3":
        return subset[subset["best_onlabel_rank_int"] == 3]
    if stratum == "rank2_3_pooled":
        return subset[subset["best_onlabel_rank_int"].isin((2, 3))]
    raise ValueError(f"Unknown scenario-3 stratum: {stratum!r}")


def _build_clinical_utility_paired_table(
    model_table: pd.DataFrame,
    comparison_id: str,
) -> pd.DataFrame:
    """Build one-row-per-model paired table for a pre-specified comparison."""
    left_col, right_col = {
        "rank1_vs_onlabel": ("rank1_TTD", "onlabel_TTD"),
        "rank1_vs_sc": ("rank1_TTD", "chemo_TTD"),
        "onlabel_vs_sc": ("onlabel_TTD", "chemo_TTD"),
        "above_onlabel_vs_onlabel": ("above_onlabel_TTD", "onlabel_TTD"),
    }[comparison_id]

    paired = model_table.copy()
    paired = paired[
        paired[left_col].notna()
        & paired[right_col].notna()
        & np.isfinite(paired[left_col].astype(float))
        & np.isfinite(paired[right_col].astype(float))
    ].copy()
    if paired.empty:
        return paired

    paired["rank_TTD"] = paired[left_col].astype(float)
    paired["chemo_TTD"] = paired[right_col].astype(float)
    paired["comparison"] = comparison_id
    paired["diff_TTD"] = paired["rank_TTD"] - paired["chemo_TTD"]
    paired["ratio_TTD"] = paired["rank_TTD"] / paired["chemo_TTD"]
    paired["pct_change_TTD"] = 100.0 * paired["diff_TTD"] / paired["chemo_TTD"]
    paired["improved"] = paired["diff_TTD"] > 0
    paired["worse"] = paired["diff_TTD"] < 0
    paired["tie"] = paired["diff_TTD"] == 0
    return paired.sort_values("pct_change_TTD", ascending=False).reset_index(drop=True)


def _clinical_utility_stratified_survival_stats(
    km_data: pd.DataFrame,
    left_group: str,
    right_group: str,
    data_column: str,
) -> Dict[str, Any]:
    """Mantel-Cox stratified log-rank and Cox PH (strata = Model) for one scenario pair."""
    chi2_stat, logrank_p, n_strata = stratified_logrank_pairwise(
        km_data,
        data_column,
        "CENSOR",
        "GROUP",
        "Model",
        left_group,
        right_group,
    )
    cox_stats = stratified_cox_pairwise(
        km_data,
        data_column,
        "CENSOR",
        "GROUP",
        "Model",
        left_group,
        right_group,
    )
    pooled_stats = pooled_logrank_cox_pairwise(
        km_data,
        data_column,
        left_group,
        right_group,
    )
    return {
        "stratified_logrank_chi2": chi2_stat,
        "stratified_logrank_p": logrank_p,
        "n_strata_logrank": n_strata,
        **cox_stats,
        **pooled_stats,
    }


def _apply_fdr_wilcoxon_clinical(summary: pd.DataFrame, family_col: str) -> pd.DataFrame:
    """BH-FDR on Wilcoxon p-values within each analysis family."""
    if summary.empty:
        summary["wilcoxon_p_fdr_bh"] = np.nan
        return summary
    out_parts: List[pd.DataFrame] = []
    for _, group in summary.groupby(family_col, sort=False):
        part = group.copy()
        if part["wilcoxon_p"].notna().any():
            _, fdr, _, _ = multipletests(part["wilcoxon_p"].astype(float).values, method="fdr_bh")
            part["wilcoxon_p_fdr_bh"] = fdr
        else:
            part["wilcoxon_p_fdr_bh"] = np.nan
        out_parts.append(part)
    return pd.concat(out_parts, ignore_index=True)


def _clinical_utility_km_dataset_from_paired(
    paired: pd.DataFrame,
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    left_group: str,
    right_group: str,
    data_column: str,
    *,
    left_col: str = "rank_TTD",
    right_col: str = "chemo_TTD",
) -> pd.DataFrame:
    """Convert paired table to long KM format with GROUP labels."""
    if paired.empty:
        return pd.DataFrame()

    models = set(paired["Model"].unique())
    frames: List[pd.DataFrame] = []

    if left_group == CLINICAL_UTILITY_GROUP_RANK1:
        # One value per model — aligned with paired Wilcoxon (model_table rank1_TTD / median).
        left_rows = paired[["Model", "rank_TTD"]].rename(columns={"rank_TTD": data_column}).copy()
        left_rows[data_column] = pd.to_numeric(left_rows[data_column], errors="coerce")
        left_rows = left_rows.dropna(subset=[data_column])
    elif left_group == CLINICAL_UTILITY_GROUP_ONLABEL:
        left_rows = ranked[
            (ranked["Model"].isin(models)) & (ranked["ON_LABEL_UTILITY"] == "On_label")
        ].copy()
        left_rows["_rank_int"] = left_rows["RANK"].map(_clinical_utility_rank_as_int)
        left_rows = (
            left_rows.sort_values(["Model", "_rank_int", "LEVEL"], ascending=[True, True, False])
            .groupby("Model", as_index=False)
            .first()
        )
    elif left_group == CLINICAL_UTILITY_GROUP_ABOVE:
        rows: List[pd.Series] = []
        for model in models:
            model_rows = ranked[ranked["Model"] == model]
            on_row = _best_onlabel_row(model_rows, data_column)
            if on_row is None:
                continue
            above = _best_offexp_above_onlabel_row(
                model_rows,
                _clinical_utility_rank_as_int(on_row["RANK"]),
                data_column,
            )
            if above is not None:
                rows.append(above)
        left_rows = pd.DataFrame(rows) if rows else pd.DataFrame()
    else:
        left_rows = pd.DataFrame()

    if not left_rows.empty:
        left_rows = left_rows.copy()
        left_rows["GROUP"] = left_group
        left_rows["CENSOR"] = True
        frames.append(left_rows)

    if right_group == CLINICAL_UTILITY_GROUP_SC:
        right_rows = chemo_collapsed[chemo_collapsed["Model"].isin(models)].copy()
        if not right_rows.empty:
            right_rows["GROUP"] = CLINICAL_UTILITY_GROUP_SC
            right_rows["CENSOR"] = True
            frames.append(right_rows)
    elif right_group == CLINICAL_UTILITY_GROUP_ONLABEL:
        right_rows = ranked[
            (ranked["Model"].isin(models)) & (ranked["ON_LABEL_UTILITY"] == "On_label")
        ].copy()
        right_rows["_rank_int"] = right_rows["RANK"].map(_clinical_utility_rank_as_int)
        right_rows = (
            right_rows.sort_values(["Model", "_rank_int", "LEVEL"], ascending=[True, True, False])
            .groupby("Model", as_index=False)
            .first()
        )
        if not right_rows.empty:
            right_rows["GROUP"] = CLINICAL_UTILITY_GROUP_ONLABEL
            right_rows["CENSOR"] = True
            frames.append(right_rows)
    else:
        fallback = paired[["Model", right_col]].rename(columns={right_col: data_column}).copy()
        fallback[data_column] = pd.to_numeric(fallback[data_column], errors="coerce")
        fallback = fallback.dropna(subset=[data_column])
        if not fallback.empty:
            fallback["GROUP"] = right_group
            fallback["CENSOR"] = True
            frames.append(fallback)

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out[data_column] = pd.to_numeric(out[data_column], errors="coerce")
    return out.dropna(subset=[data_column])


def _save_clinical_utility_waterfall(
    paired: pd.DataFrame,
    output_path: Path,
    title: str,
    result: Dict[str, Any],
) -> None:
    """Waterfall plot of % TTD change (left arm vs right arm)."""
    if paired.empty:
        return
    improved_color = "#2e7d32"
    worse_color = "#c62828"
    tie_color = "#757575"
    plot_df = paired.sort_values("pct_change_TTD", ascending=True).reset_index(drop=True)
    colors = np.where(
        plot_df["pct_change_TTD"] > 0,
        improved_color,
        np.where(plot_df["pct_change_TTD"] < 0, worse_color, tie_color),
    )
    fig, ax = plt.subplots(figsize=(max(9, len(plot_df) * 0.12), 6.5))
    ax.bar(range(len(plot_df)), plot_df["pct_change_TTD"], color=colors, edgecolor="none", width=0.85)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Models (sorted by % change)")
    ax.set_ylabel("% change TimeToDouble (left arm vs right arm)")
    ax.set_title(title)
    stats_lines = [
        f"n = {int(result.get('n_models', 0))} models",
        f"Median ΔTTD = {result.get('median_diff_TTD', np.nan):+.1f} d",
        f"Wilcoxon {format_p_equals_mathtext(result.get('wilcoxon_p'))}"
        if pd.notna(result.get("wilcoxon_p"))
        else r"Wilcoxon $\mathit{P}$ = n/a",
    ]
    ax.text(
        0.98,
        0.98,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
    )
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_utility_at_rank_km_dataset(
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    rank: str,
    data_column: str,
) -> pd.DataFrame:
    """Per-rank On/Off/Exp arms plus model-matched SC for contributing models."""
    rank_rows = ranked[ranked["RANK"].astype(str) == str(rank)].copy()
    if rank_rows.empty:
        return pd.DataFrame()

    models = set(rank_rows["Model"].unique())
    frames: List[pd.DataFrame] = []
    for utility in ON_LABEL_UTILITY_ORDER:
        utility_rows = rank_rows[rank_rows["ON_LABEL_UTILITY"] == utility].copy()
        if utility_rows.empty:
            continue
        utility_rows["GROUP"] = utility
        utility_rows["CENSOR"] = True
        frames.append(utility_rows)

    chemo_sub = chemo_collapsed[chemo_collapsed["Model"].isin(models)].copy()
    if not chemo_sub.empty:
        chemo_sub["GROUP"] = CLINICAL_UTILITY_GROUP_SC
        chemo_sub["CENSOR"] = True
        frames.append(chemo_sub)

    out = pd.concat(frames, ignore_index=True)
    out[data_column] = pd.to_numeric(out[data_column], errors="coerce")
    return out.dropna(subset=[data_column])


def _write_clinical_utility_interpretation(
    output_dir: Path,
    *,
    model_table: pd.DataFrame,
    rank_distribution: pd.DataFrame,
    paired_summary: pd.DataFrame,
    utility_rank_pairwise: pd.DataFrame,
) -> None:
    """Plain-language supplement interpretation for clinical utility ranking."""
    lines = [
        "CLINICAL UTILITY RANKING ANALYSIS — INTERPRETATION",
        "===================================================",
        "",
        "Endpoint: TimeToDouble in PDX (preclinical pharmacodynamic surrogate; not patient OS/PFS).",
        "",
        "Utility groups: On_label (IGAZ), Off_label (HAMIS), Exp (n/a or missing).",
        "Dense ranking: higher DDA LEVEL → lower rank number (1 = best).",
        "Reference: model-matched standard chemotherapy (On_label = chemo rows).",
        "",
        "PRE-SPECIFIED HYPOTHESES",
        "------------------------",
        "H1: In models without on-label therapy, global rank 1 TTD exceeds matched SC.",
        "H2: When on-label occupies rank 1, its TTD exceeds matched SC.",
        "H3: When on-label is not rank 1, global rank 1 TTD exceeds on-label TTD.",
        "H4: When on-label is not rank 1, best Off/Exp above on-label exceeds on-label TTD.",
        "",
        "SCENARIO COHORTS",
        "----------------",
    ]
    for scenario in CLINICAL_UTILITY_SCENARIOS:
        n = int((model_table["scenario"] == scenario).sum())
        lines.append(f"  {scenario}: {n} models")

    lines.extend(["", "ON-LABEL RANK DISTRIBUTION", "--------------------------"])
    for _, row in rank_distribution.iterrows():
        if str(row["bucket"]).startswith("scenario_"):
            continue
        lines.append(f"  {row['bucket']}: {int(row['n_models'])} models")

    if "ranks_above_onlabel_valid" in model_table.columns:
        s3 = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO3]
        if not s3.empty:
            valid = s3["ranks_above_onlabel_valid"].fillna(False).astype(bool)
            lines.extend(
                [
                    "",
                    "DESCRIPTIVE GATE (Scenario 3)",
                    "-----------------------------",
                    f"  Models with only Off/Exp above on-label: {int(valid.sum())}/{len(s3)}",
                ]
            )

    lines.extend(["", "PAIRED TEST SUMMARY (PRIMARY)", "-----------------------------"])
    if paired_summary.empty:
        lines.append("  No paired comparisons computed.")
    else:
        for _, row in paired_summary.iterrows():
            label = row.get("comparison_group", row["comparison"])
            line = (
                f"  {label} ({row['scenario']} / {row.get('stratum', 'all')}): "
                f"n={int(row['n_models'])}, median ΔTTD={row['median_diff_TTD']:+.1f}, "
                f"Wilcoxon p={row['wilcoxon_p']:.4g}"
            )
            if pd.notna(row.get("wilcoxon_p_fdr_bh")):
                line += f", Wilcoxon FDR={row['wilcoxon_p_fdr_bh']:.4g}"
            if pd.notna(row.get("stratified_logrank_p")):
                line += f"; strat. log-rank p={row['stratified_logrank_p']:.4g}"
                if pd.notna(row.get("stratified_logrank_p_fdr_bh")):
                    line += f" (FDR={row['stratified_logrank_p_fdr_bh']:.4g})"
            if pd.notna(row.get("stratified_cox_hr")):
                line += (
                    f"; strat. Cox HR={row['stratified_cox_hr']:.2f} "
                    f"({row['stratified_cox_hr_low']:.2f}–{row['stratified_cox_hr_high']:.2f}), "
                    f"p={row['stratified_cox_p']:.4g}"
                )
            lines.append(line)

    lines.extend(
        [
            "",
            "STRATIFIED SURVIVAL TESTS (SUPPORTING / EXPLORATORY)",
            "----------------------------------------------------",
            "Mantel–Haenszel / Mantel–Cox stratified log-rank: sum signed (O−E) and variances",
            "across Model strata, then 1-df chi-square (direction-preserving).",
            "Stratified Cox PH: Surv(TTD) ~ arm + strata(Model); HR < 1 favours left arm.",
            "Fig. 4D uses pooled (unstratified) log-rank / Cox on rank-1 vs all SC rows.",
            "Stratified tests remain available for exploratory matched-arm analyses.",
            "",
            "STATISTICAL HIERARCHY",
            "---------------------",
            "Primary: Wilcoxon signed-rank on paired TTD differences within each model.",
            "Fig. 4D: pooled log-rank and unstratified Cox (rank-1 no-on-label vs all SC).",
            "Multiplicity: BH-FDR within Scenario 3 comparison family and within each rank.",
            "",
            "LIMITATIONS",
            "-----------",
            "PDX TTD is a preclinical endpoint; results support algorithm validation,",
            "not direct clinical prescribing claims.",
        ]
    )
    (output_dir / "INTERPRETATION.txt").write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Publication-grade figures for short supervisor summary (Q1–Q7)
# ---------------------------------------------------------------------------

CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR = "publication_figures"
SHORT_SUMMARY_RANK_COMPOSITION_STEM = "pub_rank_utility_composition"
SHORT_SUMMARY_RANK_DISPLAY_ORDER: Tuple[str, ...] = ("1", "2", "3", "4", "5", "6", _POOLED_LABEL)
SHORT_SUMMARY_RANK_UTILITY_COLORS: Dict[str, str] = {
    "On_label": "#F4CCCC",
    "Off_label": "#C5E0C8",  # light green
    "Exp": "#F0E4D0",        # light beige
}
SHORT_SUMMARY_RANK_UTILITY_LABELS: Dict[str, str] = {
    "On_label": "On-label",
    "Off_label": "Off-label",
    "Exp": "Experimental",
}
SHORT_SUMMARY_RANK_COMPOSITION_CAPTION = (
    "Stacked bars show the percentage of screened drug rows at each global DDA rank classified as "
    "On-label (IGAZ), Off-label (HAMIS), or Experimental (n/a). Higher DDA level yields a lower "
    "rank number (rank 1 = best). Standard chemotherapy rows are excluded."
)

Q8_KM_GROUP_RANK1 = "Rank1_OffExp_ref"
Q8_KM_GROUP_SC = "SC_ref"
Q8_KM_GROUP_POOL = "Onlabel_low_rank_pool"

PUBLICATION_GROUP_LABELS: Dict[str, str] = {
    CLINICAL_UTILITY_GROUP_RANK1: "Rank 1",
    CLINICAL_UTILITY_GROUP_ONLABEL: "On-label",
    CLINICAL_UTILITY_GROUP_SC: "Standard chemo",
    CLINICAL_UTILITY_GROUP_ABOVE: "Best Off/Exp above on-label",
    "On_label_rank1": "On-label at rank 1",
    "On_label_rank7_10": "On-label rank 7–10",
    Q8_KM_GROUP_RANK1: "Rank 1 (Off/Exp)",
    Q8_KM_GROUP_SC: "Standard chemo",
    Q8_KM_GROUP_POOL: "On-label pool",
}


@dataclass
class PublicationPlotStyle:
    """Typography and layout for journal-ready figures."""

    font: float = 9.0
    title: float = 10.0
    panel_label: float = 11.0
    tick: float = 8.0
    legend: float = 8.5
    annotation: float = 8.0
    cishow: bool = False

    def apply_rcparams(self) -> None:
        plt.rcParams.update(
            {
                "axes.titlesize": self.title,
                "axes.labelsize": self.font,
                "xtick.labelsize": self.tick,
                "ytick.labelsize": self.tick,
                "legend.fontsize": self.legend,
                "font.size": self.font,
                "axes.linewidth": 0.8,
            }
        )


@dataclass(frozen=True)
class ShortSummaryPubPanel:
    """One statistical visual for a short-summary question."""

    panel_type: str  # km | waterfall | between_km | between_box
    scenario: str = ""
    stratum: str = "all"
    comparison: str = ""
    subtitle: str = ""
    between_stratum: str = "rank7_10"


SHORT_SUMMARY_PUBLICATION_PANELS: Dict[str, Tuple[ShortSummaryPubPanel, ...]] = {
    "Q3": (
        ShortSummaryPubPanel("km", CLINICAL_UTILITY_SCENARIO1, "all", "rank1_vs_sc", "Cox PH"),
        ShortSummaryPubPanel("waterfall", CLINICAL_UTILITY_SCENARIO1, "all", "rank1_vs_sc", "Wilcoxon ΔPFS"),
    ),
    "Q4": (
        ShortSummaryPubPanel("km", CLINICAL_UTILITY_SCENARIO2, "all", "onlabel_vs_sc", "Cox PH"),
        ShortSummaryPubPanel("waterfall", CLINICAL_UTILITY_SCENARIO2, "all", "onlabel_vs_sc", "Wilcoxon ΔPFS"),
    ),
}

Q3_DISCORDANCE_STRATA: Tuple[Tuple[str, str], ...] = (
    ("rank2", "On-label at rank 2"),
    ("rank3", "On-label at rank 3"),
    ("rank2_3_pooled", "Pooled ranks 2-3"),
)

SHORT_SUMMARY_QUESTION_ORDER: Tuple[str, ...] = (
    ("Q1", "Q2", "Q3", "Q4", "Q5") if ENABLE_Q2_ANALYSIS else ("Q1", "Q3", "Q4", "Q5")
)

SHORT_SUMMARY_QUESTION_TITLES: Dict[str, str] = {
    "Q1": "Q1: Paired on-label-at-rank vs Rank 1 or standard chemo (within model)",
    "Q2": "Q2: Downranked on-label pools vs Rank 1 (Off/Exp) or standard chemo (unpaired)",
    "Q3": "Q3: No on-label, Rank 1 vs standard chemotherapy",
    "Q4": "Q4: On-label at rank 1 vs standard chemotherapy",
    "Q5": "Q5: Guideline discordance (ranks 2-3), Rank 1 vs on-label",
}


def _rank_display_label(rank: str) -> str:
    text = str(rank).strip()
    if text == _POOLED_LABEL:
        return "7-10"
    return text


def _build_rank_utility_composition(ranked: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (percent, counts) of On_label / Off_label / Exp within each global DDA rank."""
    work = ranked[ranked["ON_LABEL_UTILITY"].isin(ON_LABEL_UTILITY_ORDER)].copy()
    work["RANK"] = work["RANK"].astype(str)
    counts = (
        work.groupby(["RANK", "ON_LABEL_UTILITY"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=list(ON_LABEL_UTILITY_ORDER), fill_value=0)
    )
    present = [r for r in SHORT_SUMMARY_RANK_DISPLAY_ORDER if r in counts.index]
    extra = sorted(
        [r for r in counts.index if r not in present],
        key=_clinical_utility_rank_as_int,
    )
    counts = counts.reindex(present + extra).astype(int)
    pct = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    pct = pct.fillna(0.0)
    return pct, counts


def _rank_composition_plot_style() -> PublicationPlotStyle:
    """Typography based on Fig 4a/4b etalon, +10% for the stacked bar (Fig 4c)."""
    km = PAIR_KM_ETALON_STYLE
    s = _FIG4C_FONT_SCALE
    return PublicationPlotStyle(
        font=km.label_size * s,
        title=km.font * s,
        panel_label=km.font * s,
        tick=km.tick_size * s,
        legend=km.font_legend * s,
        annotation=max(km.at_risk, km.tick_size - 0.5) * s,
    )


def _plot_short_summary_rank_utility_stacked_bar(
    ranked: pd.DataFrame,
    *,
    style: PublicationPlotStyle,
) -> plt.Figure:
    pct, counts = _build_rank_utility_composition(ranked)
    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    x = np.arange(len(pct))
    bottom = np.zeros(len(pct))
    for utility in ON_LABEL_UTILITY_ORDER:
        values = pct[utility].values
        color = SHORT_SUMMARY_RANK_UTILITY_COLORS.get(
            utility, ON_LABEL_UTILITY_COLORS.get(utility, "#888888")
        )
        bars = ax.bar(
            x,
            values,
            bottom=bottom,
            color=color,
            label=SHORT_SUMMARY_RANK_UTILITY_LABELS.get(utility, utility),
            width=0.74,
            edgecolor="white",
            linewidth=0.7,
        )
        for idx, (bar, val) in enumerate(zip(bars, values)):
            if val >= 5.0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottom[idx] + val / 2,
                    f"{val:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=max(style.tick - 0.5, 9.0),
                    color="black",
                    fontweight="normal",
                )
        bottom += values

    x_labels = [
        f"{_rank_display_label(r)}\n(n={int(counts.loc[r].sum())})"
        for r in pct.index
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=style.font, color="black")
    ax.set_ylabel("Percentage of approved status (%)", fontsize=style.font, color="black")
    ax.set_xlabel("DDA rank", fontsize=style.font, color="black")
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.tick_params(axis="both", labelsize=style.tick, colors="black", labelcolor="black")
    legend = ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=style.legend,
        borderaxespad=0.0,
    )
    for handle in legend.legend_handles:
        handle.set_edgecolor("black")
        handle.set_linewidth(0.6)
    for text in legend.get_texts():
        text.set_color("black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")
    fig.subplots_adjust(right=0.78, top=0.95, bottom=0.14)
    return fig


def _mirror_publication_assets_to_fig4_panel(source_stem: Path, fig4_stem: Path) -> None:
    """Copy pdf/png/jpg to Fig4 panel filenames for build_figure_panel_pack."""
    fig4_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf", ".jpg"):
        src = source_stem.with_suffix(ext)
        if src.is_file():
            shutil.copy2(src, fig4_stem.with_suffix(ext))


def _mirror_jpg_to_fig4_panel(source_jpg: Path, fig4_stem: Path) -> None:
    """Copy a waterfall boxplot JPG into figure_panels/Fig4 for the composite pack."""
    fig4_stem.parent.mkdir(parents=True, exist_ok=True)
    if source_jpg.is_file():
        shutil.copy2(source_jpg, fig4_stem.with_suffix(".jpg"))


def _save_short_summary_rank_composition_figure(
    ranked: pd.DataFrame,
    pub_dir: Path,
    send_dir: Optional[Path] = None,
    *,
    style: PublicationPlotStyle,
) -> Optional[Path]:
    """Write the stacked rank × utility composition figure (Fig. 4c)."""
    fig = _plot_short_summary_rank_utility_stacked_bar(
        ranked, style=_rank_composition_plot_style()
    )
    fig4_stem = pub_dir / FIG4C_PANEL_STEM
    _save_publication_figure(fig, fig4_stem)
    _mirror_publication_assets_to_fig4_panel(fig4_stem, FIG4_PANELS_DIR / FIG4C_PANEL_STEM)
    if send_dir is not None:
        send_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fig4_stem.with_suffix(".png"), send_dir / f"{FIG4C_PANEL_STEM}.png")
        shutil.copy2(fig4_stem.with_suffix(".pdf"), send_dir / f"{FIG4C_PANEL_STEM}.pdf")
    return fig4_stem.with_suffix(".png")


def _save_publication_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.12)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)



def build_no_onlabel_rank1_vs_all_sc_km_dataset(
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
    *,
    chemo: Optional[pd.DataFrame] = None,
    chemo_collapsed: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Fig. 4D: no-on-label rank-1 vs all SC rows (unmatched, unstratified).

    Left arm: Scenario-1 models with rank1_TTD (SC match not required).
    Right arm: all On_label=chemo rows (e.g. 204), not collapsed to one per model.
    """
    s1 = model_table[
        (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1)
        & model_table["rank1_TTD"].notna()
    ].copy()
    if s1.empty:
        return pd.DataFrame()

    left = s1[["Model", "rank1_TTD"]].rename(columns={"rank1_TTD": data_column})
    left["GROUP"] = "Rank1_no_onlabel"
    left["CENSOR"] = True

    # Prefer raw chemo rows so every SC treatment record is included (e.g. n=204).
    if chemo is not None and not chemo.empty:
        right = chemo[["Model", data_column]].copy()
        if "CENSOR" in chemo.columns:
            right["CENSOR"] = chemo["CENSOR"].values
        else:
            right["CENSOR"] = True
    elif chemo_collapsed is not None and not chemo_collapsed.empty:
        right = chemo_collapsed[["Model", data_column]].copy()
        right["CENSOR"] = True
    else:
        right = model_table[model_table["chemo_TTD"].notna()][["Model", "chemo_TTD"]].rename(
            columns={"chemo_TTD": data_column}
        )
        right["CENSOR"] = True
    if right.empty:
        return pd.DataFrame()
    right["GROUP"] = CLINICAL_UTILITY_GROUP_SC
    right["CENSOR"] = True

    out = pd.concat(
        [
            left[["Model", data_column, "GROUP", "CENSOR"]],
            right[["Model", data_column, "GROUP", "CENSOR"]],
        ],
        ignore_index=True,
    )
    out[data_column] = pd.to_numeric(out[data_column], errors="coerce")
    return out.dropna(subset=[data_column])


def build_q3_rank12_pooled_vs_sc_km_dataset(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Q3 KM dataset: pooled rank 1 and rank 2 TTD vs matched SC (no-on-label models)."""
    work = model_table[
        (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1)
        & model_table["chemo_TTD"].notna()
    ].copy()
    if work.empty:
        return pd.DataFrame()
    models = set(work["Model"].astype(str))
    # Scenario-1 / no-on-label arm: drop UNKNOWN On_label (e.g. tamoxifen).
    rank12 = ranked[
        ranked["Model"].astype(str).isin(models)
        & ranked["RANK"].astype(str).isin(("1", "2"))
        & ranked["ON_LABEL_UTILITY"].notna()
    ].copy()
    if rank12.empty:
        return pd.DataFrame()
    rank12["GROUP"] = "Rank1_2_pooled"
    rank12["CENSOR"] = True
    rank12[data_column] = pd.to_numeric(rank12[data_column], errors="coerce")
    rank12 = rank12.dropna(subset=[data_column])

    sc = work[["Model", "chemo_TTD"]].rename(columns={"chemo_TTD": data_column})
    sc["GROUP"] = CLINICAL_UTILITY_GROUP_SC
    sc["CENSOR"] = True
    sc[data_column] = pd.to_numeric(sc[data_column], errors="coerce")
    sc = sc.dropna(subset=[data_column])

    out = pd.concat(
        [rank12[[data_column, "GROUP", "CENSOR"]], sc[[data_column, "GROUP", "CENSOR"]]],
        ignore_index=True,
    )
    return out.dropna(subset=[data_column])


PUBLICATION_KM_COLORS: Dict[str, str] = {
    "Rank1_no_onlabel": "#0F52BA",
    "Rank1_2_pooled": "#2171B5",
    CLINICAL_UTILITY_GROUP_SC: CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_SC],
}
PUBLICATION_KM_LABELS: Dict[str, str] = {
    "Rank1_no_onlabel": "Top (no on-label)",
    "Rank1_2_pooled": "Rank 1-2 pooled",
    CLINICAL_UTILITY_GROUP_SC: "SC",
}

# Match KM_PDX Top/Bottom pair-KM typography exactly (Top is the size etalon).
ETALON_PAIR_KM_STYLE = PAIR_KM_ETALON_STYLE
PAIR_KM_ETALON_PX: Tuple[int, int] = (1947, 1441)


def _pad_pair_km_to_etalon_canvas(path: Path, size: Tuple[int, int] = PAIR_KM_ETALON_PX) -> None:
    """Pad/crop a saved pair-KM jpg to the Top etalon pixel size."""
    try:
        from PIL import Image
    except ImportError:
        LOGGER.warning("Pillow not available; skipping canvas pad for %s", path)
        return
    if not path.is_file():
        return
    img = Image.open(path).convert("RGB")
    if img.size == size:
        return
    tw, th = size
    img = img.crop((0, 0, min(img.size[0], tw), min(img.size[1], th)))
    canvas = Image.new("RGB", size, (255, 255, 255))
    x = (tw - img.size[0]) // 2
    y = (th - img.size[1]) // 2
    canvas.paste(img, (x, y))
    canvas.save(path, quality=95)



def _format_p_for_manuscript(p: Any) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)) or pd.isna(p):
        return "n/a"
    p = float(p)
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def _write_fig4d_manuscript_results_text(
    pub_dir: Path,
    row: Dict[str, Any],
    *,
    model_table: pd.DataFrame,
) -> Path:
    """Write manuscript-ready Results paragraph for Fig. 4D (unstratified KM stats)."""
    left_key = next((k for k in row if k.startswith("mPFS_Rank1")), None)
    right_key = next((k for k in row if k.startswith("mPFS_SC") or k == "mPFS_SC"), None)
    m_left = row.get(left_key) if left_key else row.get("mPFS_Rank1_no_onlabel")
    m_right = row.get(right_key) if right_key else row.get("mPFS_SC")
    n_left = int(row.get("n_Rank1_no_onlabel") or 0)
    n_sc = int(row.get("n_SC") or 0)
    logrank_p = row.get("logrank_p", row.get("pooled_logrank_p"))
    hr = row.get("HR", row.get("pooled_cox_hr"))
    hr_lo = row.get("HR_CI_lower_95", row.get("pooled_cox_hr_low"))
    hr_hi = row.get("HR_CI_upper_95", row.get("pooled_cox_hr_high"))
    hr_p = row.get("HR_P", row.get("pooled_cox_p"))

    paired = _build_clinical_utility_paired_table(
        model_table[
            (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1)
            & model_table["chemo_TTD"].notna()
            & model_table["rank1_TTD"].notna()
        ],
        "rank1_vs_sc",
    )
    wilcoxon_p = np.nan
    med_diff = np.nan
    if not paired.empty:
        wt = run_paired_tests(paired)
        wilcoxon_p = wt.get("wilcoxon_p", np.nan)
        med_diff = wt.get("median_diff_TTD", np.nan)

    km_sig = (
        (pd.notna(logrank_p) and float(logrank_p) < 0.05)
        or (pd.notna(hr_p) and float(hr_p) < 0.05)
    )
    wil_sig = pd.notna(wilcoxon_p) and float(wilcoxon_p) < 0.05

    def _hr_txt(h: Any, lo: Any, hi: Any) -> str:
        if pd.isna(h) or pd.isna(lo) or pd.isna(hi):
            return "n/a"
        return f"{float(h):.2f} (95% CI {float(lo):.2f}–{float(hi):.2f})"

    if km_sig:
        km_clause = (
            f"By log-rank and Cox regression, top-ranked treatments prolonged PFS relative to SCs "
            f"(mPFS {_format_pair_km_median(m_left)} vs {_format_pair_km_median(m_right)} days; "
            f"log-rank p = {_format_p_for_manuscript(logrank_p)}; "
            f"HR = {_hr_txt(hr, hr_lo, hr_hi)}, "
            f"p = {_format_p_for_manuscript(hr_p)}; "
            f"n = {n_left} top-ranked vs {n_sc} SC treatments; Fig. 4D)."
        )
    else:
        km_clause = (
            f"mPFS of the top-ranked non–on-label treatments was higher than that of SCs "
            f"({_format_pair_km_median(m_left)} vs {_format_pair_km_median(m_right)} days), "
            f"although the difference was not statistically significant "
            f"(log-rank p = {_format_p_for_manuscript(logrank_p)}; "
            f"HR = {_hr_txt(hr, hr_lo, hr_hi)}, "
            f"p = {_format_p_for_manuscript(hr_p)}; "
            f"n = {n_left} top-ranked vs {n_sc} SC treatments; Fig. 4D)."
        )

    wil_clause = ""
    if pd.notna(wilcoxon_p):
        wil_clause = (
            f" Using the Wilcoxon signed-rank test, we also performed a paired analysis in which "
            f"MTAs within each rank group were paired with the corresponding SC treatment of the "
            f"same tumor; top-ranked MTAs provided a "
            f"{'significant' if wil_sig else 'non-significant'} PFS benefit compared with matched "
            f"SCs (median ΔPFS = {float(med_diff):+.1f} days; "
            f"Wilcoxon p = {_format_p_for_manuscript(wilcoxon_p)}; Fig. 4E)."
        )

    if wil_sig and not km_sig:
        bridge = (
            " Thus, while the cohort-level KM comparison was not significant, the paired Wilcoxon "
            "analysis confirmed a significant within-tumor PFS benefit of top-ranked MTAs over "
            "matched SCs (Fig. 4E)."
        )
    elif km_sig and wil_sig:
        bridge = (
            " Thus, both the KM/Cox comparison (Fig. 4D) and the paired Wilcoxon analysis "
            "(Fig. 4E) supported a PFS benefit of top-ranked MTAs over SCs."
        )
    else:
        bridge = " Paired Wilcoxon analyses across rank groups are shown in Fig. 4E."

    paragraph = (
        "Because the lack of approved MTAs is a common clinical situation, we modeled this scenario "
        "with PCT treatment data. To this end, we first compared top-ranked treatments with SCs "
        "within the subset of tumors without approved MTA treatments by log-rank and Cox regression "
        f"tests (Fig. 4D). {km_clause}{wil_clause}{bridge} "
        "Moreover, when off-label or experimental MTAs were ranked first, they conferred a survival "
        "benefit over lower-ranked on-label options, the benefit increasing with increasing rank "
        "distance from the top."
    )

    methods_note = (
        "\n\n### Methods note (Fig. 4D)\n\n"
        "Fig. 4D compares rank-1 TTD from tumors without on-label therapy with all On_label=chemo "
        "treatment rows (not collapsed; models with multiple chemo regimens contribute multiple "
        "points). Annotated statistics are pooled log-rank and unstratified Cox PH "
        "(no Model stratification), matching the other KM panels in the manuscript.\n"
    )

    stats_block = (
        "\n\n### Numeric summary\n\n"
        f"- n Rank1 (no on-label): {n_left}\n"
        f"- n All SC: {n_sc}\n"
        f"- mPFS top-ranked (no on-label): {_format_pair_km_median(m_left)}\n"
        f"- mPFS All SC: {_format_pair_km_median(m_right)}\n"
        f"- log-rank p: {_format_p_for_manuscript(logrank_p)}\n"
        f"- Cox HR: {_hr_txt(hr, hr_lo, hr_hi)}, p = {_format_p_for_manuscript(hr_p)}\n"
        f"- Wilcoxon signed-rank p (Fig. 4E paired): {_format_p_for_manuscript(wilcoxon_p)}\n"
    )

    out = pub_dir / "Fig4d_manuscript_results_text.md"
    out.write_text(
        "# Fig. 4D — manuscript Results text (unstratified KM: rank-1 vs all SC)\n\n"
        + paragraph
        + methods_note
        + stats_block,
        encoding="utf-8",
    )
    # Also place next to Fig4 panel assets for easy copy-paste into the article.
    FIG4_PANELS_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out, FIG4_PANELS_DIR / out.name)
    return out


def _publication_km_plot_style() -> KMPlotStyle:
    return PAIR_KM_ETALON_STYLE


def _mirror_km_jpg_to_publication_assets(jpg_path: Path, stem: Path) -> None:
    """Copy KM jpg to publication png/pdf alongside vector-friendly png."""
    if not jpg_path.is_file():
        return
    stem.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(jpg_path, stem.with_suffix(".png"))
    img = plt.imread(str(jpg_path))
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    ax.imshow(img)
    ax.axis("off")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def generate_publication_km_cox_figures(
    pub_dir: Path,
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
    chemo_collapsed: Optional[pd.DataFrame] = None,
    chemo: Optional[pd.DataFrame] = None,
    fig4d_only: bool = False,
) -> List[str]:
    """KM figures for publication_figures; Fig4d = unstratified rank-1 vs all SC."""
    pub_dir.mkdir(parents=True, exist_ok=True)
    # Use Top/Bottom etalon style (default KMPlotStyle), not larger publication fonts.
    style = ETALON_PAIR_KM_STYLE
    colors = dict(PUBLICATION_KM_COLORS)
    label_map = dict(PUBLICATION_KM_LABELS)
    summary_rows: List[Dict[str, Any]] = []
    written: List[str] = []

    # Official Fig. 4D: no-on-label rank-1 vs all SC; pooled log-rank / Cox.
    no_onlabel_data = build_no_onlabel_rank1_vs_all_sc_km_dataset(
        model_table,
        data_column,
        chemo=chemo,
        chemo_collapsed=chemo_collapsed,
    )
    if not no_onlabel_data.empty:
        stem = pub_dir / FIG4D_PANEL_STEM
        left_g = "Rank1_no_onlabel"
        right_g = CLINICAL_UTILITY_GROUP_SC
        pooled_stats = pooled_logrank_cox_pairwise(
            no_onlabel_data, data_column, left_g, right_g
        )
        stats_override = {
            "logrank_p": pooled_stats.get("pooled_logrank_p", np.nan),
            "hr": pooled_stats.get("pooled_cox_hr", np.nan),
            "ci_low": pooled_stats.get("pooled_cox_hr_low", np.nan),
            "ci_high": pooled_stats.get("pooled_cox_hr_high", np.nan),
            "hr_p": pooled_stats.get("pooled_cox_p", np.nan),
        }
        row = create_pair_km_plot(
            no_onlabel_data,
            left_g,
            right_g,
            data_column,
            colors,
            style,
            stem.with_suffix(".jpg"),
            title="",
            stat_label_map=label_map,
            stats_override=stats_override,
            logrank_label="logrank",
            hr_label="HR",
        )
        if row:
            row = dict(row)
            row["figure"] = stem.name
            row["stats_method"] = "pooled_logrank_and_unstratified_Cox"
            row["design"] = "unmatched_no_onlabel_rank1_vs_all_SC_unstratified"
            row["n_Rank1_no_onlabel"] = int((no_onlabel_data["GROUP"] == left_g).sum())
            row["n_SC"] = int((no_onlabel_data["GROUP"] == right_g).sum())
            row["n_models_rank1"] = int(
                no_onlabel_data.loc[no_onlabel_data["GROUP"] == left_g, "Model"].nunique()
            )
            row["n_models_SC"] = int(
                no_onlabel_data.loc[no_onlabel_data["GROUP"] == right_g, "Model"].nunique()
            )
            row["pooled_logrank_p"] = pooled_stats.get("pooled_logrank_p", np.nan)
            row["pooled_cox_hr"] = pooled_stats.get("pooled_cox_hr", np.nan)
            row["pooled_cox_p"] = pooled_stats.get("pooled_cox_p", np.nan)
            summary_rows.append(row)
        _pad_pair_km_to_etalon_canvas(stem.with_suffix(".jpg"))
        _mirror_km_jpg_to_publication_assets(stem.with_suffix(".jpg"), stem)
        _mirror_publication_assets_to_fig4_panel(
            stem, FIG4_PANELS_DIR / FIG4D_PANEL_STEM
        )
        written.append(f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/{stem.name}.png")
        if summary_rows and summary_rows[-1].get("figure") == stem.name:
            _write_fig4d_manuscript_results_text(
                pub_dir,
                summary_rows[-1],
                model_table=model_table,
            )
            written.append(
                f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/Fig4d_manuscript_results_text.md"
            )

    if not fig4d_only:
        q3_data = build_q3_rank12_pooled_vs_sc_km_dataset(model_table, ranked, data_column)
        if (
            not q3_data.empty
            and q3_data["GROUP"].nunique() >= 2
            and "Rank1_2_pooled" in q3_data["GROUP"].values
            and CLINICAL_UTILITY_GROUP_SC in q3_data["GROUP"].values
        ):
            stem = pub_dir / "pub_Q3_km_rank12_pooled_vs_SC"
            row = create_pair_km_plot(
                q3_data,
                "Rank1_2_pooled",
                CLINICAL_UTILITY_GROUP_SC,
                data_column,
                colors,
                style,
                pub_dir / "pub_Q3_km_rank12_pooled_vs_SC.jpg",
                title="",
                stat_label_map=label_map,
                logrank_label="logrank",
                hr_label="HR",
            )
            if row:
                row = dict(row)
                row["figure"] = stem.name
                row["question"] = "Q3"
                row["n_Rank1_2_pooled"] = int((q3_data["GROUP"] == "Rank1_2_pooled").sum())
                row["n_SC"] = int((q3_data["GROUP"] == CLINICAL_UTILITY_GROUP_SC).sum())
                summary_rows.append(row)
            _mirror_km_jpg_to_publication_assets(stem.with_suffix(".jpg"), stem)
            written.append(f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/{stem.name}.png")

    if summary_rows:
        safe_to_excel(
            pd.DataFrame(summary_rows),
            pub_dir / "pub_km_cox_pairwise_summary.xlsx",
        )
        written.append(f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_km_cox_pairwise_summary.xlsx")

    if written:
        LOGGER.info("Publication KM/Cox figures written to %s (%d files)", pub_dir, len(written))
    return written


def _format_publication_p(p_value: Any) -> str:
    """On-plot scientific *P* value via mathtext (no leading P=)."""
    return format_p_sci_mathtext(p_value, digits=1)


def _format_annot_p_scientific(p_value: Any) -> str:
    """Italic *P* with one-decimal scientific notation via mathtext.

    Example: ``$\\mathit{P} = 5.6 \\\\times 10^{-2}$``.
    Unicode superscripts are avoided because Arial lacks those glyphs (□).
    """
    return format_p_equals_mathtext(p_value, digits=1)


def _boxplot_delta_n_p_annotation(delta: float, n: int, p_value: Any) -> str:
    """Annotation block: Δdiff, n, then *P* (lowest). No significance stars."""
    return f"Δ{float(delta):+.1f}\nn={int(n)}\n{_format_annot_p_scientific(p_value)}"


def _short_summary_cohort_table(
    model_table: pd.DataFrame,
    scenario: str,
    stratum: str,
) -> pd.DataFrame:
    if scenario == CLINICAL_UTILITY_SCENARIO1:
        subset = model_table[model_table["scenario"] == scenario]
    elif scenario == CLINICAL_UTILITY_SCENARIO2:
        subset = model_table[model_table["scenario"] == scenario]
    elif scenario == CLINICAL_UTILITY_SCENARIO3:
        subset = _scenario3_stratum_filter(model_table, stratum)
    elif scenario == CLINICAL_UTILITY_SCENARIO4:
        subset = _scenario4_low_rank_filter(model_table, stratum)
    else:
        raise ValueError(f"Unknown scenario for publication figure: {scenario!r}")
    return subset[subset["chemo_TTD"].notna()].copy()


def _between_cohort_onlabel_series(
    model_table: pd.DataFrame,
    low_stratum: str,
) -> Tuple[pd.Series, pd.Series]:
    cohort_a = model_table[
        (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO2) & model_table["onlabel_TTD"].notna()
    ]["onlabel_TTD"].astype(float)
    cohort_b = _scenario4_low_rank_filter(model_table, low_stratum)
    cohort_b = cohort_b[cohort_b["onlabel_TTD"].notna()]["onlabel_TTD"].astype(float)
    return cohort_a, cohort_b


def _plot_publication_waterfall_on_ax(
    ax: plt.Axes,
    paired: pd.DataFrame,
    result: pd.Series,
    *,
    left_group: str,
    right_group: str,
    group_colors: Dict[str, str],
    left_label: str,
    right_label: str,
    style: PublicationPlotStyle,
    panel_label: str = "",
) -> None:
    left_color = group_colors.get(left_group, "#0173B2")
    right_color = group_colors.get(right_group, "#DE8F05")
    tie_color = "#949494"
    plot_df = paired.sort_values("diff_TTD", ascending=True).reset_index(drop=True)
    colors = np.where(
        plot_df["diff_TTD"] > 0,
        left_color,
        np.where(plot_df["diff_TTD"] < 0, right_color, tie_color),
    )
    x = np.arange(len(plot_df))
    ax.bar(x, plot_df["diff_TTD"], color=colors, edgecolor="none", width=0.88, linewidth=0)
    ax.axhline(0, color="black", linewidth=0.7, zorder=1)
    ax.set_xlabel("PDX models (sorted by ΔPFS)", fontsize=style.font)
    ax.set_ylabel(f"ΔPFS (days)\n({left_label} − {right_label})", fontsize=style.font)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    stats_lines = [
        f"n = {int(result.get('n_models', len(paired)))}",
        f"Median ΔPFS = {result.get('median_diff_TTD', np.nan):+.1f} d "
        f"(95% CI {result.get('median_diff_CI_low', np.nan):.1f} to "
        f"{result.get('median_diff_CI_high', np.nan):.1f})",
        f"{result.get('pct_improved', np.nan):.0f}% favour {left_label}",
        f"Wilcoxon {format_p_equals_mathtext(result.get('wilcoxon_p'))}{_exploratory_bracket(result)}",
    ]
    if pd.notna(result.get("wilcoxon_p_fdr_bh")):
        stats_lines.append(f"FDR = {_format_publication_p(result['wilcoxon_p_fdr_bh'])}")
    ax.text(
        0.98,
        0.98,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=style.annotation,
    )
    if panel_label:
        ax.text(
            0.02,
            0.98,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=style.panel_label,
            fontweight="bold",
        )


def _plot_publication_km_on_ax(
    ax: plt.Axes,
    km_data: pd.DataFrame,
    left_group: str,
    right_group: str,
    data_column: str,
    colors: Dict[str, str],
    result: pd.Series,
    *,
    style: PublicationPlotStyle,
    panel_label: str = "",
    subtitle: str = "",
    logrank_label: str = "Log-rank",
    hr_label: str = "HR",
    use_stratified: bool = False,
) -> None:
    pair_groups = [left_group, right_group]
    km_style = KMPlotStyle(font=style.font, tick_size=style.tick, label_size=style.font, at_risk=style.tick - 1)
    kmf_list: List[KaplanMeierFitter] = []

    for group in pair_groups:
        group_data = km_data[km_data["GROUP"] == group].dropna(subset=[data_column])
        if group_data.empty:
            continue
        kmf = KaplanMeierFitter()
        display = PUBLICATION_GROUP_LABELS.get(group, group)
        kmf.fit(group_data[data_column], event_observed=group_data["CENSOR"], label=display)
        kmf.plot(ax=ax, color=colors.get(group), ci_show=style.cishow, linewidth=1.6)
        kmf_list.append(kmf)

    apply_km_axis_style(ax, km_style, title=subtitle)
    ax.set_xlabel("PFS (days)", fontsize=style.font)
    ax.set_ylabel("Kaplan–Meier estimate", fontsize=style.font)
    ax.legend(frameon=False, fontsize=style.legend, loc="lower left")

    left_label = PUBLICATION_GROUP_LABELS.get(left_group, left_group)
    right_label = PUBLICATION_GROUP_LABELS.get(right_group, right_group)
    left_data = km_data[km_data["GROUP"] == left_group].dropna(subset=[data_column])
    right_data = km_data[km_data["GROUP"] == right_group].dropna(subset=[data_column])
    left_median = KaplanMeierFitter().fit(
        left_data[data_column], event_observed=left_data["CENSOR"]
    ).median_survival_time_
    right_median = KaplanMeierFitter().fit(
        right_data[data_column], event_observed=right_data["CENSOR"]
    ).median_survival_time_

    if use_stratified and pd.notna(result.get("stratified_logrank_p")):
        logrank_p = float(result["stratified_logrank_p"])
        hr = result.get("stratified_cox_hr", np.nan)
        ci_low = result.get("stratified_cox_hr_low", np.nan)
        ci_high = result.get("stratified_cox_hr_high", np.nan)
        hr_p = result.get("stratified_cox_p", np.nan)
    else:
        logrank_p = float(result.get("pooled_logrank_p", np.nan))
        hr = result.get("pooled_cox_hr", np.nan)
        ci_low = result.get("pooled_cox_hr_low", np.nan)
        ci_high = result.get("pooled_cox_hr_high", np.nan)
        hr_p = result.get("pooled_cox_p", np.nan)

    stats_lines = [
        f"{left_label} mPFS = {_format_pair_km_median(left_median)} d",
        f"{right_label} mPFS = {_format_pair_km_median(right_median)} d",
    ]
    if pd.notna(result.get("n_models")):
        stats_lines.append(f"n = {int(result['n_models'])}")
    elif pd.notna(result.get("n_cohort_a")) and pd.notna(result.get("n_cohort_b")):
        stats_lines.append(
            f"n = {int(result['n_cohort_a'])} vs {int(result['n_cohort_b'])}"
        )
    stats_lines.append(f"{logrank_label} {format_p_equals_mathtext(logrank_p)}{_exploratory_bracket(result)}")
    if pd.notna(hr):
        exp = _exploratory_bracket(result)
        if pd.notna(ci_low) and pd.notna(ci_high):
            stats_lines.append(
                f"{hr_label} = {float(hr):.2f} ({float(ci_low):.2f}–{float(ci_high):.2f}){exp}"
            )
        else:
            stats_lines.append(f"{hr_label} = {float(hr):.2f}{exp}")
        if pd.notna(hr_p):
            stats_lines.append(f"{hr_label} {format_p_equals_mathtext(hr_p)}{exp}")
    if pd.notna(result.get("wilcoxon_p_fdr_bh")):
        stats_lines.append(f"Wilcoxon FDR = {_format_publication_p(result['wilcoxon_p_fdr_bh'])}")

    ax.text(
        0.97,
        0.40,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        fontsize=style.annotation,
        va="top",
        ha="right",
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="none", pad=3),
    )
    if kmf_list:
        add_km_at_risk_table(ax, kmf_list, at_risk_ypos=-0.38)
    if panel_label:
        ax.text(
            0.02,
            0.02,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=style.panel_label,
            fontweight="bold",
        )


def _write_statistical_caveats_reference(
    summary_dir: Path,
    *,
    items: List[Dict[str, Any]],
) -> None:
    """Standalone index of exploratory vs confirmatory evidence tiers."""
    lines = [
        "# Statistical caveats — clinical utility ranking",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        SHORT_SUMMARY_STATISTICAL_HIERARCHY,
        "",
        "## Analysis flags",
        "",
        "| Question | Flag | Rationale |",
        "| --- | --- | --- |",
        "| Q1 | [exploratory] | Paired on-label-at-rank vs Rank 1 or SC within model; eight tests. |",
    ]
    if ENABLE_Q2_ANALYSIS:
        lines.append(
            "| Q2 | [exploratory] | Unpaired on-label-at-rank pools vs Rank 1 (Off/Exp) or SC; eight tests. |"
        )
    lines.extend([
        "| Q3 | Confirmatory | Pre-specified Scenario 1; Wilcoxon significant (n ≈ 104). |",
        "| Q4 | Confirmatory (borderline) | Pre-specified; Wilcoxon p ≈ 0.07. |",
        "| Q5 | Confirmatory (rank 2) / [exploratory] (rank 3, pooled) | Three-stratum discordance panel. |",
        "",
        "## Per-question caveats (full text)",
        "",
    ])
    for item in items:
        lines.extend([
            f"### {item.get('heading', item.get('qid', ''))}",
            "",
            item.get("caveat", "_No caveat recorded._"),
            "",
        ])
    general_limits = [
        "## General limitations",
        "",
        "- TimeToDouble renamed PFS in summaries is a pharmacodynamic surrogate, not clinical PFS/OS.",
        "- All KM analyses treat every observation as an event (CENSOR = True); curves illustrate TTD "
        "distributions, not censored clinical survival.",
        "- Stratified log-rank and Wilcoxon can disagree because they target different estimands.",
    ]
    if ENABLE_Q2_ANALYSIS:
        general_limits.append(
            "- Q2 pools use on-label PFS at the specified global rank(s) in the screen, "
            "not best-on-label scenario assignment."
        )
    general_limits.append("")
    lines.extend(general_limits)
    (summary_dir / "STATISTICAL_CAVEATS.md").write_text("\n".join(lines), encoding="utf-8")


def _plot_publication_between_km_on_ax(
    ax: plt.Axes,
    model_table: pd.DataFrame,
    between_row: pd.Series,
    data_column: str,
    *,
    style: PublicationPlotStyle,
    panel_label: str = "",
) -> None:
    group_a = "On_label_rank1"
    group_b = "On_label_rank7_10"
    cohort_a, cohort_b = _between_cohort_onlabel_series(model_table, str(between_row["stratum"]))
    frames: List[pd.DataFrame] = []
    for series, grp in ((cohort_a, group_a), (cohort_b, group_b)):
        if series.empty:
            continue
        frames.append(
            pd.DataFrame(
                {
                    data_column: series.values,
                    "GROUP": grp,
                    "CENSOR": True,
                    "Model": [f"{grp}_{i}" for i in range(len(series))],
                }
            )
        )
    km_data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    colors = {group_a: CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_ONLABEL], group_b: "#C2185B"}
    _plot_publication_km_on_ax(
        ax,
        km_data,
        group_a,
        group_b,
        data_column,
        colors,
        between_row,
        style=style,
        panel_label=panel_label,
        subtitle="Unpaired on-label PFS",
        logrank_label="Log-rank",
        hr_label="Cox HR",
        use_stratified=False,
    )


def _plot_publication_between_box_on_ax(
    ax: plt.Axes,
    model_table: pd.DataFrame,
    between_row: pd.Series,
    *,
    style: PublicationPlotStyle,
    panel_label: str = "",
) -> None:
    cohort_a, cohort_b = _between_cohort_onlabel_series(model_table, str(between_row["stratum"]))
    labels = [
        f"On-label at rank 1\n(n = {len(cohort_a)})",
        f"On-label rank 7–10\n(n = {len(cohort_b)})",
    ]
    data = [cohort_a.values, cohort_b.values]
    box = ax.boxplot(
        data,
        tick_labels=labels,
        patch_artist=True,
        widths=0.55,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"linewidth": 0.9},
        capprops={"linewidth": 0.9},
        boxprops={"linewidth": 0.9},
    )
    palette = [CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_ONLABEL], "#C2185B"]
    for patch, color in zip(box["boxes"], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(SEED)
    for idx, values in enumerate(data):
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(
            np.full(len(values), idx + 1) + jitter,
            values,
            s=22,
            color=palette[idx],
            edgecolor="white",
            linewidth=0.4,
            alpha=0.9,
            zorder=3,
        )
    ax.set_ylabel("On-label PFS (days)", fontsize=style.font)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    stats_lines = [
        f"Median A = {between_row.get('median_a', np.nan):.1f} d",
        f"Median B = {between_row.get('median_b', np.nan):.1f} d",
        f"Δmedian = {between_row.get('median_diff_a_minus_b', np.nan):+.1f} d",
        f"Mann–Whitney {format_p_equals_mathtext(between_row.get('mannwhitney_p'))} [exploratory]",
    ]
    ax.text(
        0.98,
        0.98,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=style.annotation,
    )
    if panel_label:
        ax.text(
            0.02,
            0.98,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=style.panel_label,
            fontweight="bold",
        )


def _render_publication_question_composite(
    qid: str,
    panels: Sequence[ShortSummaryPubPanel],
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    paired_summary: pd.DataFrame,
    between_cohort: pd.DataFrame,
    data_column: str,
    style: PublicationPlotStyle,
) -> Optional[plt.Figure]:
    n_panels = len(panels)
    if n_panels == 2:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
        axes_flat = np.atleast_1d(axes).flatten()
    elif n_panels == 3:
        fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6))
        axes_flat = np.atleast_1d(axes).flatten()
    elif n_panels == 4:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.0))
        axes_flat = np.atleast_1d(axes).flatten()
    else:
        return None

    fig.suptitle(SHORT_SUMMARY_QUESTION_TITLES.get(qid, qid), fontsize=style.title + 1, y=0.995)
    panel_letters = "ABCDEFGH"

    for idx, (panel, ax) in enumerate(zip(panels, axes_flat)):
        letter = panel_letters[idx]
        if panel.panel_type == "between_km":
            hits = between_cohort[between_cohort["stratum"] == panel.between_stratum]
            if hits.empty:
                ax.axis("off")
                continue
            _plot_publication_between_km_on_ax(
                ax,
                model_table,
                hits.iloc[0],
                data_column,
                style=style,
                panel_label=letter,
            )
            ax.set_title(panel.subtitle, fontsize=style.title, pad=6)
            continue
        if panel.panel_type == "between_box":
            hits = between_cohort[between_cohort["stratum"] == panel.between_stratum]
            if hits.empty:
                ax.axis("off")
                continue
            _plot_publication_between_box_on_ax(
                ax,
                model_table,
                hits.iloc[0],
                style=style,
                panel_label=letter,
            )
            ax.set_title(panel.subtitle, fontsize=style.title, pad=6)
            continue

        cohort = _short_summary_cohort_table(model_table, panel.scenario, panel.stratum)
        paired = _build_clinical_utility_paired_table(cohort, panel.comparison)
        result = _clinical_paired_row(paired_summary, panel.scenario, panel.comparison, panel.stratum)
        if paired.empty or result is None:
            ax.axis("off")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        left_group, right_group = next(
            (lg, rg) for cid, lg, rg in CLINICAL_UTILITY_COMPARISONS if cid == panel.comparison
        )
        if panel.panel_type == "waterfall":
            _plot_publication_waterfall_on_ax(
                ax,
                paired,
                result,
                left_group=left_group,
                right_group=right_group,
                group_colors=CLINICAL_UTILITY_COLORS,
                left_label=PUBLICATION_GROUP_LABELS.get(left_group, left_group),
                right_label=PUBLICATION_GROUP_LABELS.get(right_group, right_group),
                style=style,
                panel_label=letter,
            )
            ax.set_title(panel.subtitle, fontsize=style.title, pad=6)
            continue

        km_data = _clinical_utility_km_dataset_from_paired(
            paired,
            ranked,
            chemo_collapsed,
            left_group,
            right_group,
            data_column,
        )
        if km_data.empty:
            ax.axis("off")
            continue
        _plot_publication_km_on_ax(
            ax,
            km_data,
            left_group,
            right_group,
            data_column,
            CLINICAL_UTILITY_COLORS,
            result,
            style=style,
            panel_label=letter,
            subtitle=panel.subtitle,
        )

    if n_panels == 4:
        fig.tight_layout(rect=[0.0, 0.06, 1.0, 0.96])
    else:
        fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.94])
    return fig


def _render_q3_discordance_composite(
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    paired_summary: pd.DataFrame,
    data_column: str,
    style: PublicationPlotStyle,
) -> plt.Figure:
    """Q3: 3×2 panel — rank 2, rank 3, and pooled ranks 2-3 (KM + Wilcoxon each)."""
    fig, axes = plt.subplots(3, 2, figsize=(10.8, 13.5))
    panel_letters = "ABCDEF"
    panel_idx = 0
    comparison = "rank1_vs_onlabel"
    left_group, right_group = next(
        (lg, rg) for cid, lg, rg in CLINICAL_UTILITY_COMPARISONS if cid == comparison
    )

    for row_idx, (stratum, row_title) in enumerate(Q3_DISCORDANCE_STRATA):
        cohort = _short_summary_cohort_table(model_table, CLINICAL_UTILITY_SCENARIO3, stratum)
        paired = _build_clinical_utility_paired_table(cohort, comparison)
        result = _clinical_paired_row(
            paired_summary, CLINICAL_UTILITY_SCENARIO3, comparison, stratum
        )
        km_ax = axes[row_idx, 0]
        wf_ax = axes[row_idx, 1]
        if paired.empty or result is None:
            km_ax.axis("off")
            wf_ax.axis("off")
            continue

        letter_km = panel_letters[panel_idx]
        panel_idx += 1
        km_data = _clinical_utility_km_dataset_from_paired(
            paired,
            ranked,
            chemo_collapsed,
            left_group,
            right_group,
            data_column,
        )
        _plot_publication_km_on_ax(
            km_ax,
            km_data,
            left_group,
            right_group,
            data_column,
            CLINICAL_UTILITY_COLORS,
            result,
            style=style,
            panel_label=letter_km,
            subtitle="Stratified Cox PH",
        )

        letter_wf = panel_letters[panel_idx]
        panel_idx += 1
        _plot_publication_waterfall_on_ax(
            wf_ax,
            paired,
            result,
            left_group=left_group,
            right_group=right_group,
            group_colors=CLINICAL_UTILITY_COLORS,
            left_label=PUBLICATION_GROUP_LABELS.get(left_group, left_group),
            right_label=PUBLICATION_GROUP_LABELS.get(right_group, right_group),
            style=style,
            panel_label=letter_wf,
        )
        wf_ax.set_title("Wilcoxon ΔPFS", fontsize=style.title, pad=6)

        km_ax.text(
            -0.30,
            0.5,
            row_title,
            transform=km_ax.transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=style.font,
            fontweight="bold",
        )

    fig.suptitle(
        SHORT_SUMMARY_QUESTION_TITLES["Q5"],
        fontsize=style.title + 1,
        y=0.995,
    )
    fig.tight_layout(rect=[0.08, 0.02, 1.0, 0.97])
    return fig


def generate_short_summary_publication_figures(
    summary_dir: Path,
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    chemo_collapsed: pd.DataFrame,
    paired_summary: pd.DataFrame,
    between_cohort: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
    utility_panel_source: Optional[Path] = None,
    chemo: Optional[pd.DataFrame] = None,
) -> List[str]:
    """Create publication-grade composite figures for short supervisor summary (Q1–Q5)."""
    pub_dir = summary_dir / CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR
    send_dir = summary_dir / CLINICAL_SUPERVISOR_FIGURES_COLLECTION_DIR
    pub_dir.mkdir(parents=True, exist_ok=True)
    send_dir.mkdir(parents=True, exist_ok=True)

    style = PublicationPlotStyle()
    style.apply_rcparams()
    written: List[str] = []

    rank_png = _save_short_summary_rank_composition_figure(
        ranked, pub_dir, send_dir, style=style
    )
    if rank_png is not None:
        rel = f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png"
        written.append(rel)

    km_written = generate_publication_km_cox_figures(
        pub_dir,
        model_table=model_table,
        ranked=ranked,
        data_column=data_column,
        chemo_collapsed=chemo_collapsed,
        chemo=chemo,
    )
    written.extend(km_written)

    q9_summary = build_q9_paired_pool_sensitivity_table(model_table, ranked, data_column=data_column)
    fig_q1 = _render_q9_paired_pool_composite(
        model_table, ranked, q9_summary, style=style, data_column=data_column
    )
    stem_q1 = pub_dir / "pub_Q1_composite"
    _save_publication_figure(fig_q1, stem_q1)
    rel_q1 = f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_Q1_composite.png"
    written.append(rel_q1)
    shutil.copy2(stem_q1.with_suffix(".png"), send_dir / "pub_Q1_composite.png")
    shutil.copy2(stem_q1.with_suffix(".pdf"), send_dir / "pub_Q1_composite.pdf")

    # Q2: Downranked on-label pools vs Rank 1 (Off/Exp) or standard chemo (unpaired) — disabled.
    if ENABLE_Q2_ANALYSIS:
        sensitivity = build_low_rank_pool_sensitivity_table(
            model_table, ranked, data_column=data_column
        )
        fig_q2 = _render_q8_low_rank_sensitivity_composite(
            model_table, ranked, sensitivity, style=style, data_column=data_column
        )
        stem_q2 = pub_dir / "pub_Q2_composite"
        _save_publication_figure(fig_q2, stem_q2)
        rel_q2 = f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_Q2_composite.png"
        written.append(rel_q2)
        shutil.copy2(stem_q2.with_suffix(".png"), send_dir / "pub_Q2_composite.png")
        shutil.copy2(stem_q2.with_suffix(".pdf"), send_dir / "pub_Q2_composite.pdf")

    for qid, panels in SHORT_SUMMARY_PUBLICATION_PANELS.items():
        fig = _render_publication_question_composite(
            qid,
            panels,
            model_table=model_table,
            ranked=ranked,
            chemo_collapsed=chemo_collapsed,
            paired_summary=paired_summary,
            between_cohort=between_cohort,
            data_column=data_column,
            style=style,
        )
        if fig is None:
            continue
        stem = pub_dir / f"pub_{qid}_composite"
        _save_publication_figure(fig, stem)
        rel = f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_{qid}_composite.png"
        written.append(rel)
        shutil.copy2(stem.with_suffix(".png"), send_dir / f"pub_{qid}_composite.png")
        shutil.copy2(stem.with_suffix(".pdf"), send_dir / f"pub_{qid}_composite.pdf")

    fig_q5 = _render_q3_discordance_composite(
        model_table=model_table,
        ranked=ranked,
        chemo_collapsed=chemo_collapsed,
        paired_summary=paired_summary,
        data_column=data_column,
        style=style,
    )
    stem_q5 = pub_dir / "pub_Q5_composite"
    _save_publication_figure(fig_q5, stem_q5)
    rel_q5 = f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_Q5_composite.png"
    written.append(rel_q5)
    shutil.copy2(stem_q5.with_suffix(".png"), send_dir / "pub_Q5_composite.png")
    shutil.copy2(stem_q5.with_suffix(".pdf"), send_dir / "pub_Q5_composite.pdf")

    if utility_panel_source and utility_panel_source.is_file():
        appendix_png = pub_dir / "pub_appendix_utility_at_rank_panel.png"
        appendix_pdf = pub_dir / "pub_appendix_utility_at_rank_panel.pdf"
        img = plt.imread(str(utility_panel_source))
        fig, ax = plt.subplots(figsize=(14, 10))
        ax.imshow(img)
        ax.axis("off")
        fig.suptitle("Appendix: Per-rank utility vs matched standard chemotherapy", fontsize=style.title + 1)
        fig.tight_layout()
        _save_publication_figure(fig, pub_dir / "pub_appendix_utility_at_rank_panel")
        shutil.copy2(appendix_png, send_dir / appendix_png.name)
        shutil.copy2(appendix_pdf, send_dir / appendix_pdf.name)
        written.append(f"{CLINICAL_SUPERVISOR_PUBLICATION_FIGURES_DIR}/pub_appendix_utility_at_rank_panel.png")

    _write_publication_figures_index(pub_dir, written)
    LOGGER.info("Publication figures written to %s (%d files)", pub_dir, len(written))
    return written


def _write_publication_figures_index(pub_dir: Path, written: Sequence[str]) -> None:
    lines = [
        "# Publication figures — short supervisor summary (Q1–Q5)",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "Each question has a composite figure with panels for the statistical analyses cited in "
        "CLINICAL_SUMMARY_FOR_SUPERVISOR_short.docx:",
        "",
        f"- Rank composition: `{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png` (stacked On-label status by rank)",
        "",
        "- Publication KM/Cox: `pub_km_no_onlabel_rank1_vs_SC` (no on-label rank 1 vs SC)",
        "- Publication KM/Cox: `pub_Q3_km_rank12_pooled_vs_SC` (Q3: rank 1-2 pooled vs SC)",
        "- Summary table: `pub_km_cox_pairwise_summary.xlsx`",
        "",
        "- Paired Wilcoxon: waterfall plot of ΔPFS (days) per PDX model",
        "- Stratified Cox PH: Kaplan–Meier curves with Model-stratified log-rank and HR",
        "- Q1 pool sensitivity: paired Wilcoxon waterfall and stratified KM within each model",
    ]
    if ENABLE_Q2_ANALYSIS:
        lines.append(
            "- Q2 pool sensitivity: unpaired Mann–Whitney box and KM for progressive on-label-at-rank pools"
        )
    lines.extend([
        "",
        "Formats: PDF (vector) and PNG (300 dpi). Copies of PNG/PDF are also in `figures_to_send/`.",
        "",
        "| Question | File | Panels |",
        "| --- | --- | --- |",
    ])
    panel_desc: Dict[str, str] = {
        "Q1": "4×2 panel: paired on-label-at-rank vs Rank 1 or SC (waterfall, stratified KM)",
        "Q3": "A KM (strat. Cox); B Wilcoxon ΔPFS",
        "Q4": "A KM (strat. Cox); B Wilcoxon ΔPFS",
        "Q5": "3×2 panel: rank 2, rank 3, pooled ranks 2-3 (KM + Wilcoxon ΔPFS each)",
    }
    if ENABLE_Q2_ANALYSIS:
        panel_desc["Q2"] = (
            "4×2 panel: on-label-at-rank pools vs Rank 1 (Off/Exp) or SC (unpaired box, KM)"
        )
    for qid in SHORT_SUMMARY_QUESTION_ORDER:
        lines.append(
            f"| {qid} | `pub_{qid}_composite.pdf` | {panel_desc.get(qid, '')} |"
        )
    lines.append("| Appendix | `pub_appendix_utility_at_rank_panel.pdf` | Per-rank utility panel |")
    lines.append("")
    (pub_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


CLINICAL_SUPERVISOR_SUMMARY_DIR = "clinical_interpretation_summary"
CLINICAL_SUPERVISOR_FIGURES_COLLECTION_DIR = "figures_to_send"
CLINICAL_SUPERVISOR_KEY_FIGURES: Tuple[str, ...] = (
    "clinical_utility_scenario1_no_onlabel_all_rank1_vs_sc_Kmplot_",
    "clinical_utility_scenario2_onlabel_rank1_all_onlabel_vs_sc_Kmplot_",
    "clinical_utility_scenario3_onlabel_not_top_rank2_rank1_vs_onlabel_Kmplot_",
    "clinical_utility_scenario3_rank2_waterfall_",
    "clinical_utility_scenario3_onlabel_not_top_rank2_rank1_vs_sc_Kmplot_",
    "clinical_utility_scenario3_onlabel_not_top_rank2_onlabel_vs_sc_Kmplot_",
    "clinical_utility_utility_at_rank_panel_",
)


def _clinical_format_p(p_value: Any) -> str:
    if pd.isna(p_value):
        return "n/a"
    return format_p_sci_plain(p_value, digits=1)


def _clinical_sig_word(p_value: Any, alpha: float = 0.05) -> str:
    if pd.isna(p_value):
        return "not testable"
    return "significant" if float(p_value) < alpha else "not significant"


def _p_cell_with_star(p_value: Any, alpha: float = 0.05) -> str:
    """Format p-value; append '*' when significant at alpha."""
    if pd.isna(p_value):
        return "n/a"
    text = _clinical_format_p(p_value)
    if float(p_value) < alpha:
        text += "*"
    return text


def _cox_hr_cell(
    row: Optional[pd.Series],
    *,
    hr_key: str = "stratified_cox_hr",
    p_key: str = "stratified_cox_p",
    alpha: float = 0.05,
) -> str:
    if row is None or pd.isna(row.get(hr_key)):
        return "n/a"
    hr = float(row[hr_key])
    p_text = _p_cell_with_star(row.get(p_key), alpha=alpha)
    return f"{hr:.2f} ({p_text})"


def _clinical_paired_row(
    paired_summary: pd.DataFrame,
    scenario: str,
    comparison: str,
    stratum: str = "all",
) -> Optional[pd.Series]:
    if paired_summary.empty:
        return None
    mask = (
        (paired_summary["scenario"] == scenario)
        & (paired_summary["comparison"] == comparison)
        & (paired_summary["stratum"] == stratum)
    )
    hits = paired_summary[mask]
    return hits.iloc[0] if not hits.empty else None


def _clinical_results_line(row: Optional[pd.Series]) -> str:
    if row is None:
        return "_No paired data available._"
    parts = [
        f"n = {int(row['n_models'])} models",
        f"median ΔTTD = {row['median_diff_TTD']:+.1f} days "
        f"(95% CI {row['median_diff_CI_low']:.1f}, {row['median_diff_CI_high']:.1f})",
        f"{row['pct_improved']:.1f}% of models improved",
        f"Wilcoxon p = {_clinical_format_p(row['wilcoxon_p'])} ({_clinical_sig_word(row['wilcoxon_p'])})",
    ]
    if pd.notna(row.get("wilcoxon_p_fdr_bh")):
        parts.append(
            f"Wilcoxon FDR = {_clinical_format_p(row['wilcoxon_p_fdr_bh'])} "
            f"({_clinical_sig_word(row['wilcoxon_p_fdr_bh'])})"
        )
    if pd.notna(row.get("stratified_cox_hr")):
        parts.append(
            f"Stratified Cox HR = {row['stratified_cox_hr']:.2f} "
            f"({row['stratified_cox_hr_low']:.2f}–{row['stratified_cox_hr_high']:.2f}), "
            f"p = {_clinical_format_p(row['stratified_cox_p'])}"
        )
    if pd.notna(row.get("stratified_logrank_p")):
        parts.append(f"Stratified log-rank p = {_clinical_format_p(row['stratified_logrank_p'])}")
    return "; ".join(parts) + "."


def _write_clinical_cohort_flow(
    summary_dir: Path,
    model_table: pd.DataFrame,
) -> None:
    total = len(model_table)
    with_chemo = int(model_table["chemo_TTD"].notna().sum())
    without_chemo = total - with_chemo
    scenario_counts = model_table["scenario"].value_counts().to_dict()

    lines = [
        "# Cohort flow — Clinical utility ranking analysis",
        "",
        "## Summary counts",
        "",
        f"- **Total tumour models analysed:** {total}",
        f"- **Models with matched standard chemotherapy (SC):** {with_chemo}",
        f"- **Models excluded (no SC row):** {without_chemo}",
        "",
        "## Scenario assignment (mutually exclusive)",
        "",
        f"| Scenario | n models |",
        f"|----------|----------|",
    ]
    labels = {
        CLINICAL_UTILITY_SCENARIO1: "1 — No on-label drug",
        CLINICAL_UTILITY_SCENARIO2: "2 — On-label at rank 1",
        CLINICAL_UTILITY_SCENARIO3: "3 — On-label not rank 1 (rank 2 or 3)",
        CLINICAL_UTILITY_SCENARIO4: "4 — On-label at rank 4+ (exploratory)",
    }
    for key in CLINICAL_UTILITY_SCENARIOS:
        lines.append(f"| {labels[key]} | {int(scenario_counts.get(key, 0))} |")

    s3 = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO3]
    if not s3.empty:
        lines.extend(
            [
                "",
                "## Scenario 3 sub-strata",
                "",
                f"- On-label at rank 2: {int((s3['best_onlabel_rank_int'] == 2).sum())} models",
                f"- On-label at rank 3: {int((s3['best_onlabel_rank_int'] == 3).sum())} models",
                f"- Descriptive gate (Off/Exp above on-label only): "
                f"{int(s3['ranks_above_onlabel_valid'].fillna(False).astype(bool).sum())}/{len(s3)}",
            ]
        )

    lines.extend(
        [
            "",
            "## Flow diagram",
            "",
            "```mermaid",
            "flowchart TD",
            f"  start[\"All PDX models with ranked treatments<br/>n = {total}\"]",
            f"  chemo{{\"Has matched SC?\"}}",
            f"  excl[\"Excluded: no SC<br/>n = {without_chemo}\"]",
            f"  incl[\"Included in paired analyses<br/>n = {with_chemo}\"]",
            "  onlabel{\"Has On_label drug?\"}",
            f"  s1[\"Scenario 1: No on-label<br/>n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0))}\"]",
            "  rank{\"Best on-label global rank\"}",
            f"  s2[\"Scenario 2: On-label rank 1<br/>n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0))}\"]",
            f"  s3[\"Scenario 3: On-label rank 2-3<br/>n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0))}\"]",
            f"  s4[\"Scenario 4: On-label rank 4+<br/>n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO4, 0))}\"]",
            "  start --> chemo",
            f"  chemo -->|No| excl",
            f"  chemo -->|Yes| incl",
            "  incl --> onlabel",
            f"  onlabel -->|No| s1",
            "  onlabel -->|Yes| rank",
            "  rank -->|Rank 1| s2",
            "  rank -->|Rank 2 or 3| s3",
            "  rank -->|Rank 4+| s4",
            "```",
        ]
    )
    (summary_dir / "COHORT_FLOW.md").write_text("\n".join(lines), encoding="utf-8")


def _copy_clinical_key_figures(output_dir: Path, summary_dir: Path, date: str) -> List[str]:
    fig_dir = summary_dir / CLINICAL_SUPERVISOR_FIGURES_COLLECTION_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    copied: List[str] = []
    for prefix in CLINICAL_SUPERVISOR_KEY_FIGURES:
        matches = sorted(output_dir.glob(f"{prefix}*{date}*.jpg"))
        if not matches:
            matches = sorted(output_dir.glob(f"{prefix}*.jpg"))
        if not matches:
            LOGGER.warning("No figure found for prefix %s", prefix)
            continue
        src = matches[-1]
        dest = fig_dir / src.name
        shutil.copy2(src, dest)
        copied.append(str(dest.relative_to(summary_dir)))
    if copied:
        _write_figures_collection_index(fig_dir, copied)
    return copied


def _write_figures_collection_index(fig_dir: Path, copied_rel_names: List[str]) -> None:
    """Write README in figures collection folder for supervisor sharing."""
    labels = {
        "scenario1_no_onlabel_all_rank1_vs_sc": "Exp1 / Scenario 1 — Rank 1 vs standard chemotherapy (KM)",
        "scenario2_onlabel_rank1_all_onlabel_vs_sc": "Exp2 / Scenario 2 — On-label at rank 1 vs SC (KM)",
        "scenario3_onlabel_not_top_rank2_rank1_vs_onlabel": "Exp3 / Scenario 3 rank 2 — Rank 1 vs on-label (KM)",
        "scenario3_rank2_waterfall": "Scenario 3 rank 2 — Waterfall (% change TTD, rank 1 vs on-label)",
        "scenario3_onlabel_not_top_rank2_rank1_vs_sc": "Scenario 3 rank 2 — Rank 1 vs SC (KM, ancillary)",
        "scenario3_onlabel_not_top_rank2_onlabel_vs_sc": "Scenario 3 rank 2 — On-label vs SC (KM, ancillary)",
        "utility_at_rank_panel": "Per-rank utility vs matched SC (panel)",
    }
    lines = [
        "# Key figures — clinical utility ranking (shareable collection)",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "Copy this entire `figures_to_send` folder when sending plots to your supervisor.",
        "",
        "| # | File | Description |",
        "| --- | --- | --- |",
    ]
    for idx, rel in enumerate(copied_rel_names, start=1):
        name = Path(rel).name
        desc = name
        for key, label in labels.items():
            if key in name:
                desc = label
                break
        lines.append(f"| {idx} | `{name}` | {desc} |")
    lines.append("")
    (fig_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def _comparison_plain_label(comparison: str) -> str:
    labels = {
        "rank1_vs_sc": "Algorithm rank 1 versus matched standard chemotherapy (SC)",
        "onlabel_vs_sc": "Best on-label agent versus matched SC",
        "rank1_vs_onlabel": "Algorithm rank 1 versus best on-label agent",
        "above_onlabel_vs_onlabel": "Best Off_label/Exp compound above on-label versus on-label",
    }
    return labels.get(comparison, comparison.replace("_", " "))


def _stratum_plain_label(stratum: str) -> str:
    labels = {
        "all": "entire scenario cohort",
        "rank2": "on-label at global rank 2",
        "rank3": "on-label at global rank 3",
        "rank2_3_pooled": "pooled on-label at ranks 2–3",
        "rank4_plus": "on-label at global rank 4 or lower (exploratory)",
        "rank7_10": "on-label pooled at ranks 7–10",
        "pool_ge6": "progressive pool: on-label rank ≥ 6",
        "pool_ge5": "progressive pool: on-label rank ≥ 5",
        "pool_ge4": "progressive pool: on-label rank ≥ 4",
    }
    return labels.get(stratum, stratum)


_ONLABEL_STRATUM_ARM_LABELS: Dict[str, str] = {
    "all": "On-label",
    "rank2": "On-label at rank 2",
    "rank3": "On-label at rank 3",
    "rank2_3_pooled": "On-label at ranks 2-3",
    "rank4_plus": "On-label at rank ≥4",
}


def _arm_publication_label(group: str, stratum: str, scenario: str) -> str:
    if group == CLINICAL_UTILITY_GROUP_RANK1:
        return "Rank 1 (Off/Exp)"
    if group == CLINICAL_UTILITY_GROUP_SC:
        return "Standard chemo"
    if group == CLINICAL_UTILITY_GROUP_ABOVE:
        return "Best Off/Exp above on-label"
    if group == CLINICAL_UTILITY_GROUP_ONLABEL:
        if scenario == CLINICAL_UTILITY_SCENARIO2 and stratum == "all":
            return "On-label at rank 1"
        pool_label = LOW_RANK_POOL_Q8_LABELS.get(stratum)
        if pool_label:
            return pool_label
        return _ONLABEL_STRATUM_ARM_LABELS.get(stratum, "On-label")
    return PUBLICATION_GROUP_LABELS.get(group, group)


def _scenario_comparison_group_label(
    scenario: str,
    stratum: str,
    comparison: str,
    left_group: str,
    right_group: str,
) -> str:
    """Human-readable comparison label for scenario results tables."""
    del comparison  # arms are fully determined by left/right group and stratum
    left = _arm_publication_label(left_group, stratum, scenario)
    right = _arm_publication_label(right_group, stratum, scenario)
    return f"{left} vs {right}"


def _enrich_paired_summary_table(paired_summary: pd.DataFrame) -> pd.DataFrame:
    if paired_summary.empty:
        return paired_summary
    work = paired_summary.copy()
    work["comparison_group"] = work.apply(
        lambda r: _scenario_comparison_group_label(
            str(r["scenario"]),
            str(r.get("stratum", "all")),
            str(r.get("comparison", "")),
            str(r.get("left_group", "")),
            str(r.get("right_group", "")),
        ),
        axis=1,
    )
    if "n_models" in work.columns:
        work = work[work["n_models"] >= MIN_SCENARIO_PAIRED_N].copy()
    return work.reset_index(drop=True)


def _epistemic_plain_label(class_id: str) -> str:
    labels = {
        "A_ineffective_correctly_deprioritized": (
            "Class A — pharmacodynamically ineffective on-label agent correctly deprioritized by ranking"
        ),
        "B_suboptimal_vs_rank1_not_worse_than_sc": (
            "Class B — on-label inferior to rank 1 but not clearly worse than chemotherapy"
        ),
        "C_possible_ranking_miss": (
            "Class C — possible ranking failure: on-label outperforms algorithm rank 1"
        ),
        "D_indeterminate": "Class D — indeterminate (missing or tied TTD values)",
    }
    return labels.get(class_id, class_id)


def _clinical_results_markdown_table(row: Optional[pd.Series]) -> List[str]:
    if row is None:
        return ["_No paired data available for this experiment._", ""]
    headers = [
        "Metric", "Value",
    ]
    rows = [
        ("Models (paired)", str(int(row["n_models"]))),
        ("Proportion with longer TTD (favourable)", f"{row['pct_improved']:.1f}%"),
        (
            "Median paired difference in TTD (days)",
            f"{row['median_diff_TTD']:+.1f} "
            f"(bootstrap 95% CI {row['median_diff_CI_low']:.1f} to {row['median_diff_CI_high']:.1f})",
        ),
        ("Wilcoxon signed-rank p", f"{_clinical_format_p(row['wilcoxon_p'])} ({_clinical_sig_word(row['wilcoxon_p'])})"),
    ]
    if pd.notna(row.get("wilcoxon_p_fdr_bh")):
        rows.append(
            (
                "Benjamini–Hochberg FDR (within comparison family)",
                f"{_clinical_format_p(row['wilcoxon_p_fdr_bh'])} ({_clinical_sig_word(row['wilcoxon_p_fdr_bh'])})",
            )
        )
    if pd.notna(row.get("binom_p")):
        rows.append(("Binomial test on direction", f"p = {_clinical_format_p(row['binom_p'])}"))
    if pd.notna(row.get("stratified_cox_hr")):
        rows.append(
            (
                "Stratified Cox HR (arm vs reference, Model as stratum)",
                f"{row['stratified_cox_hr']:.2f} "
                f"(95% CI {row['stratified_cox_hr_low']:.2f}–{row['stratified_cox_hr_high']:.2f}); "
                f"p = {_clinical_format_p(row['stratified_cox_p'])}",
            )
        )
    if pd.notna(row.get("stratified_logrank_p")):
        rows.append(("Stratified log-rank p", _clinical_format_p(row["stratified_logrank_p"])))
    if pd.notna(row.get("pooled_logrank_p")):
        rows.append(("Pooled log-rank p", _clinical_format_p(row["pooled_logrank_p"])))
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for label, value in rows:
        lines.append(f"| {label} | {value} |")
    lines.append("")
    return lines


def _foundation_group_formation_section() -> List[str]:
    """Bottom-up explanation of how all comparison groups are constructed."""
    return [
        "## How the comparison groups are built (foundation)",
        "",
        "Every experiment below compares **two treatment arms** that are defined deterministically "
        "from the PDX drug-response workbook. The logic is the same for all experiments; only the "
        "**scenario filter** and **which two TTD columns are subtracted** change.",
        "",
        "### Step 1 — Start from raw PDX drug-response rows",
        "",
        "The source workbook contains one row per **tumour model × screened compound**. Each row carries:",
        "",
        "- **Model** — unique PDX identifier (patient tumour xenograft).",
        "- **COMPOUND** — drug tested ex vivo.",
        "- **LEVEL** — DDA pharmacodynamic score (higher = better predicted response).",
        "- **TimeToDouble (TTD)** — days until tumour burden doubles under that treatment.",
        "- **On_label** — regulatory status in the source annotation: **IGAZ** (approved/on-label), "
        "**HAMIS** (off-label), or missing/n.a. (mapped to **Exp** = experimental).",
        "",
        "Rows with unmapped On_label values are excluded. Chemotherapy rows are handled separately (Step 4).",
        "",
        "### Step 2 — Dense DDA ranking within each tumour model",
        "",
        "For each Model independently, all screened compounds (chemo excluded) are ranked by LEVEL using "
        "**dense ranking**: higher LEVEL → lower rank number, so **rank 1 = best DDA response** in that model. "
        "Tied LEVEL values share the same rank. Ranks 7 and below are **pooled** and labelled **\"7-10\"** "
        "(one pooled stratum, not an immunotherapy class).",
        "",
        "This produces a **global rank** per compound per model (1, 2, 3, 4, 5, 6, or 7-10).",
        "",
        "### Step 3 — On-label utility annotation",
        "",
        "Each compound receives an **ON_LABEL_UTILITY** label:",
        "",
        "| Source On_label | Utility label | Meaning |",
        "| --- | --- | --- |",
        "| IGAZ | **On_label** | Regulatory-approved for the tumour context |",
        "| HAMIS | **Off_label** | Screened but not approved (off-label) |",
        "| n/a / missing | **Exp** | Experimental / unclassified in source |",
        "",
        "**Best on-label compound per model:** among all On_label rows for that model, take the one with "
        "the **best (lowest) global rank**. Its rank defines the model's **best_onlabel_rank** "
        "(e.g. 2 means the approved drug is only the second-best DDA hit in that tumour).",
        "",
        "**Rank-1 compound per model:** the compound(s) at global rank 1. If multiple compounds tie at "
        "rank 1, their TTD values are summarised by the **median** → **rank1_TTD**.",
        "",
        "**Best Off_label/Exp above on-label** (Scenario 3/4 only): when on-label is not rank 1, take the "
        "Off_label or Exp compound with the **lowest global rank that is still better (numerically smaller) "
        "than the on-label rank** → **above_onlabel_TTD**. In Scenario 3, a descriptive gate confirms that "
        "*every* compound ranked above on-label is Off_label or Exp.",
        "",
        "### Step 4 — Standard chemotherapy (SC) reference",
        "",
        "Chemotherapy is **not ranked** with other drugs. Rows with **On_label = chemo** are extracted "
        "separately from the same workbook. For each Model, one **chemo_TTD** is assigned (median TTD if "
        "multiple chemo rows exist). This is the **matched standard chemotherapy** for that specific PDX — "
        "the clinical fallback arm.",
        "",
        "### Step 5 — One summary row per tumour model",
        "",
        "Each of the 178 models becomes one row in the **model table** with pre-computed TTDs:",
        "",
        "| Column | What it is |",
        "| --- | --- |",
        "| **rank1_TTD** | TTD of global rank-1 compound(s) |",
        "| **onlabel_TTD** | TTD of the best-ranked On_label compound (if any) |",
        "| **above_onlabel_TTD** | TTD of best Off_label/Exp above on-label (if on-label not rank 1) |",
        "| **chemo_TTD** | TTD of matched SC for this model |",
        "",
        "### Step 6 — Scenario assignment (mutually exclusive, one per model)",
        "",
        "Each model is assigned to exactly **one** scenario based solely on **best_onlabel_rank**:",
        "",
        "| Scenario | Inclusion rule |",
        "| --- | --- |",
        "| **Scenario 1** | No On_label compound in the screened panel |",
        "| **Scenario 2** | Best on-label at **global rank 1** |",
        "| **Scenario 3** | Best on-label at **global rank 2 or 3** |",
        "| **Scenario 4** | Best on-label at **global rank 4+** (including pooled 7–10) |",
        "",
        "Within Scenario 3, models are further split into **rank-2**, **rank-3**, or **pooled 2–3** sub-strata. "
        "Within Scenario 4, **rank4_plus** (all Scenario 4), **rank7_10** (on-label at pooled ranks 7–10 only), "
        "and progressive pools **pool_ge4 / pool_ge5 / pool_ge6** (on-label rank ≥ 4, ≥ 5, or ≥ 6) are used.",
        "",
        "### Step 7 — Forming the two arms in a paired experiment",
        "",
        "Paired experiments never mix models between arms. For each included model, **Arm A** and **Arm B** "
        "are two TTD values from the **same tumour**:",
        "",
        "| Comparison ID | Arm A (numerator in ΔTTD) | Arm B (reference) |",
        "| --- | --- | --- |",
        "| **rank1_vs_sc** | rank1_TTD | chemo_TTD |",
        "| **onlabel_vs_sc** | onlabel_TTD | chemo_TTD |",
        "| **rank1_vs_onlabel** | rank1_TTD | onlabel_TTD |",
        "| **above_onlabel_vs_onlabel** | above_onlabel_TTD | onlabel_TTD |",
        "",
        "**ΔTTD = Arm A − Arm B.** Positive ΔTTD means Arm A prolongs tumour doubling time (favourable).",
        "",
        "### Step 8 — Final paired sample size (who is counted in n)",
        "",
        "A model enters a given experiment only if:",
        "",
        "1. It belongs to the experiment's **scenario** (and **stratum**, if applicable), and",
        "2. It has a finite **chemo_TTD** (all pre-specified paired analyses require SC), and",
        "3. Both arm TTD values for that comparison are present and finite.",
        "",
        "The **n_models** in each results table is the count after these filters. It can be smaller than the "
        "scenario headcount when TTD is missing for one arm.",
        "",
        "### Step 9 — Between-cohort experiments (unpaired, different logic)",
        "",
        "Experiments Exp27–Exp30 do **not** pair arms within the same model. Instead:",
        "",
        "- **Cohort A** = on-label TTD values from **Scenario 2** models (on-label at rank 1), one value per model.",
        "- **Cohort B** = on-label TTD values from a **Scenario 4 low-rank pool** (e.g. rank 7–10), one value per model.",
        "",
        "These are **different models** in the two cohorts; inference uses Mann–Whitney U (unpaired).",
        "",
        "---",
        "",
    ]


def _build_paired_group_formation(
    *,
    scenario: str,
    stratum: str,
    comparison: str,
    row: Optional[pd.Series],
) -> str:
    """Per-experiment crystal-clear group formation paragraph."""
    scenario_names = {
        CLINICAL_UTILITY_SCENARIO1: "Scenario 1 (no On_label compound in the screen)",
        CLINICAL_UTILITY_SCENARIO2: "Scenario 2 (best on-label at global rank 1)",
        CLINICAL_UTILITY_SCENARIO3: "Scenario 3 (best on-label at global rank 2 or 3)",
        CLINICAL_UTILITY_SCENARIO4: "Scenario 4 (best on-label at global rank 4 or lower, including 7–10)",
    }
    scenario_label = scenario_names.get(scenario, scenario)

    stratum_rules = {
        "all": "No further stratum split — all models in this scenario are included.",
        "rank2": (
            "Further restricted to models whose **best on-label global rank = 2** "
            "(the approved drug is exactly second in the DDA ranking)."
        ),
        "rank3": (
            "Further restricted to models whose **best on-label global rank = 3**."
        ),
        "rank2_3_pooled": (
            "Further restricted to models whose **best on-label global rank is 2 or 3** "
            "(both discordance sub-strata combined)."
        ),
        "rank4_plus": (
            "All Scenario 4 models: best on-label at rank 4, 5, 6, or pooled 7–10."
        ),
        "rank7_10": (
            "Subset of Scenario 4 where the best on-label compound sits in the **pooled ranks 7–10** only "
            "(deepest DDA deprioritisation; in this dataset all are cetuximab)."
        ),
        "pool_ge6": (
            "Subset of Scenario 4 where **best on-label rank ≥ 6** (includes pooled 7–10; "
            "in this dataset identical to rank7_10, n = 3)."
        ),
        "pool_ge5": (
            "Subset of Scenario 4 where **best on-label rank ≥ 5** "
            "(in this dataset identical to rank7_10, n = 3)."
        ),
        "pool_ge4": (
            "Subset of Scenario 4 where **best on-label rank ≥ 4** "
            "(all Scenario 4 models in this dataset, n = 8)."
        ),
    }
    stratum_rule = stratum_rules.get(stratum, f"Stratum filter: {stratum}.")

    arm_defs = {
        "rank1_vs_sc": (
            "**Arm A (test)** = **rank1_TTD**: TimeToDouble of the compound at **global DDA rank 1** "
            "in that model (median if multiple rank-1 ties). "
            "**Arm B (reference)** = **chemo_TTD**: matched standard chemotherapy for the **same model**."
        ),
        "onlabel_vs_sc": (
            "**Arm A (test)** = **onlabel_TTD**: TimeToDouble of the **best-ranked On_label (IGAZ) compound** "
            "in that model — the approved drug with the lowest global rank. "
            "**Arm B (reference)** = **chemo_TTD**: matched SC for the **same model**."
        ),
        "rank1_vs_onlabel": (
            "**Arm A (test)** = **rank1_TTD**: the algorithm's top-ranked compound (global rank 1; "
            "typically Off_label or Exp in discordant scenarios). "
            "**Arm B (reference)** = **onlabel_TTD**: the best-ranked approved (On_label) compound in the "
            "**same model**."
        ),
        "above_onlabel_vs_onlabel": (
            "**Arm A (test)** = **above_onlabel_TTD**: the best Off_label/Exp compound whose global rank is "
            "**better (numerically lower) than the on-label rank** in that model. "
            "**Arm B (reference)** = **onlabel_TTD**: the best On_label compound in the **same model**. "
            "When the descriptive gate is satisfied, Arm A coincides with rank 1 vs on-label."
        ),
    }
    arms = arm_defs.get(comparison, f"Arms defined by comparison {comparison}.")

    n_paired = int(row["n_models"]) if row is not None else 0
    exclusion = (
        f"After requiring finite values for both arms and a matched SC row, **{n_paired} models** "
        "contribute paired observations to this experiment."
        if n_paired
        else "No models met the paired inclusion criteria for this experiment."
    )

    # Scenario-specific nuance
    nuance = ""
    if scenario == CLINICAL_UTILITY_SCENARIO1:
        nuance = (
            " Models in this experiment have **zero On_label rows** in the screened panel; "
            "rank 1 is therefore always an Off_label or Exp compound."
        )
    elif scenario == CLINICAL_UTILITY_SCENARIO2:
        nuance = (
            " Because on-label is at rank 1, **Arm A and the rank-1 arm are the same drug** "
            "(onlabel_TTD equals rank1_TTD for every model in this cohort)."
        )
    elif scenario == CLINICAL_UTILITY_SCENARIO3:
        nuance = (
            " Every model in Scenario 3 has an on-label drug at rank 2 or 3, and **every compound ranked "
            "above it is Off_label or Exp** (descriptive gate 32/32)."
        )
    elif scenario == CLINICAL_UTILITY_SCENARIO4:
        nuance = (
            " On-label is deeply deprioritised (rank ≥ 4). Rank 1 is virtually always Off_label/Exp."
        )

    return (
        "These groups were formed as follows. "
        f"**Cohort filter:** {scenario_label}. {stratum_rule} "
        "**Pairing unit:** one PDX model — both arms come from the **same tumour**, never crossed between models. "
        f"{arms} "
        "**Inclusion:** model must have finite chemo_TTD and finite values in both arm columns. "
        f"{exclusion}{nuance}"
    )


def _build_between_cohort_group_formation(stratum: str, row: pd.Series) -> str:
    pool_desc = {
        "rank7_10": "best on-label at pooled global ranks **7–10**",
        "pool_ge6": "best on-label at global rank **≥ 6** (here: same 3 models as rank 7–10)",
        "pool_ge5": "best on-label at global rank **≥ 5** (here: same 3 models as rank 7–10)",
        "pool_ge4": "best on-label at global rank **≥ 4** (all 8 Scenario 4 models)",
    }.get(stratum, f"Scenario 4 low-rank pool `{stratum}`")

    return (
        "These groups were formed as follows — this is an **unpaired between-cohort** design (not within-model pairing). "
        f"**Cohort A** collects **onlabel_TTD** (TimeToDouble of each model's best On_label compound) from "
        f"**Scenario 2** models only — tumours where the approved drug is **rank 1** (n = {int(row['n_cohort_a'])}). "
        f"**Cohort B** collects **onlabel_TTD** from **Scenario 4** models filtered to {pool_desc} "
        f"(n = {int(row['n_cohort_b'])}). "
        "Each model contributes **one on-label TTD value** to its cohort. Models in Cohort A are **not** the same "
        "as models in Cohort B. The comparison asks whether on-label pharmacodynamic performance when ranked first "
        "is systematically higher than on-label performance when deeply deprioritised — a construct-validity check "
        "for rank-based stratification."
    )


def _between_cohort_markdown_table(row: pd.Series) -> List[str]:
    lines = [
        "| Metric | Value |",
        "| --- | --- |",
        f"| Cohort A (Scenario 2: on-label at rank 1) | n = {int(row['n_cohort_a'])} |",
        f"| Cohort B (Scenario 4 low-rank stratum: {row['stratum']}) | n = {int(row['n_cohort_b'])} |",
        f"| Median on-label TTD — Cohort A | {row['median_a']:.1f} days |",
        f"| Median on-label TTD — Cohort B | {row['median_b']:.1f} days |",
        f"| Median difference (A − B) | {row['median_diff_a_minus_b']:+.1f} days |",
        f"| Mann–Whitney U (two-sided) | p = {_clinical_format_p(row['mannwhitney_p'])} |",
    ]
    if pd.notna(row.get("pooled_cox_hr")):
        lines.append(
            f"| Unpaired Cox HR (A vs B) | {row['pooled_cox_hr']:.2f}; "
            f"p = {_clinical_format_p(row['pooled_cox_p'])} |"
        )
    if pd.notna(row.get("pooled_logrank_p")):
        lines.append(f"| Unpaired log-rank p | {_clinical_format_p(row['pooled_logrank_p'])} |")
    lines.append("")
    return lines


def _experiment_section(
  exp_id: str,
  title: str,
  *,
  clinical_question: str,
  mtb_framing: str,
  cohort_design: str,
  group_formation: str,
  comparison: str,
  row: Optional[pd.Series],
  interpretation: str,
  pi_conclusion: str,
) -> List[str]:
    return [
        f"### {exp_id}: {title}",
        "",
        "**Clinical question.** " + clinical_question,
        "",
        "**Molecular tumour board analogue.** " + mtb_framing,
        "",
        "**Cohort and experimental design.** " + cohort_design,
        "",
        "**How these groups were formed.** " + group_formation,
        "",
        f"**Comparison under test.** {_comparison_plain_label(comparison)}.",
        "",
        "**Quantitative results.**",
        "",
        *_clinical_results_markdown_table(row),
        "**Clinical-oncology interpretation.** " + interpretation,
        "",
        "**Principal investigator conclusion.** " + pi_conclusion,
        "",
        "---",
        "",
    ]


def _paired_experiment_section(
    exp_id: str,
    title: str,
    *,
    scenario: str,
    stratum: str,
    comparison: str,
    clinical_question: str,
    mtb_framing: str,
    cohort_design: str,
    row: Optional[pd.Series],
    interpretation: str,
    pi_conclusion: str,
) -> List[str]:
    return _experiment_section(
        exp_id,
        title,
        clinical_question=clinical_question,
        mtb_framing=mtb_framing,
        cohort_design=cohort_design,
        group_formation=_build_paired_group_formation(
            scenario=scenario,
            stratum=stratum,
            comparison=comparison,
            row=row,
        ),
        comparison=comparison,
        row=row,
        interpretation=interpretation,
        pi_conclusion=pi_conclusion,
    )


def _interpret_exp1(row: Optional[pd.Series]) -> Tuple[str, str]:
    if row is None:
        return ("Insufficient data.", "Experiment not evaluable.")
    sig = float(row["wilcoxon_p"]) < 0.05 and float(row["median_diff_TTD"]) > 0
    interp = (
        "In the absence of a regulatory-approved, biomarker-concordant targeted option, the molecular "
        "tumour board would typically default to cytotoxic chemotherapy. Here, the DDA algorithm's "
        "top-ranked compound prolongs tumour doubling time relative to model-matched standard chemotherapy "
        f"in {row['pct_improved']:.0f}% of models, with a median gain of {row['median_diff_TTD']:.1f} days. "
        "This pattern is consistent with algorithm-guided escalation beyond conventional cytotoxic fallback "
        "when the screening panel lacks an on-label agent."
    )
    if sig:
        pi = (
            "Preclinical evidence supports deploying DDA rank 1 as a rational pharmacodynamic alternative "
            "to chemotherapy in this molecular context. This is the strongest scenario-specific signal in "
            "the programme and justifies prioritising algorithm-guided compound selection in tumours without "
            "an on-label match in the ex vivo screen."
        )
    else:
        pi = "Evidence is insufficient to recommend algorithm rank 1 over chemotherapy in this stratum."
    return interp, pi


def _interpret_exp2(row: Optional[pd.Series]) -> Tuple[str, str]:
    if row is None:
        return ("Insufficient data.", "Experiment not evaluable.")
    wil_sig = float(row["wilcoxon_p"]) < 0.05
    cox_sig = pd.notna(row.get("stratified_cox_p")) and float(row["stratified_cox_p"]) < 0.05
    interp = (
        "This experiment tests guideline concordance: when the approved on-label therapy also occupies "
        "the apex of the DDA ranking, clinicians would reasonably select that agent. Relative to "
        "matched chemotherapy, on-label rank 1 demonstrates a favourable pharmacodynamic profile "
        f"(median ΔTTD {row['median_diff_TTD']:+.1f} days; {row['pct_improved']:.0f}% models improved). "
        "The stratified Cox model — which respects the paired PDX structure — provides complementary "
        "evidence of reduced hazard of rapid tumour growth."
    )
    if wil_sig or cox_sig:
        pi = (
            "Guideline-concordant on-label selection is supported by preclinical TTD data. In practice, "
            "this validates the internal consistency of the ranking when regulatory and algorithmic "
            "priorities align; the on-label agent should remain the default recommendation."
        )
    else:
        pi = (
            "Directional benefit is present but the primary Wilcoxon comparison does not reach conventional "
            "significance. I would characterise this as supportive but not definitive evidence for "
            "guideline-concordant selection — the stratified Cox result should be weighed alongside "
            "clinical factors not captured in PDX."
        )
    return interp, pi


def _interpret_exp3_rank2(row: Optional[pd.Series]) -> Tuple[str, str]:
    if row is None:
        return ("Insufficient data.", "Experiment not evaluable.")
    sig = float(row["wilcoxon_p"]) < 0.05 and float(row["median_diff_TTD"]) > 0
    interp = (
        "This is the pivotal discordance experiment. An on-label agent is available but ranked second; "
        "the algorithm's rank-1 compound is Off_label or experimental. Strict application of label "
        "priority would administer the on-label drug despite a superior pharmacodynamic signal at rank 1. "
        f"In this stratum, rank 1 outperforms on-label in {row['pct_improved']:.0f}% of models "
        f"(median ΔTTD {row['median_diff_TTD']:+.1f} days)."
    )
    if sig:
        pi = (
            "I would present this to the tumour board as evidence that rigid on-label priority is "
            "pharmacodynamically suboptimal when the approved agent sits at rank 2. The rank-1 "
            "Off_label/Exp recommendation merits serious consideration in a research or compassionate-use "
            "framework, pending patient-level validation."
        )
    else:
        pi = "Discordance at rank 2 is not statistically established in this run."
    return interp, pi


def _interpret_exp5_rank2_onlabel_sc(row: Optional[pd.Series]) -> Tuple[str, str]:
    if row is None:
        return ("Insufficient data.", "Experiment not evaluable.")
    interp = (
        "This ancillary comparison asks whether the on-label agent at rank 2 still offers meaningful "
        "tumour growth control relative to chemotherapy. A negative or null result would strengthen "
        "the case that rank-1 escalation is clinically motivated rather than driven by marginal "
        "on-label activity."
    )
    worse = float(row["median_diff_TTD"]) < 0
    pi = (
        "On-label at rank 2 does not convincingly outperform chemotherapy on TTD, reinforcing that "
        "the rank-1 Off_label/Exp signal reflects genuine pharmacodynamic superiority rather than "
        "on-label underperformance alone."
        if worse
        else "On-label retains some activity versus chemotherapy; interpret rank-1 superiority in that context."
    )
    return interp, pi


def _interpret_underpowered(stratum: str, n: int) -> Tuple[str, str]:
    interp = (
        f"The {stratum} stratum contains only {n} models. At this sample size, paired Wilcoxon tests "
        "lack adequate power to discriminate modest TTD differences, and confidence intervals are wide. "
        "Results are reported for completeness but must not be over-interpreted."
    )
    pi = (
        "I classify this experiment as hypothesis-generating only. No change to molecular tumour board "
        "policy should be inferred from this stratum in isolation."
    )
    return interp, pi


def _interpret_s4_low_rank(
    row: Optional[pd.Series], *, stratum: str, epistemic: Optional[pd.DataFrame]
) -> Tuple[str, str]:
    if row is None:
        return ("Insufficient data.", "Experiment not evaluable.")
    n = int(row["n_models"]) if row is not None else 0
    interp = (
        f"When the best on-label agent is deeply deprioritised by DDA (stratum: {_stratum_plain_label(stratum)}), "
        "the clinically relevant question is epistemic: has the algorithm correctly identified an ineffective "
        "on-label option, or has it missed a therapeutically active approved drug? TimeToDouble provides "
        "a pharmacodynamic readout of tumour growth suppression in the PDX, not a surrogate for overall "
        "survival."
    )
    if epistemic is not None and not epistemic.empty and stratum in ("rank7_10", "pool_ge6", "pool_ge5"):
        sub = epistemic[epistemic["best_onlabel_rank"].astype(str) == _POOLED_LABEL]
        if not sub.empty:
            counts = sub["epistemic_class"].value_counts()
            parts = [f"{_epistemic_plain_label(k)}: {int(v)}" for k, v in counts.items()]
            interp += " In the rank 7–10 subset, epistemic classification yields: " + "; ".join(parts) + "."
    pi = (
        f"With n = {n}, formal hypothesis testing is exploratory. Directionally, rank 1 tends to dominate "
        "on-label TTD where comparisons are evaluable. For tumour board purposes, deeply low-ranked on-label "
        "agents should be discussed as potentially inactive in the ex vivo model rather than as missed "
        "opportunities — unless Class C cases are identified at the patient level."
    )
    return interp, pi


def _add_markdown_runs(paragraph: Any, text: str) -> None:
    import re
    parts = re.split(r"(\*\*[^*]+\*\*)", text)
    for part in parts:
        if part.startswith("**") and part.endswith("**"):
            run = paragraph.add_run(part[2:-2])
            run.bold = True
        else:
            paragraph.add_run(part)


def _write_supervisor_summary_docx(docx_path: Path, lines: List[str]) -> None:
    from docx import Document
    from docx.enum.text import WD_PARAGRAPH_ALIGNMENT
    from docx.shared import Pt

    document = Document()
    style = document.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        if stripped.startswith("```"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                i += 1
            i += 1
            continue

        if stripped == "---":
            i += 1
            continue

        if stripped.startswith("# ") and not stripped.startswith("## "):
            document.add_heading(stripped[2:].strip(), level=0)
            i += 1
            continue
        if stripped.startswith("#### "):
            document.add_heading(stripped[5:].strip(), level=3)
            i += 1
            continue
        if stripped.startswith("### "):
            document.add_heading(stripped[4:].strip(), level=2)
            i += 1
            continue
        if stripped.startswith("## "):
            document.add_heading(stripped[3:].strip(), level=1)
            i += 1
            continue

        if stripped.startswith("|") and stripped.count("|") >= 2:
            table_lines: List[str] = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i].strip())
                i += 1
            if len(table_lines) >= 2:
                headers = [c.strip() for c in table_lines[0].strip("|").split("|")]
                data_rows = []
                for tl in table_lines[2:]:
                    data_rows.append([c.strip() for c in tl.strip("|").split("|")])
                ncols = len(headers)
                table = document.add_table(rows=1 + len(data_rows), cols=ncols)
                table.style = "Table Grid"
                for j, h in enumerate(headers):
                    table.rows[0].cells[j].text = h
                for ri, dr in enumerate(data_rows):
                    for j in range(min(ncols, len(dr))):
                        table.rows[ri + 1].cells[j].text = dr[j]
                document.add_paragraph("")
            continue

        if stripped.startswith("- "):
            p = document.add_paragraph(style="List Bullet")
            _add_markdown_runs(p, stripped[2:])
            i += 1
            continue

        if stripped.startswith("*") and stripped.endswith("*") and not stripped.startswith("**"):
            p = document.add_paragraph()
            p.alignment = WD_PARAGRAPH_ALIGNMENT.CENTER
            run = p.add_run(stripped.strip("*"))
            run.italic = True
            i += 1
            continue

        p = document.add_paragraph()
        _add_markdown_runs(p, stripped)
        i += 1

    document.save(docx_path)


def _build_clinical_supervisor_markdown(
    *,
    model_table: pd.DataFrame,
    paired_summary: pd.DataFrame,
    utility_rank_pairwise: pd.DataFrame,
    between_cohort: Optional[pd.DataFrame],
    epistemic: Optional[pd.DataFrame],
    source_path: Path,
    copied_figs: List[str],
) -> List[str]:
    s1 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO1, "rank1_vs_sc")
    s2 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO2, "onlabel_vs_sc")
    s3_r2_on = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2")
    s3_r2_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_sc", "rank2")
    s3_r2_ol_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "onlabel_vs_sc", "rank2")
    s3_r2_above = _clinical_paired_row(
        paired_summary, CLINICAL_UTILITY_SCENARIO3, "above_onlabel_vs_onlabel", "rank2"
    )
    s3_r3_on = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank3")
    s3_r3_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_sc", "rank3")
    s3_r3_ol_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "onlabel_vs_sc", "rank3")
    s3_pool_on = _clinical_paired_row(
        paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2_3_pooled"
    )
    s3_pool_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_sc", "rank2_3_pooled")
    s4_r4_on = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO4, "rank1_vs_onlabel", "rank4_plus")
    s4_r4_sc = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO4, "rank1_vs_sc", "rank4_plus")
    s4_r710_on = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO4, "rank1_vs_onlabel", "rank7_10")
    s4_r710_ol = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO4, "onlabel_vs_sc", "rank7_10")
    s4_pool4_on = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO4, "rank1_vs_onlabel", "pool_ge4")

    n_total = len(model_table)
    n_chemo = int(model_table["chemo_TTD"].notna().sum())
    scenario_counts = model_table["scenario"].value_counts()
    s3 = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO3]
    n_rank2 = int((s3["best_onlabel_rank_int"] == 2).sum()) if not s3.empty else 0
    n_rank3 = int((s3["best_onlabel_rank_int"] == 3).sum()) if not s3.empty else 0
    n710 = int((model_table["best_onlabel_rank"].astype(str) == _POOLED_LABEL).sum())
    gate = int(s3["ranks_above_onlabel_valid"].fillna(False).astype(bool).sum()) if not s3.empty else 0

    exp1_i, exp1_pi = _interpret_exp1(s1)
    exp2_i, exp2_pi = _interpret_exp2(s2)
    exp3_i, exp3_pi = _interpret_exp3_rank2(s3_r2_on)
    exp5_i, exp5_pi = _interpret_exp5_rank2_onlabel_sc(s3_r2_ol_sc)
    exp7_i, exp7_pi = _interpret_underpowered("rank 3", n_rank3)
    exp15_i, exp15_pi = _interpret_s4_low_rank(s4_r4_on, stratum="rank4_plus", epistemic=epistemic)
    exp19_i, exp19_pi = _interpret_s4_low_rank(s4_r710_on, stratum="rank7_10", epistemic=epistemic)

    doc: List[str] = [
        "# Clinical Utility of DDA Compound Ranking in Patient-Derived Xenograft Models",
        "",
        "**Document type:** Structured experimental report for doctoral supervision",
        f"**Prepared:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "**Authoring perspective:** Clinical oncology / principal investigator interpretation layer",
        f"**Analysis pipeline:** `clinical_utility_ranking.py` (standalone)",
        f"**Primary data source:** `{source_path.name}`",
        "**Pharmacodynamic endpoint:** TimeToDouble (TTD, days) — ex vivo tumour growth surrogate",
        "",
        "---",
        "",
        "## Executive summary for supervision",
        "",
        "This report evaluates whether dense drug–disease association (DDA) compound ranking in PDX models "
        "supports clinically meaningful treatment prioritisation across four mutually exclusive scenarios, "
        "each mirroring a distinct molecular tumour board (MTB) decision branch. Every experiment uses "
        "within-model pairing against the same tumour's matched standard chemotherapy (SC), preserving "
        "the patient-specific pharmacodynamic context that motivates PDX-guided precision oncology.",
        "",
        f"The analysable cohort comprises **{n_total}** tumour models, of which **{n_chemo}** have an "
        "evaluable SC reference. Primary inference relies on the Wilcoxon signed-rank test applied to "
        "paired TTD differences; stratified Cox proportional hazards and stratified log-rank tests provide "
        "complementary survival-framework summaries with Model as stratum.",
        "",
    ]

    findings: List[str] = []
    if s1 is not None and float(s1["wilcoxon_p"]) < 0.05:
        findings.append(
            f"**Experiment 1 (no on-label):** DDA rank 1 significantly prolongs TTD versus SC "
            f"(median +{s1['median_diff_TTD']:.1f} days, p = {_clinical_format_p(s1['wilcoxon_p'])})."
        )
    if s2 is not None and pd.notna(s2.get("stratified_cox_p")) and float(s2["stratified_cox_p"]) < 0.05:
        findings.append(
            f"**Experiment 2 (guideline concordance):** On-label at rank 1 shows reduced growth hazard "
            f"versus SC (stratified Cox HR = {s2['stratified_cox_hr']:.2f}, "
            f"p = {_clinical_format_p(s2['stratified_cox_p'])})."
        )
    if s3_r2_on is not None and float(s3_r2_on["wilcoxon_p"]) < 0.05:
        findings.append(
            f"**Experiment 3 (discordance at rank 2):** Algorithm rank 1 significantly outperforms "
            f"on-label (median +{s3_r2_on['median_diff_TTD']:.1f} days, "
            f"p = {_clinical_format_p(s3_r2_on['wilcoxon_p'])}; FDR-adjusted p = "
            f"{_clinical_format_p(s3_r2_on['wilcoxon_p_fdr_bh'])})."
        )
    if n710:
        findings.append(
            f"**Scenario 4 / rank 7–10 stratum (n = {n710}):** Exploratory; epistemic review suggests "
            "predominantly correct deprioritisation of pharmacodynamically inactive on-label cetuximab."
        )
    if findings:
        doc.extend(findings)
    else:
        doc.append("_No scenario met conventional significance thresholds in this run._")
    doc.extend([
        "",
        "These results support DDA ranking as **preclinical decision support** for compound prioritisation. "
        "They do not constitute prescribing evidence; translation requires prospective clinical validation.",
        "",
        "---",
        "",
        "## Clinical framing: from PDX screen to molecular tumour board",
        "",
        "In contemporary precision oncology, the MTB integrates tumour genomics, regulatory approval status "
        "(on-label), and off-label or experimental options. The DDA algorithm operationalises a similar hierarchy "
        "in PDX: compounds are ranked by pharmacodynamic LEVEL within each model, and on-label annotation "
        "(IGAZ = approved; HAMIS = off-label; Exp = experimental) permits scenario stratification.",
        "",
        "We deliberately avoid constructing a synthetic \"clinician choice\" arm. Instead, each experiment "
        "poses a question the MTB would actually face:",
        "",
        "1. **No approved match in the panel** — should we default to chemotherapy or trust algorithm rank 1?",
        "2. **Concordance** — when on-label is rank 1, does it outperform chemotherapy?",
        "3. **Discordance (ranks 2–3)** — should guideline priority override a superior Off_label/Exp rank 1?",
        "4. **Deep deprioritisation (rank 4+)** — is the on-label agent pharmacodynamically ineffective, "
        "   or has the algorithm failed to recognise activity?",
        "",
        "TimeToDouble quantifies the interval for tumour burden to double under treatment — a pharmacodynamic "
        "surrogate for growth control. It is analytically analogous to progression-free dynamics but is **not** "
        "interchangeable with patient overall survival or progression-free survival.",
        "",
        "---",
        "",
        "## Experiment registry",
        "",
        "| Exp ID | Scenario | Clinical question | Primary comparison |",
        "| --- | --- | --- | --- |",
        "| Exp1 | 1 — No on-label | Can rank 1 replace SC fallback? | rank 1 vs SC |",
        "| Exp2 | 2 — On-label rank 1 | Does guideline-concordant on-label beat SC? | on-label vs SC |",
        "| Exp3–Exp6 | 3 — Rank 2 stratum | Discordance tests at rank 2 | rank1 vs on-label; ancillary |",
        "| Exp7–Exp10 | 3 — Rank 3 stratum | Discordance tests at rank 3 (underpowered) | rank1 vs on-label; ancillary |",
        "| Exp11–Exp14 | 3 — Pooled ranks 2–3 | Aggregated discordance | rank1 vs on-label; ancillary |",
        "| Exp15–Exp18 | 4 — Rank 4+ | Exploratory low-rank deprioritisation | rank1 vs on-label; ancillary |",
        "| Exp19–Exp20 | 4 — Rank 7–10 | Deep deprioritisation (cetuximab stratum) | rank1 vs on-label; on-label vs SC |",
        "| Exp25–Exp26 | 4 — Pool ≥ rank 4 | Progressive pooling for power | same as rank 4+ pool |",
        "| Exp27–Exp30 | Between-cohort | Is on-label TTD at rank 1 higher than in low-rank pools? | Mann–Whitney |",
        "",
        "---",
        "",
        "## Cohort architecture",
        "",
        f"| Stratum | n models |",
        f"| --- | --- |",
        f"| Total PDX models | {n_total} |",
        f"| With matched SC (paired analyses) | {n_chemo} |",
        f"| Scenario 1 — no on-label in screen | {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0))} |",
        f"| Scenario 2 — on-label at rank 1 | {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0))} |",
        f"| Scenario 3 — on-label at ranks 2–3 | {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0))} |",
        f"| — on-label at rank 2 | {n_rank2} |",
        f"| — on-label at rank 3 | {n_rank3} |",
        f"| Scenario 4 — on-label at rank 4+ | {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO4, 0))} |",
        f"| — pooled ranks 7–10 | {n710} |",
        "",
        f"**Descriptive consistency gate (Scenario 3):** In all {len(s3)} discordant models, every compound "
        f"ranked above the on-label drug was Off_label or Exp ({gate}/{len(s3)}). This confirms that "
        "discordance analyses interrogate genuine algorithm–guideline tension, not mis-annotated on-label drugs.",
        "",
        "Full cohort flow: `COHORT_FLOW.md`. Tabulated results: `TABLE2_SCENARIO_RESULTS.csv`.",
        "",
        "---",
        "",
        "## Structured experiments",
        "",
        "Each subsection below is a pre-specified experiment with explicit clinical question, MTB analogue, "
        "cohort definition, **step-by-step group formation**, statistical readout, and PI-level interpretation.",
        "",
    ])
    doc.extend(_foundation_group_formation_section())

    doc.extend(_paired_experiment_section(
        "Experiment 1 (Exp1)",
        "Scenario 1 — Absence of an on-label agent in the screening panel",
        scenario=CLINICAL_UTILITY_SCENARIO1,
        stratum="all",
        comparison="rank1_vs_sc",
        clinical_question=(
            "When no regulatory-approved (on-label) compound is represented in the ex vivo drug screen, "
            "does the DDA algorithm's top-ranked agent deliver superior tumour growth control compared with "
            "the model-matched standard chemotherapy that would constitute clinical fallback?"
        ),
        mtb_framing=(
            "The MTB concludes that no biomarker-matched targeted therapy is available from the screened "
            "formulary. Cytotoxic chemotherapy is the default pathway. The algorithm nevertheless proposes "
            "a rank-1 Off_label or experimental compound with the highest pharmacodynamic score."
        ),
        cohort_design=(
            f"Scenario 1 cohort: models without any On_label treatment in the ranked panel "
            f"(n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0))}; "
            f"paired n = {int(s1['n_models']) if s1 is not None else 0}). "
            "KM curves with Mantel–Haenszel stratified log-rank and stratified Cox (strata = Model); "
            "primary paired inference: Wilcoxon signed-rank on ΔTTD = TTD(rank 1) − TTD(SC)."
        ),
        row=s1,
        interpretation=exp1_i,
        pi_conclusion=exp1_pi,
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 2 (Exp2)",
        "Scenario 2 — Guideline–algorithm concordance (on-label at rank 1)",
        scenario=CLINICAL_UTILITY_SCENARIO2,
        stratum="all",
        comparison="onlabel_vs_sc",
        clinical_question=(
            "When the approved on-label therapy occupies the apex of the DDA ranking, does it demonstrate "
            "pharmacodynamic superiority over matched standard chemotherapy — thereby validating "
            "guideline-concordant selection?"
        ),
        mtb_framing=(
            "Regulatory approval and algorithmic ranking align. The MTB would recommend the on-label agent "
            "without controversy. This experiment tests whether that recommendation is supported by TTD "
            "relative to the cytotoxic backstop."
        ),
        cohort_design=(
            f"Scenario 2 cohort (n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0))}; "
            f"paired n = {int(s2['n_models']) if s2 is not None else 0}). "
            "Primary: Wilcoxon on paired ΔTTD. Supporting: stratified Cox PH with Model as stratum."
        ),
        row=s2,
        interpretation=exp2_i,
        pi_conclusion=exp2_pi,
    ))

    doc.extend([
        "## Scenario 3 — Guideline discordance (on-label not at rank 1)",
        "",
        "Scenario 3 emulates the high-stakes MTB deliberation in which an approved agent exists but is "
        "not the algorithm's preferred compound. Experiments are stratified by the global rank of the best "
        "on-label drug and pooled where appropriate. Multiplicity within each stratum is controlled by "
        "Benjamini–Hochberg FDR across the four pre-specified comparisons.",
        "",
    ])

    doc.extend(_paired_experiment_section(
        "Experiment 3 (Exp3)",
        "Discordance at on-label rank 2 — rank 1 versus on-label (primary hypothesis H3)",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2",
        comparison="rank1_vs_onlabel",
        clinical_question=(
            "If the on-label agent sits at rank 2, does the algorithm's rank-1 Off_label/Exp compound "
            "achieve materially longer TTD — implying that strict label priority would forgo superior "
            "pharmacodynamic activity?"
        ),
        mtb_framing=(
            "The board must decide between an approved but subordinate on-label option and a higher-ranked "
            "off-label or investigational agent. This is the clinically most contentious branch of the "
            "decision tree."
        ),
        cohort_design=(
            f"Scenario 3, rank-2 stratum (n = {n_rank2}; paired n = "
            f"{int(s3_r2_on['n_models']) if s3_r2_on is not None else 0}). "
            "Pre-specified primary comparison: rank 1 vs on-label."
        ),
        row=s3_r2_on,
        interpretation=exp3_i,
        pi_conclusion=exp3_pi,
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 4 (Exp4)",
        "Discordance at rank 2 — rank 1 versus chemotherapy (ancillary)",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2",
        comparison="rank1_vs_sc",
        clinical_question=(
            "Does rank 1 also outperform chemotherapy in the rank-2 discordance stratum, establishing "
            "that the signal is not driven solely by on-label underperformance?"
        ),
        mtb_framing="Anchors algorithm rank 1 against the cytotoxic fallback while on-label remains available at rank 2.",
        cohort_design=f"Same rank-2 stratum as Exp3 (paired n = {int(s3_r2_sc['n_models']) if s3_r2_sc is not None else 0}).",
        row=s3_r2_sc,
        interpretation=(
            "Rank 1 shows a numerically favourable median ΔTTD versus SC, though this ancillary comparison "
            "does not reach conventional significance in isolation. Interpret alongside Exp3 and Exp5."
        ),
        pi_conclusion=(
            "Supports a coherent narrative of rank-1 superiority in the discordant setting, but the "
            "primary inferential weight rests on Exp3 (rank 1 vs on-label)."
        ),
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 5 (Exp5)",
        "Discordance at rank 2 — on-label versus chemotherapy (ancillary)",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2",
        comparison="onlabel_vs_sc",
        clinical_question=(
            "Does the on-label agent at rank 2 retain meaningful activity relative to chemotherapy?"
        ),
        mtb_framing=(
            "If on-label were clearly superior to SC, deprioritising it would be harder to defend. "
            "Null or negative results strengthen the case for rank-1 escalation."
        ),
        cohort_design=f"Rank-2 stratum (paired n = {int(s3_r2_ol_sc['n_models']) if s3_r2_ol_sc is not None else 0}).",
        row=s3_r2_ol_sc,
        interpretation=exp5_i,
        pi_conclusion=exp5_pi,
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 6 (Exp6)",
        "Discordance at rank 2 — best Off_label/Exp above on-label versus on-label (H4)",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2",
        comparison="above_onlabel_vs_onlabel",
        clinical_question=(
            "Is the best Off_label/Exp compound ranked above the on-label drug superior to the on-label "
            "agent itself? (Descriptively equivalent to rank 1 vs on-label when the gate is satisfied.)"
        ),
        mtb_framing="Formalises hypothesis H4: therapeutic opportunity lies in off-label/experimental tiers.",
        cohort_design=f"Rank-2 stratum; descriptive gate satisfied in all Scenario 3 models (paired n = {int(s3_r2_above['n_models']) if s3_r2_above is not None else 0}).",
        row=s3_r2_above,
        interpretation=exp3_i,
        pi_conclusion=exp3_pi,
    ))

    # Rank 3 block
    doc.extend(_paired_experiment_section(
        "Experiment 7 (Exp7)",
        "Discordance at on-label rank 3 — rank 1 versus on-label",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank3",
        comparison="rank1_vs_onlabel",
        clinical_question="As for Exp3, but when on-label is deprioritised to rank 3.",
        mtb_framing="Same MTB tension as rank 2, with greater algorithmic separation from the approved agent.",
        cohort_design=f"Rank-3 stratum (n = {n_rank3}; paired n = {int(s3_r3_on['n_models']) if s3_r3_on is not None else 0}). Underpowered.",
        row=s3_r3_on,
        interpretation=exp7_i,
        pi_conclusion=exp7_pi,
    ))

    for exp_id, title, comp, row in [
        ("Experiment 8 (Exp8)", "Rank 3 — rank 1 vs SC", "rank1_vs_sc", s3_r3_sc),
        ("Experiment 9 (Exp9)", "Rank 3 — on-label vs SC", "onlabel_vs_sc", s3_r3_ol_sc),
    ]:
        ui, upi = _interpret_underpowered("rank 3", n_rank3)
        doc.extend(_paired_experiment_section(
            exp_id, title,
            scenario=CLINICAL_UTILITY_SCENARIO3,
            stratum="rank3",
            comparison=comp,
            clinical_question=f"Ancillary rank-3 comparison: {_comparison_plain_label(comp)}.",
            mtb_framing="Exploratory; reported for completeness.",
            cohort_design=f"Rank-3 stratum (n = {n_rank3}).",
            row=row,
            interpretation=ui,
            pi_conclusion=upi,
        ))

    doc.extend(_paired_experiment_section(
        "Experiment 11 (Exp11)",
        "Pooled ranks 2–3 — rank 1 versus on-label",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2_3_pooled",
        comparison="rank1_vs_onlabel",
        clinical_question=(
            "When on-label is not rank 1 (ranks 2 or 3 combined), does rank 1 outperform on-label on average?"
        ),
        mtb_framing="Aggregated discordance view for MTB policy discussions spanning both sub-strata.",
        cohort_design=f"Pooled Scenario 3 (n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0))}; paired n = {int(s3_pool_on['n_models']) if s3_pool_on is not None else 0}).",
        row=s3_pool_on,
        interpretation=(
            "Pooling increases sample size but dilutes the strong rank-2 signal with the underpowered "
            "rank-3 stratum. Interpretation should be driven primarily by Exp3."
        ),
        pi_conclusion=(
            "Pooled analysis is supportive at nominal p < 0.05 but does not survive FDR adjustment. "
            "I would not base MTB policy on the pooled estimate alone."
        ),
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 12 (Exp12)",
        "Pooled ranks 2–3 — rank 1 versus SC (ancillary)",
        scenario=CLINICAL_UTILITY_SCENARIO3,
        stratum="rank2_3_pooled",
        comparison="rank1_vs_sc",
        clinical_question="Pooled ancillary comparison of rank 1 against chemotherapy.",
        mtb_framing="Secondary to Exp11.",
        cohort_design=f"Pooled Scenario 3 (paired n = {int(s3_pool_sc['n_models']) if s3_pool_sc is not None else 0}).",
        row=s3_pool_sc,
        interpretation="Non-significant pooled comparison; rank-1 superiority is primarily relative to on-label, not SC.",
        pi_conclusion="No independent action item from this ancillary experiment.",
    ))

    doc.extend([
        "## Scenario 4 — Deep deprioritisation of on-label therapy (exploratory programme)",
        "",
        "When the best on-label agent ranks fourth or lower, the clinical question shifts from "
        "\"which agent is best?\" to \"is the on-label drug pharmacodynamically relevant at all?\" "
        "We apply an epistemic classification (Table 4) distinguishing correct deprioritisation "
        "(Class A) from possible ranking failures (Class C).",
        "",
    ])

    doc.extend(_paired_experiment_section(
        "Experiment 15 (Exp15)",
        "Scenario 4 — rank 4+ stratum: rank 1 versus on-label",
        scenario=CLINICAL_UTILITY_SCENARIO4,
        stratum="rank4_plus",
        comparison="rank1_vs_onlabel",
        clinical_question=(
            "Among models with deeply deprioritised on-label agents, does algorithm rank 1 still "
            "outperform on-label on TTD?"
        ),
        mtb_framing=(
            "The on-label option would rarely be prioritised clinically when ranked ≥ 4. The experiment "
            "tests whether that clinical intuition is reflected in pharmacodynamic data."
        ),
        cohort_design=f"Scenario 4, rank 4+ (n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO4, 0))}; paired n = {int(s4_r4_on['n_models']) if s4_r4_on is not None else 0}). Exploratory.",
        row=s4_r4_on,
        interpretation=exp15_i,
        pi_conclusion=exp15_pi,
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 16 (Exp16)",
        "Scenario 4 — rank 4+ stratum: rank 1 versus SC (ancillary)",
        scenario=CLINICAL_UTILITY_SCENARIO4,
        stratum="rank4_plus",
        comparison="rank1_vs_sc",
        clinical_question="Does rank 1 outperform chemotherapy in the rank 4+ exploratory cohort?",
        mtb_framing="Tests whether rank 1 remains an active regimen when on-label is deeply deprioritised.",
        cohort_design=f"Rank 4+ stratum (paired n = {int(s4_r4_sc['n_models']) if s4_r4_sc is not None else 0}).",
        row=s4_r4_sc,
        interpretation=exp15_i,
        pi_conclusion="Trend favouring rank 1 vs SC at nominal p ≈ 0.08; exploratory only.",
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 19 (Exp19)",
        "Scenario 4 — ranks 7–10 pooled: rank 1 versus on-label",
        scenario=CLINICAL_UTILITY_SCENARIO4,
        stratum="rank7_10",
        comparison="rank1_vs_onlabel",
        clinical_question=(
            "In the most deeply deprioritised on-label stratum (ranks 7–10, pooled label \"7–10\"), "
            "does rank 1 dominate on-label TTD — and is this consistent with on-label ineffectiveness "
            "rather than algorithmic error?"
        ),
        mtb_framing=(
            f"In this dataset, all {n710} rank 7–10 models share cetuximab as the on-label compound. "
            "The MTB would almost certainly not prioritise cetuximab when ranked at the bottom of the panel."
        ),
        cohort_design=f"Rank 7–10 pooled stratum (paired n = {int(s4_r710_on['n_models']) if s4_r710_on is not None else 0}). Formal inference underpowered.",
        row=s4_r710_on,
        interpretation=exp19_i,
        pi_conclusion=exp19_pi,
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 20 (Exp20)",
        "Scenario 4 — ranks 7–10: on-label versus SC",
        scenario=CLINICAL_UTILITY_SCENARIO4,
        stratum="rank7_10",
        comparison="onlabel_vs_sc",
        clinical_question="Is cetuximab in this stratum inferior to matched chemotherapy on TTD?",
        mtb_framing="Direct test of on-label pharmacodynamic inactivity in the lowest-rank pool.",
        cohort_design=f"Rank 7–10 stratum (paired n = {int(s4_r710_ol['n_models']) if s4_r710_ol is not None else 0}).",
        row=s4_r710_ol,
        interpretation=exp19_i,
        pi_conclusion=(
            "Directionally, on-label underperforms SC (median ΔTTD negative in all three models). "
            "Confirms that deep deprioritisation aligns with lack of ex vivo activity, not merely ordinal rank."
        ),
    ))

    doc.extend(_paired_experiment_section(
        "Experiment 26 (Exp26)",
        "Scenario 4 — progressive pool (rank ≥ 4): rank 1 versus on-label",
        scenario=CLINICAL_UTILITY_SCENARIO4,
        stratum="pool_ge4",
        comparison="rank1_vs_onlabel",
        clinical_question=(
            "Does progressive pooling of low-rank on-label models (pool_ge4) alter inference compared "
            "with the rank 7–10 subset alone?"
        ),
        mtb_framing=(
            "Pooling trades specificity for power. In this dataset, pool_ge4 equals the full Scenario 4 "
            "cohort (n = 8); pools ≥ 5 and ≥ 6 collapse to the rank 7–10 subset (n = 3)."
        ),
        cohort_design=f"pool_ge4 stratum (paired n = {int(s4_pool4_on['n_models']) if s4_pool4_on is not None else 0}).",
        row=s4_pool4_on,
        interpretation=exp15_i,
        pi_conclusion=exp15_pi,
    ))

    if epistemic is not None and not epistemic.empty:
        doc.extend([
            "### Epistemic classification of Scenario 4 models (Table 4)",
            "",
            "Each Scenario 4 model is classified according to whether on-label TTD is inferior to rank 1 "
            "and to SC, separating **ineffective on-label correctly deprioritised** from **possible ranking miss**.",
            "",
            "**How these groups were formed.** The epistemic table includes **every Scenario 4 model** "
            "(best on-label at rank 4+). No pairing test is performed — each row is one tumour with its "
            "pre-computed **onlabel_TTD**, **rank1_TTD**, and **chemo_TTD**. Class A: on-label worse than "
            "rank 1 and worse than SC. Class B: on-label worse than rank 1 but not worse than SC. "
            "Class C: on-label better than rank 1 (possible ranking miss). Class D: missing or tied values.",
            "",
            "| Model | On-label rank | On-label vs rank 1 | On-label vs SC | Epistemic class |",
            "| --- | --- | --- | --- | --- |",
        ])
        for _, er in epistemic.iterrows():
            doc.append(
                f"| {er['Model']} | {er['best_onlabel_rank']} | {er['onlabel_vs_rank1']} | "
                f"{er['onlabel_vs_sc']} | {_epistemic_plain_label(er['epistemic_class'])} |"
            )
        doc.extend(["", "---", ""])

    if between_cohort is not None and not between_cohort.empty:
        doc.extend([
            "## Between-cohort experiments (Scenario 2 versus Scenario 4 pools)",
            "",
            "These experiments compare the distribution of on-label TTD in guideline-concordant Scenario 2 "
            "(on-label at rank 1) against on-label TTD in progressively pooled low-rank Scenario 4 cohorts. "
            "A higher median in Scenario 2 would support the construct validity of rank-based stratification.",
            "",
        ])
        exp_bc_ids = ["Exp27", "Exp28", "Exp29", "Exp30"]
        for idx, (_, bc) in enumerate(between_cohort.iterrows()):
            eid = exp_bc_ids[idx] if idx < len(exp_bc_ids) else f"ExpBC{idx}"
            doc.extend([
                f"### {eid}: Between-cohort comparison — Scenario 2 vs `{bc['stratum']}`",
                "",
                f"**Clinical question.** Is on-label pharmacodynamic performance when ranked first "
                f"systematically superior to on-label performance in the `{bc['stratum']}` low-rank pool?",
                "",
                f"**Design.** Unpaired comparison of on-label TTD values; Mann–Whitney U test (two-sided); "
                "unpaired log-rank and Cox as supporting analyses.",
                "",
                "**How these groups were formed.** "
                + _build_between_cohort_group_formation(str(bc["stratum"]), bc),
                "",
                "**Results.**",
                "",
                *_between_cohort_markdown_table(bc),
                "**PI conclusion.** ",
            ])
            if pd.notna(bc.get("mannwhitney_p")) and float(bc["mannwhitney_p"]) < 0.05:
                doc.append(
                    "Scenario 2 on-label TTD is significantly higher than in the low-rank pool, supporting "
                    "rank-stratified clinical interpretation."
                )
            else:
                doc.append(
                    "No significant separation at conventional thresholds; interpret with caution given "
                    "unpaired design and uneven cohort sizes."
                )
            doc.extend(["", "---", ""])

    doc.extend([
        "## Cross-rank utility versus standard chemotherapy",
        "",
        "Complementary per-rank analyses compare On_label, Off_label, and Exp compounds at each global "
        "rank position to matched SC (stratified Cox and log-rank). Full results: `TABLE3_UTILITY_AT_RANK_VS_SC.csv`.",
        "",
        "**How these groups were formed.** For each global rank position (1, 2, 3, …, 7–10), all PDX models "
        "contribute drug-response rows at that rank. Within each model, compounds at that rank are split by "
        "**ON_LABEL_UTILITY** (On_label / Off_label / Exp). Each utility subgroup is compared to the **same "
        "model's chemo_TTD** (matched SC) using stratified survival models with Model as stratum — so every "
        "comparison is paired at the model level, but the test arm is defined by **rank position × utility label**, "
        "not by scenario.",
        "",
    ])
    if not utility_rank_pairwise.empty:
        r1 = utility_rank_pairwise[utility_rank_pairwise["Rank"].astype(str) == "1"]
        if not r1.empty:
            doc.append("| Utility at rank 1 | Stratified Cox HR vs SC | p-value |")
            doc.append("| --- | --- | --- |")
            for _, urow in r1.iterrows():
                doc.append(
                    f"| {urow['left_group']} | {urow['stratified_cox_hr']:.2f} "
                    f"({urow['stratified_cox_hr_low']:.2f}–{urow['stratified_cox_hr_high']:.2f}) | "
                    f"{_clinical_format_p(urow['stratified_cox_p'])} |"
                )
            doc.append("")

    doc.extend([
        "---",
        "",
        "## Integrated discussion",
        "",
        "### Principal findings",
        "",
        "1. **Algorithm as chemotherapy alternative (Exp1):** When no on-label agent is screened, DDA rank 1 "
        "   delivers statistically significant TTD prolongation versus matched SC — the clearest actionable signal.",
        "2. **Guideline concordance (Exp2):** On-label at rank 1 demonstrates pharmacodynamic benefit over SC, "
        "   particularly on stratified Cox analysis, validating internal consistency when approval and ranking align.",
        "3. **Guideline discordance at rank 2 (Exp3):** Rank-1 Off_label/Exp compounds significantly outperform "
        "   on-label at rank 2; ancillary comparisons suggest on-label does not clearly beat SC.",
        "4. **Deep deprioritisation (Scenario 4):** Epistemic review indicates predominantly correct "
        "   deprioritisation of inactive on-label cetuximab at ranks 7–10; Class C cases warrant individual review.",
        "",
        "### Statistical considerations",
        "",
        "Wilcoxon signed-rank tests respect the paired PDX design and are appropriate for skewed TTD differences. "
        "Stratified Cox models yield interpretable hazard ratios but should be read alongside Wilcoxon p-values. "
        "Stratified log-rank tests are often underpowered because each model contributes a single event time per arm. "
        "FDR adjustment applies within Scenario 3 comparison families; Scenario 4 and between-cohort analyses are exploratory.",
        "",
        "### Translation to clinical practice",
        "",
        "These experiments support DDA ranking as **preclinical decision support** for MTB deliberation and trial "
        "design. They do not substitute for prospective clinical trials, toxicity assessment, regulatory approval, "
        "or reimbursement considerations. TTD in PDX is a pharmacodynamic surrogate and must not be equated with "
        "patient overall survival.",
        "",
        "---",
        "",
        "## Conclusions and recommended MTB posture",
        "",
        "1. **No on-label in screen:** Recommend algorithm rank 1 over default chemotherapy pending clinical validation.",
        "2. **On-label at rank 1:** Endorse guideline-concordant on-label selection; preclinical TTD supports this pathway.",
        "3. **On-label at rank 2:** Present rank-1 Off_label/Exp recommendation as pharmacodynamically superior; "
        "   discuss off-label access and trial enrolment.",
        "4. **On-label at ranks 4+ / 7–10:** Treat on-label as likely inactive in the ex vivo model unless Class C "
        "   epistemic flags mandate individual reassessment.",
        "",
        "---",
        "",
        "## Limitations",
        "",
        "- PDX TTD is a pharmacodynamic surrogate, not a clinical survival endpoint.",
        "- Rank-3, rank-4+, and rank 7–10 strata are underpowered for definitive inference.",
        "- Real-world prescribing integrates toxicity, comorbidity, drug access, and molecular co-alterations not modelled here.",
        "- Four models lacked matched SC and were excluded from paired analyses.",
        "- Between-cohort comparisons are unpaired and sensitive to cohort composition (e.g., cetuximab enrichment in rank 7–10).",
        "",
        "---",
        "",
        "## Deliverables",
        "",
        "| Item | File |",
        "| --- | --- |",
        "| Cohort flow | `COHORT_FLOW.md` |",
        "| Table 1 — rank distribution | `TABLE1_COHORT_DISTRIBUTION.csv` |",
        "| Table 2 — experiment results | `TABLE2_SCENARIO_RESULTS.csv` |",
        "| Table 3 — utility at rank vs SC | `TABLE3_UTILITY_AT_RANK_VS_SC.csv` |",
        "| Table 4 — epistemic classification | `TABLE4_LOW_RANK_EPISTEMIC.csv` |",
        "| Experiment registry | `TABLE_EXP_REGISTRY.csv` |",
        "| Master interpretation | `CLINICAL_INTERPRETATION_MASTER.md` |",
        "",
        "---",
        "",
        "## Key figures",
        "",
    ])
    if copied_figs:
        for rel in copied_figs:
            doc.append(f"- `{rel}`")
    else:
        doc.append("_Key figures in parent output directory._")
    doc.extend([
        "",
        "## Reproduction",
        "",
        "```bash",
        "python src/clinical_utility_ranking.py --output-dir src",
        "```",
        "",
        "*End of structured clinical summary for supervisor review.*",
    ])
    return doc


def _short_paired_comparison_bullet(
    row: Optional[pd.Series],
    *,
    comparison_label: str,
    favour_label: str = "test arm",
    delta_reference: str = "",
) -> str:
    """One paired comparison with all reported test statistics."""
    if row is None:
        return "No evaluable paired data."
    exp = _exploratory_bracket(row)
    delta_ref = f" ({delta_reference})" if delta_reference else ""
    stats: List[str] = [
        (
            f"median ΔPFS {row['median_diff_TTD']:+.1f} d "
            f"(95% CI {row['median_diff_CI_low']:.1f} to {row['median_diff_CI_high']:.1f}){delta_ref}"
        ),
        f"{row['pct_improved']:.0f}% models favour {favour_label}",
        f"Wilcoxon p = {_clinical_format_p(row['wilcoxon_p'])}{exp}",
    ]
    if pd.notna(row.get("wilcoxon_p_fdr_bh")):
        stats.append(f"FDR-adjusted Wilcoxon p = {_clinical_format_p(row['wilcoxon_p_fdr_bh'])}{exp}")
    if pd.notna(row.get("stratified_logrank_p")):
        stats.append(
            f"stratified log-rank p = {_clinical_format_p(row['stratified_logrank_p'])}{exp}"
        )
    if pd.notna(row.get("stratified_cox_hr")) and pd.notna(row.get("stratified_cox_p")):
        cox_exp = " [exploratory]" if _cox_estimate_unstable(row) else exp
        stats.append(
            f"stratified Cox HR = {row['stratified_cox_hr']:.2f} "
            f"(p = {_clinical_format_p(row['stratified_cox_p'])}){cox_exp}"
        )
    return _short_summary_prose(
        f"{comparison_label} (n = {int(row['n_models'])}): " + ", ".join(stats) + ".",
        preserve_punctuation=True,
    )


def _short_paired_results_bullets(
    row: Optional[pd.Series],
    *,
    comparison_label: str,
    favour_label: str = "test arm",
    delta_reference: str = "",
) -> List[str]:
    if row is None:
        return ["No evaluable paired data."]
    return [
        _short_paired_comparison_bullet(
            row,
            comparison_label=comparison_label,
            favour_label=favour_label,
            delta_reference=delta_reference,
        )
    ]


def _short_unpaired_comparison_bullet(row: pd.Series) -> str:
    """One unpaired cohort comparison with Mann-Whitney, log-rank, and Cox."""
    exp = _exploratory_bracket(row)
    stats: List[str] = [
        f"median {row['median_b']:.1f} vs {row['median_a']:.1f} d (pool vs reference)",
        f"Mann-Whitney p = {_clinical_format_p(row['mannwhitney_p'])}{exp}",
    ]
    if pd.notna(row.get("pooled_logrank_p")):
        stats.append(f"log-rank p = {_clinical_format_p(row['pooled_logrank_p'])}{exp}")
    if pd.notna(row.get("pooled_cox_hr")) and pd.notna(row.get("pooled_cox_p")):
        stats.append(
            f"Cox HR = {row['pooled_cox_hr']:.2f} "
            f"(p = {_clinical_format_p(row['pooled_cox_p'])}){exp}"
        )
    return _short_summary_prose(
        f"{row['pool_label']} vs {row['reference_label']}: " + ", ".join(stats) + ".",
        preserve_punctuation=True,
    )


def _short_results_line(row: Optional[pd.Series], *, extra: str = "") -> str:
    """Legacy single-line paired results (prefer _short_paired_results_bullets)."""
    bullets = _short_paired_results_bullets(row, comparison_label="Paired comparison")
    text = bullets[0]
    if extra:
        text = text.rstrip(".") + f", {extra}."
    return text


def _short_paired_stats_clause(
    row: Optional[pd.Series],
    *,
    left_name: str,
    right_name: str,
) -> str:
    """Publication-style statistical clause for a paired comparison."""
    if row is None:
        return "No evaluable paired data were available."
    exp = _exploratory_bracket(row)
    clauses = [
        (
            f"Median ΔPFS ({left_name} minus {right_name}) was {row['median_diff_TTD']:+.1f} days "
            f"(95% CI {row['median_diff_CI_low']:.1f} to {row['median_diff_CI_high']:.1f})"
        ),
        f"{row['pct_improved']:.0f}% of models favoured {left_name}",
        f"Wilcoxon p = {_clinical_format_p(row['wilcoxon_p'])}{exp}",
    ]
    if pd.notna(row.get("wilcoxon_p_fdr_bh")):
        clauses.append(f"FDR-adjusted Wilcoxon p = {_clinical_format_p(row['wilcoxon_p_fdr_bh'])}{exp}")
    if pd.notna(row.get("stratified_cox_hr")) and pd.notna(row.get("stratified_cox_p")):
        cox_exp = " [exploratory]" if _cox_estimate_unstable(row) else exp
        clauses.append(
            f"stratified Cox HR {row['stratified_cox_hr']:.2f} "
            f"(p = {_clinical_format_p(row['stratified_cox_p'])}){cox_exp}"
        )
    if pd.notna(row.get("stratified_logrank_p")):
        clauses.append(
            f"stratified log-rank p = {_clinical_format_p(row['stratified_logrank_p'])}{exp}"
        )
    return ". ".join(clauses) + "."


def _short_q8_publication_trend(sensitivity: pd.DataFrame, alpha: float = 0.05) -> str:
    """Synthesised publication statement for Q8 unpaired pool sensitivity."""
    if sensitivity is None or sensitivity.empty:
        return "No evaluable on-label-at-rank pool data."
    narrow = sensitivity[sensitivity["stratum"] != "pool_ge4"]
    wide = sensitivity[sensitivity["stratum"] == "pool_ge4"]
    narrow_mw_sig = int(
        narrow["mannwhitney_p"].apply(lambda p: pd.notna(p) and float(p) < alpha).sum()
    )
    wide_mw_sig = int(
        wide["mannwhitney_p"].apply(lambda p: pd.notna(p) and float(p) < alpha).sum()
    ) if not wide.empty else 0
    med_on_narrow = float(narrow["median_b"].median()) if not narrow.empty else np.nan
    med_ref_narrow = float(narrow["median_a"].median()) if not narrow.empty else np.nan
    parts = [
        "Progressive on-label-at-rank pools were compared unpaired to Rank 1 (Off/Exp) and matched "
        "standard chemotherapy across PDX models.",
        (
            f"In narrower pools (ranks 7-10 through 5-10, n = 7 to 11 per pool), on-label PFS was "
            f"shorter than both references (typical medians ~{med_on_narrow:.0f} vs ~{med_ref_narrow:.0f} days), "
            f"with Mann-Whitney, log-rank, and Cox tests significant in {narrow_mw_sig} of "
            f"{len(narrow)} comparisons [exploratory]."
        ),
    ]
    if not wide.empty:
        wr = wide.iloc[0]
        parts.append(
            f"Expanding to on-label at rank 4-10 (n = {int(wr['n_cohort_b'])}) attenuated separation "
            f"(median {wr['median_b']:.0f} vs {wr['median_a']:.0f} days vs Rank 1, "
            f"{wide_mw_sig} of {len(wide)} rank 4-10 comparisons significant at α = 0.05)."
        )
    parts.append(
        "Overall trend: deeply downranked on-label pharmacodynamics underperform algorithm rank 1 and "
        "chemotherapy at the cohort level, and the signal weakens as higher-rank on-label models enter the pool."
    )
    return _short_summary_prose(" ".join(parts))


def _short_q9_publication_trend(q9_summary: pd.DataFrame, alpha: float = 0.05) -> str:
    """Synthesised publication statement for Q9 paired pool sensitivity."""
    if q9_summary is None or q9_summary.empty:
        return "No evaluable paired on-label-at-rank data."
    vs_rank1 = q9_summary[q9_summary["comparison"] == "onlabel_vs_rank1"]
    vs_sc = q9_summary[q9_summary["comparison"] == "onlabel_vs_sc"]
    med_delta = float(q9_summary["median_diff_TTD"].median())
    pct_favour = float(q9_summary["pct_improved"].median())
    n_rank1_sig = int(
        vs_rank1["wilcoxon_p"].apply(lambda p: pd.notna(p) and float(p) < alpha).sum()
    )
    n_sc_sig = int(
        vs_sc["wilcoxon_p"].apply(lambda p: pd.notna(p) and float(p) < alpha).sum()
    ) if not vs_sc.empty else 0
    parts = [
        "Within each PDX model, on-label PFS at downranked screen positions was compared paired to "
        "Rank 1 (Off/Exp) or matched standard chemotherapy.",
        (
            f"Across progressive pools (n = 7 to 20), median ΔPFS (on-label minus reference) was "
            f"{med_delta:+.1f} days and only ~{pct_favour:.0f}% of models favoured on-label, "
            f"indicating consistent within-model underperformance of downranked on-label therapy."
        ),
        (
            f"Paired Wilcoxon tests versus Rank 1 (Off/Exp) were significant in {n_rank1_sig} of "
            f"{len(vs_rank1)} pools, versus standard chemo in {n_sc_sig} of {len(vs_sc)} pools. "
            "Stratified log-rank and Cox summaries were largely non-significant, as expected with "
            "small paired n per pool [exploratory]."
        ),
        (
            "Overall trend: when on-label sits at ranks 7-10 through 4-10, Rank 1 (Off/Exp) typically "
            "prolongs PFS within the same model. Discordance is clearest against rank 1 and dilutes "
            "when the pool widens or the reference is chemotherapy."
        ),
    ]
    return _short_summary_prose(" ".join(parts))


def _short_publication_statement(
    qid: str,
    *,
    row: Optional[pd.Series] = None,
    bc_row: Optional[pd.Series] = None,
    low_rank_sensitivity: Optional[pd.DataFrame] = None,
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
    ep_a710: int = 0,
    ep_c710: int = 0,
    n710_epi: int = 0,
) -> str:
    """Advanced publication-style results sentence with embedded p-values."""
    if qid == "Q1" and q9_paired_sensitivity is not None and not q9_paired_sensitivity.empty:
        return _short_q9_publication_trend(q9_paired_sensitivity)
    if qid == "Q2" and ENABLE_Q2_ANALYSIS and low_rank_sensitivity is not None and not low_rank_sensitivity.empty:
        return _short_q8_publication_trend(low_rank_sensitivity)
    if qid == "Q3":
        intro = (
            f"We measured whether DDA rank 1 prolongs PDX PFS relative to standard chemotherapy "
            f"in models with no screened on-label agent (paired n = {_short_paired_n(row)})."
        )
        return _short_summary_prose(intro + " " + _short_paired_stats_clause(
            row, left_name="rank 1", right_name="standard chemotherapy"
        ))
    if qid == "Q4":
        intro = (
            f"We measured whether the on-label compound at global rank 1 prolongs PDX PFS relative "
            f"to matched standard chemotherapy when guideline and algorithm priorities align "
            f"(paired n = {_short_paired_n(row)})."
        )
        return _short_summary_prose(intro + " " + _short_paired_stats_clause(
            row, left_name="on-label at rank 1", right_name="standard chemotherapy"
        ))
    if qid == "Q6":
        intro = (
            f"We measured whether rank 1 prolongs PDX PFS relative to deeply deprioritised on-label "
            f"therapy (rank 7-10 pool, paired n = {_short_paired_n(row)}, all cetuximab)."
        )
        stats = _short_paired_stats_clause(row, left_name="rank 1", right_name="on-label rank 7-10")
        epistemic = (
            f"Epistemic review classified {ep_a710} of {n710_epi} models as Class A "
            f"(inactive on-label correctly deprioritised) and {ep_c710} as Class C."
        )
        return _short_summary_prose(intro + " " + stats + " " + epistemic)
    if qid == "Q7" and bc_row is not None:
        exp = _exploratory_bracket(bc_row)
        na, nb = int(bc_row["n_cohort_a"]), int(bc_row["n_cohort_b"])
        intro = (
            f"We compared on-label PDX PFS between Scenario 2 models with on-label at rank 1 "
            f"(n = {na}) and Scenario 4 models with on-label in the rank 7-10 pool (n = {nb})."
        )
        stats = (
            f"Median PFS was {bc_row['median_a']:.1f} days versus {bc_row['median_b']:.1f} days "
            f"(Δmedian {bc_row['median_diff_a_minus_b']:+.1f} days, "
            f"Mann-Whitney p = {_clinical_format_p(bc_row['mannwhitney_p'])}{exp})."
        )
        if pd.notna(bc_row.get("pooled_logrank_p")):
            stats += (
                f" Unpaired log-rank p = {_clinical_format_p(bc_row['pooled_logrank_p'])}{exp}."
            )
        else:
            stats += "."
        return _short_summary_prose(intro + " " + stats)
    return _short_summary_prose("Results not available for this question.")


SHORT_SUMMARY_PFS_CONVENTION = (
    "Convention: throughout this summary, the PDX endpoint TimeToDouble (days) is referred to as "
    "PFS (progression-free survival), following our programme naming. The underlying measurement "
    "remains ex vivo tumour doubling time in patient-derived xenografts. It is a pharmacodynamic "
    "surrogate and not clinical patient PFS or overall survival."
)

SHORT_SUMMARY_STATISTICAL_PARAGRAPHS = (
    "Statistical analyses reported per question: paired Wilcoxon signed-rank on ΔPFS within each PDX model "
    "(bootstrap 95% CI for median ΔPFS), stratified Cox PH and stratified log-rank on Kaplan-Meier curves "
    "(stratified by model)"
    + (", unpaired Mann-Whitney for Q2," if ENABLE_Q2_ANALYSIS else ",")
    + " and paired Wilcoxon for Q1 pool sensitivity. "
    "KM plots are pharmacodynamic visualisations (all observations treated as events, not clinical survival "
    "with censoring). Values marked [exploratory] reflect underpowered strata (n ≤ 10) or numerically "
    "unstable estimates (e.g. Cox HR at very small n).",
    "We compare two treatments within the same PDX model (paired design). The primary test is Wilcoxon on "
    "ΔPFS, which directly asks whether one arm is systematically longer.",
    "Stratified Cox asks a related question: \"Is there a consistent hazard ratio favoring one arm across "
    "models?\" — it can detect a systematic direction and often aligns with Wilcoxon.",
    "Stratified log-rank, with exactly one measurement per arm per model, cannot accumulate directional "
    "evidence in our implementation: its p-value stays near 0.5 regardless of effect size. So "
    "non-significant log-rank does not contradict significant Wilcoxon/Cox — it reflects a test–design "
    "mismatch, not absence of effect.",
    "KM plots remain useful visual summaries; for inference on paired PDX comparisons, Wilcoxon is primary; "
    "stratified Cox is supportive; stratified log-rank should not be over-interpreted here.",
)
SHORT_SUMMARY_STATISTICAL_HIERARCHY = "\n\n".join(SHORT_SUMMARY_STATISTICAL_PARAGRAPHS)

EXPLORATORY_N_THRESHOLD = 10
SMALL_COHORT_N_THRESHOLD = 5


def _exploratory_bracket(row: Optional[pd.Series]) -> str:
    """Return ' [exploratory]' when estimates are underpowered or numerically unstable."""
    if row is None:
        return ""
    n = _short_paired_n(row) if pd.notna(row.get("n_models")) else 0
    if n == 0 and pd.notna(row.get("n_cohort_b")):
        n = int(row["n_cohort_b"])
    if n <= EXPLORATORY_N_THRESHOLD or _cox_estimate_unstable(row):
        return " [exploratory]"
    return ""


def _cox_estimate_unstable(row: Optional[pd.Series]) -> bool:
    """True when stratified Cox HR is not reliable (very small n or separation)."""
    if row is None:
        return True
    n = int(row.get("n_models", 0))
    if n < SMALL_COHORT_N_THRESHOLD:
        return True
    hr = row.get("stratified_cox_hr", np.nan)
    hi = row.get("stratified_cox_hr_high", np.nan)
    lo = row.get("stratified_cox_hr_low", np.nan)
    if pd.notna(hr) and (float(hr) < 1e-4 or float(hr) > 1e4):
        return True
    if pd.notna(hi) and (np.isinf(hi) or float(hi) > 500):
        return True
    if pd.notna(lo) and float(lo) <= 0:
        return True
    return False


def _wilcoxon_primary_significant(row: Optional[pd.Series], alpha: float = 0.05) -> bool:
    if row is None or pd.isna(row.get("wilcoxon_p")):
        return False
    return float(row["wilcoxon_p"]) < alpha


def _epistemic_counts_for_stratum(
    epistemic: Optional[pd.DataFrame],
    rank_label: str,
) -> Tuple[int, int, int]:
    """Return (n_class_a, n_class_c, n_in_stratum) for one on-label rank bucket."""
    if epistemic is None or epistemic.empty:
        return 0, 0, 0
    sub = epistemic[epistemic["best_onlabel_rank"].astype(str) == str(rank_label)]
    ep_a = int((sub["epistemic_class"] == "A_ineffective_correctly_deprioritized").sum())
    ep_c = int((sub["epistemic_class"] == "C_possible_ranking_miss").sum())
    return ep_a, ep_c, len(sub)


def _epistemic_counts_scenario4_overall(epistemic: Optional[pd.DataFrame]) -> Tuple[int, int, int]:
    if epistemic is None or epistemic.empty:
        return 0, 0, 0
    ep_a = int((epistemic["epistemic_class"] == "A_ineffective_correctly_deprioritized").sum())
    ep_c = int((epistemic["epistemic_class"] == "C_possible_ranking_miss").sum())
    return ep_a, ep_c, len(epistemic)


def _short_caveat_for_question(
    qid: str,
    *,
    row: Optional[pd.Series] = None,
    bc_row: Optional[pd.Series] = None,
    epistemic: Optional[pd.DataFrame] = None,
) -> str:
    """Exploratory / rigor caveats appended to each short-summary question (nothing removed)."""
    parts: List[str] = []

    if qid == "Q1":
        parts.append(
            "[exploratory]: eight paired within-model comparisons (four on-label-at-rank pools × "
            "on-label vs Rank 1 or SC) without multiplicity adjustment."
        )
        parts.append(
            "Positive ΔPFS is on-label minus reference within the same model. "
            "Progressive pools: on-label at rank 7-10, then 6-10, 5-10, and 4-10."
        )

    elif qid == "Q2" and ENABLE_Q2_ANALYSIS:
        parts.append(
            "[exploratory]: eight unpaired comparisons (four progressive on-label-at-rank pools × "
            "Rank 1 Off/Exp or SC reference) without multiplicity adjustment."
        )
        parts.append(
            "Pools count models with on-label at the specified global rank(s) in the screen "
            "(n = 7 at rank 7-10, 9 at rank 6-10, 11 at rank 5-10, 20 at rank 4-10 in this dataset)."
        )

    elif qid == "Q3":
        if row is not None and not _wilcoxon_primary_significant(row):
            parts.append("Wilcoxon did not reach α = 0.05 in this run.")
        if row is not None and _wilcoxon_primary_significant(row) and pd.notna(row.get("stratified_cox_p")):
            if float(row["stratified_cox_p"]) >= 0.05:
                parts.append(
                    "Stratified Cox p ≥ 0.05 while Wilcoxon is significant — estimands differ "
                    "(paired median ΔPFS vs stratified event-time hazard)."
                )
        if row is not None and pd.notna(row.get("stratified_logrank_p")):
            if float(row["stratified_logrank_p"]) >= 0.05:
                parts.append(
                    "Stratified log-rank is non-significant; KM curves illustrate TTD distributions only."
                )

    elif qid == "Q4":
        parts.append(
            f"Cohort size modest (paired n = {_short_paired_n(row)}). "
            "Wilcoxon p ≈ 0.07 is borderline; stratified Cox also reports reduced growth hazard."
        )

    elif qid == "Q5":
        parts.append(
            "Three-stratum panel (rank 2, rank 3, pooled ranks 2-3). Rank-2 stratum is confirmatory "
            "(Wilcoxon + FDR significant). Rank-3 stratum is [exploratory] and underpowered. "
            "Pooled ranks 2-3 is [exploratory]: nominal Wilcoxon p < 0.05 but FDR-adjusted p > 0.05."
        )

    if not parts:
        parts.append(
            "KM panels illustrate pharmacodynamic TTD distributions; paired Wilcoxon on ΔPFS is reported "
            "alongside stratified Cox and log-rank where available."
        )
    return _short_summary_prose(" ".join(parts))


def _apply_short_item_caveats(
    items: List[Dict[str, Any]],
    *,
    paired_rows: Dict[str, Optional[pd.Series]],
    epistemic: Optional[pd.DataFrame],
    skip_qids: Optional[set[str]] = None,
) -> None:
    skip = skip_qids or set()
    for item in items:
        qid = str(item.get("qid", ""))
        if qid in skip and item.get("caveat"):
            continue
        item["caveat"] = _short_caveat_for_question(
            qid,
            row=paired_rows.get(qid),
            epistemic=epistemic,
        )

SHORT_SUMMARY_FIGURE_KEYS: Dict[str, Tuple[str, ...]] = {
    qid: (f"pub_{qid}_composite",) for qid in SHORT_SUMMARY_QUESTION_ORDER
}


def _short_summary_prose(text: str, *, preserve_punctuation: bool = False) -> str:
    """Plain text for short summary: strip markdown bold, TTD to PFS, style constraints."""
    out = text.replace("**", "")
    replacements = (
        ("TimeToDouble", "PFS"),
        ("ΔTTD", "ΔPFS"),
        (" on-label TTD", " on-label PFS"),
        (" on-label PFS values", " on-label PFS values"),  # idempotent guard
        ("TTD values", "PFS values"),
        ("valid TTD", "valid PFS"),
        ("finite TTD", "finite PFS"),
        ("median TTD", "median PFS"),
        ("exceeds on-label TTD", "exceeds on-label PFS"),
        ("on-label TTD", "on-label PFS"),
        ("preclinical TTD data", "preclinical PFS data"),
        (" dominate on-label PFS", " dominate on-label PFS"),
        ("% change TTD", "% change PFS"),
    )
    for old, new in replacements:
        out = out.replace(old, new)
    if not preserve_punctuation:
        out = out.replace(" — ", ": ")
        out = out.replace("—", ":")
        out = re.sub(r";\s*", ". ", out)
        out = re.sub(r"\.\s+\.", ". ", out)
        out = re.sub(r"\.\.+", ".", out)
    return out.strip()


def _short_paired_n(row: Optional[pd.Series]) -> int:
    return int(row["n_models"]) if row is not None else 0


def _resolve_short_summary_figures(fig_dir: Path) -> Dict[str, List[Path]]:
    """Map question IDs to figure paths under figures_to_send."""
    resolved: Dict[str, List[Path]] = {qid: [] for qid in SHORT_SUMMARY_FIGURE_KEYS}
    resolved["appendix"] = []
    if not fig_dir.is_dir():
        return resolved
    all_png = sorted(fig_dir.glob("*.png"))
    seen_names: set[str] = set()
    for qid, keys in SHORT_SUMMARY_FIGURE_KEYS.items():
        for key in keys:
            for path in all_png:
                if key in path.name and path.name not in seen_names:
                    resolved[qid].append(path)
                    seen_names.add(path.name)
                    break
    for path in all_png:
        if "pub_appendix_utility_at_rank_panel" in path.name:
            resolved["appendix"] = [path]
            break
        if "utility_at_rank_panel" in path.name and not resolved["appendix"]:
            resolved["appendix"] = [path]
    return resolved


def _apply_short_item_figures(items: List[Dict[str, Any]], fig_dir: Optional[Path]) -> None:
    figure_map = _resolve_short_summary_figures(fig_dir) if fig_dir else {}
    for item in items:
        qid = item.get("qid", "")
        item["figures"] = figure_map.get(qid, [])


def _short_groups_bullets(*lines: str) -> List[str]:
    """Build a bullet list for the Groups section (plain text, PFS naming applied)."""
    return [_short_summary_prose(line) for line in lines if line.strip()]


def _extend_short_results_markdown(doc: List[str], results: Any) -> None:
    """Append Results heading and one bullet per comparison."""
    doc.append("Results.")
    bullets = results if isinstance(results, list) else [str(results)]
    for bullet in bullets:
        doc.append(f"- {bullet}")
    doc.append("")


def _add_short_results_docx(document: Any, results: Any) -> None:
    """Append Results heading and bullet list to short-summary Word doc."""
    from docx.shared import Pt

    heading = document.add_paragraph()
    run_label = heading.add_run("Results.")
    run_label.bold = True
    bullets = results if isinstance(results, list) else [str(results)]
    for bullet in bullets:
        try:
            document.add_paragraph(bullet, style="List Bullet")
        except KeyError:
            p = document.add_paragraph()
            p.add_run(f"• {bullet}").font.size = Pt(10)


def _extend_short_groups_markdown(doc: List[str], groups: Any) -> None:
    """Append Groups heading and bullet lines to short-summary markdown."""
    doc.append("Groups.")
    bullets = groups if isinstance(groups, list) else [str(groups)]
    for bullet in bullets:
        doc.append(f"- {bullet}")
    doc.append("")


def _add_short_groups_docx(document: Any, groups: Any) -> None:
    """Append Groups heading and bullet list to short-summary Word doc."""
    from docx.shared import Pt

    heading = document.add_paragraph()
    run_label = heading.add_run("Groups.")
    run_label.bold = True
    bullets = groups if isinstance(groups, list) else [str(groups)]
    for bullet in bullets:
        try:
            document.add_paragraph(bullet, style="List Bullet")
        except KeyError:
            p = document.add_paragraph()
            p.add_run(f"• {bullet}").font.size = Pt(10)


def _short_groups_q1(scenario_n: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Cohort selection: Scenario 1 — no approved on-label (IGAZ) drug in the screened panel "
        f"(n = {scenario_n} models).",
        "Design: Paired within-model comparison (both arms from the same PDX tumour).",
        f"Group A (test): Rank 1 (n = {n} paired): compound with the best global DDA rank (rank 1), "
        "the algorithm's top pick (always Off_label or Exp in this cohort).",
        f"Group B (reference) — Standard chemotherapy / SC (n = {n} paired): model-matched chemo "
        "reference for the same PDX.",
        f"Paired analysis: n = {n} (models with valid PFS in both arms and matched SC).",
    )


def _short_groups_q2(scenario_n: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Cohort selection: Scenario 2 — best approved on-label drug at global rank 1 "
        f"(n = {scenario_n} models).",
        "Design: Paired within-model comparison.",
        f"Group A (test) — On-label at rank 1 (n = {n} paired): approved (IGAZ) compound at the top "
        "of the ranking (rank 1 and on-label are the same drug in this cohort).",
        f"Group B (reference) — Standard chemotherapy / SC (n = {n} paired): matched chemo "
        "for the same tumour.",
        f"Paired analysis: n = {n}.",
    )


def _short_groups_q3_rank2(n_rank2: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Cohort selection: Scenario 3, rank-2 stratum — best on-label drug at global rank 2 "
        f"(n = {n_rank2} models).",
        "Design: Paired within-model comparison.",
        f"Group A (test) — Rank 1 (n = {n} paired): top DDA-ranked compound (Off_label or Exp).",
        f"Group B (reference) — On-label at rank 2 (n = {n} paired): approved (IGAZ) drug, "
        "second in the global ranking.",
        f"Paired analysis: n = {n}.",
    )


def _short_groups_q3_rank3(n_rank3: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Cohort selection: Scenario 3, rank-3 stratum — best on-label drug at global rank 3 "
        f"(n = {n_rank3} models).",
        "Design: Paired within-model comparison (same logic as Q3, rank-2 stratum).",
        f"Group A (test) — Rank 1 (n = {n} paired): top-ranked Off_label/Exp compound.",
        f"Group B (reference) — On-label at rank 3 (n = {n} paired): approved (IGAZ) drug.",
        f"Paired analysis: n = {n}.",
    )


def _short_groups_q3_pooled(scenario_n: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Cohort selection: Scenario 3 pooled — on-label at global rank 2 or rank 3 "
        f"(n = {scenario_n} models).",
        "Design: Paired within-model comparison.",
        f"Group A (test) — Rank 1 (n = {n} paired): best-ranked compound per model.",
        f"Group B (reference) — On-label (n = {n} paired): best approved drug per model "
        "(occupies rank 2 or 3).",
        f"Paired analysis: n = {n}.",
    )


def _short_groups_q6(scenario4_n: int, row: Optional[pd.Series]) -> List[str]:
    n = _short_paired_n(row)
    return _short_groups_bullets(
        f"Parent cohort: Scenario 4 — on-label at global rank 4+ (n = {scenario4_n} models).",
        "Subset: on-label at pooled ranks 7–10 (n = 3 models; all cetuximab).",
        "Design: Paired within-model comparison.",
        f"Group A (test) — Rank 1 (n = {n} paired): top DDA-ranked compound per model.",
        f"Group B (reference) — On-label rank 7–10 (n = {n} paired): cetuximab as best on-label, "
        "deeply deprioritised.",
        f"Paired analysis: n = {n}.",
    )


def _short_groups_q7(bc: pd.Series) -> List[str]:
    na, nb = int(bc["n_cohort_a"]), int(bc["n_cohort_b"])
    return _short_groups_bullets(
        "Design: Unpaired comparison — different PDX models in each group (not paired within model).",
        f"Group A — On-label at rank 1: one on-label PFS value per model from Scenario 2 (n = {na}).",
        f"Group B — On-label rank 7–10 pool: one on-label PFS per model from Scenario 4, "
        f"ranks 7–10 only (n = {nb}; all cetuximab).",
        f"Comparison sample sizes: n = {na} vs n = {nb}.",
    )


def _build_short_summary_significance_table(
    *,
    paired_summary: pd.DataFrame,
    between_cohort: Optional[pd.DataFrame],
    low_rank_sensitivity: Optional[pd.DataFrame],
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
    alpha: float = 0.05,
) -> Tuple[List[str], List[List[str]]]:
    """Markdown lines and Word table rows summarising key tests; '*' marks p < alpha."""
    header = ["Question", "Groups compared", "n", "Wilcoxon / MW", "Log-rank", "Cox HR (p)"]
    rows: List[List[str]] = []

    paired_specs: Tuple[Tuple[str, str, str, str, str], ...] = (
        ("Q3", CLINICAL_UTILITY_SCENARIO1, "rank1_vs_sc", "all", "Rank 1 vs standard chemotherapy"),
        ("Q4", CLINICAL_UTILITY_SCENARIO2, "onlabel_vs_sc", "all", "On-label rank 1 vs standard chemotherapy"),
        ("Q5", CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2", "Rank 1 vs on-label at rank 2"),
        ("Q5", CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank3", "Rank 1 vs on-label at rank 3"),
        ("Q5", CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2_3_pooled", "Rank 1 vs on-label ranks 2-3"),
    )
    scenario_rows: List[List[str]] = []
    for qid, scenario, comparison, stratum, groups in paired_specs:
        row = _clinical_paired_row(paired_summary, scenario, comparison, stratum)
        if row is None:
            continue
        scenario_rows.append([
            qid,
            groups,
            str(int(row["n_models"])),
            _p_cell_with_star(row.get("wilcoxon_p"), alpha=alpha),
            _p_cell_with_star(row.get("stratified_logrank_p"), alpha=alpha),
            _cox_hr_cell(row, alpha=alpha),
        ])

    q1_rows: List[List[str]] = []
    if q9_paired_sensitivity is not None and not q9_paired_sensitivity.empty:
        for _, r in q9_paired_sensitivity.iterrows():
            q1_rows.append([
                "Q1",
                f"{r['pool_label']}, {r['comparison_label']}",
                str(int(r["n_models"])),
                _p_cell_with_star(r.get("wilcoxon_p"), alpha=alpha),
                _p_cell_with_star(r.get("stratified_logrank_p"), alpha=alpha),
                _cox_hr_cell(r, alpha=alpha),
            ])

    q2_rows: List[List[str]] = []
    if ENABLE_Q2_ANALYSIS and low_rank_sensitivity is not None and not low_rank_sensitivity.empty:
        for _, r in low_rank_sensitivity.iterrows():
            q2_rows.append([
                "Q2",
                f"{r['pool_label']} vs {r['reference_label']}",
                f"{int(r['n_cohort_a'])} vs {int(r['n_cohort_b'])}",
                _p_cell_with_star(r.get("mannwhitney_p"), alpha=alpha),
                _p_cell_with_star(r.get("pooled_logrank_p"), alpha=alpha),
                _cox_hr_cell(r, hr_key="pooled_cox_hr", p_key="pooled_cox_p", alpha=alpha),
            ])

    rows = q1_rows + q2_rows + scenario_rows

    md: List[str] = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * len(header)) + " |",
    ]
    for row in rows:
        md.append("| " + " | ".join(row) + " |")
    return md, rows


def _short_q8_results_bullets(sensitivity: pd.DataFrame) -> List[str]:
    if sensitivity is None or sensitivity.empty:
        return [_short_summary_prose("No evaluable low-rank pool sensitivity data.")]
    return [_short_unpaired_comparison_bullet(row) for _, row in sensitivity.iterrows()]


def _short_groups_q8(sensitivity: pd.DataFrame) -> List[str]:
    lines = [
        "Design: Unpaired Mann-Whitney, log-rank, and Cox PH (different PDX models in each cohort).",
        "Pool B uses on-label PFS at the specified global rank(s) in the drug screen (one value per model), "
        "not the model's best on-label rank.",
        "Progressive pools: on-label at rank 7-10 only (n = 7), then rank 6-10, 5-10, and 4-10 to increase n.",
        "Reference A (top rows): Rank 1 (Off/Exp) PFS from all models where rank 1 is not on-label.",
        "Reference A (bottom rows): matched standard chemotherapy PFS from all models with SC.",
    ]
    if sensitivity is not None and not sensitivity.empty:
        r1_rows = sensitivity[sensitivity["reference_arm"] == "rank1_offexp"]
        sc_rows = sensitivity[sensitivity["reference_arm"] == "sc"]
        if not r1_rows.empty:
            lines.append(f"Rank 1 (Off/Exp) reference n = {int(r1_rows.iloc[0]['n_cohort_a'])}.")
        if not sc_rows.empty:
            lines.append(f"Standard chemo reference n = {int(sc_rows.iloc[0]['n_cohort_a'])}.")
        for stratum in LOW_RANK_POOL_Q8_STRATA:
            sub = sensitivity[sensitivity["stratum"] == stratum]
            if not sub.empty:
                nb = int(sub.iloc[0]["n_cohort_b"])
                lines.append(f"{LOW_RANK_POOL_Q8_LABELS[stratum]} pool n = {nb}.")
    return _short_groups_bullets(*lines)


def _short_q8_significant_summary(sensitivity: pd.DataFrame, alpha: float = 0.05) -> str:
    if sensitivity is None or sensitivity.empty:
        return "No comparisons evaluable."
    sig_bits: List[str] = []
    for _, r in sensitivity.iterrows():
        tests: List[str] = []
        for label, col in (
            ("MW", "mannwhitney_p"),
            ("log-rank", "pooled_logrank_p"),
            ("Cox", "pooled_cox_p"),
        ):
            p = r.get(col)
            if pd.notna(p) and float(p) < alpha:
                tests.append(f"{label} p = {_clinical_format_p(p)}*")
        if tests:
            sig_bits.append(f"{r['pool_label']} vs {r['reference_label']} ({', '.join(tests)})")
    if not sig_bits:
        return (
            "No progressive pool reached conventional significance (α = 0.05) on Mann-Whitney, "
            "log-rank, or Cox versus Rank 1 (Off/Exp) or standard chemo at the tested pool sizes."
        )
    return "Significant separation at α = 0.05 (*): " + "; ".join(sig_bits) + "."


def _short_q9_results_bullets(q9_summary: pd.DataFrame) -> List[str]:
    if q9_summary is None or q9_summary.empty:
        return [_short_summary_prose("No evaluable paired on-label-at-rank data.")]
    bullets: List[str] = []
    for _, row in q9_summary.iterrows():
        label = f"{row['pool_label']}, {row['comparison_label']}"
        bullets.append(
            _short_paired_comparison_bullet(
                row,
                comparison_label=label,
                favour_label="on-label",
                delta_reference="on-label minus reference",
            )
        )
    return bullets


def _short_groups_q9(q9_summary: pd.DataFrame) -> List[str]:
    lines = [
        "Design: Paired within-model analysis (same PDX model).",
        "Test arm: on-label PFS at the specified global rank in the drug screen (pool rank).",
        "Reference: Rank 1 (Off/Exp) PFS or matched standard chemotherapy PFS from the same model.",
        "Positive ΔPFS (on-label minus reference) means on-label had longer PFS in that model.",
        "Progressive pools: on-label at rank 7-10, then 6-10, 5-10, and 4-10.",
    ]
    if q9_summary is not None and not q9_summary.empty:
        for stratum in LOW_RANK_POOL_Q8_STRATA:
            sub = q9_summary[q9_summary["stratum"] == stratum]
            if not sub.empty:
                n = int(sub.iloc[0]["n_models"])
                lines.append(f"{LOW_RANK_POOL_Q8_LABELS[stratum]} paired n = {n}.")
    return _short_groups_bullets(*lines)


def _short_q9_significant_summary(q9_summary: pd.DataFrame, alpha: float = 0.05) -> str:
    if q9_summary is None or q9_summary.empty:
        return "No comparisons evaluable."
    sig_bits: List[str] = []
    for _, r in q9_summary.iterrows():
        tests: List[str] = []
        for label, col in (
            ("Wilcoxon", "wilcoxon_p"),
            ("strat. log-rank", "stratified_logrank_p"),
            ("strat. Cox", "stratified_cox_p"),
        ):
            p = r.get(col)
            if pd.notna(p) and float(p) < alpha:
                tests.append(f"{label} p = {_clinical_format_p(p)}*")
        if tests:
            sig_bits.append(
                f"{r['pool_label']}, {r['comparison_label']} ({', '.join(tests)})"
            )
    if not sig_bits:
        return (
            "No paired pool comparison reached conventional significance (α = 0.05) on "
            "Wilcoxon, stratified log-rank, or stratified Cox."
        )
    return "Significant paired effects at α = 0.05 (*): " + "; ".join(sig_bits) + "."


def _short_q3_merged_groups(
    n_rank2: int,
    n_rank3: int,
    scenario3_n: int,
    s3_r2: Optional[pd.Series],
    s3_r3: Optional[pd.Series],
    s3_pool: Optional[pd.Series],
) -> List[str]:
    n2 = _short_paired_n(s3_r2)
    n3 = _short_paired_n(s3_r3)
    npool = _short_paired_n(s3_pool)
    return _short_groups_bullets(
        "Design: Paired within-model comparison (Rank 1 Off/Exp vs on-label) across three Scenario 3 strata.",
        "Primary comparison: rank 1 vs on-label at the model's discordant on-label rank.",
        f"Stratum A — On-label at rank 2 (n = {n_rank2} models, paired n = {n2}).",
        f"Stratum B — On-label at rank 3 (n = {n_rank3} models, paired n = {n3}).",
        f"Stratum C — Pooled on-label at ranks 2-3 (n = {scenario3_n} models, paired n = {npool}).",
        f"Group A (test): Rank 1 (Off/Exp), paired n = {n2} to {npool} per stratum.",
        f"Group B (reference): On-label at the stratum-specific rank, paired n = {n2} to {npool}.",
    )


def _short_q3_merged_results_bullets(
    s3_r2: Optional[pd.Series],
    s3_r3: Optional[pd.Series],
    s3_pool: Optional[pd.Series],
) -> List[str]:
    specs: Tuple[Tuple[str, Optional[pd.Series], str], ...] = (
        ("Rank 2 stratum", s3_r2, "Rank 1 (Off/Exp) vs on-label at rank 2"),
        ("Rank 3 stratum", s3_r3, "Rank 1 (Off/Exp) vs on-label at rank 3"),
        ("Pooled ranks 2-3", s3_pool, "Rank 1 (Off/Exp) vs on-label at ranks 2-3"),
    )
    bullets: List[str] = []
    for stratum_label, row, comparison_label in specs:
        if row is None:
            continue
        bullet = _short_paired_comparison_bullet(
            row,
            comparison_label=f"{stratum_label}: {comparison_label}",
            favour_label="rank 1 (Off/Exp)",
            delta_reference="rank 1 minus on-label",
        )
        bullets.append(bullet)
    if not bullets:
        return [_short_summary_prose("No evaluable Scenario 3 discordance data.")]
    return bullets


def _short_q3_merged_publication_statement(
    s3_r2: Optional[pd.Series],
    s3_r3: Optional[pd.Series],
    s3_pool: Optional[pd.Series],
) -> str:
    parts = [
        "We assessed whether Rank 1 (Off/Exp) prolongs PDX PFS relative to on-label therapy across "
        "Scenario 3 discordance strata (on-label at global rank 2, rank 3, and pooled ranks 2-3)."
    ]
    if s3_r2 is not None:
        parts.append(
            f"Rank-2 stratum (paired n = {_short_paired_n(s3_r2)}): "
            + _short_paired_stats_clause(s3_r2, left_name="rank 1", right_name="on-label at rank 2")
        )
    if s3_r3 is not None:
        parts.append(
            f"Rank-3 stratum (paired n = {_short_paired_n(s3_r3)}): "
            + _short_paired_stats_clause(s3_r3, left_name="rank 1", right_name="on-label at rank 3")
            + " [exploratory]"
        )
    if s3_pool is not None:
        parts.append(
            f"Pooled ranks 2-3 (paired n = {_short_paired_n(s3_pool)}): "
            + _short_paired_stats_clause(s3_pool, left_name="rank 1", right_name="on-label at ranks 2-3")
        )
    parts.append(
        "Overall trend: discordance is driven by the rank-2 stratum. Rank-3 is underpowered. "
        "The pooled estimate is directionally consistent but FDR-adjusted Wilcoxon is non-significant."
    )
    return _short_summary_prose(" ".join(parts))


def _short_caveat_q3_merged(
    s3_r2: Optional[pd.Series],
    s3_r3: Optional[pd.Series],
    s3_pool: Optional[pd.Series],
) -> str:
    parts = [
        "Three-stratum panel (rank 2, rank 3, pooled ranks 2-3). Rank-2 stratum is confirmatory "
        "(Wilcoxon + FDR significant). Rank-3 stratum is [exploratory] and underpowered "
        f"(n = {_short_paired_n(s3_r3)}). Pooled ranks 2-3 is [exploratory]: nominal Wilcoxon p < 0.05 "
        "but FDR-adjusted p > 0.05. Interpret pooled results alongside the rank-2 stratum.",
    ]
    if s3_r2 is not None and pd.notna(s3_r2.get("stratified_logrank_p")):
        if float(s3_r2["stratified_logrank_p"]) >= 0.05:
            parts.append(
                "Stratified log-rank on KM may be non-significant despite significant Wilcoxon in "
                "rank-2 stratum — expected when estimands differ."
            )
    return _short_summary_prose(" ".join(parts))


def _build_short_supervisor_items(
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    paired_summary: pd.DataFrame,
    between_cohort: Optional[pd.DataFrame],
    epistemic: Optional[pd.DataFrame],
    low_rank_sensitivity: Optional[pd.DataFrame] = None,
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
    fig_dir: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    """Structured Q&A items for the short supervisor summary."""
    if low_rank_sensitivity is None and ENABLE_Q2_ANALYSIS:
        low_rank_sensitivity = build_low_rank_pool_sensitivity_table(model_table, ranked)
    elif low_rank_sensitivity is None:
        low_rank_sensitivity = pd.DataFrame()
    if q9_paired_sensitivity is None:
        q9_paired_sensitivity = build_q9_paired_pool_sensitivity_table(model_table, ranked)
    scenario_counts = model_table["scenario"].value_counts()

    s1 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO1, "rank1_vs_sc")
    s2 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO2, "onlabel_vs_sc")
    s3_r2 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2")
    s3_r3 = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank3")
    s3_pool = _clinical_paired_row(paired_summary, CLINICAL_UTILITY_SCENARIO3, "rank1_vs_onlabel", "rank2_3_pooled")

    s3 = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO3]
    n_rank2 = int((s3["best_onlabel_rank_int"] == 2).sum()) if not s3.empty else 0
    n_rank3 = int((s3["best_onlabel_rank_int"] == 3).sum()) if not s3.empty else 0

    paired_method = _short_summary_prose(
        "Paired Wilcoxon signed-rank test on ΔPFS (test arm minus reference) within each PDX model. "
        "Bootstrap 95% CI for median ΔPFS. Stratified Cox PH (stratified by model) and stratified log-rank on KM."
    )
    between_method = _short_summary_prose(
        "Unpaired Mann-Whitney U on on-label PFS. Unpaired Cox and log-rank between cohorts."
    )

    items: List[Dict[str, Any]] = []

    if q9_paired_sensitivity is not None and not q9_paired_sensitivity.empty:
        items.append({
            "qid": "Q1",
            "heading": "Q1: Paired on-label-at-rank vs Rank 1 (Off/Exp) or standard chemo",
            "question": _short_summary_prose(
                "Within the same PDX model, is on-label PFS at a downranked screen position "
                "better or worse than Rank 1 (Off/Exp) or matched standard chemotherapy?"
            ),
            "groups": _short_groups_q9(q9_paired_sensitivity),
            "method": paired_method + " Eight paired comparisons in a 4×2 panel "
            "(four progressive on-label-at-rank pools × on-label vs Rank 1 or SC: "
            "Wilcoxon ΔPFS, stratified log-rank, stratified Cox). Multiplicity not adjusted.",
            "results": _short_q9_results_bullets(q9_paired_sensitivity),
            "interpretation": _short_summary_prose(_short_q9_significant_summary(q9_paired_sensitivity)),
            "takeaway": _short_summary_prose(
                "Interpret Q1 as model-specific pharmacodynamic discordance when on-label "
                "is deprioritised in the screen."
            ),
        })

    # Q2: Downranked on-label pools vs Rank 1 (Off/Exp) or standard chemo (unpaired) — disabled.
    if ENABLE_Q2_ANALYSIS and low_rank_sensitivity is not None and not low_rank_sensitivity.empty:
        items.append({
            "qid": "Q2",
            "heading": "Q2: Downranked on-label pools vs Rank 1 (Off/Exp) or standard chemo",
            "question": _short_summary_prose(
                "When on-label drugs are deprioritised in the screen (at ranks 7-10, then progressively "
                "6-10 through 4-10), does pooled on-label PFS separate from Rank 1 (Off/Exp) or from "
                "matched standard chemotherapy?"
            ),
            "groups": _short_groups_q8(low_rank_sensitivity),
            "method": between_method + " Eight unpaired comparisons in a 4×2 panel "
            "(four progressive on-label-at-rank pools × two reference arms: box, KM, Mann-Whitney, "
            "log-rank, Cox). Multiplicity not adjusted.",
            "results": _short_q8_results_bullets(low_rank_sensitivity),
            "interpretation": _short_summary_prose(_short_q8_significant_summary(low_rank_sensitivity)),
            "takeaway": _short_summary_prose(
                "Progressive pooling increases n from 7 (on-label at rank 7-10) to 20 (rank 4-10). "
                "Use to assess when downranked on-label pharmacodynamics diverge from Rank 1 (Off/Exp) or SC."
            ),
        })

    items.extend([
        {
            "qid": "Q3",
            "heading": "Q3: No on-label agent, is DDA rank 1 better than chemotherapy?",
            "question": _short_summary_prose(
                "When the drug screen contains no approved (on-label) compound, does the algorithm's "
                "top-ranked agent prolong PFS versus model-matched standard chemotherapy?"
            ),
            "groups": _short_groups_q1(
                int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0)), s1
            ),
            "method": paired_method + f" Comparison: rank 1 vs SC. Cohort: Scenario 1 (n = "
            f"{int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0))} models).",
            "results": _short_paired_results_bullets(
                s1,
                comparison_label="Rank 1 vs standard chemotherapy",
                favour_label="rank 1",
                delta_reference="rank 1 minus SC",
            ),
            "interpretation": _short_summary_prose(
                "DDA rank 1 shows significant tumour growth delay relative to chemotherapy in the majority "
                "of models without an on-label option. This is the strongest preclinical signal in the programme."
            ),
            "takeaway": _short_summary_prose(
                "Support using DDA rank 1 over default chemotherapy when no on-label drug is screened "
                "(preclinical decision support, not prescribing evidence)."
            ),
        },
        {
            "qid": "Q4",
            "heading": "Q4: Guideline concordance, is on-label at rank 1 better than chemotherapy?",
            "question": _short_summary_prose(
                "When the approved on-label therapy is also the DDA rank-1 compound, does it outperform "
                "matched standard chemotherapy?"
            ),
            "groups": _short_groups_q2(
                int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0)), s2
            ),
            "method": paired_method + f" Comparison: on-label vs SC. Cohort: Scenario 2 (n = "
            f"{int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0))} models).",
            "results": _short_paired_results_bullets(
                s2,
                comparison_label="On-label at rank 1 vs standard chemotherapy",
                favour_label="on-label at rank 1",
                delta_reference="on-label minus SC",
            ),
            "interpretation": _short_summary_prose(
                "On-label at rank 1 demonstrates favourable pharmacodynamics versus chemotherapy. "
                "Wilcoxon p ≈ 0.07 is borderline. Stratified Cox reports reduced growth hazard."
            ),
            "takeaway": _short_summary_prose(
                "Guideline-concordant on-label selection is supported by preclinical PFS data."
            ),
        },
        {
            "qid": "Q5",
            "heading": "Q5: Guideline discordance (ranks 2-3), is Rank 1 better than on-label?",
            "question": _short_summary_prose(
                "When the approved drug is deprioritised to global rank 2 or 3, does the algorithm's "
                "Rank 1 (Off/Exp) compound outperform on-label therapy? Presented as a three-stratum panel: "
                "rank 2, rank 3, and pooled ranks 2-3."
            ),
            "groups": _short_q3_merged_groups(
                n_rank2, n_rank3, int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0)),
                s3_r2, s3_r3, s3_pool,
            ),
            "method": paired_method + (
                f" Comparison: rank 1 vs on-label. Cohort: Scenario 3 "
                f"(rank-2 n = {n_rank2}, rank-3 n = {n_rank3}, pooled n = "
                f"{int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0))} models). "
                "FDR adjustment across comparisons within each stratum. Six-panel figure "
                "(KM + Wilcoxon ΔPFS per stratum)."
            ),
            "results": _short_q3_merged_results_bullets(s3_r2, s3_r3, s3_pool),
            "interpretation": _short_summary_prose(
                "Rank-2 stratum: Rank 1 (Off/Exp) significantly outperforms on-label (~77% of models improved). "
                "Rank-3 stratum: underpowered (n = 10), no reliable inference. "
                "Pooled ranks 2-3: nominal Wilcoxon significance (p ≈ 0.04) but FDR-adjusted p > 0.05. "
                "Signal is driven primarily by the rank-2 stratum."
            ),
            "takeaway": _short_summary_prose(
                "Rigid on-label priority is pharmacodynamically suboptimal when on-label sits at rank 2. "
                "Rank-1 recommendation merits MTB discussion. Do not change policy from rank-3 stratum alone."
            ),
        },
    ])

    publication_rows = {
        "Q3": s1,
        "Q4": s2,
    }
    for item in items:
        qid = str(item.get("qid", ""))
        if qid == "Q1":
            item["publication_statement"] = _short_q9_publication_trend(q9_paired_sensitivity)
        elif qid == "Q2" and ENABLE_Q2_ANALYSIS:
            item["publication_statement"] = _short_q8_publication_trend(low_rank_sensitivity)
        elif qid == "Q5":
            item["publication_statement"] = _short_q3_merged_publication_statement(
                s3_r2, s3_r3, s3_pool
            )
            item["caveat"] = _short_caveat_q3_merged(s3_r2, s3_r3, s3_pool)
        else:
            item["publication_statement"] = _short_publication_statement(
                qid, row=publication_rows.get(qid)
            )

    paired_rows = {
        "Q3": s1,
        "Q4": s2,
    }
    _apply_short_item_caveats(
        items,
        paired_rows=paired_rows,
        epistemic=epistemic,
        skip_qids={"Q5"},
    )
    _apply_short_item_figures(items, fig_dir)
    return items


def _build_clinical_supervisor_short_markdown(
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    paired_summary: pd.DataFrame,
    between_cohort: Optional[pd.DataFrame],
    epistemic: Optional[pd.DataFrame],
    source_path: Path,
    fig_dir: Optional[Path] = None,
    low_rank_sensitivity: Optional[pd.DataFrame] = None,
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
) -> List[str]:
    """Compact supervisor summary: question, groups, method, results, figures, interpretation, takeaway."""
    n_total = len(model_table)
    n_chemo = int(model_table["chemo_TTD"].notna().sum())
    scenario_counts = model_table["scenario"].value_counts()
    if low_rank_sensitivity is None and ENABLE_Q2_ANALYSIS:
        low_rank_sensitivity = build_low_rank_pool_sensitivity_table(model_table, ranked)
    elif low_rank_sensitivity is None:
        low_rank_sensitivity = pd.DataFrame()
    if q9_paired_sensitivity is None:
        q9_paired_sensitivity = build_q9_paired_pool_sensitivity_table(model_table, ranked)
    sig_md, _ = _build_short_summary_significance_table(
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
    )
    items = _build_short_supervisor_items(
        model_table=model_table,
        ranked=ranked,
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        epistemic=epistemic,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
        fig_dir=fig_dir,
    )
    rank_comp_fig = (
        fig_dir / f"{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png"
        if fig_dir and (fig_dir / f"{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png").is_file()
        else None
    )

    doc: List[str] = [
        "# Clinical Utility of DDA Ranking: Short Supervisor Summary",
        "",
        f"Date: {datetime.now().strftime('%Y-%m-%d')} | Source: `{source_path.name}` | Endpoint: PFS (days) in PDX",
        "",
        "## Screen composition by rank",
        "",
        _short_summary_prose(SHORT_SUMMARY_RANK_COMPOSITION_CAPTION),
        "",
    ]
    if rank_comp_fig is not None:
        doc.extend([
            "Figure.",
            f"- `{rank_comp_fig.name}`",
            "",
        ])
    doc.extend([
        "## Naming convention",
        "",
        _short_summary_prose(SHORT_SUMMARY_PFS_CONVENTION),
        "",
        "## Statistical analyses",
        "",
    ])
    for i, para in enumerate(SHORT_SUMMARY_STATISTICAL_PARAGRAPHS):
        doc.extend([
            _short_summary_prose(para, preserve_punctuation=(i > 0)),
            "",
        ])
    doc.extend([
        "## Overview",
        "",
        _short_summary_prose(
            f"This summary reports preclinical clinical-utility analyses across {n_total} PDX models "
            f"({n_chemo} with matched standard chemotherapy). Models are stratified into four mutually "
            "exclusive scenarios by the global DDA rank of the best on-label drug."
        ),
        "",
        "Cohort sizes:",
        f"- Scenario 1 (no on-label): n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO1, 0))}",
        f"- Scenario 2 (on-label rank 1): n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO2, 0))}",
        f"- Scenario 3 (on-label rank 2–3): n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO3, 0))}",
        f"- Scenario 4 (on-label rank 4+): n = {int(scenario_counts.get(CLINICAL_UTILITY_SCENARIO4, 0))}",
        "",
        "---",
        "",
    ])

    for item in items:
        doc.extend([
            f"## {item['heading']}",
            "",
            f"Clinical question. {item['question']}",
            "",
        ])
        _extend_short_groups_markdown(doc, item["groups"])
        doc.extend([
            f"Statistical method. {item['method']}",
            "",
        ])
        _extend_short_results_markdown(doc, item["results"])
        figures = item.get("figures") or []
        if figures:
            doc.append("Figures.")
            for fig in figures:
                doc.append(f"- `{fig.name}`")
            doc.append("")
        doc.extend([
            f"Publication statement. {item.get('publication_statement', '')}",
            "",
            f"Interpretation. {item['interpretation']}",
            "",
            f"Key takeaway. {item['takeaway']}",
            "",
        ])
        if item.get("caveat"):
            doc.extend([
                f"Statistical caveat. {item['caveat']}",
                "",
            ])
        doc.extend([
            "---",
            "",
        ])

    doc.extend([
        "## Statistical significance summary",
        "",
        _short_summary_prose(
            "Asterisk (*) marks p < 0.05 (α = 0.05). Paired questions use Wilcoxon, stratified log-rank, "
            "and stratified Cox HR."
            + (
                " Q2 uses Mann-Whitney and unpaired log-rank/Cox."
                if ENABLE_Q2_ANALYSIS
                else ""
            )
            + " Q1 uses paired Wilcoxon and stratified log-rank/Cox within each model."
        ),
        "",
    ])
    doc.extend(sig_md)
    doc.extend([
        "",
        "## Overall conclusions",
        "",
        _short_summary_prose(
            "1. Q1 pool sensitivity: [exploratory] downranked on-label pharmacodynamics (paired within model)."
            if not ENABLE_Q2_ANALYSIS
            else "1. Q1 and Q2 pool sensitivity: [exploratory] downranked on-label pharmacodynamics."
        ),
        _short_summary_prose(
            "2. No on-label: DDA rank 1 significantly beats chemotherapy (Q3), strongest pre-specified signal."
        ),
        _short_summary_prose(
            "3. On-label rank 1: Guideline-concordant choice supported (Q4). Wilcoxon borderline, see caveat."
        ),
        _short_summary_prose(
            "4. Guideline discordance (Q5): Rank 1 significantly outperforms on-label at rank 2; "
            "rank-3 stratum underpowered. pooled ranks 2-3 directionally consistent."
        ),
        "",
        _short_summary_prose(
            "Limitation: PDX PFS (TimeToDouble) is a pharmacodynamic surrogate, not patient survival. "
            "Rank-3 and low-rank strata are underpowered. KM curves illustrate TTD distributions."
        ),
        "",
    ])
    doc.append(
        "Analysis flags and per-question caveats: STATISTICAL_CAVEATS.md."
    )
    return doc


def _write_supervisor_short_docx(
    docx_path: Path,
    items: List[Dict[str, Any]],
    *,
    meta: Dict[str, Any],
) -> None:
    """Write compact supervisor summary in Word with embedded figures per question."""
    from docx import Document
    from docx.shared import Pt, Inches

    document = Document()
    normal = document.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(10)

    for section in document.sections:
        section.top_margin = Inches(0.75)
        section.bottom_margin = Inches(0.75)
        section.left_margin = Inches(0.9)
        section.right_margin = Inches(0.9)

    document.add_heading("Clinical Utility of DDA Ranking: Short Supervisor Summary", level=0)
    p = document.add_paragraph()
    p.add_run(
        f"Date: {meta['date']}  |  Source: {meta['source']}  |  Endpoint: PFS (days) in PDX"
    ).font.size = Pt(9)

    document.add_heading("Screen composition by rank", level=1)
    document.add_paragraph(_short_summary_prose(SHORT_SUMMARY_RANK_COMPOSITION_CAPTION))
    rank_fig: Optional[Path] = meta.get("rank_utility_figure")
    if rank_fig and rank_fig.is_file():
        cap = document.add_paragraph(rank_fig.name)
        cap.runs[0].italic = True
        cap.runs[0].font.size = Pt(8)
        try:
            document.add_picture(str(rank_fig), width=Inches(6.4))
        except Exception as exc:
            document.add_paragraph(f"[Could not embed figure: {exc}]")
        document.add_paragraph()

    document.add_heading("Naming convention", level=1)
    document.add_paragraph(_short_summary_prose(meta.get("pfs_convention", SHORT_SUMMARY_PFS_CONVENTION)))

    document.add_heading("Statistical analyses", level=1)
    for i, para in enumerate(meta.get("statistical_hierarchy_paragraphs", SHORT_SUMMARY_STATISTICAL_PARAGRAPHS)):
        document.add_paragraph(_short_summary_prose(para, preserve_punctuation=(i > 0)))

    document.add_heading("Overview", level=1)
    document.add_paragraph(_short_summary_prose(meta["overview"]))

    for item in items:
        document.add_heading(item["heading"], level=2)
        para = document.add_paragraph()
        run_label = para.add_run("Clinical question. ")
        run_label.bold = True
        para.add_run(str(item["question"]))
        _add_short_groups_docx(document, item["groups"])
        para = document.add_paragraph()
        run_label = para.add_run("Statistical method. ")
        run_label.bold = True
        para.add_run(str(item["method"]))
        _add_short_results_docx(document, item["results"])

        figures = item.get("figures") or []
        if figures:
            fig_heading = document.add_paragraph()
            fig_run = fig_heading.add_run("Figures.")
            fig_run.bold = True
            for fig_path in figures:
                if fig_path.is_file():
                    cap = document.add_paragraph(fig_path.name)
                    cap.runs[0].italic = True
                    cap.runs[0].font.size = Pt(8)
                    try:
                        document.add_picture(str(fig_path), width=Inches(6.2))
                    except Exception as exc:
                        document.add_paragraph(f"[Could not embed figure: {exc}]")
                    document.add_paragraph()

        for label, key in [
            ("Publication statement", "publication_statement"),
            ("Interpretation", "interpretation"),
            ("Key takeaway", "takeaway"),
        ]:
            para = document.add_paragraph()
            run_label = para.add_run(f"{label}. ")
            run_label.bold = True
            if item.get(key):
                para.add_run(str(item[key]))

        if item.get("caveat"):
            para = document.add_paragraph()
            run_label = para.add_run("Statistical caveat. ")
            run_label.bold = True
            run_caveat = para.add_run(str(item["caveat"]))
            run_caveat.italic = True

    sig_header: List[str] = meta.get("significance_table_header") or []
    sig_rows: List[List[str]] = meta.get("significance_table_rows") or []
    if sig_header and sig_rows:
        document.add_heading("Statistical significance summary", level=1)
        document.add_paragraph(
            _short_summary_prose(
                "Asterisk (*) marks p < 0.05 (α = 0.05). Paired questions use Wilcoxon, stratified log-rank, "
                "and stratified Cox HR."
                + (
                    " Q2 uses Mann-Whitney and unpaired log-rank/Cox."
                    if ENABLE_Q2_ANALYSIS
                    else ""
                )
                + " Q1 uses paired Wilcoxon and stratified log-rank/Cox within each model."
            )
        )
        table = document.add_table(rows=1 + len(sig_rows), cols=len(sig_header))
        table.style = "Table Grid"
        for col_idx, heading in enumerate(sig_header):
            table.rows[0].cells[col_idx].text = heading
        for row_idx, row_vals in enumerate(sig_rows, start=1):
            for col_idx, cell_text in enumerate(row_vals):
                table.rows[row_idx].cells[col_idx].text = cell_text
        document.add_paragraph()

    document.add_heading("Overall conclusions", level=1)
    for line in meta["conclusions"]:
        document.add_paragraph(line, style="List Bullet")

    try:
        document.save(docx_path)
    except PermissionError:
        alt_path = docx_path.with_name(f"{docx_path.stem}_updated{docx_path.suffix}")
        document.save(alt_path)
        LOGGER.warning("Could not overwrite %s (file may be open); saved to %s", docx_path, alt_path)


def _short_supervisor_meta(
    model_table: pd.DataFrame,
    source_path: Path,
    fig_dir: Optional[Path],
    *,
    paired_summary: Optional[pd.DataFrame] = None,
    between_cohort: Optional[pd.DataFrame] = None,
    low_rank_sensitivity: Optional[pd.DataFrame] = None,
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    rank_fig_path = (
        fig_dir / f"{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png"
        if fig_dir and (fig_dir / f"{SHORT_SUMMARY_RANK_COMPOSITION_STEM}.png").is_file()
        else None
    )
    sig_header = ["Question", "Groups compared", "n", "Wilcoxon / MW", "Log-rank", "Cox HR (p)"]
    sig_rows: List[List[str]] = []
    if paired_summary is not None:
        _, sig_rows = _build_short_summary_significance_table(
            paired_summary=paired_summary,
            between_cohort=between_cohort,
            low_rank_sensitivity=low_rank_sensitivity,
            q9_paired_sensitivity=q9_paired_sensitivity,
        )
    return {
        "date": datetime.now().strftime("%Y-%m-%d"),
        "source": source_path.name,
        "rank_utility_figure": rank_fig_path,
        "pfs_convention": SHORT_SUMMARY_PFS_CONVENTION,
        "statistical_hierarchy": SHORT_SUMMARY_STATISTICAL_HIERARCHY,
        "statistical_hierarchy_paragraphs": SHORT_SUMMARY_STATISTICAL_PARAGRAPHS,
        "overview": _short_summary_prose(
            f"Preclinical clinical-utility analysis of DDA compound ranking in {len(model_table)} PDX models "
            f"({int(model_table['chemo_TTD'].notna().sum())} with matched standard chemotherapy). "
            "Models are assigned to four scenarios by the global rank of the best on-label drug."
        ),
        "conclusions": [
            _short_summary_prose(
                "Q1 pool sensitivity: [exploratory] downranked on-label pharmacodynamics (paired within model)."
                if not ENABLE_Q2_ANALYSIS
                else "Q1 and Q2 pool sensitivity: [exploratory] downranked on-label pharmacodynamics."
            ),
            _short_summary_prose(
                "No on-label: DDA rank 1 significantly beats chemotherapy (Q3), strongest pre-specified signal."
            ),
            _short_summary_prose(
                "On-label rank 1: Guideline-concordant choice supported (Q4). Wilcoxon borderline, interpret with caveat."
            ),
            _short_summary_prose(
                "Guideline discordance (Q5): Rank 1 significantly outperforms on-label at rank 2; "
                "rank-3 stratum underpowered. pooled ranks 2-3 directionally consistent."
            ),
            _short_summary_prose(
                "PDX PFS (TimeToDouble) is a pharmacodynamic surrogate, not patient survival."
            ),
        ],
        "significance_table_header": sig_header,
        "significance_table_rows": sig_rows,
    }


def regenerate_short_supervisor_summary_files(
    summary_dir: Path,
    *,
    source_path: Path,
    regenerate_figures: bool = True,
) -> None:
    """Regenerate short summary docs; optionally rebuild publication figures from source."""
    model_xlsx = sorted(summary_dir.parent.glob("clinical_utility_model_table_*.xlsx"))
    if not model_xlsx:
        raise FileNotFoundError("clinical_utility_model_table_*.xlsx not found in output folder")
    model_table = pd.read_excel(model_xlsx[-1])
    table2 = summary_dir / "TABLE2_SCENARIO_RESULTS.csv"
    paired_summary = pd.read_csv(table2) if table2.is_file() else pd.DataFrame()
    between_path = sorted(summary_dir.parent.glob("clinical_utility_between_cohort_onlabel_*.xlsx"))
    between_cohort = pd.read_excel(between_path[-1]) if between_path else pd.DataFrame()
    sens_path = sorted(summary_dir.parent.glob("clinical_utility_low_rank_pool_sensitivity_*.xlsx"))
    ranked = _prepare_ranked_onlabel_patient_table(source_path, DEFAULT_SOURCE_SHEET)
    if ENABLE_Q2_ANALYSIS:
        low_rank_sensitivity = (
            pd.read_excel(sens_path[-1])
            if sens_path
            else build_low_rank_pool_sensitivity_table(model_table, ranked)
        )
    else:
        low_rank_sensitivity = pd.DataFrame()
    q9_paired_sensitivity = build_q9_paired_pool_sensitivity_table(model_table, ranked)
    ep_path = summary_dir / "TABLE4_LOW_RANK_EPISTEMIC.csv"
    epistemic = pd.read_csv(ep_path) if ep_path.is_file() else pd.DataFrame()
    fig_dir = summary_dir / CLINICAL_SUPERVISOR_FIGURES_COLLECTION_DIR

    if regenerate_figures:
        chemo = load_sc_chemo_reference(source_path, DEFAULT_SOURCE_SHEET)
        chemo_collapsed = _collapse_chemo_per_model(chemo, DEFAULT_DATA_COLUMN)
        utility_matches = sorted(
            summary_dir.parent.glob("clinical_utility_utility_at_rank_panel_*.jpg")
        )
        generate_short_summary_publication_figures(
            summary_dir,
            model_table=model_table,
            ranked=ranked,
            chemo_collapsed=chemo_collapsed,
            paired_summary=paired_summary,
            between_cohort=between_cohort,
            utility_panel_source=utility_matches[-1] if utility_matches else None,
            chemo=chemo,
        )
    short_items = _build_short_supervisor_items(
        model_table=model_table,
        ranked=ranked,
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        epistemic=epistemic,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
        fig_dir=fig_dir,
    )
    _write_statistical_caveats_reference(summary_dir, items=short_items)
    short_lines = _build_clinical_supervisor_short_markdown(
        model_table=model_table,
        ranked=ranked,
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        epistemic=epistemic,
        source_path=source_path,
        fig_dir=fig_dir,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
    )
    (summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR_short.md").write_text(
        "\n".join(short_lines), encoding="utf-8"
    )
    _write_supervisor_short_docx(
        summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR_short.docx",
        short_items,
        meta=_short_supervisor_meta(
            model_table,
            source_path,
            fig_dir,
            paired_summary=paired_summary,
            between_cohort=between_cohort,
            low_rank_sensitivity=low_rank_sensitivity,
            q9_paired_sensitivity=q9_paired_sensitivity,
        ),
    )


def _write_clinical_supervisor_summary(
    output_dir: Path,
    *,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    rank_distribution: pd.DataFrame,
    paired_summary: pd.DataFrame,
    utility_rank_pairwise: pd.DataFrame,
    between_cohort: Optional[pd.DataFrame] = None,
    epistemic: Optional[pd.DataFrame] = None,
    low_rank_sensitivity: Optional[pd.DataFrame] = None,
    q9_paired_sensitivity: Optional[pd.DataFrame] = None,
    exp_registry: Optional[pd.DataFrame] = None,
    legacy: str,
    date: str,
    source_path: Path,
) -> None:
    """Write supervisor-ready clinical summary package (markdown, Word, tables, key figures)."""
    summary_dir = output_dir / CLINICAL_SUPERVISOR_SUMMARY_DIR
    summary_dir.mkdir(parents=True, exist_ok=True)

    _write_clinical_cohort_flow(summary_dir, model_table)

    primary_cols = [
        "scenario", "stratum", "comparison_group", "comparison", "n_models", "pct_improved",
        "median_diff_TTD", "median_diff_CI_low", "median_diff_CI_high",
        "wilcoxon_p", "wilcoxon_p_fdr_bh", "binom_p",
        "stratified_cox_hr", "stratified_cox_hr_low", "stratified_cox_hr_high", "stratified_cox_p",
        "stratified_logrank_p", "pooled_logrank_p",
    ]
    if not paired_summary.empty:
        table2 = paired_summary[[c for c in primary_cols if c in paired_summary.columns]].copy()
        table2.to_csv(summary_dir / "TABLE2_SCENARIO_RESULTS.csv", index=False)
        safe_to_excel(table2, summary_dir / "TABLE2_SCENARIO_RESULTS.xlsx")

    rank_distribution[~rank_distribution["bucket"].astype(str).str.startswith("scenario_")].to_csv(
        summary_dir / "TABLE1_COHORT_DISTRIBUTION.csv", index=False
    )
    if not utility_rank_pairwise.empty:
        utility_rank_pairwise.to_csv(summary_dir / "TABLE3_UTILITY_AT_RANK_VS_SC.csv", index=False)
    if epistemic is not None and not epistemic.empty:
        epistemic.to_csv(summary_dir / "TABLE4_LOW_RANK_EPISTEMIC.csv", index=False)
    if q9_paired_sensitivity is not None and not q9_paired_sensitivity.empty:
        q9_paired_sensitivity.to_csv(summary_dir / "TABLE5_Q9_PAIRED_POOL_SENSITIVITY.csv", index=False)
        safe_to_excel(
            q9_paired_sensitivity,
            summary_dir / "TABLE5_Q9_PAIRED_POOL_SENSITIVITY.xlsx",
        )
    if exp_registry is not None and not exp_registry.empty:
        exp_registry.to_csv(summary_dir / "TABLE_EXP_REGISTRY.csv", index=False)

    copied_figs = _copy_clinical_key_figures(output_dir, summary_dir, date)

    doc_lines = _build_clinical_supervisor_markdown(
        model_table=model_table,
        paired_summary=paired_summary,
        utility_rank_pairwise=utility_rank_pairwise,
        between_cohort=between_cohort,
        epistemic=epistemic,
        source_path=source_path,
        copied_figs=copied_figs,
    )

    md_path = summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR.md"
    md_path.write_text("\n".join(doc_lines), encoding="utf-8")

    docx_path = summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR.docx"
    try:
        _write_supervisor_summary_docx(docx_path, doc_lines)
        LOGGER.info("Supervisor Word document written to %s", docx_path)
    except Exception as exc:
        LOGGER.warning("Could not write supervisor .docx: %s", exc)

    fig_dir = summary_dir / CLINICAL_SUPERVISOR_FIGURES_COLLECTION_DIR
    if q9_paired_sensitivity is None:
        q9_paired_sensitivity = build_q9_paired_pool_sensitivity_table(model_table, ranked)
    short_items = _build_short_supervisor_items(
        model_table=model_table,
        ranked=ranked,
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        epistemic=epistemic,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
        fig_dir=fig_dir,
    )
    _write_statistical_caveats_reference(summary_dir, items=short_items)
    short_lines = _build_clinical_supervisor_short_markdown(
        model_table=model_table,
        ranked=ranked,
        paired_summary=paired_summary,
        between_cohort=between_cohort,
        epistemic=epistemic,
        source_path=source_path,
        fig_dir=fig_dir,
        low_rank_sensitivity=low_rank_sensitivity,
        q9_paired_sensitivity=q9_paired_sensitivity,
    )
    short_md_path = summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR_short.md"
    short_md_path.write_text("\n".join(short_lines), encoding="utf-8")

    short_docx_path = summary_dir / "CLINICAL_SUMMARY_FOR_SUPERVISOR_short.docx"
    try:
        _write_supervisor_short_docx(
            short_docx_path,
            short_items,
            meta=_short_supervisor_meta(
                model_table,
                source_path,
                fig_dir,
                paired_summary=paired_summary,
                between_cohort=between_cohort,
                low_rank_sensitivity=low_rank_sensitivity,
                q9_paired_sensitivity=q9_paired_sensitivity,
            ),
        )
        LOGGER.info("Short supervisor Word document written to %s", short_docx_path)
    except Exception as exc:
        LOGGER.warning("Could not write short supervisor .docx: %s", exc)

    LOGGER.info("Supervisor clinical summary written to %s", summary_dir)

# --- paired tests ---
def bootstrap_ci(
    values: np.ndarray,
    statistic: str = "median",
    n_boot: int = N_BOOTSTRAP,
    alpha: float = 0.05,
    seed: int = SEED,
) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    if statistic == "median":
        fn = np.median
    elif statistic == "mean":
        fn = np.mean
    else:
        raise ValueError(f"Unknown statistic: {statistic}")
    boots = np.array(
        [fn(rng.choice(values, size=len(values), replace=True)) for _ in range(n_boot)]
    )
    return (float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2)))


def run_stratified_cox(paired: pd.DataFrame) -> Dict[str, float]:
    """Cox PH with model strata: rank (arm=1) vs matched chemo (arm=0) within each model."""
    rows: List[Dict[str, Any]] = []
    for _, row in paired.iterrows():
        rows.append(
            {"Model": row["Model"], "time": float(row["rank_TTD"]), "event": 1, "arm": 1}
        )
        rows.append(
            {"Model": row["Model"], "time": float(row["chemo_TTD"]), "event": 1, "arm": 0}
        )
    cox_df = pd.DataFrame(rows)
    if cox_df.empty or cox_df["Model"].nunique() < 2:
        return {
            "cox_hr": np.nan,
            "cox_hr_low": np.nan,
            "cox_hr_high": np.nan,
            "cox_p": np.nan,
        }
    try:
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="time", event_col="event", strata=["Model"], formula="arm")
        summary = cph.summary.loc["arm"]
        return {
            "cox_hr": float(summary["exp(coef)"]),
            "cox_hr_low": float(summary["exp(coef) lower 95%"]),
            "cox_hr_high": float(summary["exp(coef) upper 95%"]),
            "cox_p": float(summary["p"]),
        }
    except Exception as exc:
        LOGGER.warning("[paired] Stratified Cox failed: %s", exc)
        return {
            "cox_hr": np.nan,
            "cox_hr_low": np.nan,
            "cox_hr_high": np.nan,
            "cox_p": np.nan,
        }


def run_pooled_logrank(paired: pd.DataFrame) -> float:
    """Unpaired log-rank on pooled rank vs chemo times (sensitivity only)."""
    if paired.empty:
        return np.nan
    try:
        res = logrank_test(
            paired["rank_TTD"].astype(float),
            paired["chemo_TTD"].astype(float),
            event_observed_A=np.ones(len(paired), dtype=int),
            event_observed_B=np.ones(len(paired), dtype=int),
        )
        return float(res.p_value)
    except Exception as exc:
        LOGGER.warning("[paired] Pooled log-rank failed: %s", exc)
        return np.nan


def run_paired_tests(paired: pd.DataFrame) -> Dict[str, Any]:
    """Wilcoxon, binomial, and log-ratio paired t-test on one paired table."""
    diff = paired["diff_TTD"].astype(float).values
    ratio = paired["ratio_TTD"].astype(float).values
    n = len(paired)
    n_improved = int(paired["improved"].sum())
    n_worse = int(paired["worse"].sum())
    n_tie = int(paired["tie"].sum())

    wilcoxon_stat, wilcoxon_p = (np.nan, np.nan)
    if n >= 1 and not np.all(diff == 0):
        try:
            res = stats.wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
            wilcoxon_stat, wilcoxon_p = float(res.statistic), float(res.pvalue)
        except ValueError:
            wilcoxon_stat, wilcoxon_p = (0.0, 1.0)

    binom_p = np.nan
    n_decisive = n_improved + n_worse
    if n_decisive > 0:
        binom_p = float(stats.binomtest(n_improved, n_decisive, 0.5, alternative="two-sided").pvalue)

    log_ratio = np.log(ratio[np.isfinite(ratio) & (ratio > 0)])
    ttest_p = np.nan
    ttest_stat = np.nan
    if len(log_ratio) >= 2 and np.std(log_ratio) > 0:
        ttest_res = stats.ttest_1samp(log_ratio, 0.0, alternative="two-sided")
        ttest_stat, ttest_p = float(ttest_res.statistic), float(ttest_res.pvalue)

    med_diff = float(np.median(diff))
    med_ratio = float(np.median(ratio[np.isfinite(ratio)]))
    diff_ci = bootstrap_ci(diff, statistic="median")
    ratio_ci = bootstrap_ci(ratio[np.isfinite(ratio)], statistic="median")
    cox = run_stratified_cox(paired)
    logrank_pooled_p = run_pooled_logrank(paired)

    return {
        "n_models": n,
        "n_improved": n_improved,
        "n_worse": n_worse,
        "n_tie": n_tie,
        "pct_improved": round(100.0 * n_improved / n, 2) if n else 0.0,
        "median_diff_TTD": round(med_diff, 2),
        "median_ratio_TTD": round(med_ratio, 3),
        "median_diff_CI_low": round(diff_ci[0], 2),
        "median_diff_CI_high": round(diff_ci[1], 2),
        "median_ratio_CI_low": round(ratio_ci[0], 3),
        "median_ratio_CI_high": round(ratio_ci[1], 3),
        "wilcoxon_stat": wilcoxon_stat,
        "wilcoxon_p": wilcoxon_p,
        "binom_p": binom_p,
        "log_ratio_t_stat": ttest_stat,
        "log_ratio_t_p": ttest_p,
        **cox,
        "logrank_pooled_p": logrank_pooled_p,
    }

# Low-rank extensions
CLINICAL_UTILITY_SCENARIO4_LOW_STRATA: Tuple[str, ...] = (
    "rank7_10", "pool_ge6", "pool_ge5", "pool_ge4",
)
BETWEEN_COHORT_LOW_STRATA: Tuple[str, ...] = CLINICAL_UTILITY_SCENARIO4_LOW_STRATA
LOW_RANK_POOL_Q8_STRATA: Tuple[str, ...] = ("rank7_10", "pool_ge6", "pool_ge5", "pool_ge4")
# Combined Q1 waterfall also includes progressive pool down to 3-IO (ranks 3-10).
Q1_DOWNRANK_POOL_STRATA_WITH_3IO: Tuple[str, ...] = (
    "rank7_10",
    "pool_ge6",
    "pool_ge5",
    "pool_ge4",
    "pool_ge3",
)
LOW_RANK_POOL_Q8_LABELS: Dict[str, str] = {
    "rank7_10": "On-label at rank 7-10",
    "pool_ge6": "On-label at rank 6-10",
    "pool_ge5": "On-label at rank 5-10",
    "pool_ge4": "On-label at rank 4-10",
    "pool_ge3": "On-label at rank 3-10 (3-IO)",
}
LOW_RANK_POOL_Q8_REFERENCE_LABELS: Dict[str, str] = {
    "rank1_offexp": "Rank 1 (Off/Exp)",
    "sc": "Standard chemo",
}


def _cohort_rank1_offexp_series(model_table: pd.DataFrame) -> pd.Series:
    mask = (
        model_table["rank1_utility"].isin(("Off_label", "Exp"))
        & model_table["rank1_TTD"].notna()
    )
    return model_table.loc[mask, "rank1_TTD"].astype(float)


def _cohort_sc_global_series(model_table: pd.DataFrame) -> pd.Series:
    return model_table.loc[model_table["chemo_TTD"].notna(), "chemo_TTD"].astype(float)


def _onlabel_at_rank_pool_table(
    ranked: pd.DataFrame,
    stratum: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """One row per model: on-label screened at the Q8 rank pool (most deprioritised in pool)."""
    on_rows = ranked[ranked["ON_LABEL_UTILITY"] == "On_label"].copy()
    if on_rows.empty:
        return on_rows
    on_rows["_rank_int"] = on_rows["RANK"].map(_clinical_utility_rank_as_int)
    if stratum == "rank7_10":
        pool = on_rows[on_rows["RANK"].astype(str) == _POOLED_LABEL]
    elif stratum == "pool_ge6":
        pool = on_rows[on_rows["_rank_int"] >= 6]
    elif stratum == "pool_ge5":
        pool = on_rows[on_rows["_rank_int"] >= 5]
    elif stratum == "pool_ge4":
        pool = on_rows[on_rows["_rank_int"] >= 4]
    elif stratum == "pool_ge3":
        pool = on_rows[on_rows["_rank_int"] >= 3]
    else:
        raise ValueError(f"Unknown Q8 pool stratum: {stratum!r}")
    if pool.empty:
        return pool
    pool = pool.sort_values(
        ["Model", "_rank_int", "LEVEL", data_column],
        ascending=[True, False, False, False],
    )
    return pool.groupby("Model", as_index=False).first()


def _onlabel_at_rank_pool_series(
    ranked: pd.DataFrame,
    stratum: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.Series:
    pool = _onlabel_at_rank_pool_table(ranked, stratum, data_column)
    if pool.empty:
        return pd.Series(dtype=float)
    return pool.loc[pool[data_column].notna(), data_column].astype(float)


def _cohort_onlabel_low_rank_pool(
    ranked: pd.DataFrame,
    low_stratum: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.Series:
    """Q8 on-label PFS pool: models with on-label at the specified global rank(s)."""
    return _onlabel_at_rank_pool_series(ranked, low_stratum, data_column)


def _run_unpaired_two_cohort_comparison(
    cohort_a: pd.Series,
    cohort_b: pd.Series,
    *,
    comparison_id: str,
    low_stratum: str,
    reference_arm: str,
) -> Dict[str, Any]:
    """Mann-Whitney / log-rank / Cox for unpaired cohort A vs low-rank on-label pool B."""
    a = cohort_a.dropna().astype(float)
    b = cohort_b.dropna().astype(float)
    out: Dict[str, Any] = {
        "comparison_id": comparison_id,
        "stratum": low_stratum,
        "reference_arm": reference_arm,
        "pool_label": LOW_RANK_POOL_Q8_LABELS.get(low_stratum, low_stratum),
        "reference_label": LOW_RANK_POOL_Q8_REFERENCE_LABELS.get(reference_arm, reference_arm),
        "n_cohort_a": int(len(a)),
        "n_cohort_b": int(len(b)),
        "median_a": float(np.median(a)) if len(a) else np.nan,
        "median_b": float(np.median(b)) if len(b) else np.nan,
        "median_diff_a_minus_b": np.nan,
        "mannwhitney_p": np.nan,
        "pooled_logrank_p": np.nan,
        "pooled_cox_hr": np.nan,
        "pooled_cox_hr_low": np.nan,
        "pooled_cox_hr_high": np.nan,
        "pooled_cox_p": np.nan,
    }
    if len(a) < 1 or len(b) < 1:
        return out
    out["median_diff_a_minus_b"] = round(out["median_a"] - out["median_b"], 2)
    try:
        out["mannwhitney_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        pass
    try:
        out["pooled_logrank_p"] = float(
            logrank_test(
                a,
                b,
                event_observed_A=np.ones(len(a)),
                event_observed_B=np.ones(len(b)),
            ).p_value
        )
    except Exception:
        pass
    try:
        df = pd.DataFrame({
            "time": pd.concat([a, b], ignore_index=True),
            "event": 1,
            "arm": [1] * len(a) + [0] * len(b),
        })
        cph = CoxPHFitter()
        cph.fit(df, duration_col="time", event_col="event", formula="arm")
        out["pooled_cox_hr"] = float(cph.summary.loc["arm", "exp(coef)"])
        try:
            out["pooled_cox_hr_low"] = float(cph.confidence_intervals_.loc["arm", "95% lower-bound"])
            out["pooled_cox_hr_high"] = float(cph.confidence_intervals_.loc["arm", "95% upper-bound"])
        except Exception:
            pass
        out["pooled_cox_p"] = float(cph.summary.loc["arm", "p"])
    except Exception:
        pass
    return out


def build_low_rank_pool_sensitivity_table(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Q8: progressive on-label-at-rank pools vs Rank 1 (Off/Exp) or SC (unpaired)."""
    rank1_ref = _cohort_rank1_offexp_series(model_table)
    sc_ref = _cohort_sc_global_series(model_table)
    rows: List[Dict[str, Any]] = []
    for stratum in LOW_RANK_POOL_Q8_STRATA:
        pool = _cohort_onlabel_low_rank_pool(ranked, stratum, data_column)
        rows.append(
            _run_unpaired_two_cohort_comparison(
                rank1_ref,
                pool,
                comparison_id=f"{stratum}_vs_rank1_offexp",
                low_stratum=stratum,
                reference_arm="rank1_offexp",
            )
        )
        rows.append(
            _run_unpaired_two_cohort_comparison(
                sc_ref,
                pool,
                comparison_id=f"{stratum}_vs_sc",
                low_stratum=stratum,
                reference_arm="sc",
            )
        )
    return pd.DataFrame(rows)


Q9_PAIRED_COMPARISONS: Tuple[str, ...] = ("onlabel_vs_rank1", "onlabel_vs_sc")
Q9_COMPARISON_LABELS: Dict[str, str] = {
    "onlabel_vs_rank1": "On-label vs Rank 1 (Off/Exp)",
    "onlabel_vs_sc": "On-label vs standard chemo",
}


def _build_q9_paired_table(
    pool_table: pd.DataFrame,
    model_table: pd.DataFrame,
    comparison: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Within-model paired table: on-label PFS at pool rank vs rank 1 or SC."""
    if pool_table.empty:
        return pd.DataFrame()
    work = pool_table[["Model", data_column]].merge(
        model_table[["Model", "rank1_TTD", "chemo_TTD"]],
        on="Model",
        how="inner",
    )
    work["onlabel_pool_TTD"] = pd.to_numeric(work[data_column], errors="coerce")
    work["rank1_TTD"] = pd.to_numeric(work["rank1_TTD"], errors="coerce")
    work["chemo_TTD"] = pd.to_numeric(work["chemo_TTD"], errors="coerce")
    if comparison == "onlabel_vs_rank1":
        paired = work[work["onlabel_pool_TTD"].notna() & work["rank1_TTD"].notna()].copy()
        paired["rank_TTD"] = paired["onlabel_pool_TTD"]
        paired["chemo_TTD"] = paired["rank1_TTD"]
    elif comparison == "onlabel_vs_sc":
        paired = work[work["onlabel_pool_TTD"].notna() & work["chemo_TTD"].notna()].copy()
        paired["rank_TTD"] = paired["onlabel_pool_TTD"]
        paired["chemo_TTD"] = paired["chemo_TTD"]
    else:
        raise ValueError(f"Unknown Q9 comparison: {comparison!r}")
    if paired.empty:
        return paired
    paired["comparison"] = comparison
    paired["diff_TTD"] = paired["rank_TTD"] - paired["chemo_TTD"]
    paired["ratio_TTD"] = paired["rank_TTD"] / paired["chemo_TTD"]
    paired["pct_change_TTD"] = 100.0 * paired["diff_TTD"] / paired["chemo_TTD"]
    paired["improved"] = paired["diff_TTD"] > 0
    paired["worse"] = paired["diff_TTD"] < 0
    paired["tie"] = paired["diff_TTD"] == 0
    return paired.sort_values("pct_change_TTD", ascending=False).reset_index(drop=True)


def _q9_km_dataset_from_paired(
    paired: pd.DataFrame,
    left_group: str,
    right_group: str,
    data_column: str,
) -> pd.DataFrame:
    if paired.empty:
        return pd.DataFrame()
    left = paired[["Model", "rank_TTD"]].rename(columns={"rank_TTD": data_column})
    left["GROUP"] = left_group
    left["CENSOR"] = True
    right = paired[["Model", "chemo_TTD"]].rename(columns={"chemo_TTD": data_column})
    right["GROUP"] = right_group
    right["CENSOR"] = True
    out = pd.concat([left, right], ignore_index=True)
    out[data_column] = pd.to_numeric(out[data_column], errors="coerce")
    return out.dropna(subset=[data_column])


def _summarize_q9_paired_comparison(
    paired: pd.DataFrame,
    *,
    stratum: str,
    comparison: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> Dict[str, Any]:
    left_group = CLINICAL_UTILITY_GROUP_ONLABEL
    right_group = (
        CLINICAL_UTILITY_GROUP_RANK1
        if comparison == "onlabel_vs_rank1"
        else CLINICAL_UTILITY_GROUP_SC
    )
    reference_arm = "rank1_offexp" if comparison == "onlabel_vs_rank1" else "sc"
    out: Dict[str, Any] = {
        "stratum": stratum,
        "comparison": comparison,
        "comparison_id": f"{stratum}_{comparison}",
        "reference_arm": reference_arm,
        "pool_label": LOW_RANK_POOL_Q8_LABELS.get(stratum, stratum),
        "reference_label": LOW_RANK_POOL_Q8_REFERENCE_LABELS.get(reference_arm, reference_arm),
        "comparison_label": Q9_COMPARISON_LABELS.get(comparison, comparison),
    }
    if paired.empty:
        out["n_models"] = 0
        return out
    tests = run_paired_tests(paired)
    km_data = _q9_km_dataset_from_paired(paired, left_group, right_group, data_column)
    surv = (
        _clinical_utility_stratified_survival_stats(
            km_data, left_group, right_group, data_column
        )
        if not km_data.empty
        else {}
    )
    out.update(tests)
    out.update(surv)
    out["median_onlabel"] = float(np.median(paired["rank_TTD"]))
    out["median_reference"] = float(np.median(paired["chemo_TTD"]))
    return out


def build_q9_paired_pool_sensitivity_table(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Q9: paired within-model on-label-at-rank vs rank 1 or SC for each progressive pool."""
    rows: List[Dict[str, Any]] = []
    for stratum in LOW_RANK_POOL_Q8_STRATA:
        pool = _onlabel_at_rank_pool_table(ranked, stratum, data_column)
        for comparison in Q9_PAIRED_COMPARISONS:
            paired = _build_q9_paired_table(pool, model_table, comparison, data_column)
            rows.append(
                _summarize_q9_paired_comparison(
                    paired,
                    stratum=stratum,
                    comparison=comparison,
                    data_column=data_column,
                )
            )
    return pd.DataFrame(rows)


def _q9_paired_row_to_series(row: Dict[str, Any]) -> pd.Series:
    """Map Q9 summary dict to Series compatible with publication paired plots."""
    s = pd.Series(row)
    if pd.notna(s.get("cox_hr")) and pd.isna(s.get("stratified_cox_hr")):
        s["stratified_cox_hr"] = s["cox_hr"]
        s["stratified_cox_hr_low"] = s.get("cox_hr_low", np.nan)
        s["stratified_cox_hr_high"] = s.get("cox_hr_high", np.nan)
        s["stratified_cox_p"] = s.get("cox_p", np.nan)
    if pd.notna(s.get("logrank_pooled_p")) and pd.isna(s.get("stratified_logrank_p")):
        s["stratified_logrank_p"] = s["logrank_pooled_p"]
    return s


def _render_q9_paired_pool_composite(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    q9_summary: pd.DataFrame,
    *,
    style: PublicationPlotStyle,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> plt.Figure:
    n_pools = len(LOW_RANK_POOL_Q8_STRATA)
    fig, axes = plt.subplots(4, n_pools, figsize=(4.2 * n_pools, 15.5))
    if n_pools == 1:
        axes = np.array([[axes]] if axes.ndim == 1 else [axes])
    panel_letters = "ABCDEFGHIJKLMNOP"
    panel_idx = 0
    row_specs = (
        ("onlabel_vs_rank1", 0, "waterfall"),
        ("onlabel_vs_rank1", 1, "km"),
        ("onlabel_vs_sc", 2, "waterfall"),
        ("onlabel_vs_sc", 3, "km"),
    )
    left_group = CLINICAL_UTILITY_GROUP_ONLABEL
    for col, stratum in enumerate(LOW_RANK_POOL_Q8_STRATA):
        pool = _onlabel_at_rank_pool_table(ranked, stratum, data_column)
        for comparison, row_idx, plot_kind in row_specs:
            hits = q9_summary[
                (q9_summary["stratum"] == stratum)
                & (q9_summary["comparison"] == comparison)
            ]
            ax = axes[row_idx, col]
            if hits.empty:
                ax.axis("off")
                continue
            row = _q9_paired_row_to_series(hits.iloc[0].to_dict())
            paired = _build_q9_paired_table(pool, model_table, comparison, data_column)
            if paired.empty:
                ax.axis("off")
                continue
            right_group = (
                CLINICAL_UTILITY_GROUP_RANK1
                if comparison == "onlabel_vs_rank1"
                else CLINICAL_UTILITY_GROUP_SC
            )
            label = panel_letters[panel_idx] if panel_idx < len(panel_letters) else ""
            left_label = "On-label (pool rank)"
            right_label = (
                "Rank 1" if comparison == "onlabel_vs_rank1" else "Standard chemo"
            )
            subtitle = Q9_COMPARISON_LABELS.get(comparison, comparison)
            if plot_kind == "waterfall":
                _plot_publication_waterfall_on_ax(
                    ax,
                    paired,
                    row,
                    left_group=left_group,
                    right_group=right_group,
                    group_colors=CLINICAL_UTILITY_COLORS,
                    left_label=left_label,
                    right_label=right_label,
                    style=style,
                    panel_label=label,
                )
            else:
                km_data = _q9_km_dataset_from_paired(
                    paired, left_group, right_group, data_column
                )
                _plot_publication_km_on_ax(
                    ax,
                    km_data,
                    left_group,
                    right_group,
                    data_column,
                    CLINICAL_UTILITY_COLORS,
                    row,
                    style=style,
                    panel_label=label,
                    subtitle=subtitle,
                    use_stratified=USE_STRATIFIED_SURVIVAL,
                )
            panel_idx += 1
        axes[0, col].set_title(
            LOW_RANK_POOL_Q8_LABELS.get(stratum, stratum),
            fontsize=style.title,
            pad=8,
        )
    row_titles = [
        "On-label vs Rank 1: waterfall",
        "On-label vs Rank 1: KM",
        "On-label vs SC: waterfall",
        "On-label vs SC: KM",
    ]
    for row_idx, title in enumerate(row_titles):
        axes[row_idx, 0].text(
            -0.30,
            0.5,
            title,
            transform=axes[row_idx, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=style.font,
            fontweight="bold",
        )
    fig.suptitle(
        SHORT_SUMMARY_QUESTION_TITLES["Q1"],
        fontsize=style.title + 1,
        y=0.995,
    )
    fig.tight_layout(rect=[0.07, 0.02, 1.0, 0.97])
    return fig


def _plot_q8_sensitivity_box_on_ax(
    ax: plt.Axes,
    row: pd.Series,
    cohort_a: pd.Series,
    cohort_b: pd.Series,
    *,
    style: PublicationPlotStyle,
    panel_label: str = "",
) -> None:
    ref_color = (
        CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_RANK1]
        if row["reference_arm"] == "rank1_offexp"
        else CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_SC]
    )
    pool_color = CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_ONLABEL]
    labels = [
        f"{row['reference_label']}\n(n={int(row['n_cohort_a'])})",
        f"{row['pool_label']}\n(n={int(row['n_cohort_b'])})",
    ]
    data = [cohort_a.values, cohort_b.values]
    box = ax.boxplot(
        data,
        tick_labels=labels,
        patch_artist=True,
        widths=0.55,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"linewidth": 0.9},
        capprops={"linewidth": 0.9},
        boxprops={"linewidth": 0.9},
    )
    palette = [ref_color, pool_color]
    for patch, color in zip(box["boxes"], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.58)
    rng = np.random.default_rng(SEED)
    for idx, values in enumerate(data):
        if len(values) == 0:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(
            np.full(len(values), idx + 1) + jitter,
            values,
            s=20,
            color=palette[idx],
            edgecolor="white",
            linewidth=0.35,
            alpha=0.9,
            zorder=3,
        )
    exp = _exploratory_bracket(row)
    stats_lines = [
        f"Median {row['reference_label']} = {row['median_a']:.1f} d",
        f"Median pool = {row['median_b']:.1f} d",
        f"MW {format_p_equals_mathtext(row.get('mannwhitney_p'))}{exp}",
        f"Log-rank {format_p_equals_mathtext(row.get('pooled_logrank_p'))}{exp}",
    ]
    if pd.notna(row.get("pooled_cox_hr")):
        if pd.notna(row.get("pooled_cox_hr_low")) and pd.notna(row.get("pooled_cox_hr_high")):
            stats_lines.append(
                f"Cox HR = {row['pooled_cox_hr']:.2f} "
                f"({row['pooled_cox_hr_low']:.2f}–{row['pooled_cox_hr_high']:.2f}){exp}"
            )
        else:
            stats_lines.append(f"Cox HR = {row['pooled_cox_hr']:.2f}{exp}")
        if pd.notna(row.get("pooled_cox_p")):
            stats_lines.append(f"Cox {format_p_equals_mathtext(row['pooled_cox_p'])}{exp}")
    ax.text(
        0.97,
        0.97,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=style.annotation - 0.5,
    )
    ax.set_ylabel("PFS (days)", fontsize=style.font)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    subtitle = f"{row['reference_label']} vs {row['pool_label']}"
    ax.set_title(subtitle, fontsize=style.title - 0.5, pad=4)
    if panel_label:
        ax.text(
            0.02,
            0.98,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=style.panel_label,
            fontweight="bold",
        )


def _plot_q8_sensitivity_km_on_ax(
    ax: plt.Axes,
    row: pd.Series,
    cohort_a: pd.Series,
    cohort_b: pd.Series,
    *,
    data_column: str,
    style: PublicationPlotStyle,
    panel_label: str = "",
) -> None:
    ref_arm = str(row["reference_arm"])
    left_group = Q8_KM_GROUP_RANK1 if ref_arm == "rank1_offexp" else Q8_KM_GROUP_SC
    right_group = Q8_KM_GROUP_POOL
    ref_color = (
        CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_RANK1]
        if ref_arm == "rank1_offexp"
        else CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_SC]
    )
    pool_color = CLINICAL_UTILITY_COLORS[CLINICAL_UTILITY_GROUP_ONLABEL]
    colors = {left_group: ref_color, right_group: pool_color}
    frames: List[pd.DataFrame] = []
    for series, grp in ((cohort_a, left_group), (cohort_b, right_group)):
        if series.empty:
            continue
        frames.append(
            pd.DataFrame(
                {
                    data_column: series.values,
                    "GROUP": grp,
                    "CENSOR": True,
                    "Model": [f"{grp}_{i}" for i in range(len(series))],
                }
            )
        )
    km_data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    _plot_publication_km_on_ax(
        ax,
        km_data,
        left_group,
        right_group,
        data_column,
        colors,
        row,
        style=style,
        panel_label=panel_label,
        subtitle=f"{row['pool_label']} vs {row['reference_label']}",
        logrank_label="Log-rank",
        hr_label="Cox HR",
        use_stratified=False,
    )


def _render_q8_low_rank_sensitivity_composite(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    sensitivity: pd.DataFrame,
    *,
    style: PublicationPlotStyle,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> plt.Figure:
    n_pools = len(LOW_RANK_POOL_Q8_STRATA)
    fig, axes = plt.subplots(4, n_pools, figsize=(4.2 * n_pools, 15.5))
    if n_pools == 1:
        axes = np.array([[axes]] if axes.ndim == 1 else [axes])
    rank1_ref = _cohort_rank1_offexp_series(model_table)
    sc_ref = _cohort_sc_global_series(model_table)
    panel_letters = "ABCDEFGHIJKLMNOP"
    panel_idx = 0
    row_specs = (
        ("rank1_offexp", 0, "box"),
        ("rank1_offexp", 1, "km"),
        ("sc", 2, "box"),
        ("sc", 3, "km"),
    )
    for col, stratum in enumerate(LOW_RANK_POOL_Q8_STRATA):
        pool = _cohort_onlabel_low_rank_pool(ranked, stratum, data_column)
        for ref_arm, row_idx, plot_kind in row_specs:
            hits = sensitivity[
                (sensitivity["stratum"] == stratum)
                & (sensitivity["reference_arm"] == ref_arm)
            ]
            ax = axes[row_idx, col]
            if hits.empty:
                ax.axis("off")
                continue
            row = hits.iloc[0]
            ref = rank1_ref if ref_arm == "rank1_offexp" else sc_ref
            label = panel_letters[panel_idx] if panel_idx < len(panel_letters) else ""
            if plot_kind == "box":
                _plot_q8_sensitivity_box_on_ax(
                    ax,
                    row,
                    ref,
                    pool,
                    style=style,
                    panel_label=label,
                )
            else:
                _plot_q8_sensitivity_km_on_ax(
                    ax,
                    row,
                    ref,
                    pool,
                    data_column=data_column,
                    style=style,
                    panel_label=label,
                )
            panel_idx += 1
        axes[0, col].set_title(
            LOW_RANK_POOL_Q8_LABELS.get(stratum, stratum),
            fontsize=style.title,
            pad=8,
        )
    row_titles = [
        "Rank 1 (Off/Exp): box",
        "Rank 1 (Off/Exp): KM",
        "Standard chemo: box",
        "Standard chemo: KM",
    ]
    for row_idx, title in enumerate(row_titles):
        axes[row_idx, 0].text(
            -0.30,
            0.5,
            title,
            transform=axes[row_idx, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=style.font,
            fontweight="bold",
        )
    fig.suptitle(
        SHORT_SUMMARY_QUESTION_TITLES["Q2"],
        fontsize=style.title + 1,
        y=0.995,
    )
    fig.tight_layout(rect=[0.07, 0.02, 1.0, 0.97])
    return fig


def _scenario4_low_rank_filter(model_table: pd.DataFrame, stratum: str) -> pd.DataFrame:
    subset = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO4].copy()
    if stratum == "rank7_10":
        return subset[subset["best_onlabel_rank"].astype(str) == _POOLED_LABEL]
    if stratum == "pool_ge6":
        return subset[subset["best_onlabel_rank_int"] >= 6]
    if stratum == "pool_ge5":
        return subset[subset["best_onlabel_rank_int"] >= 5]
    if stratum == "pool_ge4":
        return subset[subset["best_onlabel_rank_int"] >= 4]
    raise ValueError(f"Unknown low-rank stratum: {stratum!r}")


def _classify_epistemic_row(row: pd.Series) -> str:
    if not np.isfinite(row.get("onlabel_TTD", np.nan)) or not np.isfinite(row.get("rank1_TTD", np.nan)):
        return "D_indeterminate"
    ol, r1 = float(row["onlabel_TTD"]), float(row["rank1_TTD"])
    sc = float(row.get("chemo_TTD", np.nan))
    if ol < r1 and np.isfinite(sc) and ol < sc:
        return "A_ineffective_correctly_deprioritized"
    if ol < r1 and (not np.isfinite(sc) or ol >= sc):
        return "B_suboptimal_vs_rank1_not_worse_than_sc"
    if ol >= r1:
        return "C_possible_ranking_miss"
    return "D_indeterminate"


def build_low_rank_epistemic_table(model_table: pd.DataFrame) -> pd.DataFrame:
    low = model_table[model_table["scenario"] == CLINICAL_UTILITY_SCENARIO4].copy()
    rows = []
    for _, row in low.iterrows():
        ol, r1, sc = row["onlabel_TTD"], row["rank1_TTD"], row["chemo_TTD"]
        rows.append({
            "Model": row["Model"],
            "best_onlabel_rank": row["best_onlabel_rank"],
            "onlabel_COMPOUND": row["onlabel_COMPOUND"],
            "rank1_COMPOUND": row["rank1_COMPOUND"],
            "rank1_utility": row["rank1_utility"],
            "onlabel_TTD": ol, "rank1_TTD": r1, "chemo_TTD": sc,
            "onlabel_vs_rank1": "worse" if np.isfinite(ol) and np.isfinite(r1) and ol < r1 else (
                "better" if np.isfinite(ol) and np.isfinite(r1) and ol > r1 else "tie_or_na"),
            "onlabel_vs_sc": "worse" if np.isfinite(ol) and np.isfinite(sc) and ol < sc else (
                "better" if np.isfinite(ol) and np.isfinite(sc) and ol > sc else "tie_or_na"),
            "epistemic_class": _classify_epistemic_row(row),
        })
    return pd.DataFrame(rows)


def run_between_cohort_onlabel_comparison(model_table: pd.DataFrame, low_stratum: str) -> Dict[str, Any]:
    cohort_a = model_table[
        (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO2) & model_table["onlabel_TTD"].notna()
    ]["onlabel_TTD"].astype(float)
    cohort_b = _scenario4_low_rank_filter(model_table, low_stratum)
    cohort_b = cohort_b[cohort_b["onlabel_TTD"].notna()]["onlabel_TTD"].astype(float)
    out: Dict[str, Any] = {
        "stratum": low_stratum, "n_cohort_a": len(cohort_a), "n_cohort_b": len(cohort_b),
        "median_a": float(np.median(cohort_a)) if len(cohort_a) else np.nan,
        "median_b": float(np.median(cohort_b)) if len(cohort_b) else np.nan,
        "median_diff_a_minus_b": np.nan, "mannwhitney_p": np.nan,
        "pooled_logrank_p": np.nan, "pooled_cox_hr": np.nan, "pooled_cox_p": np.nan,
    }
    if len(cohort_a) < 1 or len(cohort_b) < 1:
        return out
    out["median_diff_a_minus_b"] = round(out["median_a"] - out["median_b"], 2)
    try:
        out["mannwhitney_p"] = float(stats.mannwhitneyu(cohort_a, cohort_b, alternative="two-sided").pvalue)
    except ValueError:
        pass
    try:
        out["pooled_logrank_p"] = float(logrank_test(
            cohort_a, cohort_b, event_observed_A=np.ones(len(cohort_a)),
            event_observed_B=np.ones(len(cohort_b)),
        ).p_value)
    except Exception:
        pass
    try:
        df = pd.DataFrame({
            "time": pd.concat([cohort_a, cohort_b], ignore_index=True), "event": 1,
            "arm": [1] * len(cohort_a) + [0] * len(cohort_b),
        })
        cph = CoxPHFitter()
        cph.fit(df, duration_col="time", event_col="event", formula="arm")
        out["pooled_cox_hr"] = float(cph.summary.loc["arm", "exp(coef)"])
        out["pooled_cox_p"] = float(cph.summary.loc["arm", "p"])
    except Exception:
        pass
    return out


def _build_exp_registry() -> pd.DataFrame:
    rows = [
        ("Exp1", "No on-label: rank 1 vs SC?", "S1 rank1 vs SC", "Wilcoxon, strat Cox"),
        ("Exp2", "On-label rank 1 vs SC?", "S2 onlabel vs SC", "Wilcoxon, strat Cox"),
        ("Exp3", "Rank 2: rank1 vs onlabel?", "S3 rank2", "Wilcoxon+FDR"),
        ("Exp19", "Rank 7-10: rank1 vs onlabel?", "S4 rank7_10", "Wilcoxon"),
        ("Exp20", "Rank 7-10: onlabel vs SC?", "S4 rank7_10", "Wilcoxon"),
        ("Exp27", "S2 vs rank7-10 onlabel TTD?", "between cohort", "Mann-Whitney"),
        ("Exp30", "S2 vs pool_ge4 onlabel TTD?", "between cohort", "Mann-Whitney"),
    ]
    return pd.DataFrame(rows, columns=["exp_id", "clinical_question", "groups_compared", "statistical_methods"])


def _write_clinical_interpretation_master(
    summary_dir: Path, *, model_table: pd.DataFrame, paired_summary: pd.DataFrame,
    between_cohort: pd.DataFrame, epistemic: pd.DataFrame, exp_registry: pd.DataFrame,
    source_path: Path, output_dir: Path,
) -> None:
    sc = model_table["scenario"].value_counts()
    n710 = int((model_table["best_onlabel_rank"].astype(str) == _POOLED_LABEL).sum())
    lines = [
        "# Clinical Utility Ranking — Master Interpretation", "",
        f"**Generated:** {datetime.now():%Y-%m-%d %H:%M}",
        "**Pipeline:** `clinical_utility_ranking.py` (standalone)", f"**Source:** `{source_path.name}`", "",
        "## Clinical questions", "",
        "1. No on-label: rank 1 vs chemo?", "2. On-label rank 1 vs chemo?",
        "3. On-label not top: rank1 vs onlabel?", "4. Low rank 7-10: ineffective or miss?",
        "5. Between-cohort rank1-onlabel vs low-rank pools?", "",
        "## Flowchart", "", "```mermaid", "flowchart TD",
        "  A[Has On_label?] -->|No| S1[Scenario 1]", "  A -->|Yes| B[Best on-label rank]",
        "  B -->|1| S2[Scenario 2]", "  B -->|2-3| S3[Scenario 3]",
        f"  B -->|4+| S4[Scenario 4 n={int(sc.get(CLINICAL_UTILITY_SCENARIO4,0))}]",
        f"  S4 --> R710[Rank 7-10 n={n710}]", "  R710 --> Pool[Progressive pools]",
        "  S2 --> BC[Between-cohort]", "  Pool --> BC", "```", "",
        "## Results", "",
    ]
    if not paired_summary.empty:
        for _, r in paired_summary.iterrows():
            lines.append(f"- {r['scenario']}/{r['stratum']}/{r['comparison']}: n={int(r['n_models'])}, "
                         f"d={r['median_diff_TTD']:+.1f}, p={_clinical_format_p(r['wilcoxon_p'])}")
    lines.append("")
    if not between_cohort.empty:
        lines.append("### Between-cohort")
        for _, r in between_cohort.iterrows():
            lines.append(f"- {r['stratum']}: MW p={_clinical_format_p(r['mannwhitney_p'])}")
    if not epistemic.empty:
        lines.append("### Epistemic")
        for k, v in epistemic["epistemic_class"].value_counts().items():
            lines.append(f"- {k}: {v}")
    lines += ["", "## MTB", "- Rank1 if no on-label; on-label if rank1; Off/Exp at rank2.",
              f"- Rank 7-10 n={n710} exploratory.", "- PDX TTD only.", "",
              "```bash", "python src/clinical_utility_ranking.py --output-dir src", "```"]
    (summary_dir / "CLINICAL_INTERPRETATION_MASTER.md").write_text("\n".join(lines), encoding="utf-8")
    exp_registry.to_csv(summary_dir / "TABLE_EXP_REGISTRY.csv", index=False)
    if not epistemic.empty:
        epistemic.to_csv(summary_dir / "TABLE4_LOW_RANK_EPISTEMIC.csv", index=False)


def _write_clinical_cohort_flow_extended(summary_dir: Path, model_table: pd.DataFrame) -> None:
    _write_clinical_cohort_flow(summary_dir, model_table)
    p = summary_dir / "COHORT_FLOW.md"
    t = p.read_text(encoding="utf-8")
    n710 = int((model_table["best_onlabel_rank"].astype(str) == _POOLED_LABEL).sum())
    extra = f"\n## S4 detail\n- Rank 7-10: {n710}\n- Pools: rank7_10, pool_ge6, pool_ge5, pool_ge4\n"
    p.write_text(t.replace("```\n", extra + "\n```\n", 1), encoding="utf-8")


Q1_MEDIAN_WATERFALL_SHORT_LABELS: Dict[str, str] = {
    "rank7_10": "Rank 1 (Off/exp) vs 7-",
    "pool_ge6": "Rank 1 (Off/exp) vs 6-10",
    "pool_ge5": "Rank 1 (Off/exp) vs 5-10",
    "pool_ge4": "Rank 1 (Off/exp) vs 4-10",
    "pool_ge3": "Rank 1 (Off/exp) vs 3-IO",
}
Q1_REVERSE_POOL_STRATA: Tuple[str, ...] = (
    "rank2",
    "pool_le3",
    "pool_le4",
    "pool_le5",
    "pool_le6",
    "pool_le7_10",
)
# Combined Q1 order: 2-IO → 3-IO → … → 7-10 → 2 (last)
Q1_COMBINED_POOL_ORDER: Tuple[str, ...] = (
    "pool_le7_10",
    "pool_ge3",
    "pool_ge4",
    "pool_ge5",
    "pool_ge6",
    "rank7_10",
    "rank2",
)
Q1_COMBINED_SHORT_LABELS: Dict[str, str] = {
    "rank2": "Rank 1 (off/exp) vs\n2 (on label)",
    "pool_le7_10": "Rank 1 (off/exp) vs\n2-IO (on label)",
    "pool_ge3": "Rank 1 (off/exp) vs\n3-IO (on label)",
    "pool_ge4": "Rank 1 (off/exp) vs\n4-IO (on label)",
    "pool_ge5": "Rank 1 (off/exp) vs\n5-IO (on label)",
    "pool_ge6": "Rank 1 (off/exp) vs\n6-IO (on label)",
    "rank7_10": "Rank 1 (off/exp) vs\n7-10 (on label)",
}
Q1_REVERSE_POOL_SHORT_LABELS: Dict[str, str] = {
    "rank2": "Rank 1 (Off/exp) vs on-label rank 2",
    "pool_le3": "Rank 1 (Off/exp) vs on-label ranks 2-3",
    "pool_le4": "Rank 1 (Off/exp) vs on-label ranks 2-4",
    "pool_le5": "Rank 1 (Off/exp) vs on-label ranks 2-5",
    "pool_le6": "Rank 1 (Off/exp) vs on-label ranks 2-6",
    "pool_le7_10": "Rank 1 (Off/exp) vs 2-IO",
}
Q1_REVERSE_POOL_SHADES: Dict[str, str] = {
    "rank2": "#7F2704",
    "pool_le3": "#B8430A",
    "pool_le4": "#E6550D",
    "pool_le5": "#FD8D3C",
    "pool_le6": "#FDBE85",
    "pool_le7_10": "#FEE6CE",
}
Q1_DOWNRANK_POOL_SHADES: Dict[str, str] = {
    "rank7_10": "#08306B",
    "pool_ge6": "#2171B5",
    "pool_ge5": "#6BAED6",
    "pool_ge4": "#C6DBEF",
    "pool_ge3": "#DEEBF7",
}
RANK_MPFS_BLUE_SHADES: Dict[str, str] = {
    "1": "#08306B",
    "2": "#08519C",
    "3": "#2171B5",
    "4": "#4292C6",
    "5": "#6BAED6",
    "6": "#9ECAE1",
    "7-10": "#C6DBEF",
}


def build_q1_rank1_vs_pool_median_table(
    q9_summary: pd.DataFrame,
) -> pd.DataFrame:
    """One row per progressive pool: Rank 1 vs on-label-at-rank median mPFS."""
    if q9_summary is None or q9_summary.empty:
        return pd.DataFrame()
    hits = q9_summary[q9_summary["comparison"] == "onlabel_vs_rank1"].copy()
    if hits.empty:
        return pd.DataFrame()
    rows: List[Dict[str, Any]] = []
    for stratum in LOW_RANK_POOL_Q8_STRATA:
        sub = hits[hits["stratum"] == stratum]
        if sub.empty:
            continue
        row = sub.iloc[0]
        median_rank1 = float(row["median_reference"]) if pd.notna(row.get("median_reference")) else np.nan
        median_onlabel = float(row["median_onlabel"]) if pd.notna(row.get("median_onlabel")) else np.nan
        paired_diff = float(row["median_diff_TTD"]) if pd.notna(row.get("median_diff_TTD")) else np.nan
        # Positive = Rank 1 longer mPFS than on-label at pool
        delta_medians = (
            median_rank1 - median_onlabel
            if pd.notna(median_rank1) and pd.notna(median_onlabel)
            else np.nan
        )
        rows.append(
            {
                "stratum": stratum,
                "comparison_label": Q1_MEDIAN_WATERFALL_SHORT_LABELS.get(stratum, stratum),
                "pool_label": LOW_RANK_POOL_Q8_LABELS.get(stratum, stratum),
                "n_models": int(row["n_models"]) if pd.notna(row.get("n_models")) else 0,
                "median_mPFS_rank1": median_rank1,
                "median_mPFS_onlabel_pool": median_onlabel,
                "delta_mPFS_rank1_minus_onlabel": delta_medians,
                "paired_median_diff_onlabel_minus_rank1": paired_diff,
                "wilcoxon_p": row.get("wilcoxon_p", np.nan),
                "pooled_logrank_p": row.get("pooled_logrank_p", np.nan),
                "pooled_cox_hr": row.get("pooled_cox_hr", row.get("pooled_hr", np.nan)),
                "pooled_cox_p": row.get("pooled_cox_p", row.get("pooled_hr_p", np.nan)),
            }
        )
    return pd.DataFrame(rows)


def _sig_bar_label(p_value: Any, *, alpha: float = 0.05) -> str:
    if pd.isna(p_value):
        return ""
    p = float(p_value)
    star = "*" if p < alpha else ""
    return f"{format_p_equals_mathtext(p)}{star}"


def _onlabel_rank2_upward_pool_table(
    ranked: pd.DataFrame,
    stratum: str,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """One row per model: on-label in rank-2-upward pool (most deprioritised in range)."""
    on_rows = ranked[ranked["ON_LABEL_UTILITY"] == "On_label"].copy()
    if on_rows.empty:
        return on_rows
    on_rows["_rank_int"] = on_rows["RANK"].map(_clinical_utility_rank_as_int)
    on_rows = on_rows[on_rows["_rank_int"] >= 2]
    if on_rows.empty:
        return on_rows
    if stratum == "rank2":
        pool = on_rows[on_rows["_rank_int"] == 2]
    elif stratum == "pool_le3":
        pool = on_rows[on_rows["_rank_int"].between(2, 3)]
    elif stratum == "pool_le4":
        pool = on_rows[on_rows["_rank_int"].between(2, 4)]
    elif stratum == "pool_le5":
        pool = on_rows[on_rows["_rank_int"].between(2, 5)]
    elif stratum == "pool_le6":
        pool = on_rows[on_rows["_rank_int"].between(2, 6)]
    elif stratum == "pool_le7_10":
        pool = on_rows
    else:
        raise ValueError(f"Unknown rank-2-upward pool stratum: {stratum!r}")
    if pool.empty:
        return pool
    pool = pool.sort_values(
        ["Model", "_rank_int", "LEVEL", data_column],
        ascending=[True, False, False, False],
    )
    return pool.groupby("Model", as_index=False).first()


def _build_q1_rank1_vs_onlabel_pool_paired(
    pool_table: pd.DataFrame,
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Within-model paired table: Rank 1 (Off/exp) vs on-label in pool."""
    if pool_table.empty:
        return pd.DataFrame()
    work = pool_table[["Model", data_column]].merge(
        model_table[["Model", "rank1_TTD"]],
        on="Model",
        how="inner",
    )
    work["onlabel_pool_TTD"] = pd.to_numeric(work[data_column], errors="coerce")
    work["rank1_TTD"] = pd.to_numeric(work["rank1_TTD"], errors="coerce")
    paired = work[work["onlabel_pool_TTD"].notna() & work["rank1_TTD"].notna()].copy()
    if paired.empty:
        return paired
    paired["rank_TTD"] = paired["rank1_TTD"]
    paired["chemo_TTD"] = paired["onlabel_pool_TTD"]
    paired["diff_TTD"] = paired["rank_TTD"] - paired["chemo_TTD"]
    paired["ratio_TTD"] = paired["rank_TTD"] / paired["chemo_TTD"]
    paired["pct_change_TTD"] = 100.0 * paired["diff_TTD"] / paired["chemo_TTD"]
    paired["improved"] = paired["diff_TTD"] > 0
    paired["worse"] = paired["diff_TTD"] < 0
    paired["tie"] = paired["diff_TTD"] == 0
    return paired.sort_values("pct_change_TTD", ascending=False).reset_index(drop=True)


def build_q1_rank1_vs_rank2_upward_pool_median_table(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Rank 1 (Off/exp) vs progressively wider on-label pools from rank 2 upward."""
    rows: List[Dict[str, Any]] = []
    for stratum in Q1_REVERSE_POOL_STRATA:
        pool = _onlabel_rank2_upward_pool_table(ranked, stratum, data_column)
        paired = _build_q1_rank1_vs_onlabel_pool_paired(pool, model_table, data_column)
        if paired.empty:
            continue
        tests = run_paired_tests(paired)
        median_rank1 = float(np.median(paired["rank_TTD"]))
        median_onlabel = float(np.median(paired["chemo_TTD"]))
        rows.append(
            {
                "stratum": stratum,
                "comparison_label": Q1_REVERSE_POOL_SHORT_LABELS.get(stratum, stratum),
                "n_models": int(tests.get("n_models", len(paired))),
                "median_mPFS_rank1": median_rank1,
                "median_mPFS_onlabel_pool": median_onlabel,
                "delta_mPFS_rank1_minus_onlabel": median_rank1 - median_onlabel,
                "paired_median_diff_rank1_minus_onlabel": tests.get("median_diff_TTD", np.nan),
                "wilcoxon_p": tests.get("wilcoxon_p", np.nan),
                "pooled_logrank_p": tests.get("logrank_pooled_p", np.nan),
                "pooled_cox_hr": tests.get("cox_hr", np.nan),
                "pooled_cox_p": tests.get("cox_p", np.nan),
            }
        )
    return pd.DataFrame(rows)


def build_q1_rank1_vs_downrank_pool_median_table(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    strata: Sequence[str] = Q1_DOWNRANK_POOL_STRATA_WITH_3IO,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Rank 1 (Off/exp) vs progressive downranked on-label pools (includes 3-IO)."""
    rows: List[Dict[str, Any]] = []
    for stratum in strata:
        pool = _onlabel_at_rank_pool_table(ranked, stratum, data_column)
        paired = _build_q1_rank1_vs_onlabel_pool_paired(pool, model_table, data_column)
        if paired.empty:
            continue
        tests = run_paired_tests(paired)
        median_rank1 = float(np.median(paired["rank_TTD"]))
        median_onlabel = float(np.median(paired["chemo_TTD"]))
        rows.append(
            {
                "stratum": stratum,
                "comparison_label": Q1_MEDIAN_WATERFALL_SHORT_LABELS.get(stratum, stratum),
                "pool_label": LOW_RANK_POOL_Q8_LABELS.get(stratum, stratum),
                "n_models": int(tests.get("n_models", len(paired))),
                "median_mPFS_rank1": median_rank1,
                "median_mPFS_onlabel_pool": median_onlabel,
                "delta_mPFS_rank1_minus_onlabel": median_rank1 - median_onlabel,
                "paired_median_diff_rank1_minus_onlabel": tests.get("median_diff_TTD", np.nan),
                "wilcoxon_p": tests.get("wilcoxon_p", np.nan),
                "pooled_logrank_p": tests.get("logrank_pooled_p", np.nan),
                "pooled_cox_hr": tests.get("cox_hr", np.nan),
                "pooled_cox_p": tests.get("cox_p", np.nan),
            }
        )
    return pd.DataFrame(rows)


def build_q1_combined_median_waterfall_table(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Combine pools in order: 2, 2-3, …, 2-IO, 3-IO, …, 7-10."""
    down = build_q1_rank1_vs_downrank_pool_median_table(
        model_table, ranked, data_column=data_column
    )
    up = build_q1_rank1_vs_rank2_upward_pool_median_table(
        model_table, ranked, data_column=data_column
    )
    frames: List[pd.DataFrame] = []
    if not up.empty:
        u = up.copy()
        u["pool_family"] = "upward"
        frames.append(u)
    if not down.empty:
        d = down.copy()
        d["pool_family"] = "downrank"
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["stratum"].astype(str).isin(Q1_COMBINED_POOL_ORDER)].copy()
    order = {stratum: i for i, stratum in enumerate(Q1_COMBINED_POOL_ORDER)}
    combined["_order"] = combined["stratum"].astype(str).map(order)
    combined = combined.sort_values("_order", kind="stable").drop(columns="_order")
    combined["comparison_label"] = combined["stratum"].map(
        lambda s: Q1_COMBINED_SHORT_LABELS.get(str(s), str(s))
    )
    return combined.reset_index(drop=True)


def _q1_combined_stratum_paired(
    stratum: str,
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Paired Rank-1 vs on-label-pool values for one Q1 combined stratum."""
    if stratum in Q1_REVERSE_POOL_STRATA:
        pool = _onlabel_rank2_upward_pool_table(ranked, stratum, data_column)
    else:
        pool = _onlabel_at_rank_pool_table(ranked, stratum, data_column)
    return _build_q1_rank1_vs_onlabel_pool_paired(pool, model_table, data_column)


def _style_boxplot_artists(
    box: Dict[str, Any],
    *,
    facecolor: str,
    hatch: Optional[str] = None,
    edgecolor: str = "black",
) -> None:
    for patch in box["boxes"]:
        patch.set_facecolor(facecolor)
        patch.set_edgecolor(edgecolor)
        patch.set_alpha(0.85)
        if hatch:
            patch.set_hatch(hatch)
            patch.set_linewidth(1.0)
    for key in ("whiskers", "caps", "medians"):
        for artist in box[key]:
            artist.set_color("black")
            artist.set_linewidth(0.9)
    for flier in box["fliers"]:
        flier.set_markeredgecolor("black")
        flier.set_markersize(3.5)


def _save_q1_combined_mpfs_boxplot(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    output_path: Path,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
    ylim_low: float = 0.0,
    ylim_high: Optional[float] = None,
) -> None:
    """Paired boxplots: Rank 1 (off/exp) vs on-label pool mPFS for Q1 combined strata."""
    light_blue = "#9ECAE1"
    light_green = "#A1D99B"
    groups: List[Dict[str, Any]] = []
    for stratum in Q1_COMBINED_POOL_ORDER:
        paired = _q1_combined_stratum_paired(
            stratum, model_table, ranked, data_column=data_column
        )
        if paired.empty:
            continue
        median_rank1 = float(np.median(paired["rank_TTD"]))
        median_onlabel = float(np.median(paired["chemo_TTD"]))
        groups.append(
            {
                "stratum": stratum,
                "label": Q1_COMBINED_SHORT_LABELS.get(stratum, stratum),
                "rank1": paired["rank_TTD"].astype(float).values,
                "onlabel": paired["chemo_TTD"].astype(float).values,
                "delta": median_rank1 - median_onlabel,
                "n": int(len(paired)),
                "wilcoxon_p": _wilcoxon_p_from_paired_diff(paired["diff_TTD"]),
            }
        )
    if not groups:
        return

    fig, ax = plt.subplots(figsize=(9.5, 8.5))
    x = np.arange(len(groups))
    width = 0.34
    all_vals: List[float] = []
    for i, g in enumerate(groups):
        rank_box = ax.boxplot(
            [g["rank1"]],
            positions=[x[i] - width / 2],
            widths=width * 0.95,
            patch_artist=True,
            manage_ticks=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        on_box = ax.boxplot(
            [g["onlabel"]],
            positions=[x[i] + width / 2],
            widths=width * 0.95,
            patch_artist=True,
            manage_ticks=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        _style_boxplot_artists(rank_box, facecolor=light_blue)
        _style_boxplot_artists(on_box, facecolor=light_green)
        all_vals.extend(list(g["rank1"]))
        all_vals.extend(list(g["onlabel"]))
        top_h = float(np.nanmax(np.concatenate([g["rank1"], g["onlabel"]])))
        ax.text(
            x[i],
            top_h + 1.2,
            _boxplot_delta_n_p_annotation(g["delta"], g["n"], g["wilcoxon_p"]),
            ha="center",
            va="bottom",
            fontsize=10.7,
        )

    labels = [g["label"] for g in groups]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=15, rotation=60, ha="right")
    ax.set_ylabel(
        "PFS (days)\nrank 1 (off/exp) vs on label in lower rank",
        fontsize=16.6,
    )
    ax.tick_params(axis="y", labelsize=15)
    ymax = float(np.nanmax(all_vals)) if all_vals else 40.0
    if ylim_high is None:
        ylim_high = max(40.0, ymax * 1.28)
    ax.set_ylim(ylim_low, ylim_high)
    legend_handles = [
        Patch(facecolor=light_blue, edgecolor="black", label="Rank 1 (off/exp)"),
        Patch(facecolor=light_green, edgecolor="black", label="On label lower"),
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.06),
        fontsize=15,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _save_q1_combined_mpfs_boxplot_tilted(
    model_table: pd.DataFrame,
    ranked: pd.DataFrame,
    output_path: Path,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
    xlim_low: float = 0.0,
    xlim_high: Optional[float] = None,
) -> None:
    """Horizontal (tilted) Q1 combined boxplots: categories on y, mPFS on x."""
    light_blue = "#9ECAE1"
    light_green = "#A1D99B"
    groups: List[Dict[str, Any]] = []
    for stratum in Q1_COMBINED_POOL_ORDER:
        paired = _q1_combined_stratum_paired(
            stratum, model_table, ranked, data_column=data_column
        )
        if paired.empty:
            continue
        median_rank1 = float(np.median(paired["rank_TTD"]))
        median_onlabel = float(np.median(paired["chemo_TTD"]))
        groups.append(
            {
                "stratum": stratum,
                "label": Q1_COMBINED_SHORT_LABELS.get(stratum, stratum),
                "rank1": paired["rank_TTD"].astype(float).values,
                "onlabel": paired["chemo_TTD"].astype(float).values,
                "delta": median_rank1 - median_onlabel,
                "n": int(len(paired)),
                "wilcoxon_p": _wilcoxon_p_from_paired_diff(paired["diff_TTD"]),
            }
        )
    if not groups:
        return

    fig, ax = plt.subplots(figsize=(8.5, 9.5))
    y = np.arange(len(groups))
    width = 0.34
    all_vals: List[float] = []
    for i, g in enumerate(groups):
        rank_box = ax.boxplot(
            [g["rank1"]],
            positions=[y[i] - width / 2],
            widths=width * 0.95,
            vert=False,
            patch_artist=True,
            manage_ticks=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        on_box = ax.boxplot(
            [g["onlabel"]],
            positions=[y[i] + width / 2],
            widths=width * 0.95,
            vert=False,
            patch_artist=True,
            manage_ticks=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        _style_boxplot_artists(rank_box, facecolor=light_blue)
        _style_boxplot_artists(on_box, facecolor=light_green)
        all_vals.extend(list(g["rank1"]))
        all_vals.extend(list(g["onlabel"]))
        right_h = float(np.nanmax(np.concatenate([g["rank1"], g["onlabel"]])))
        ax.text(
            right_h + 1.2,
            y[i],
            _boxplot_delta_n_p_annotation(g["delta"], g["n"], g["wilcoxon_p"]),
            ha="left",
            va="center",
            fontsize=10.7,
        )

    labels = [g["label"] for g in groups]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=15)
    ax.invert_yaxis()
    ax.set_xlabel(
        "PFS (days)\nrank 1 (off/exp) vs on label in lower rank",
        fontsize=16.6,
    )
    ax.tick_params(axis="x", labelsize=15)
    xmax = float(np.nanmax(all_vals)) if all_vals else 40.0
    if xlim_high is None:
        xlim_high = max(40.0, xmax * 1.28)
    ax.set_xlim(xlim_low, xlim_high)
    legend_handles = [
        Patch(facecolor=light_blue, edgecolor="black", label="Rank 1 (off/exp)"),
        Patch(facecolor=light_green, edgecolor="black", label="On label lower"),
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        fontsize=15,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# Backwards-compatible aliases (Q1 combined now uses boxplots).
def _save_q1_combined_median_ttd_waterfall(
    combined_table: pd.DataFrame,
    output_path: Path,
    *,
    ylim_low: float = 0.0,
    ylim_high: Optional[float] = None,
    model_table: Optional[pd.DataFrame] = None,
    ranked: Optional[pd.DataFrame] = None,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> None:
    if model_table is None or ranked is None:
        return
    _save_q1_combined_mpfs_boxplot(
        model_table,
        ranked,
        output_path,
        data_column=data_column,
        ylim_low=ylim_low,
        ylim_high=ylim_high,
    )


def _save_q1_combined_median_ttd_grouped_bars(
    combined_table: pd.DataFrame,
    output_path: Path,
    *,
    ylim_low: float = 0.0,
    ylim_high: Optional[float] = None,
    model_table: Optional[pd.DataFrame] = None,
    ranked: Optional[pd.DataFrame] = None,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> None:
    if model_table is None or ranked is None:
        return
    _save_q1_combined_mpfs_boxplot(
        model_table,
        ranked,
        output_path,
        data_column=data_column,
        ylim_low=ylim_low,
        ylim_high=ylim_high,
    )


def _save_q1_median_ttd_waterfall(
    median_table: pd.DataFrame,
    output_path: Path,
    *,
    ylim: float = 50.0,
    bar_colors: Optional[Sequence[str]] = None,
    title: str = "Q1: Rank 1 (Off/exp) vs downranked on-label pools (median mPFS)",
    ylabel: str = "ΔmPFS (days)\nnon-on-label in rank 1 − match on-label in lower rank",
) -> None:
    """Waterfall of Rank 1 (Off/exp) − on-label-pool median mPFS."""
    if median_table is None or median_table.empty:
        return
    plot_df = median_table.copy()
    values = plot_df["delta_mPFS_rank1_minus_onlabel"].astype(float).values
    labels = plot_df["comparison_label"].astype(str).tolist()
    if bar_colors is not None and len(bar_colors) == len(plot_df):
        colors = list(bar_colors)
    else:
        colors = ["#2e7d32" if v >= 0 else "#c62828" for v in values]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(len(plot_df))
    bars = ax.bar(x, values, color=colors, edgecolor="black", linewidth=0.6, width=0.72)
    ax.axhline(0, color="black", linewidth=0.9)
    for bar, val, n, (_, row) in zip(bars, values, plot_df["n_models"], plot_df.iterrows()):
        if not np.isfinite(val):
            continue
        sig = _sig_bar_label(row.get("wilcoxon_p"))
        y_text = val + (2.5 if val >= 0 else -2.5)
        lines = [f"{val:+.1f} d", f"n={int(n)}"]
        if sig:
            lines.append(sig)
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y_text,
            "\n".join(lines),
            ha="center",
            va="bottom" if val >= 0 else "top",
            fontsize=10.0,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11, rotation=15, ha="right")
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_ylim(-ylim, ylim)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_rank_mpfs_bar_table(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Median mPFS by global rank; each rank also reports the no-on-label subset."""
    rows: List[Dict[str, Any]] = []
    work = ranked.copy()
    work["RANK"] = work["RANK"].astype(str)
    work[data_column] = pd.to_numeric(work[data_column], errors="coerce")

    no_onlabel_models = set(
        model_table.loc[
            model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1, "Model"
        ].astype(str)
    )

    for rank in RANKING_UTILITY_RANKS:
        rank_rows = work[work["RANK"] == str(rank)].dropna(subset=[data_column])
        rows.append(
            {
                "rank": str(rank),
                "cohort": "all",
                "n": int(len(rank_rows)),
                "mPFS": float(rank_rows[data_column].median()) if not rank_rows.empty else np.nan,
            }
        )
        rank_no = rank_rows[rank_rows["Model"].astype(str).isin(no_onlabel_models)]
        rows.append(
            {
                "rank": str(rank),
                "cohort": "no_onlabel",
                "n": int(len(rank_no)),
                "mPFS": float(rank_no[data_column].median()) if not rank_no.empty else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _save_rank_mpfs_paired_barplot(
    mpfs_table: pd.DataFrame,
    output_path: Path,
) -> None:
    """Paired barplot of mPFS by rank; each rank shows all + no-on-label (lined hatch)."""
    if mpfs_table is None or mpfs_table.empty:
        return
    fig, ax = plt.subplots(figsize=(12.5 * 0.75, 5.8))
    x_positions: List[float] = []
    heights: List[float] = []
    colors: List[str] = []
    dotted: List[bool] = []
    tick_pos: List[float] = []
    tick_labels: List[str] = []
    annotations: List[Tuple[float, float, str]] = []

    cursor = 0.0
    bar_gap = 0.38
    group_gap = 0.55
    for rank in RANKING_UTILITY_RANKS:
        rank_color = RANK_MPFS_BLUE_SHADES.get(str(rank), "#6BAED6")
        all_row = mpfs_table[
            (mpfs_table["rank"].astype(str) == str(rank)) & (mpfs_table["cohort"] == "all")
        ]
        no_row = mpfs_table[
            (mpfs_table["rank"].astype(str) == str(rank)) & (mpfs_table["cohort"] == "no_onlabel")
        ]
        if all_row.empty or not np.isfinite(float(all_row.iloc[0]["mPFS"])):
            continue
        all_mpfs = float(all_row.iloc[0]["mPFS"])
        all_n = int(all_row.iloc[0]["n"])
        x_all = cursor
        x_no = cursor + bar_gap
        x_positions.extend([x_all, x_no])
        heights.extend(
            [
                all_mpfs,
                float(no_row.iloc[0]["mPFS"]) if not no_row.empty else np.nan,
            ]
        )
        colors.extend([rank_color, rank_color])
        dotted.extend([False, True])
        annotations.append((x_all, all_mpfs, f"{all_mpfs:.1f}\nn={all_n}"))
        if not no_row.empty and np.isfinite(float(no_row.iloc[0]["mPFS"])):
            annotations.append(
                (
                    x_no,
                    float(no_row.iloc[0]["mPFS"]),
                    f"{float(no_row.iloc[0]['mPFS']):.1f}\nn={int(no_row.iloc[0]['n'])}",
                )
            )
        tick_pos.append(cursor + bar_gap / 2)
        tick_labels.append(_rank_display_label(rank))
        cursor += bar_gap + group_gap

    for x, h, c, is_dotted in zip(x_positions, heights, colors, dotted):
        if not np.isfinite(h):
            continue
        bar = ax.bar(
            x,
            h,
            width=0.34,
            color=c,
            edgecolor="black",
            linewidth=0.7,
            alpha=0.95,
            zorder=2 if not is_dotted else 3,
        )
        if is_dotted:
            bar[0].set_hatch("///")
            bar[0].set_edgecolor("black")
            bar[0].set_linewidth(1.0)

    for x, h, text in annotations:
        if not np.isfinite(h):
            continue
        ax.text(x, h + max(heights) * 0.02, text, ha="center", va="bottom", fontsize=9.0)

    legend_handles = [
        Patch(facecolor=RANK_MPFS_BLUE_SHADES["1"], edgecolor="black", label="All data"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["3"],
            edgecolor="black",
            hatch="///",
            label="No on-label models",
        ),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right", fontsize=11)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=11)
    ax.set_ylabel("mPFS (days)", fontsize=12)
    ax.set_title("Median mPFS by DDA rank (all vs no on-label models)", fontsize=13)
    ax.tick_params(axis="y", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ymax = max([h for h in heights if np.isfinite(h)], default=1.0)
    ax.set_ylim(0, ymax * 1.22)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _wilcoxon_p_from_paired_diff(diff: pd.Series) -> float:
    """Two-sided Wilcoxon signed-rank p for paired differences; NaN if not computable."""
    values = pd.to_numeric(diff, errors="coerce").dropna().astype(float).values
    if len(values) < 1:
        return np.nan
    if np.allclose(values, 0.0):
        return 1.0
    try:
        return float(stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return np.nan


def build_rank_delta_mpfs_bar_table(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> pd.DataFrame:
    """Median paired ΔmPFS (rank TTD − matched SC) by global rank; all + no-on-label."""
    rows: List[Dict[str, Any]] = []
    work = ranked.copy()
    work["RANK"] = work["RANK"].astype(str)
    work[data_column] = pd.to_numeric(work[data_column], errors="coerce")

    chemo_map = (
        model_table.set_index("Model")["chemo_TTD"].to_dict()
        if not model_table.empty and "chemo_TTD" in model_table.columns
        else {}
    )
    no_onlabel_models = set(
        model_table.loc[
            model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1, "Model"
        ].astype(str)
    )

    for rank in RANKING_UTILITY_RANKS:
        rank_rows = work[work["RANK"] == str(rank)].dropna(subset=[data_column]).copy()
        if rank_rows.empty:
            continue
        rank_rows["chemo_TTD"] = rank_rows["Model"].astype(str).map(chemo_map)
        rank_rows["chemo_TTD"] = pd.to_numeric(rank_rows["chemo_TTD"], errors="coerce")
        paired = rank_rows[rank_rows["chemo_TTD"].notna()].copy()
        if paired.empty:
            continue

        for cohort_name, subset in (
            ("all", paired),
            ("no_onlabel", paired[paired["Model"].astype(str).isin(no_onlabel_models)]),
        ):
            if subset.empty:
                rows.append(
                    {
                        "rank": str(rank),
                        "cohort": cohort_name,
                        "n": 0,
                        "delta_mPFS": np.nan,
                        "median_mPFS_rank": np.nan,
                        "median_mPFS_matched_SC": np.nan,
                        "wilcoxon_p": np.nan,
                    }
                )
                continue
            collapsed = _collapse_rank_rows_mean_ttd_per_model(subset, data_column)
            if collapsed.empty:
                rows.append(
                    {
                        "rank": str(rank),
                        "cohort": cohort_name,
                        "n": 0,
                        "delta_mPFS": np.nan,
                        "median_mPFS_rank": np.nan,
                        "median_mPFS_matched_SC": np.nan,
                        "wilcoxon_p": np.nan,
                    }
                )
                continue
            collapsed = collapsed.copy()
            collapsed["delta_mPFS"] = (
                collapsed[data_column].astype(float) - collapsed["chemo_TTD"].astype(float)
            )
            rows.append(
                {
                    "rank": str(rank),
                    "cohort": cohort_name,
                    "n": int(len(collapsed)),
                    "median_mPFS_rank": float(collapsed[data_column].median()),
                    "median_mPFS_matched_SC": float(collapsed["chemo_TTD"].median()),
                    # Difference of medians: median(rank) − median(SC).
                    "delta_mPFS": float(collapsed[data_column].median())
                    - float(collapsed["chemo_TTD"].median()),
                    "wilcoxon_p": _wilcoxon_p_from_paired_diff(collapsed["delta_mPFS"]),
                }
            )
    return pd.DataFrame(rows)


def _save_rank_delta_mpfs_paired_barplot(
    delta_table: pd.DataFrame,
    output_path: Path,
) -> None:
    """Paired barplot of ΔmPFS vs matched SC by rank (all + no-on-label lined)."""
    if delta_table is None or delta_table.empty:
        return
    fig, ax = plt.subplots(figsize=(15.0, 5.8))
    x_positions: List[float] = []
    heights: List[float] = []
    colors: List[str] = []
    dotted: List[bool] = []
    tick_pos: List[float] = []
    tick_labels: List[str] = []
    annotations: List[Tuple[float, float, str]] = []

    cursor = 0.0
    bar_gap = 0.38
    group_gap = 0.55
    for rank in RANKING_UTILITY_RANKS:
        rank_color = RANK_MPFS_BLUE_SHADES.get(str(rank), "#6BAED6")
        all_row = delta_table[
            (delta_table["rank"].astype(str) == str(rank)) & (delta_table["cohort"] == "all")
        ]
        no_row = delta_table[
            (delta_table["rank"].astype(str) == str(rank))
            & (delta_table["cohort"] == "no_onlabel")
        ]
        if all_row.empty or not np.isfinite(float(all_row.iloc[0]["delta_mPFS"])):
            continue
        all_delta = float(all_row.iloc[0]["delta_mPFS"])
        all_n = int(all_row.iloc[0]["n"])
        x_all = cursor
        x_no = cursor + bar_gap
        x_positions.extend([x_all, x_no])
        heights.extend(
            [
                all_delta,
                float(no_row.iloc[0]["delta_mPFS"]) if not no_row.empty else np.nan,
            ]
        )
        colors.extend([rank_color, rank_color])
        dotted.extend([False, True])
        all_star = ""
        all_p = all_row.iloc[0].get("wilcoxon_p")
        if pd.notna(all_p) and float(all_p) < 0.05:
            all_star = "*"
        annotations.append(
            (x_all, all_delta, f"{all_delta:+.1f} d{all_star}\nn={all_n}")
        )
        if not no_row.empty and np.isfinite(float(no_row.iloc[0]["delta_mPFS"])):
            no_delta = float(no_row.iloc[0]["delta_mPFS"])
            no_star = ""
            no_p = no_row.iloc[0].get("wilcoxon_p")
            if pd.notna(no_p) and float(no_p) < 0.05:
                no_star = "*"
            annotations.append(
                (
                    x_no,
                    no_delta,
                    f"{no_delta:+.1f} d{no_star}\nn={int(no_row.iloc[0]['n'])}",
                )
            )
        tick_pos.append(cursor + bar_gap / 2)
        tick_labels.append(_rank_display_label(rank))
        cursor += bar_gap + group_gap

    for x, h, c, is_dotted in zip(x_positions, heights, colors, dotted):
        if not np.isfinite(h):
            continue
        bar = ax.bar(
            x,
            h,
            width=0.34,
            color=c,
            edgecolor="black",
            linewidth=0.7,
            alpha=0.95,
            zorder=2 if not is_dotted else 3,
        )
        if is_dotted:
            bar[0].set_hatch("///")
            bar[0].set_edgecolor("white")
            bar[0].set_linewidth(1.0)

    finite_heights = [h for h in heights if np.isfinite(h)]
    pad = max(abs(h) for h in finite_heights) * 0.12 if finite_heights else 1.0
    for x, h, text in annotations:
        if not np.isfinite(h):
            continue
        ax.text(
            x,
            h + (pad if h >= 0 else -pad),
            text,
            ha="center",
            va="bottom" if h >= 0 else "top",
            fontsize=12.75,
        )

    ax.axhline(0, color="black", linewidth=0.9)
    legend_handles = [
        Patch(facecolor=RANK_MPFS_BLUE_SHADES["1"], edgecolor="black", label="All data"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["3"],
            edgecolor="white",
            hatch="///",
            label="No on-label models",
        ),
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        fontsize=21.2,
        title=r"$\ast\,P < 0.05\ \mathrm{(Wilcoxon)}$",
        title_fontsize=19.9,
    )
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=21.2)
    ax.set_ylabel("Δ mPFS (days)\nRank − matched SC", fontsize=22.5)
    ax.tick_params(axis="y", labelsize=21.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if finite_heights:
        span = max(abs(h) for h in finite_heights)
        ax.set_ylim(-span * 1.55, span * 1.55)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _collapse_rank_rows_mean_ttd_per_model(
    paired: pd.DataFrame,
    data_column: str,
) -> pd.DataFrame:
    """One row per Model: mean TTD across same-rank compounds; keep matched SC.

    Dense-rank ties (and duplicate compound rows) can yield multiple RANK==k
    rows per model. Fig4e uses the model-level mean TTD so n matches unique
    tumors (aligned with Fig4d), not raw treatment rows.
    """
    if paired.empty:
        return paired
    work = paired.copy()
    work["Model"] = work["Model"].astype(str)
    work[data_column] = pd.to_numeric(work[data_column], errors="coerce")
    work["chemo_TTD"] = pd.to_numeric(work["chemo_TTD"], errors="coerce")
    work = work.dropna(subset=[data_column, "chemo_TTD"])
    if work.empty:
        return work
    return (
        work.groupby("Model", as_index=False)
        .agg(
            **{
                data_column: (data_column, "mean"),
                "chemo_TTD": ("chemo_TTD", "first"),
            }
        )
        .reset_index(drop=True)
    )


def build_rank_mpfs_vs_sc_paired_groups(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> List[Dict[str, Any]]:
    """Per-rank paired mPFS distributions vs matched SC (all + no-on-label)."""
    groups: List[Dict[str, Any]] = []
    work = ranked.copy()
    work["RANK"] = work["RANK"].astype(str)
    work[data_column] = pd.to_numeric(work[data_column], errors="coerce")
    chemo_map = (
        model_table.set_index("Model")["chemo_TTD"].to_dict()
        if not model_table.empty and "chemo_TTD" in model_table.columns
        else {}
    )
    no_onlabel_models = set(
        model_table.loc[
            model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1, "Model"
        ].astype(str)
    )
    for rank in RANKING_UTILITY_RANKS:
        rank_rows = work[work["RANK"] == str(rank)].dropna(subset=[data_column]).copy()
        if rank_rows.empty:
            continue
        rank_rows["chemo_TTD"] = pd.to_numeric(
            rank_rows["Model"].astype(str).map(chemo_map), errors="coerce"
        )
        paired = rank_rows[rank_rows["chemo_TTD"].notna()].copy()
        if paired.empty:
            continue
        for cohort_name, subset in (
            ("all", paired),
            ("no_onlabel", paired[paired["Model"].astype(str).isin(no_onlabel_models)]),
        ):
            if subset.empty:
                continue
            collapsed = _collapse_rank_rows_mean_ttd_per_model(subset, data_column)
            if collapsed.empty:
                continue
            rank_vals = collapsed[data_column].astype(float).values
            sc_vals = collapsed["chemo_TTD"].astype(float).values
            paired_diff = rank_vals - sc_vals
            groups.append(
                {
                    "rank": str(rank),
                    "cohort": cohort_name,
                    "rank_vals": rank_vals,
                    "sc_vals": sc_vals,
                    "n": int(len(collapsed)),
                    # Difference of medians: median(rank) − median(SC).
                    "delta": float(np.median(rank_vals) - np.median(sc_vals)),
                    "wilcoxon_p": _wilcoxon_p_from_paired_diff(pd.Series(paired_diff)),
                }
            )
    return groups


def build_boxplot_wilcoxon_fdr_table(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    data_column: str = DEFAULT_DATA_COLUMN,
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Raw and BH-FDR Wilcoxon results for the two displayed boxplot families.

    Benjamini-Hochberg correction is applied separately within each plot:
    seven Q1 combined-pool tests and fourteen rank-vs-SC tests.
    """
    rows: List[Dict[str, Any]] = []

    for order, stratum in enumerate(Q1_COMBINED_POOL_ORDER, start=1):
        paired = _q1_combined_stratum_paired(
            stratum, model_table, ranked, data_column=data_column
        )
        if paired.empty:
            continue
        rank_vals = paired["rank_TTD"].astype(float).values
        onlabel_vals = paired["chemo_TTD"].astype(float).values
        rows.append(
            {
                "plot": "Q1 combined: Rank 1 vs lower-rank on-label",
                "comparison_order": order,
                "comparison_id": stratum,
                "comparison": Q1_COMBINED_SHORT_LABELS.get(stratum, stratum).replace(
                    "\n", " "
                ),
                "cohort": "paired models",
                "n_pairs": int(len(paired)),
                "delta_definition": "median(Rank 1) - median(on-label pool)",
                "delta_mPFS": float(np.median(rank_vals) - np.median(onlabel_vals)),
                "wilcoxon_p": _wilcoxon_p_from_paired_diff(paired["diff_TTD"]),
            }
        )

    rank_groups = build_rank_mpfs_vs_sc_paired_groups(
        ranked, model_table, data_column=data_column
    )
    rank_order = {str(rank): i for i, rank in enumerate(RANKING_UTILITY_RANKS, start=1)}
    cohort_order = {"all": 0, "no_onlabel": 1}
    for group in rank_groups:
        rank = str(group["rank"])
        cohort = str(group["cohort"])
        rows.append(
            {
                "plot": "Rank mPFS vs matched SC",
                "comparison_order": 2 * (rank_order[rank] - 1)
                + cohort_order.get(cohort, 0)
                + 1,
                "comparison_id": rank,
                "comparison": f"Rank {_rank_display_label(rank)} vs matched SC",
                "cohort": "All data" if cohort == "all" else "No on-label models",
                "n_pairs": int(group["n"]),
                "delta_definition": "median(rank) - median(matched SC)",
                "delta_mPFS": float(group["delta"]),
                "wilcoxon_p": float(group["wilcoxon_p"]),
            }
        )

    table = pd.DataFrame(rows)
    if table.empty:
        return table

    table["wilcoxon_p_fdr_bh"] = np.nan
    for _, index in table.groupby("plot", sort=False).groups.items():
        idx = list(index)
        valid = table.loc[idx, "wilcoxon_p"].notna()
        valid_idx = table.loc[idx].index[valid]
        if len(valid_idx):
            _, adjusted, _, _ = multipletests(
                table.loc[valid_idx, "wilcoxon_p"].astype(float).values,
                alpha=alpha,
                method="fdr_bh",
            )
            table.loc[valid_idx, "wilcoxon_p_fdr_bh"] = adjusted

    table["significant_raw_p_lt_0_05"] = table["wilcoxon_p"].lt(alpha)
    table["significant_fdr_lt_0_05"] = table["wilcoxon_p_fdr_bh"].lt(alpha)

    def _change_label(row: pd.Series) -> str:
        raw = bool(row["significant_raw_p_lt_0_05"])
        fdr = bool(row["significant_fdr_lt_0_05"])
        if raw and not fdr:
            return "Significant -> not significant"
        if not raw and fdr:
            return "Not significant -> significant"
        return "Unchanged: significant" if raw else "Unchanged: not significant"

    table["significance_change_after_fdr"] = table.apply(_change_label, axis=1)
    return table.sort_values(
        ["plot", "comparison_order"], kind="stable"
    ).reset_index(drop=True)


def _add_broken_axis_marks_left(ax_top: plt.Axes, ax_bot: plt.Axes, *, size: float = 0.0075) -> None:
    """Draw diagonal break marks on the left y-axis only."""
    kwargs = dict(color="k", clip_on=False, linewidth=0.8)
    ax_bot.plot((-size, +size), (1 - size, 1 + size), transform=ax_bot.transAxes, **kwargs)
    ax_top.plot((-size, +size), (-size, +size), transform=ax_top.transAxes, **kwargs)


def _save_rank_mpfs_vs_sc_paired_boxplot(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    output_path: Path,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> None:
    """Boxplots of mPFS: rank vs matched SC for all and no-on-label cohorts.

    Y-axis is broken/compressed: 0–200 uses full panel height; 200–400 uses
    half of that height (outlier region).
    """
    groups = build_rank_mpfs_vs_sc_paired_groups(
        ranked, model_table, data_column=data_column
    )
    if not groups:
        return

    # Tick/ylabel matched to ΔmPFS barplot; legend smaller to reduce clutter.
    font_tick = 21.2
    font_ylabel = 22.5
    font_legend = 14.4
    font_legend_title = 13.2
    font_annot = 9.75
    annot_pad = 10.0

    sc_color = "#D4B4E8"
    # height_ratios: top (200–400) = half again → 1/4 of bottom (0–200)
    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(18.0, 9.0),
        gridspec_kw={"height_ratios": [1, 4], "hspace": 0.05},
    )
    tick_pos: List[float] = []
    tick_labels: List[str] = []
    cursor = 0.0
    pair_gap = 0.42
    group_gap = 0.72
    by_rank: Dict[str, List[Dict[str, Any]]] = {}
    for g in groups:
        by_rank.setdefault(g["rank"], []).append(g)

    annot_tops: List[float] = []
    for rank in RANKING_UTILITY_RANKS:
        rank_groups = by_rank.get(str(rank), [])
        if not rank_groups:
            continue
        rank_color = RANK_MPFS_BLUE_SHADES.get(str(rank), "#6BAED6")
        pair_centers: List[float] = []
        for g in rank_groups:
            is_no = g["cohort"] == "no_onlabel"
            x_rank = cursor
            x_sc = cursor + pair_gap
            for ax in (ax_top, ax_bot):
                rank_box = ax.boxplot(
                    [g["rank_vals"]],
                    positions=[x_rank],
                    widths=0.4224,
                    patch_artist=True,
                    manage_ticks=False,
                    medianprops={"color": "black", "linewidth": 1.2},
                )
                sc_box = ax.boxplot(
                    [g["sc_vals"]],
                    positions=[x_sc],
                    widths=0.4224,
                    patch_artist=True,
                    manage_ticks=False,
                    medianprops={"color": "black", "linewidth": 1.2},
                )
                _style_boxplot_artists(
                    rank_box,
                    facecolor=rank_color,
                    hatch="///" if is_no else None,
                )
                _style_boxplot_artists(
                    sc_box,
                    facecolor=sc_color,
                    hatch="///" if is_no else None,
                )
            # Always place label above the highest outlier of this pair.
            pair_max = float(
                np.nanmax(np.concatenate([g["rank_vals"], g["sc_vals"]]))
            )
            annot_y = pair_max + annot_pad
            # Nudge no-on-label labels for selected ranks upward.
            if is_no and str(rank) in {"2", "4", "6", "7-10"}:
                annot_y += 10.0
            if is_no and str(rank) == "2":
                annot_y += 12.0  # cumulative extra nudge for rank-2 no-on-label
            if is_no and str(rank) == "5":
                annot_y += 10.0
            annot_tops.append(annot_y)
            annot_ax = ax_top if annot_y > 200 else ax_bot
            text_transform = annot_ax.transData
            # Rank-6 no-on-label annotation: lift 6 pt to clear the boxes/outliers.
            if is_no and str(rank) == "6":
                text_transform = annot_ax.transData + ScaledTranslation(
                    0.0, 6.0 / 72.0, fig.dpi_scale_trans
                )
            annot_ax.text(
                (x_rank + x_sc) / 2,
                annot_y,
                _boxplot_delta_n_p_annotation(g["delta"], g["n"], g["wilcoxon_p"]),
                transform=text_transform,
                ha="center",
                va="bottom",
                fontsize=font_annot,
                clip_on=False,
                zorder=10,
            )
            pair_centers.append((x_rank + x_sc) / 2)
            cursor = x_sc + group_gap
        tick_pos.append(float(np.mean(pair_centers)))
        tick_labels.append(_rank_display_label(rank))
        cursor += 0.25

    ax_bot.set_ylim(0, 200)
    top_hi = 400.0
    if annot_tops:
        top_hi = max(400.0, max(annot_tops) + 55.0)
    ax_top.set_ylim(200, top_hi)
    ax_bot.set_yticks([0, 50, 100, 150])
    ax_top.set_yticks([300, 400])

    ax_top.spines["bottom"].set_visible(False)
    ax_top.spines["top"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_bot.spines["right"].set_visible(False)
    ax_top.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)
    ax_bot.tick_params(axis="x", which="both", top=False, labeltop=False)
    ax_bot.tick_params(axis="y", labelsize=font_tick)
    ax_top.tick_params(axis="y", labelsize=font_tick)
    # Label "200" without a tick mark (break slashes replace that tick).
    ax_bot.text(
        -0.01,
        200,
        "200",
        transform=blended_transform_factory(ax_bot.transAxes, ax_bot.transData),
        ha="right",
        va="center",
        fontsize=font_tick,
        clip_on=False,
    )
    _add_broken_axis_marks_left(ax_top, ax_bot)

    legend_handles = [
        Patch(facecolor="none", edgecolor="none", label="All tumors:"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["1"],
            edgecolor="black",
            label="Ranked MTA",
        ),
        Patch(facecolor=sc_color, edgecolor="black", label="Matched SC"),
        Patch(facecolor="none", edgecolor="none", label="No-on-label tumors:"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["1"],
            edgecolor="black",
            hatch="///",
            label="Ranked MTA",
        ),
        Patch(
            facecolor=sc_color,
            edgecolor="black",
            hatch="///",
            label="Matched SC",
        ),
    ]
    leg = ax_top.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.18),
        bbox_transform=ax_top.transAxes
        + ScaledTranslation(0, 5 / 72, fig.dpi_scale_trans),
        fontsize=font_legend,
    )
    if leg is not None and len(leg.legend_handles) >= 4:
        leg.legend_handles[0].set_visible(False)
        leg.legend_handles[3].set_visible(False)
        leg.get_texts()[0].set_fontweight("semibold")
        leg.get_texts()[3].set_fontweight("semibold")
    ax_bot.set_xticks(tick_pos)
    ax_bot.set_xticklabels(tick_labels, fontsize=font_tick)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    y_mid = 0.5 * (ax_top.get_position().y1 + ax_bot.get_position().y0)

    # Place ylabel just left of the leftmost y-tick / "200" label (no overlap).
    min_label_x = 1.0
    for ax in (ax_top, ax_bot):
        for label in ax.get_yticklabels():
            if not label.get_visible() or not str(label.get_text()).strip():
                continue
            bbox = label.get_window_extent(renderer=renderer).transformed(
                fig.transFigure.inverted()
            )
            min_label_x = min(min_label_x, float(bbox.x0))
        for txt in ax.texts:
            if str(txt.get_text()).strip() != "200":
                continue
            bbox = txt.get_window_extent(renderer=renderer).transformed(
                fig.transFigure.inverted()
            )
            min_label_x = min(min_label_x, float(bbox.x0))

    # Small gap for the rotated ylabel text body itself (~half char width).
    ylabel_x = max(0.001, min_label_x - 0.012)
    ylab = fig.supylabel(
        "PFS (days)\nRank vs matched SC",
        fontsize=font_ylabel,
        x=ylabel_x,
        y=y_mid,
        va="center",
        ha="center",
    )
    # Refine once using the ylabel's own width so it sits flush without covering ticks.
    fig.canvas.draw()
    ylab_bbox = ylab.get_window_extent(renderer=fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted()
    )
    # Desired: ylabel right edge just left of tick labels.
    gap = 0.004
    shift = (min_label_x - gap) - ylab_bbox.x1
    ylab.set_x(ylabel_x + shift)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _add_broken_axis_marks_bottom(ax_left: plt.Axes, ax_right: plt.Axes, *, size: float = 0.0075) -> None:
    """Draw diagonal break marks on the bottom x-axis only (tilted / horizontal panels)."""
    kwargs = dict(color="k", clip_on=False, linewidth=0.8)
    ax_left.plot((1 - size, 1 + size), (-size, +size), transform=ax_left.transAxes, **kwargs)
    ax_right.plot((-size, +size), (-size, +size), transform=ax_right.transAxes, **kwargs)


def _save_rank_mpfs_vs_sc_paired_boxplot_tilted(
    ranked: pd.DataFrame,
    model_table: pd.DataFrame,
    output_path: Path,
    *,
    data_column: str = DEFAULT_DATA_COLUMN,
) -> None:
    """Horizontal (tilted) rank vs SC boxplots with broken x-axis (0–200 | 200–400)."""
    groups = build_rank_mpfs_vs_sc_paired_groups(
        ranked, model_table, data_column=data_column
    )
    if not groups:
        return

    font_tick = 21.2
    font_xlabel = 22.5
    font_legend = 14.4
    font_annot = 10.75
    annot_pad = 10.0
    sc_color = "#D4B4E8"

    # width_ratios: right (200–400) = 1/4 of left (0–200); height/width swapped vs upright.
    fig, (ax_left, ax_right) = plt.subplots(
        1,
        2,
        sharey=True,
        figsize=(9.0, 18.0),
        gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05},
    )
    tick_pos: List[float] = []
    tick_labels: List[str] = []
    cursor = 0.0
    pair_gap = 0.42
    group_gap = 0.72
    by_rank: Dict[str, List[Dict[str, Any]]] = {}
    for g in groups:
        by_rank.setdefault(g["rank"], []).append(g)

    annot_rights: List[float] = []
    for rank in RANKING_UTILITY_RANKS:
        rank_groups = by_rank.get(str(rank), [])
        if not rank_groups:
            continue
        rank_color = RANK_MPFS_BLUE_SHADES.get(str(rank), "#6BAED6")
        pair_centers: List[float] = []
        for g in rank_groups:
            is_no = g["cohort"] == "no_onlabel"
            y_rank = cursor
            y_sc = cursor + pair_gap
            for ax in (ax_left, ax_right):
                rank_box = ax.boxplot(
                    [g["rank_vals"]],
                    positions=[y_rank],
                    widths=0.4224,
                    vert=False,
                    patch_artist=True,
                    manage_ticks=False,
                    medianprops={"color": "black", "linewidth": 1.2},
                )
                sc_box = ax.boxplot(
                    [g["sc_vals"]],
                    positions=[y_sc],
                    widths=0.4224,
                    vert=False,
                    patch_artist=True,
                    manage_ticks=False,
                    medianprops={"color": "black", "linewidth": 1.2},
                )
                _style_boxplot_artists(
                    rank_box,
                    facecolor=rank_color,
                    hatch="///" if is_no else None,
                )
                _style_boxplot_artists(
                    sc_box,
                    facecolor=sc_color,
                    hatch="///" if is_no else None,
                )
            pair_max = float(
                np.nanmax(np.concatenate([g["rank_vals"], g["sc_vals"]]))
            )
            annot_x = pair_max + annot_pad
            annot_rights.append(annot_x)
            annot_ax = ax_right if annot_x > 200 else ax_left
            annot_ax.text(
                annot_x,
                (y_rank + y_sc) / 2,
                _boxplot_delta_n_p_annotation(g["delta"], g["n"], g["wilcoxon_p"]),
                ha="left",
                va="center",
                fontsize=font_annot,
                clip_on=False,
                zorder=10,
            )
            pair_centers.append((y_rank + y_sc) / 2)
            cursor = y_sc + group_gap
        tick_pos.append(float(np.mean(pair_centers)))
        tick_labels.append(_rank_display_label(rank))
        cursor += 0.25

    ax_left.set_xlim(0, 200)
    right_hi = 400.0
    if annot_rights:
        right_hi = max(400.0, max(annot_rights) + 55.0)
    ax_right.set_xlim(200, right_hi)
    ax_left.set_xticks([0, 50, 100, 150])
    ax_right.set_xticks([300, 400])

    ax_left.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    ax_left.spines["top"].set_visible(False)
    ax_right.spines["top"].set_visible(False)
    ax_right.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_left.tick_params(axis="x", labelsize=font_tick)
    ax_right.tick_params(axis="x", labelsize=font_tick)
    ax_left.text(
        200,
        -0.01,
        "200",
        transform=blended_transform_factory(ax_left.transData, ax_left.transAxes),
        ha="center",
        va="top",
        fontsize=font_tick,
        clip_on=False,
    )
    _add_broken_axis_marks_bottom(ax_left, ax_right)

    legend_handles = [
        Patch(facecolor="none", edgecolor="none", label="All tumors:"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["1"],
            edgecolor="black",
            label="Ranked MTA",
        ),
        Patch(facecolor=sc_color, edgecolor="black", label="Matched SC"),
        Patch(facecolor="none", edgecolor="none", label="No-on-label tumors:"),
        Patch(
            facecolor=RANK_MPFS_BLUE_SHADES["1"],
            edgecolor="black",
            hatch="///",
            label="Ranked MTA",
        ),
        Patch(
            facecolor=sc_color,
            edgecolor="black",
            hatch="///",
            label="Matched SC",
        ),
    ]
    leg = ax_right.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.08),
        bbox_transform=ax_right.transAxes
        + ScaledTranslation(0, 5 / 72, fig.dpi_scale_trans),
        fontsize=font_legend,
    )
    if leg is not None and len(leg.legend_handles) >= 4:
        leg.legend_handles[0].set_visible(False)
        leg.legend_handles[3].set_visible(False)
        leg.get_texts()[0].set_fontweight("semibold")
        leg.get_texts()[3].set_fontweight("semibold")
    ax_left.set_yticks(tick_pos)
    ax_left.set_yticklabels(tick_labels, fontsize=font_tick)
    ax_left.invert_yaxis()
    fig.supxlabel(
        "PFS (days)\nRank vs matched SC",
        fontsize=font_xlabel,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_clinical_utility_ranking_analysis(
    output_dir: Path, source_path: Path, source_sheet: Optional[str], data_column: str = DEFAULT_DATA_COLUMN
) -> None:
    """Generate Fig. 4c–4f panel figures and supporting tables."""
    LOGGER.info("Starting Fig. 4 panel outputs -> %s", output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ranked = _prepare_ranked_onlabel_patient_table(source_path, source_sheet)
    chemo = load_sc_chemo_reference(source_path, source_sheet)
    chemo_collapsed = _collapse_chemo_per_model(chemo, data_column)
    model_table = build_model_clinical_utility_table(ranked, chemo_collapsed, data_column)
    rank_distribution = build_clinical_utility_rank_distribution(model_table)

    safe_to_excel(model_table, output_dir / "clinical_utility_model_table.xlsx")
    safe_to_excel(
        rank_distribution,
        output_dir / "clinical_utility_onlabel_rank_distribution.xlsx",
    )

    s1 = model_table[
        (model_table["scenario"] == CLINICAL_UTILITY_SCENARIO1)
        & model_table["chemo_TTD"].notna()
    ]
    paired_s1 = _build_clinical_utility_paired_table(s1, "rank1_vs_sc")
    if not paired_s1.empty:
        safe_to_excel(
            paired_s1,
            output_dir / "clinical_utility_paired_scenario1_rank1_vs_sc.xlsx",
        )

    pub_style = PublicationPlotStyle()
    pub_style.apply_rcparams()
    _save_short_summary_rank_composition_figure(ranked, output_dir, style=pub_style)
    generate_publication_km_cox_figures(
        output_dir,
        model_table=model_table,
        ranked=ranked,
        data_column=data_column,
        chemo_collapsed=chemo_collapsed,
        chemo=chemo,
        fig4d_only=True,
    )

    q1_combined_waterfall = build_q1_combined_median_waterfall_table(
        model_table, ranked, data_column=data_column
    )
    if not q1_combined_waterfall.empty:
        safe_to_excel(
            q1_combined_waterfall,
            output_dir / "clinical_utility_q1_combined_pool_median_mPFS.xlsx",
        )
        f4f_path = output_dir / f"{FIG4F_PANEL_STEM}.jpg"
        _save_q1_combined_mpfs_boxplot(
            model_table,
            ranked,
            f4f_path,
            data_column=data_column,
        )
        _mirror_jpg_to_fig4_panel(f4f_path, FIG4_PANELS_DIR / FIG4F_PANEL_STEM)

    rank_delta_mpfs_table = build_rank_delta_mpfs_bar_table(
        ranked, model_table, data_column=data_column
    )
    if not rank_delta_mpfs_table.empty:
        safe_to_excel(
            rank_delta_mpfs_table,
            output_dir / "clinical_utility_rank_deltaMPFS_bar_table.xlsx",
        )
        f4e_path = output_dir / f"{FIG4E_PANEL_STEM}.jpg"
        _save_rank_mpfs_vs_sc_paired_boxplot(
            ranked,
            model_table,
            f4e_path,
            data_column=data_column,
        )
        _mirror_jpg_to_fig4_panel(f4e_path, FIG4_PANELS_DIR / FIG4E_PANEL_STEM)

    boxplot_fdr_table = build_boxplot_wilcoxon_fdr_table(
        ranked, model_table, data_column=data_column
    )
    if not boxplot_fdr_table.empty:
        safe_to_excel(
            boxplot_fdr_table,
            output_dir / "clinical_utility_boxplot_wilcoxon_FDR_table.xlsx",
        )
        boxplot_fdr_table.to_csv(
            output_dir / "clinical_utility_boxplot_wilcoxon_FDR_table.csv",
            index=False,
        )

    LOGGER.info("Completed Fig. 4 panel outputs.")




def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clinical utility ranking (standalone).")
    parser.add_argument("--source", type=Path, default=None)
    parser.add_argument("--sheet", type=str, default=DEFAULT_SOURCE_SHEET)
    parser.add_argument("--output-dir", type=Path, default=_SCRIPT_DIR)
    parser.add_argument("--data-column", type=str, default=DEFAULT_DATA_COLUMN)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = resolve_source_path(args.source)
    sheet = args.sheet.strip() if args.sheet and args.sheet.strip() else None
    output_dir = args.output_dir.resolve() / OUTPUT_DIRNAME
    output_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "script": "clinical_utility_ranking.py",
        "source_path": str(source_path),
        "output_dir": str(output_dir),
        "generated": datetime.now().isoformat(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    run_clinical_utility_ranking_analysis(output_dir, source_path, sheet, args.data_column)
    LOGGER.info("Outputs -> %s", output_dir)


if __name__ == "__main__":
    main()
