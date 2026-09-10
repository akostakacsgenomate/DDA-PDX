"""Kaplan-Meier analysis for PDX treatment cohorts.

Per-patient compound ranking (dense rank 1–6, pooled 7-10 vs. SC), plus a
Non_assigned vs 7-10 comparison. Pairwise ranking KMs use the same layout as
``all_Kmplot_1_vs_SC`` / ``all_Kmplot_7-10_vs_SC``.

The analysis reads the primary input workbook in data/
(Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx, sheet PCT_DRUG_RESPONSE).

Usage
-----
python src/KM_PDX.py [--source PATH] [--output-dir DIR]

Outputs are written to <output-dir>/kaplan_meier/ranking/.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import logging
import re
import shutil
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from plot_typography import (
    apply_arial_font,
    format_p_equals_mathtext,
    format_p_sci_plain,
    panel_font_sizes,
)

apply_arial_font()
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, multivariate_logrank_test, pairwise_logrank_test
from openpyxl import load_workbook
from openpyxl.styles import Alignment
from openpyxl.styles import PatternFill
from scipy.stats import chi2 as chi2_distribution
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

LOGGER = logging.getLogger("KM_PDX")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# ---------------------------------------------------------------------------
# Repository-relative paths (no hard-coded absolute paths)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

DEFAULT_DATA_COLUMN = "TimeToDouble"
KM_XLIM = (0, 200)
KM_XTICKS = np.arange(0, 225, 25)

DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_SOURCE_SHEET = "PCT_DRUG_RESPONSE"
KM_OUTPUT_DIRNAME = "kaplan_meier"


# =============================================================================
# Pairwise log-rank Excel export with FDR highlighting
# =============================================================================

_FILL_CHANGED = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")
_FILL_UNCHANGED = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
_FILL_SIG = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
_FILL_NS = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")
_FILL_CHANGED_TRUE = PatternFill(start_color="FFD966", end_color="FFD966", fill_type="solid")
_FILL_CHANGED_FALSE = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")

_SIG_COLS = ("raw_significance", "fdr_significance")
def _prepare_pairwise_export_df(df: pd.DataFrame) -> pd.DataFrame:
    return _reorder_columns(_add_fdr_significance_changed(df))


def _apply_pairwise_worksheet_colors(ws) -> None:
    """Apply significance colour-coding to one pairwise log-rank worksheet."""
    headers = [cell.value for cell in ws[1]]
    center = Alignment(horizontal="center", vertical="center")
    for row in ws.iter_rows():
        for cell in row:
            cell.alignment = center
    for col_name in _SIG_COLS:
        if col_name not in headers:
            continue
        col_idx = headers.index(col_name) + 1
        for row_idx in range(2, ws.max_row + 1):
            cell = ws.cell(row=row_idx, column=col_idx)
            cell.fill = _FILL_SIG if _sig_is_significant(cell.value) else _FILL_NS
    if "fdr_significance_changed" in headers:
        col_idx = headers.index("fdr_significance_changed") + 1
        for row_idx in range(2, ws.max_row + 1):
            cell = ws.cell(row=row_idx, column=col_idx)
            if cell.value == "changed":
                cell.fill = _FILL_CHANGED
            elif cell.value == "unchanged":
                cell.fill = _FILL_UNCHANGED


def _sig_is_significant(label: Any) -> bool:
    """Return True if the significance label represents any star level."""
    return str(label).strip().lower() not in ("", "ns")


def _format_p_value(p_value: Any) -> str:
    return format_p_sci_plain(p_value, digits=1)


def _reorder_p_adj_columns(df: pd.DataFrame) -> pd.DataFrame:
    _p_adj_cols = ("p_fdr_bh", "raw_significance", "fdr_significance")
    out = df.copy()
    _cols = list(out.columns)
    for _c in _p_adj_cols:
        if _c in _cols:
            _cols.remove(_c)
    if "p" in _cols:
        _insert_at = _cols.index("p") + 1
        for _c in _p_adj_cols:
            if _c in out.columns:
                _cols.insert(_insert_at, _c)
                _insert_at += 1
        out = out[_cols]
    return out


def _add_fdr_significance_changed(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate rows where FDR correction changes significance at alpha=0.05."""
    out = df.copy()
    if "raw_significance" not in out.columns or "fdr_significance" not in out.columns:
        return out
    out["fdr_significance_changed"] = [
        "changed" if _sig_is_significant(r) != _sig_is_significant(f) else "unchanged"
        for r, f in zip(out["raw_significance"], out["fdr_significance"])
    ]
    return out


def _reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns)
    if "fdr_significance_changed" in cols:
        cols.remove("fdr_significance_changed")
    if "fdr_significance" in cols:
        cols.insert(cols.index("fdr_significance") + 1, "fdr_significance_changed")
    return df[cols]


def export_pairwise_logrank_excel(df: pd.DataFrame, path: str) -> None:
    """Write pairwise log-rank table with colour-coded significance columns."""
    out = _prepare_pairwise_export_df(df)
    out.to_excel(path, index=False, engine="openpyxl")
    wb = load_workbook(path)
    _apply_pairwise_worksheet_colors(wb.active)
    wb.save(path)


def apply_excel_alignment_and_significance(path: Path) -> None:
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
    wb.save(path)


# =============================================================================
# Embedded helper: dense per-patient LEVEL ranking
# (originally in codes/ranking_dense_level.py)
# =============================================================================

_RANK_POOL_FROM = 7
_POOLED_LABEL = "7-10"

# Ordinal legend / at-risk order for the primary ranking KM (Fig. 3f).
RANKING_GROUP_ORDER: Tuple[str, ...] = (
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    _POOLED_LABEL,
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
# Display labels for ranking-vs-SC pair KMs (internal GROUP keys unchanged).
RANKING_VS_SC_STAT_LABELS: Dict[str, str] = {
    "1": "Top",
    "7-10": "Bottom",
    "SC": "SC",
}
NONASSIGNED_VS_710_STAT_LABELS: Dict[str, str] = {
    "Non_assigned": "Non assigned",
    "7-10": "Bottom",
}
# Fixed output canvas matching the Top pair-KM etalon (dpi=300, tight bbox).
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


# =============================================================================
# Configuration containers
# =============================================================================


@dataclass(frozen=True)
class Paths:
    """Filesystem paths resolved relative to the repository data/ directory."""

    project_root: Path
    script_dir: Path
    source_primary_xlsx: Path
    source_primary_sheet: Optional[str]
    non_assigned_xlsx: Path


@dataclass(frozen=True)
class KMAnalysisSpec:
    """Metadata for one Kaplan-Meier analysis branch."""

    run_key: str
    source_script_stem: str
    legacy_name: str
    date: str
    data_column: str = DEFAULT_DATA_COLUMN
    primary_tumor: str = ""
    extra: Dict[str, Any] = field(default_factory=dict)


# Multi-group KMs are drawn at MULTIGROUP_KM_FIGSIZE for a half-width cell
# of the 4x6 Fig. 3 composite grid.
_MULTIGROUP_KM_SCALE = 8.0 / 5.0
_MULTIGROUP_KM_BASE = (5.0 * _MULTIGROUP_KM_SCALE, 5.5 * _MULTIGROUP_KM_SCALE)
MULTIGROUP_KM_FIGSIZE = (_MULTIGROUP_KM_BASE[0] * 1.20, _MULTIGROUP_KM_BASE[1])
_KM_FONTS = panel_font_sizes(*MULTIGROUP_KM_FIGSIZE, nrows=4, ncols=6, colspan=3)
# Clears the enlarged tick labels and x-axis title above the at-risk table.
MULTIGROUP_KM_AT_RISK_YPOS = -0.52


@dataclass
class KMPlotStyle:
    """Shared typography parameters (sized for the A4 composite panels)."""

    font: float = _KM_FONTS["tick"]
    tick_size: float = _KM_FONTS["tick"]
    label_size: float = _KM_FONTS["label"]
    at_risk: float = _KM_FONTS["annotation"]
    font_legend: float = _KM_FONTS["legend"]
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


# Top / Bottom / Non_assigned pair KMs are drawn at a dedicated figure size
# for the Fig. 4 composite grid.
_PAIR_KM_FONT_SCALE = 1.20
PAIR_KM_ETALON_STYLE = KMPlotStyle(
    font=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    tick_size=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    label_size=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
    at_risk=5.5 * 1.6 * _PAIR_KM_FONT_SCALE,
    font_legend=6.4 * 1.6 * _PAIR_KM_FONT_SCALE,
)


# =============================================================================
# CLI and path resolution
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description="Kaplan-Meier analysis for PDX treatment cohorts (ranking)."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help=(
            f"Primary input workbook override "
            f"(default: data/{DEFAULT_SOURCE_XLSX.name})."
        ),
    )
    parser.add_argument(
        "--sheet",
        type=str,
        default=DEFAULT_SOURCE_SHEET,
        help=(
            f"Excel worksheet for the primary source "
            f"(default: {DEFAULT_SOURCE_SHEET!r}; use empty string for first sheet)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help=f"Parent output directory; results in <output-dir>/{KM_OUTPUT_DIRNAME}/.",
    )
    parser.add_argument(
        "--analyses",
        nargs="+",
        choices=("ranking", "all"),
        default=["ranking"],
        help="Analyses to execute (default: ranking).",
    )
    return parser.parse_args()


def filesystem_safe_name(text: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_")
    safe = "".join(ch if ch in allowed else "_" for ch in text.strip().replace(" ", "_"))
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_meta(path: Path) -> Dict[str, str]:
    meta: Dict[str, str] = {"path": str(path), "exists": str(path.exists())}
    if path.exists():
        meta["sha256"] = sha256_of_file(path)
        meta["size_bytes"] = str(path.stat().st_size)
    return meta


def analysis_subfolder_name(spec: KMAnalysisSpec) -> str:
    fixed_names = {
        "ranking": "ranking",
    }
    return fixed_names.get(spec.run_key, filesystem_safe_name(spec.run_key))


def resolve_analysis_output_dir(output_root: Path, spec: KMAnalysisSpec) -> Path:
    nested = output_root / analysis_subfolder_name(spec)
    nested.mkdir(parents=True, exist_ok=True)
    return nested


# =============================================================================
# Shared data preparation
# =============================================================================


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


def resolve_on_label_columns(patient_df: pd.DataFrame) -> Tuple[str, str, str]:
    """Detect on-label column naming convention used in the workbook."""
    on_col = "On_label" if "On_label" in patient_df.columns else "On_off_label"
    on_yes = "IGAZ" if on_col == "On_label" else "YES_on"
    on_no = "HAMIS" if on_col == "On_label" else "no_off"
    return on_col, on_yes, on_no


def is_chemo_on_label(series: pd.Series) -> pd.Series:
    """True for rows tagged On_label = chemo (standard-of-care reference arm)."""
    return series.astype(str).str.strip().str.lower() == "chemo"


def load_sc_chemo_km_group(
    source_path: Path, source_sheet: Optional[str]
) -> pd.DataFrame:
    """SC reference rows labelled GROUP=SC for Kaplan-Meier analyses."""
    chemo = load_sc_chemo_reference(source_path, source_sheet)
    if chemo.empty:
        return chemo
    chemo = chemo.copy()
    chemo["GROUP"] = "SC"
    return chemo


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


def load_non_assigned_table(paths: Paths) -> pd.DataFrame:
    """Load Non_assigned rows; return empty frame with expected columns if file absent."""
    if not paths.non_assigned_xlsx.exists():
        LOGGER.warning(
            "Non_assigned workbook not found at %s; SHIVA subset plot will omit Non_assigned.",
            paths.non_assigned_xlsx,
        )
        return pd.DataFrame(columns=["Model", "COMPOUND", "TimeToDouble", "CENSOR", "GROUP", "LEVEL"])
    non_assigned = pd.read_excel(paths.non_assigned_xlsx)
    non_assigned["CENSOR"] = True
    non_assigned["GROUP"] = "Non_assigned"
    return non_assigned


def validate_configured_sources(paths: Paths) -> Dict[str, Any]:
    primary_meta = file_meta(paths.source_primary_xlsx)
    return {
        "shared_primary_source": str(paths.source_primary_xlsx),
        "shared_primary_source_sheet": paths.source_primary_sheet,
        "shared_primary_source_meta": primary_meta,
        "all_analyses_use_same_source": True,
    }


# =============================================================================
# Statistics and plotting utilities
# =============================================================================


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
    label_extra_spaces: int = 2,
) -> None:
    """Add lifelines at-risk count table below a KM axes.

    Cohort/group names are given trailing spaces before lifelines builds the
    condensed tick labels, which moves the name column left of the counts
    (clears the ``7-10`` / first-count collision).
    """
    if not kmf_list:
        return
    # Non-breaking spaces survive matplotlib tick-label whitespace stripping and
    # push cohort names left of the first at-risk count (~5+ pt clearance).
    pad = "\u00A0" * max(int(label_extra_spaces), 0)
    labels = [f"{kmf._label}{pad}" for kmf in kmf_list]
    add_at_risk_counts(
        *kmf_list,
        ypos=at_risk_ypos,
        ax=ax,
        labels=labels,
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


def fit_multigroup_km_plot(
    patient_data: pd.DataFrame,
    data_column: str,
    colors: Dict[str, str],
    style: KMPlotStyle,
    *,
    legend_labels: Optional[Sequence[str]] = None,
    legend_title: Optional[str] = None,
    line_styles: Optional[Dict[str, str]] = None,
    line_widths: Optional[Dict[str, float]] = None,
    group_order: Optional[Sequence[str]] = None,
    label_suffix_n: bool = False,
    at_risk_ypos: float = -0.35,
    title: str = "",
) -> Tuple[plt.Figure, plt.Axes, List[KaplanMeierFitter], pd.DataFrame]:
    """Fit and plot Kaplan-Meier curves for all requested groups."""
    style.apply_rcparams()
    fig, ax = plt.subplots(figsize=MULTIGROUP_KM_FIGSIZE)
    kmf_list: List[KaplanMeierFitter] = []
    rows: List[Dict[str, Any]] = []

    groups = list(group_order) if group_order else list(patient_data["GROUP"].unique())
    default_linewidth = 1.0
    line_styles = line_styles or {}
    line_widths = line_widths or {}

    for group in groups:
        group_data = patient_data[patient_data["GROUP"] == group].copy()
        if group_data.empty:
            LOGGER.warning("Skipping empty group: %s", group)
            continue

        group_data.loc[:, data_column] = pd.to_numeric(group_data[data_column], errors="coerce")
        group_data = group_data.dropna(subset=[data_column])
        if group_data.empty:
            continue

        label = f"{group} (n={len(group_data)})" if label_suffix_n else group
        kmf = KaplanMeierFitter()
        kmf.fit(group_data[data_column], event_observed=group_data["CENSOR"], label=label)
        kmf.plot(
            ax=ax,
            color=colors.get(group),
            ci_show=style.cishow,
            linestyle=line_styles.get(group, "-"),
            linewidth=line_widths.get(group, default_linewidth),
        )
        kmf_list.append(kmf)
        stat = collect_km_group_statistics(kmf, label, group, group_data)
        stat["Group"] = group
        rows.append(stat)

    apply_km_axis_style(ax, style, title=title)
    if legend_labels is not None:
        ax.legend(labels=list(legend_labels), frameon=False, fontsize=style.font_legend, title=legend_title)
    else:
        ax.legend(frameon=False, fontsize=style.font_legend, title=legend_title)

    if kmf_list:
        add_km_at_risk_table(ax, kmf_list, at_risk_ypos=at_risk_ypos)

    plt.subplots_adjust(bottom=0.25, top=0.9, hspace=0.5)
    plt.tight_layout()
    summary_df = pd.DataFrame(rows)
    return fig, ax, kmf_list, summary_df


def append_fdr_to_pairwise_hr(pairwise_hr_df: pd.DataFrame) -> pd.DataFrame:
    """Add BH-FDR adjusted p-values and star labels to a pairwise HR export table."""
    if pairwise_hr_df.empty or not pairwise_hr_df["logrank_p"].notna().any():
        return pairwise_hr_df

    out = pairwise_hr_df.copy()
    raw = out["logrank_p"].astype(float).values
    _, fdr_lr, _, _ = multipletests(raw, method="fdr_bh")
    out["logrank_p_fdr_bh"] = fdr_lr
    out["raw_significance"] = out["logrank_p"].apply(significance_stars)
    out["fdr_significance"] = out["logrank_p_fdr_bh"].apply(significance_stars)

    adj_cols = ("logrank_p_fdr_bh", "raw_significance", "fdr_significance")
    cols = [c for c in out.columns if c not in adj_cols]
    if "logrank_p" in cols:
        insert_at = cols.index("logrank_p") + 1
        for c in adj_cols:
            if c in out.columns:
                cols.insert(insert_at, c)
                insert_at += 1
        out = out[cols]
    return out


def format_pairwise_logrank_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Apply BH-FDR, star labels, and p-value formatting for log-rank export."""
    out = df.reset_index().rename(columns={"index": "GROUP"}) if "GROUP" not in df.columns else df.copy()
    out["p"] = out["p"].astype(float)
    raw = out["p"].values
    _, fdr_p, _, _ = multipletests(raw, method="fdr_bh")
    out["p_fdr_bh"] = fdr_p
    out["raw_significance"] = out["p"].apply(significance_stars)
    out["fdr_significance"] = out["p_fdr_bh"].apply(significance_stars)
    out["p_fdr_bh"] = out["p_fdr_bh"].apply(lambda x: format_p_sci_plain(x, digits=1))
    out["p"] = out["p"].apply(lambda x: format_p_sci_plain(x, digits=1))
    out = out.drop(columns=["significance", "test_statistic", "-log2(p)"], errors="ignore")

    adj_cols = ("p_fdr_bh", "raw_significance", "fdr_significance")
    cols = [c for c in out.columns if c not in adj_cols]
    if "p" in cols:
        insert_at = cols.index("p") + 1
        for c in adj_cols:
            if c in out.columns:
                cols.insert(insert_at, c)
                insert_at += 1
        out = out[cols]
    return out


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
        0.97,
        0.40,
        space_plot_operators("\n".join(stats_lines)),
        transform=ax.transAxes,
        fontsize=style.font * 0.9,
        va="top",
        ha="right",
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


def build_global_logrank_summary_row(
    patient_data: pd.DataFrame,
    data_column: str,
    *,
    robust: bool = False,
) -> Optional[pd.DataFrame]:
    """Return one-row omnibus multivariate log-rank summary (level_0/1 = 'global')."""
    data = patient_data.dropna(subset=[data_column]).copy()
    if data["GROUP"].nunique() < 2:
        return None

    kwargs = {"specify_robust": True} if robust else {}
    result = multivariate_logrank_test(
        data[data_column], data["GROUP"], data["CENSOR"], **kwargs
    )
    p = float(result.p_value)
    p_str = format_p_sci_plain(p, digits=1)
    stars = significance_stars(p)
    return pd.DataFrame(
        [
            {
                "level_0": "global",
                "level_1": "global",
                "p": p_str,
                "p_fdr_bh": p_str,
                "raw_significance": stars,
                "fdr_significance": stars,
            }
        ]
    )


def export_global_pairwise_logrank(
    patient_data: pd.DataFrame,
    data_column: str,
    output_path: Path,
    *,
    robust: bool = False,
    include_global: bool = False,
) -> pd.DataFrame:
    """Run all-pairs log-rank test and write colour-coded Excel export."""
    kwargs = {"specify_robust": True} if robust else {}
    results = pairwise_logrank_test(
        patient_data[data_column], patient_data["GROUP"], patient_data["CENSOR"], **kwargs
    )
    summary = format_pairwise_logrank_summary(results.summary)
    if include_global:
        global_row = build_global_logrank_summary_row(
            patient_data, data_column, robust=robust
        )
        if global_row is not None:
            summary = pd.concat([summary, global_row], ignore_index=True)
    export_pairwise_logrank_excel(summary, str(output_path))
    return summary



# =============================================================================
# Analysis 2: Per-patient ranking KM
# =============================================================================


def build_ranking_dataset(
    paths: Paths, source_path: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Assign dense rank groups per patient and append SC reference rows.

    Returns (data_without_sc, data_with_sc) for main and pairwise plots.
    """
    patient = load_primary_patient_table(source_path, paths.source_primary_sheet)
    patient = harmonize_column_names(patient)
    patient = normalize_id_columns(patient)
    patient = clean_base_patient_table(patient)
    patient["CENSOR"] = True

    chemo = load_sc_chemo_km_group(source_path, paths.source_primary_sheet)
    ranked = assign_dense_level_groups(patient, model_col="Model", level_col="LEVEL", group_col="GROUP")
    combined = pd.concat([ranked, chemo], ignore_index=True)

    with_sc = normalize_survival_columns(combined.copy(), DEFAULT_DATA_COLUMN)
    with_sc = with_sc.dropna(subset=[DEFAULT_DATA_COLUMN])

    without_sc = with_sc[with_sc["GROUP"] != "SC"].copy()
    return without_sc, with_sc


def run_ranking_analysis(
    paths: Paths, spec: KMAnalysisSpec, output_dir: Path, source_path: Path
) -> None:
    LOGGER.info("Starting ranking KM analysis -> %s", output_dir)
    patient_data, patient_data_with_sc = build_ranking_dataset(paths, source_path)
    style = KMPlotStyle()
    pair_style = PAIR_KM_ETALON_STYLE
    colors = dict(RANKING_KM_COLORS)

    fig, _, _, summary_df = fit_multigroup_km_plot(
        patient_data,
        spec.data_column,
        colors,
        style,
        group_order=list(RANKING_GROUP_ORDER),
        legend_title="Ranking",
        at_risk_ypos=MULTIGROUP_KM_AT_RISK_YPOS,
    )
    fig.subplots_adjust(left=0.14)
    primary = safe_filename(spec.primary_tumor) or "all"
    fig.savefig(
        output_dir / f"{primary}_Kmplot.jpg",
        dpi=300,
    )
    plt.close(fig)
    safe_to_excel(
        summary_df,
        output_dir / f"{primary}_median_survival_times.xlsx",
    )

    pairwise_rows: List[Dict[str, Any]] = []
    etalon_canvas: Optional[Tuple[int, int]] = None
    for rank_group in spec.extra.get("ranking_groups_vs_sc", ()):
        out_jpg = (
            output_dir
            / f"{primary}_Kmplot_{rank_group}_vs_SC.jpg"
        )
        row = create_pair_km_plot(
            patient_data_with_sc,
            rank_group,
            "SC",
            spec.data_column,
            colors,
            pair_style,
            out_jpg,
            color_overrides={(rank_group, "SC"): "#404040"} if rank_group == "7-10" else None,
            # Rank 1 → "Top", pooled 7-10 → "Bottom" on legend / mPFS / at-risk.
            stat_label_map=RANKING_VS_SC_STAT_LABELS,
        )
        if etalon_canvas is None and out_jpg.is_file():
            try:
                from PIL import Image

                etalon_canvas = Image.open(out_jpg).size
            except Exception:
                etalon_canvas = PAIR_KM_ETALON_PX
        _pad_pair_km_to_etalon_canvas(out_jpg, size=etalon_canvas or PAIR_KM_ETALON_PX)
        if row:
            pairwise_rows.append(row)

    # Non_assigned vs pooled rank 7-10 (same folder / canvas as Top & Bottom).
    non_assigned = load_non_assigned_table(paths)
    if not non_assigned.empty:
        non_assigned = harmonize_column_names(non_assigned)
        non_assigned = normalize_id_columns(non_assigned)
        non_assigned = normalize_survival_columns(non_assigned, spec.data_column)
        non_assigned = non_assigned.dropna(subset=[spec.data_column])
        rank710 = patient_data[patient_data["GROUP"].astype(str) == _POOLED_LABEL].copy()
        if not rank710.empty and not non_assigned.empty:
            pair_df = pd.concat(
                [
                    rank710[[spec.data_column, "GROUP", "CENSOR"]],
                    non_assigned[[spec.data_column, "GROUP", "CENSOR"]],
                ],
                ignore_index=True,
            )
            non_jpg = (
                output_dir
                / f"{primary}_Kmplot_Non_assigned_vs_7-10.jpg"
            )
            row = create_pair_km_plot(
                pair_df,
                "Non_assigned",
                _POOLED_LABEL,
                spec.data_column,
                colors,
                pair_style,
                non_jpg,
                color_overrides={("Non_assigned", _POOLED_LABEL): colors["Non_assigned"]},
                stat_label_map=NONASSIGNED_VS_710_STAT_LABELS,
            )
            _pad_pair_km_to_etalon_canvas(non_jpg, size=etalon_canvas or PAIR_KM_ETALON_PX)
            if row:
                pairwise_rows.append(row)
    else:
        LOGGER.warning(
            "Skipping Non_assigned vs 7-10 KM: non_assigned workbook not found."
        )

    pairwise_hr_df = append_fdr_to_pairwise_hr(pd.DataFrame(pairwise_rows))
    export_pairwise_logrank_excel(
        pairwise_hr_df,
        str(output_dir / f"{primary}_pairwise_HR_values.xlsx"),
    )
    export_global_pairwise_logrank(
        patient_data,
        spec.data_column,
        output_dir / f"{primary}_pairwise_logrank.xlsx",
    )
    LOGGER.info("Completed ranking KM analysis.")


# =============================================================================
def _default_run_date_prefix() -> str:
    return datetime.now().strftime("%Y_%m_%d_")


def build_analysis_specs(run_date: Optional[str] = None) -> List[KMAnalysisSpec]:
    stamp = run_date if run_date is not None else _default_run_date_prefix()
    return [
        KMAnalysisSpec(
            run_key="ranking",
            source_script_stem="KM_ranking",
            legacy_name="ranking",
            date=stamp,
            extra={"ranking_groups_vs_sc": ("1", "2", "3", "4", "5", "6", "7-10")},
        ),
    ]


def resolve_source_path(explicit: Optional[Path] = None) -> Path:
    """Resolve the primary drug-response workbook."""
    if explicit is not None:
        candidate = explicit.resolve()
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"Source file not found: {candidate}")

    if DEFAULT_SOURCE_XLSX.is_file():
        return DEFAULT_SOURCE_XLSX.resolve()

    raise FileNotFoundError(f"Source file not found: {DEFAULT_SOURCE_XLSX}")


def resolve_paths(args: argparse.Namespace) -> Paths:
    source_primary = resolve_source_path(args.source)
    sheet = args.sheet.strip() if args.sheet and args.sheet.strip() else None
    return Paths(
        project_root=REPO_ROOT,
        script_dir=_SCRIPT_DIR,
        source_primary_xlsx=source_primary,
        source_primary_sheet=sheet,
        non_assigned_xlsx=DATA_DIR / "non_assigned_merged.xlsx",
    )


def source_for_analysis(paths: Paths, _spec: KMAnalysisSpec, _args: argparse.Namespace) -> Path:
    """All analyses use the same primary workbook."""
    return paths.source_primary_xlsx


def write_analysis_info(output_dir: Path, spec: KMAnalysisSpec, source_path: Path, source_sheet: Optional[str]) -> None:
    info = {
        "run_key": spec.run_key,
        "source_script_stem": spec.source_script_stem,
        "legacy_name": spec.legacy_name,
        "date": spec.date,
        "data_column": spec.data_column,
        "output_folder": analysis_subfolder_name(spec),
        "source_path": str(source_path),
        "source_sheet": source_sheet,
        "source_meta": file_meta(source_path),
    }
    (output_dir.parent / "analysis_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")


ANALYSIS_ALIASES: Dict[str, Tuple[str, ...]] = {}


def analysis_is_requested(run_key: str, requested: Sequence[str]) -> bool:
    """True when run_key or a registered alias appears in requested analyses."""
    if "all" in requested:
        return True
    if run_key in requested:
        return True
    return any(alias in requested for alias in ANALYSIS_ALIASES.get(run_key, ()))


def main() -> None:
    """Run requested KM analyses and write outputs with reproducibility metadata."""
    args = parse_args()
    paths = resolve_paths(args)
    output_root = args.output_dir.resolve() / KM_OUTPUT_DIRNAME
    output_root.mkdir(parents=True, exist_ok=True)

    if not paths.source_primary_xlsx.exists():
        raise FileNotFoundError(
            f"Primary source not found: {paths.source_primary_xlsx}\n"
            f"Expected: data/{DEFAULT_SOURCE_XLSX.name}"
        )

    source_meta = validate_configured_sources(paths)
    LOGGER.info(
        "All analyses using primary source: %s (sheet=%s)",
        paths.source_primary_xlsx,
        paths.source_primary_sheet or "<first sheet>",
    )
    selected = "all" if "all" in args.analyses else ",".join(args.analyses)

    runners: Dict[str, Callable[..., None]] = {
        "ranking": run_ranking_analysis,
    }

    run_date = _default_run_date_prefix()
    analysis_folders: Dict[str, str] = {}
    for spec in build_analysis_specs(run_date):
        if spec.run_key not in runners:
            continue
        if not analysis_is_requested(spec.run_key, args.analyses):
            LOGGER.info("Skipping analysis '%s' (not in --analyses).", spec.run_key)
            continue

        source_path = source_for_analysis(paths, spec, args)
        out_dir = resolve_analysis_output_dir(output_root, spec)
        analysis_folders[spec.run_key] = str(out_dir)
        write_analysis_info(out_dir, spec, source_path, paths.source_primary_sheet)
        runners[spec.run_key](paths, spec, out_dir, source_path)

    metadata = {
        "seed": SEED,
        "generated_at": datetime.now().isoformat(),
        "output_root": str(output_root),
        "analyses_requested": selected,
        "analysis_output_folders": analysis_folders,
        **source_meta,
    }
    (output_root / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    LOGGER.info("Kaplan-Meier analysis completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
