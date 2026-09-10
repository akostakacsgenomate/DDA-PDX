"""Hazard-ratio analysis for PDX treatment cohorts.

Computes Cox proportional-hazards forest plots and pairwise HR tables for
per-patient compound rank (dense rank 1–6, pooled 7-10) — Fig. 3g.

Cox models use an unpenalized CoxPHFitter (no L2 penalizer).

All analyses read the supplementary PCT_DRUG_RESPONSE workbook (data/).
Rows with On_label = chemo are excluded from patient cohorts and retained only
as the SC reference arm where applicable.

Usage
-----
python src/HR_PDX.py [--source PATH] [--sheet SHEET] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from openpyxl import load_workbook
from openpyxl.styles import Alignment, PatternFill
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from plot_typography import (
    apply_arial_font,
    format_p_sci_mathtext,
    panel_font_sizes,
)

apply_arial_font()

# The ranking forest plot (Fig. 3g) is sized for a half-width cell of the 4x6
# Fig. 3 composite grid.
_RANKING_FOREST_SCALE = 1.6
RANKING_FOREST_WIDTH_IN = 8.0 * _RANKING_FOREST_SCALE
RANKING_FOREST_ROW_IN = 0.52 * _RANKING_FOREST_SCALE
_RANKING_FOREST_MARK_SCALE = 1.4
_RANKING_FOREST_HEIGHT_SCALE = 1.10
_RANKING_FOREST_CI_LINEWIDTH = 2.0 * _RANKING_FOREST_MARK_SCALE
_RANKING_FOREST_POINT_SIZE = 35.0 * _RANKING_FOREST_MARK_SCALE
_RANKING_FOREST_HEIGHT_IN = max(
    2.2 * _RANKING_FOREST_SCALE,
    RANKING_FOREST_ROW_IN * 6 + 1.2 * _RANKING_FOREST_SCALE,
) * _RANKING_FOREST_HEIGHT_SCALE
RANKING_FOREST_FONTS = panel_font_sizes(
    RANKING_FOREST_WIDTH_IN, _RANKING_FOREST_HEIGHT_IN, nrows=4, ncols=6, colspan=3
)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

LOGGER = logging.getLogger("HR_PDX")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

# ---------------------------------------------------------------------------
# Repository-relative paths (no hard-coded absolute paths)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_SOURCE_SHEET = "PCT_DRUG_RESPONSE"
HR_OUTPUT_DIRNAME = "hazard_ratios"
TREAT_TYPE = "single"
_RANK_POOL_FROM = 7
_POOLED_LABEL = "7-10"
RANKING_PAIRWISE_ORDER: List[str] = [
    "1 vs. 2",
    "1 vs. 3",
    "1 vs. 4",
    "1 vs. 5",
    "1 vs. 6",
    "1 vs. 7-10",
]


@dataclass
class Paths:
    project_root: Path
    source_primary: Path
    source_sheet: str


@dataclass
class AnalysisSpec:
    run_label: str
    order_list: List[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cox proportional-hazards analysis for PDX rank cohorts."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help=(
            "Input workbook (supplementary drug-response table). "
            f"Default: data/{DEFAULT_SOURCE_XLSX.name}"
        ),
    )
    parser.add_argument(
        "--sheet",
        type=str,
        default=DEFAULT_SOURCE_SHEET,
        help=f"Excel sheet name in --source (default: {DEFAULT_SOURCE_SHEET}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help=f"Parent output directory; results in <output-dir>/{HR_OUTPUT_DIRNAME}/.",
    )
    return parser.parse_args()


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_source_path(explicit: Path | None) -> Path:
    """Resolve the supplementary drug-response workbook path."""
    if explicit is not None:
        candidate = explicit.resolve()
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"Source file not found: {candidate}")

    if DEFAULT_SOURCE_XLSX.is_file():
        return DEFAULT_SOURCE_XLSX.resolve()

    raise FileNotFoundError(f"Source file not found: {DEFAULT_SOURCE_XLSX}")


def is_chemo_on_label(series: pd.Series) -> pd.Series:
    """True for rows tagged On_label = chemo (standard-of-care reference arm)."""
    return series.astype(str).str.strip().str.lower() == "chemo"


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


def harmonize_source_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename supplementary columns to the schema expected downstream."""
    out = df.copy()
    renames = {
        "PCT_ID": "Model",
        "DDA score": "LEVEL",
        "TimeToDouble (days)": "TimeToDouble",
    }
    rename = {source: target for source, target in renames.items() if source in out.columns and target not in out.columns}
    if "Treatment" in out.columns and "COMPOUND" not in out.columns:
        out = out.rename(columns={"Treatment": "COMPOUND"})
    if rename:
        out = out.rename(columns=rename)
    if "Approved" in out.columns and "APPROVED" not in out.columns:
        out["APPROVED"] = out["Approved"].map(normalize_approved_workbook_value)
    elif "APPROVED" in out.columns:
        out["APPROVED"] = out["APPROVED"].map(normalize_approved_workbook_value)
    if "On_label" in out.columns:
        out["On_label"] = out["On_label"].map(normalize_on_label_workbook_value)
    return out


def read_source_workbook(paths: Paths) -> pd.DataFrame:
    """Load the shared supplementary workbook (single-agent rows, harmonised)."""
    df = pd.read_excel(paths.source_primary, sheet_name=paths.source_sheet)
    df = harmonize_source_columns(df)

    required = {"Model", "COMPOUND", "LEVEL", "TimeToDouble", "APPROVED", "On_label"}
    missing = required.difference(df.columns)
    if missing:
        raise KeyError(
            f"Source table is missing required columns: {sorted(missing)} "
            f"(sheet '{paths.source_sheet}' of {paths.source_primary.name})."
        )

    if "Treatment type" in df.columns:
        df = df[df["Treatment type"] == TREAT_TYPE].copy()

    df = normalize_id_columns(df, model_col="Model", compound_col="COMPOUND")
    return clean_base_patient_table(df)


def load_patient_table(paths: Paths) -> pd.DataFrame:
    """Load patient cohort excluding On_label = chemo rows."""
    df = read_source_workbook(paths)
    chemo_mask = is_chemo_on_label(df["On_label"])
    if chemo_mask.any():
        LOGGER.info(
            "Excluded %d row(s) with On_label = chemo from patient cohort.",
            int(chemo_mask.sum()),
        )
    return df[~chemo_mask].copy()


def normalize_id_columns(df: pd.DataFrame, model_col: str, compound_col: str) -> pd.DataFrame:
    out = df.copy()
    out[model_col] = out[model_col].astype(str).str.replace("-", "", regex=False).str.lower()
    out[compound_col] = (
        out[compound_col].astype(str).str.replace("-", "", regex=False).str.lower()
    )
    return out


def clean_base_patient_table(patient_df: pd.DataFrame) -> pd.DataFrame:
    out = patient_df.copy()
    out = out[out["LEVEL"] != "#HIÁNYZIK"].copy()
    out["LEVEL"] = out["LEVEL"].astype(float)
    out["TimeToDouble"] = out["TimeToDouble"].astype(float)
    if "Day_Last" in out.columns:
        out["Day_Last"] = pd.to_numeric(out["Day_Last"], errors="coerce")
    out["CENSOR"] = True
    return out


def prepare_rank_groups(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    dense_rank = out.groupby("Model")["LEVEL"].rank(method="dense", ascending=False)
    dense_rank = dense_rank.fillna(0).astype(int)
    # Dense rank: 1 = highest LEVEL within patient; ranks >=7 pooled as "7-10"
    out["GROUP"] = dense_rank.map(lambda x: str(x) if 1 <= x <= 6 else "7-10")
    return out


def compute_pairwise_hr(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, j in combinations(summary_df.index, 2):
        g1, g2 = summary_df.loc[i, "covariate"], summary_df.loc[j, "covariate"]
        hr1, hr2 = summary_df.loc[i, "exp(coef)"], summary_df.loc[j, "exp(coef)"]
        se1, se2 = summary_df.loc[i, "se(coef)"], summary_df.loc[j, "se(coef)"]
        se_ratio = np.sqrt(se1**2 + se2**2)

        for left, right, ratio in ((g1, g2, hr1 / hr2), (g2, g1, hr2 / hr1)):
            # se(coef) is on the natural-log scale; CI must use ln(HR), not log2.
            log_ratio = np.log(ratio)
            low = np.exp(log_ratio - 1.96 * se_ratio)
            high = np.exp(log_ratio + 1.96 * se_ratio)
            z = log_ratio / se_ratio if se_ratio > 0 else np.nan
            p = 2 * (1 - norm.cdf(abs(z))) if np.isfinite(z) else np.nan
            rows.append(
                {
                    "Group Comparison": f"{left} vs. {right}",
                    "HR Ratio": ratio,
                    "95% CI Lower": low,
                    "95% CI Upper": high,
                    "P_Value": p,
                }
            )

    pairwise = pd.DataFrame(rows)
    first_orientation = pairwise.iloc[::2].copy()
    pvals = first_orientation["P_Value"].to_numpy(dtype=float)
    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    first_orientation["P_Value_FDR"] = fdr

    fdr_lookup = dict(zip(first_orientation["Group Comparison"], first_orientation["P_Value_FDR"]))

    def lookup_adjusted(comp: str, table: Dict[str, float]) -> float:
        if comp in table:
            return table[comp]
        if " vs. " not in comp:
            return np.nan
        left, right = comp.split(" vs. ", 1)
        return table.get(f"{right} vs. {left}", np.nan)

    pairwise["P_Value_FDR"] = pairwise["Group Comparison"].map(
        lambda c: lookup_adjusted(c, fdr_lookup)
    )
    return pairwise


def significance_label(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.0005:
        return "***"
    if p < 0.005:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def enrich_significance_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add raw/FDR significance and change flags to a pairwise HR table."""
    out = df.copy()
    out["raw_significance"] = out["P_Value"].apply(significance_label)
    out["FDR_values"] = out["P_Value_FDR"].apply(significance_label)
    out["fdr_significance_changed"] = out["raw_significance"].isin(["*", "**", "***"]) != out[
        "FDR_values"
    ].isin(["*", "**", "***"])
    return out


def apply_excel_formats(path: Path) -> None:
    """Center workbook cells and color significance labels."""
    green = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
    light_red = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")
    center = Alignment(horizontal="center", vertical="center")
    wb = load_workbook(path)
    for ws in wb.worksheets:
        headers = [cell.value for cell in ws[1]]
        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = center
        for col_name in ("raw_significance", "FDR_values", "fdr_significance", "Sig_raw", "Sig_fdr"):
            if col_name not in headers:
                continue
            col_idx = headers.index(col_name) + 1
            for row_idx in range(2, ws.max_row + 1):
                cell = ws.cell(row=row_idx, column=col_idx)
                value = str(cell.value).strip().lower()
                is_sig = value not in ("", "ns", "non-significant", "none", "nan")
                cell.fill = green if is_sig else light_red
    wb.save(path)


def save_excel(df: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(path, index=index, engine="openpyxl")
    apply_excel_formats(path)


def _prepare_plot_dataframe(df: pd.DataFrame, order_list: List[str]) -> pd.DataFrame:
    """Filter, order, and cast pairwise results for plotting."""
    plot_df = df[df["Group Comparison"].isin(order_list)].copy()
    plot_df["Group Comparison"] = pd.Categorical(
        plot_df["Group Comparison"],
        categories=order_list,
        ordered=True,
    )
    plot_df = plot_df.sort_values(by="Group Comparison").reset_index(drop=True)
    cols_to_convert = ["HR Ratio", "95% CI Upper", "95% CI Lower", "P_Value"]
    plot_df[cols_to_convert] = plot_df[cols_to_convert].apply(pd.to_numeric, errors="coerce")
    plot_df = plot_df.dropna(subset=cols_to_convert).reset_index(drop=True)
    plot_df["Group Comparison"] = plot_df["Group Comparison"].astype(str)
    return plot_df


def plot_ranking_original_style(
    df: pd.DataFrame,
    output_path: Path,
    order_list: List[str],
    *,
    title: str = "",
) -> None:
    intergroups_hr_df = _prepare_plot_dataframe(df, order_list)
    intergroups_hr_df["(95% CI)"] = intergroups_hr_df.apply(
        lambda row: f"({row['95% CI Lower']:.2f}-{row['95% CI Upper']:.2f})", axis=1
    )
    intergroups_hr_df["upper_error"] = intergroups_hr_df["95% CI Upper"] - intergroups_hr_df["HR Ratio"]
    intergroups_hr_df["lower_error"] = intergroups_hr_df["HR Ratio"] - intergroups_hr_df["95% CI Lower"]

    n_rows = len(intergroups_hr_df)
    if n_rows == 0:
        LOGGER.warning("No rows available for ranking forest plot.")
        return

    row_step = 2.0
    shared_y_positions = (np.arange(n_rows)[::-1]) * row_step
    bottom_padding = 1.15 * row_step
    top_padding = 0.9 * row_step

    # Wider overall figure; keep enough right-panel share so HR-axis ticks
    # do not crush, while the left table still has room for mathtext *P*/FDR.
    _fig_w = RANKING_FOREST_WIDTH_IN + 0.85
    _fig_h = max(
        2.2 * _RANKING_FOREST_SCALE,
        RANKING_FOREST_ROW_IN * n_rows + 1.2 * _RANKING_FOREST_SCALE,
    ) * _RANKING_FOREST_HEIGHT_SCALE
    _width_ratios = (1.55, 1.20)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(_fig_w, _fig_h),
        gridspec_kw={"wspace": 0.18, "width_ratios": list(_width_ratios)},
    )
    if title:
        fig.suptitle(
            title,
            fontsize=RANKING_FOREST_FONTS["title"],
            fontweight="bold",
            y=0.98,
        )

    ax = axes[0]
    # Gutters sized for "1 vs. 7-10", CI strings, and wide mathtext *P* values.
    # Extra 6 pt between HR (95%CI) and *P*.
    _cohort_shift_pt = 5.0
    _hr_p_gap_pt = 6.0
    _left_frac = _width_ratios[0] / (_width_ratios[0] + _width_ratios[1])
    _left_panel_in = _fig_w * _left_frac
    _data_span = 16.0
    _data_per_in = _data_span / max(_left_panel_in, 1e-6)
    _cohort_x = -0.7 - (_cohort_shift_pt / 72.0) * _data_per_in
    _hr_p_extra = (_hr_p_gap_pt / 72.0) * _data_per_in
    col_x_positions = [
        _cohort_x,
        3.6,
        9.2 + _hr_p_extra,
        13.6 + _hr_p_extra,
    ]
    _data_span = 16.0 + _hr_p_extra
    y_positions = shared_y_positions
    fontsize = RANKING_FOREST_FONTS["label"]
    donsitze_table = RANKING_FOREST_FONTS["annotation"]

    for j, (_, row) in enumerate(intergroups_hr_df.iterrows()):
        ax.text(col_x_positions[0], y_positions[j], f"{row['Group Comparison']}", va="center", ha="left", fontsize=donsitze_table)
        ax.text(
            col_x_positions[1],
            y_positions[j],
            f"{row['HR Ratio']:.2f} ({row['95% CI Lower']:.2f}-{row['95% CI Upper']:.2f})",
            va="center",
            ha="left",
            fontsize=donsitze_table,
        )
        ax.text(col_x_positions[2], y_positions[j], format_p_sci_mathtext(row["P_Value"]), va="center", ha="left", fontsize=donsitze_table)
        ax.text(col_x_positions[3], y_positions[j], format_p_sci_mathtext(row["P_Value_FDR"]), va="center", ha="left", fontsize=donsitze_table)

    y_position = y_positions.max() + top_padding
    ax.text(col_x_positions[0], y_position, "Cohorts", ha="left", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[1], y_position, "HR (95%CI)", ha="left", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[2], y_position, r"$\mathbf{\mathit{P}}$", ha="left", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[3], y_position, "FDR", ha="left", va="bottom", fontsize=fontsize, fontweight="bold")

    ax.set_xlim(col_x_positions[0] - 0.2, col_x_positions[0] + _data_span)
    ax.set_ylim(-bottom_padding, (n_rows - 1) * row_step + top_padding)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax = axes[1]
    ax.axvline(1, color="black", linestyle="--", lw=1)
    y_positions_point_plot = shared_y_positions
    sns.despine(ax=ax, left=True, right=True, top=False, bottom=False)
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True, labeltop=False)

    for j, (_, row) in enumerate(intergroups_hr_df.iterrows()):
        ax.plot(
            [row["95% CI Lower"], row["95% CI Upper"]],
            [y_positions_point_plot[j], y_positions_point_plot[j]],
            color="black",
            linewidth=_RANKING_FOREST_CI_LINEWIDTH,
        )
        ax.scatter(
            row["HR Ratio"],
            y_positions_point_plot[j],
            color="black",
            s=_RANKING_FOREST_POINT_SIZE,
            zorder=3,
        )

    ax.set_ylabel("", fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.set_xlim(0.25, 1.51)
    ax.set_xticks([0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    ax.set_ylim(-bottom_padding, (n_rows - 1) * row_step + top_padding)
    ax.set_xlabel("hazard ratio", fontsize=fontsize, labelpad=8)
    fig.subplots_adjust(bottom=0.24, top=0.88 if title else 0.95)
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def fit_cox_pairwise(input_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fit Cox PH model and compute pairwise hazard ratios.

    Uses an unpenalized CoxPHFitter (no L2 penalizer; Fig. 3g). Full
    one-hot coding is singular without a penalty, so one GROUP level is
    dropped as the reference (coef = 0, HR = 1) and re-inserted into the
    summary so pairwise contrasts still cover every group pair.
    """
    model_df = input_df[["TimeToDouble", "CENSOR", "GROUP"]].dropna().copy()
    # Preserve appearance order of groups in the input (rank labels are sorted
    # upstream); drop the first level as the unpenalized reference.
    group_levels = list(dict.fromkeys(model_df["GROUP"].astype(str).tolist()))
    dummy = pd.get_dummies(
        model_df["GROUP"].astype(str),
        drop_first=True,
    ).astype(int)
    cox_input = pd.concat([model_df[["TimeToDouble", "CENSOR"]], dummy], axis=1)

    cox = CoxPHFitter()
    cox.fit(cox_input, duration_col="TimeToDouble", event_col="CENSOR")
    summary = cox.summary.reset_index()
    if "covariate" not in summary.columns and "index" in summary.columns:
        summary = summary.rename(columns={"index": "covariate"})

    fitted = set(summary["covariate"].astype(str))
    for level in group_levels:
        if level in fitted:
            continue
        # Reference arm: HR = 1, SE = 0 so pairwise tests vs this arm use the
        # other group's SE alone (standard treatment coding).
        ref_row = {col: np.nan for col in summary.columns}
        ref_row["covariate"] = level
        ref_row["coef"] = 0.0
        ref_row["exp(coef)"] = 1.0
        ref_row["se(coef)"] = 0.0
        summary = pd.concat([pd.DataFrame([ref_row]), summary], ignore_index=True)

    pairwise = enrich_significance_columns(compute_pairwise_hr(summary))
    return summary, pairwise


def run_cox_pipeline(
    analysis_spec: AnalysisSpec,
    input_df: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    summary, pairwise = fit_cox_pairwise(input_df)

    summary_output = output_dir / "cox_summary.xlsx"
    save_excel(summary, summary_output)

    pairwise_output = output_dir / "pairwise_hr.xlsx"
    save_excel(pairwise, pairwise_output)

    forest_output = output_dir / "forest_plot.jpg"
    plot_ranking_original_style(
        pairwise,
        forest_output,
        analysis_spec.order_list,
        title="",
    )

    LOGGER.info(
        "Completed %s | rows=%d | output=%s",
        analysis_spec.run_label,
        len(input_df),
        output_dir,
    )
    return pairwise


def build_ranking_dataset(patient_df: pd.DataFrame) -> pd.DataFrame:
    """Assign dense per-patient rank groups for ranking analysis.

    Dense rank: 1 = highest LEVEL within patient; ranks >=7 pooled as '7-10'.
    """
    return prepare_rank_groups(patient_df.copy())


def main() -> None:
    args = parse_args()
    output_root = args.output_dir.resolve() / HR_OUTPUT_DIRNAME
    output_root.mkdir(parents=True, exist_ok=True)

    source_path = resolve_source_path(
        args.source.resolve() if args.source is not None else None
    )
    paths = Paths(
        project_root=REPO_ROOT,
        source_primary=source_path,
        source_sheet=args.sheet,
    )
    source_meta = {"selected_source": str(source_path), "sha256": sha256_of_file(source_path)}
    patient_source = load_patient_table(paths)
    ranking_spec = AnalysisSpec(
        run_label="ranking",
        order_list=list(RANKING_PAIRWISE_ORDER),
    )
    run_cox_pipeline(ranking_spec, build_ranking_dataset(patient_source), output_root)

    metadata = {
        "seed": SEED,
        "generated_at": datetime.now().isoformat(),
        "output_root": str(output_root),
        "source_sheet": paths.source_sheet,
        **source_meta,
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    LOGGER.info("Hazard-ratio analysis completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
