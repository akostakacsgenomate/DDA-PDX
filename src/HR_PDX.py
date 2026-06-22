"""Hazard-ratio analysis pipeline for PDX treatment cohorts.

Computes Cox proportional-hazards forest plots and pairwise HR tables for
three analysis branches:
  1. On-label / off-label DDA tier comparison (Low, Int, High, High On,
     High Off, High Off+Experimental, Standard-of-care chemotherapy).
  2. Per-patient compound rank (dense rank 1–6, pooled 7-10 vs SC).
  3. SHIVA-tier comparison (Low, Int, High, All, Non-assigned, SC).

All analyses share a single input workbook (data/), ensuring full
reproducibility from one source file.

Usage
-----
python HR_PDX.py [--source PATH] [--source-alternate PATH] [--output-dir DIR]

Random seed: SEED = 42 (set at module level via np.random.seed).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from math import log2
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


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

DEFAULT_SOURCE_XLSX = (
    DATA_DIR
    / "pdx_curve_metrics_single_treatments_dda_scores.xlsx"
)


@dataclass
class Paths:
    project_root: Path
    source_primary: Path
    chemo_csv: Path
    non_assigned_xlsx: Path


@dataclass
class AnalysisSpec:
    run_label: str
    source_script_stem: str
    legacy_name: str
    order_list: List[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cox proportional-hazards analysis for PDX treatment cohorts "
            "(on-label/off-label, rank, SHIVA tier)."
        )
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE_XLSX,
        help=(
            "Primary input workbook (recategorised PDX curve metrics). "
            f"Default: data/{DEFAULT_SOURCE_XLSX.name}"
        ),
    )
    parser.add_argument(
        "--source-alternate",
        type=Path,
        default=None,
        help="Optional alternate workbook for source-consistency verification.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help="Parent output directory. Timestamped subfolder created inside.",
    )
    return parser.parse_args()


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_and_validate_source(primary: Path, alternate: Path | None) -> Tuple[Path, Dict[str, str]]:
    if primary.exists():
        selected = primary
    elif alternate is not None and alternate.exists():
        selected = alternate
        LOGGER.warning("Primary source not found. Using alternate source: %s", alternate)
    else:
        alt_msg = f"\n- {alternate}" if alternate is not None else ""
        raise FileNotFoundError(f"Source file not found:\n- {primary}{alt_msg}")

    meta = {"selected_source": str(selected)}
    if primary.exists():
        meta["sha256_primary"] = sha256_of_file(primary)
    if alternate is not None and alternate.exists():
        meta["sha256_alternate"] = sha256_of_file(alternate)
    if "sha256_primary" in meta and "sha256_alternate" in meta:
        meta["source_files_identical"] = str(
            meta["sha256_primary"] == meta["sha256_alternate"]
        )
        if meta["source_files_identical"] == "False":
            LOGGER.warning("Primary and alternate source files differ. Using: %s", selected)
    return selected, meta


def normalize_id_columns(df: pd.DataFrame, model_col: str, compound_col: str) -> pd.DataFrame:
    out = df.copy()
    out[model_col] = out[model_col].astype(str).str.replace("-", "", regex=False).str.lower()
    out[compound_col] = (
        out[compound_col].astype(str).str.replace("-", "", regex=False).str.lower()
    )
    return out


def filesystem_safe_name(text: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_")
    normalized = [
        ch if ch in allowed else "_"
        for ch in text.strip().replace(" ", "_")
    ]
    safe = "".join(normalized)
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def clean_base_patient_table(patient_df: pd.DataFrame) -> pd.DataFrame:
    out = patient_df.copy()
    out = out[out["LEVEL"] != "#HIÁNYZIK"].copy()
    out["LEVEL"] = out["LEVEL"].astype(float)
    out["TimeToDouble"] = out["TimeToDouble"].astype(float)
    if "Day_Last" in out.columns:
        out["Day_Last"] = pd.to_numeric(out["Day_Last"], errors="coerce")
    out["CENSOR"] = True
    return out


def load_patient_table(source_xlsx: Path) -> pd.DataFrame:
    """Load and validate the shared recategorised source workbook."""
    source_df = pd.read_excel(source_xlsx)

    required_source_cols = {
        "Model",
        "COMPOUND",
        "LEVEL",
        "TimeToDouble",
        "APPROVED",
        "On_label",
    }
    missing_source = required_source_cols.difference(source_df.columns)
    if missing_source:
        raise KeyError(
            f"Source table is missing required columns: {sorted(missing_source)}"
        )

    source_df = normalize_id_columns(source_df, model_col="Model", compound_col="COMPOUND")
    return clean_base_patient_table(source_df)


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
            log_ratio = log2(ratio)
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


def plot_on_label_original_style(df: pd.DataFrame, output_path: Path, order_list: List[str]) -> None:
    intergroups_hr_df = _prepare_plot_dataframe(df, order_list)
    intergroups_hr_df["upper_error"] = (
        intergroups_hr_df["95% CI Upper"] - intergroups_hr_df["HR Ratio"]
    ).clip(lower=0)
    intergroups_hr_df["lower_error"] = (
        intergroups_hr_df["HR Ratio"] - intergroups_hr_df["95% CI Lower"]
    ).clip(lower=0)
    intergroups_hr_df_hr_table = intergroups_hr_df.copy()

    n_rows = len(intergroups_hr_df_hr_table)
    if n_rows == 0:
        LOGGER.warning("No rows available for on-label forest plot.")
        return

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.0, max(2.2, 0.52 * n_rows + 1.2)),
        gridspec_kw={"wspace": 0.21, "hspace": 8},
    )
    plt.subplots_adjust(left=0.08, right=0.98, wspace=0.8, hspace=15)

    ax = axes[0]
    col_x_positions = [-1.95, 2.0, 5.6, 7.8]
    space = 1.11
    y_positions = np.arange(n_rows) * -space + (n_rows - 0.49)
    fontsize = 10.5
    donsitze_table = 10.5

    for j, (_, row) in enumerate(intergroups_hr_df_hr_table.iterrows()):
        ax.text(col_x_positions[0], y_positions[j], f"{row['Group Comparison']}", va="center", fontsize=donsitze_table)
        ax.text(
            col_x_positions[1],
            y_positions[j],
            f"{row['HR Ratio']:.2f} ({row['95% CI Lower']:.2f}-{row['95% CI Upper']:.2f})",
            va="center",
            fontsize=donsitze_table,
        )
        ax.text(col_x_positions[2], y_positions[j], f"{float(row['P_Value']):.2e}", va="center", fontsize=donsitze_table)
        ax.text(col_x_positions[3], y_positions[j], f"{float(row['P_Value_FDR']):.2e}", va="center", fontsize=donsitze_table)

    y_position = y_positions[min(1, n_rows - 1)] + 2.1
    ax.text(col_x_positions[0] + 0.7, y_position, "Cohorts", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[1] + 1.65, y_position, "HR (95%CI)", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[2] + 0.75, y_position, r"$\mathbf{\mathit{p}}$", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[3] + 0.65, y_position, "FDR", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")

    plt.tight_layout()
    ax.set_xlim(-0.5, 8.5)
    ax.set_ylim(-1, n_rows)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax = axes[1]
    intergroups_hr_df["95% CI Lower"] = intergroups_hr_df["95% CI Lower"].astype(float)
    intergroups_hr_df["95% CI Upper"] = intergroups_hr_df["95% CI Upper"].astype(float)
    sns.pointplot(
        order=order_list,
        data=intergroups_hr_df,
        x="HR Ratio",
        y="Group Comparison",
        join=False,
        scale=1.0,
        color="black",
        capsize=0.9,
        ax=ax,
    )
    ax.axvline(1, color="gray", linestyle="--", lw=1)
    y_positions_point_plot = np.arange(len(intergroups_hr_df))
    sns.despine(ax=ax, left=True, right=True, top=False, bottom=False)
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True, labeltop=False)

    for j, (_, row) in enumerate(intergroups_hr_df.iterrows()):
        ax.plot(
            [row["95% CI Lower"], row["95% CI Upper"]],
            [y_positions_point_plot[j], y_positions_point_plot[j]],
            color="black",
            linewidth=2.0,
        )

    ax.set_ylabel("", fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize, pad=0)
    ax.set_xlabel("hazard ratio", fontsize=11, labelpad=8)
    ax.set_xlim(0.25, 1.51)
    ax.set_xticks([0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_ranking_original_style(df: pd.DataFrame, output_path: Path, order_list: List[str]) -> None:
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

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8, max(2.2, 0.52 * n_rows + 1.2)),
        gridspec_kw={"wspace": 0.21},
    )

    ax = axes[0]
    col_x_positions = [-0.6, 2.2, 6.25, 8.6]
    y_positions = shared_y_positions
    fontsize = 9.5
    donsitze_table = 9.5

    for j, (_, row) in enumerate(intergroups_hr_df.iterrows()):
        ax.text(col_x_positions[0], y_positions[j], f"{row['Group Comparison']}", va="center", fontsize=donsitze_table)
        ax.text(
            col_x_positions[1],
            y_positions[j],
            f"{row['HR Ratio']:.2f} ({row['95% CI Lower']:.2f}-{row['95% CI Upper']:.2f})",
            va="center",
            fontsize=donsitze_table,
        )
        ax.text(col_x_positions[2], y_positions[j], f"{row['P_Value']:.2e}", va="center", fontsize=donsitze_table)
        ax.text(col_x_positions[3], y_positions[j], f"{float(row['P_Value_FDR']):.2e}", va="center", fontsize=donsitze_table)

    y_position = y_positions.max() + top_padding
    ax.text(col_x_positions[0] + 0.77, y_position, "Cohorts", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[1] + 1.65, y_position, "HR (95%CI)", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[2] + 1.0, y_position, r"$\mathbf{\mathit{p}}$", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[3] + 0.7, y_position, "FDR", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")

    ax.set_xlim(-0.5, 9.0)
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
            linewidth=2,
        )
        ax.scatter(row["HR Ratio"], y_positions_point_plot[j], color="black", s=35, zorder=3)

    ax.set_ylabel("", fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.set_xlim(0.25, 1.51)
    ax.set_xticks([0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    ax.set_ylim(-bottom_padding, (n_rows - 1) * row_step + top_padding)
    ax.set_xlabel("", fontsize=40)
    fig.subplots_adjust(bottom=0.24)
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def _plot_shiva_single(
    plot_df: pd.DataFrame,
    output_path: Path,
    primary_title: str = "",
    xlim: Tuple[float, float] = (0.25, 2.7),
    xticks: Tuple[float, ...] = (0.5, 1.0, 1.5, 2.0, 2.5),
) -> None:
    """Render a single SHIVA-tier forest plot panel."""
    n_rows = len(plot_df)
    if n_rows == 0:
        LOGGER.warning("No rows available for SHIVA forest plot: %s", output_path.name)
        return

    row_step = 2.0
    shared_y_positions = (np.arange(n_rows)[::-1]) * row_step
    bottom_padding = 1.15 * row_step
    top_padding = 0.9 * row_step

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8, max(2.2, 0.52 * n_rows + 1.2)),
        gridspec_kw={"wspace": 0.21},
    )

    ax = axes[0]
    col_x_positions = [-2.7, 1.5, 6.0, 8.6]
    y_positions = shared_y_positions
    fontsize = 10.5
    donsitze_table = 10.5

    for j, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(col_x_positions[0], y_positions[j], f"{row['Group Comparison']}", va="center", fontsize=donsitze_table)
        ax.text(
            col_x_positions[1],
            y_positions[j],
            f"{row['HR Ratio']:.2f} ({row['95% CI Lower']:.2f}-{row['95% CI Upper']:.2f})",
            va="center",
            fontsize=donsitze_table,
        )
        ax.text(col_x_positions[2], y_positions[j], f"{float(row['P_Value']):.2e}", va="center", fontsize=donsitze_table)
        ax.text(col_x_positions[3], y_positions[j], f"{float(row['P_Value_FDR']):.2e}", va="center", fontsize=donsitze_table)

    y_position = y_positions.max() + top_padding
    ax.text(col_x_positions[0] + 1.3, y_position, "Cohorts", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[1] + 1.85, y_position, "HR (95%CI)", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[2] + 0.9, y_position, r"$\mathbf{\mathit{p}}$", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    ax.text(col_x_positions[3] + 0.6, y_position, "FDR", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")

    ax.set_xlim(-0.5, 9.0)
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
    ax.tick_params(axis="x", which="both", bottom=True, top=False, labelbottom=True, labeltop=False, pad=3)

    for j, (_, row) in enumerate(plot_df.iterrows()):
        ax.plot(
            [row["95% CI Lower"], row["95% CI Upper"]],
            [y_positions_point_plot[j], y_positions_point_plot[j]],
            color="black",
            linewidth=2,
        )
        ax.scatter(row["HR Ratio"], y_positions_point_plot[j], color="black", s=35, zorder=3)

    ax.set_ylabel("", fontsize=fontsize)
    ax.set_xlim(*xlim)
    ax.set_xticks(list(xticks))
    ax.set_ylim(-bottom_padding, (n_rows - 1) * row_step + top_padding)
    ax.set_title(primary_title, fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.set_xlabel("hazard ratio", fontsize=fontsize)
    fig.subplots_adjust(bottom=0.24)
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_shiva_original_style(df: pd.DataFrame, output_dir: Path) -> None:
    """Render SHIVA-tier forest plots (with and without Non-assigned group)."""
    shiva_plot_df = df.copy()
    shiva_plot_df["Group Comparison"] = (
        shiva_plot_df["Group Comparison"].astype(str).str.replace("Non_assigned", "NonA", regex=False)
    )

    order_list_nona = [
        "High vs. Int",
        "High vs. Low",
        "Int vs. Low",
        "NonA vs. High",
        "NonA vs. SC",
        "NonA vs. Int",
        "NonA vs. Low",
    ]
    order_list_levels = [
        "High vs. Int",
        "High vs. All",
        "High vs. Low",
        "Int vs. All",
        "Int vs. Low",
        "All vs. Low",
    ]

    plot_nona = _prepare_plot_dataframe(shiva_plot_df, order_list_nona)
    plot_levels = _prepare_plot_dataframe(shiva_plot_df, order_list_levels)

    _plot_shiva_single(
        plot_nona,
        output_dir / "forest_plot_nona.jpg",
        primary_title="",
    )
    # Canonical output filename for downstream traceability
    _plot_shiva_single(
        plot_levels,
        output_dir / "forest_plot.jpg",
        primary_title="",
    )


def run_cox_pipeline(
    analysis_spec: AnalysisSpec,
    input_df: pd.DataFrame,
    output_dir: Path,
) -> None:
    model_df = input_df[["TimeToDouble", "CENSOR", "GROUP"]].dropna().copy()
    dummy = pd.get_dummies(model_df["GROUP"], drop_first=False).astype(int)
    cox_input = pd.concat([model_df[["TimeToDouble", "CENSOR"]], dummy], axis=1)

    cox = CoxPHFitter(penalizer=0.1)
    cox.fit(cox_input, duration_col="TimeToDouble", event_col="CENSOR")
    summary = cox.summary.reset_index()

    summary_output = output_dir / "cox_summary.xlsx"
    summary.to_excel(summary_output, index=False)

    pairwise = compute_pairwise_hr(summary)
    pairwise = enrich_significance_columns(pairwise)
    pairwise_output = output_dir / "pairwise_hr.xlsx"
    pairwise.to_excel(pairwise_output, index=False)

    forest_output = output_dir / "forest_plot.jpg"
    if analysis_spec.run_label == "on_label":
        plot_on_label_original_style(pairwise, forest_output, analysis_spec.order_list)
    elif analysis_spec.run_label == "ranking":
        plot_ranking_original_style(pairwise, forest_output, analysis_spec.order_list)
    elif analysis_spec.run_label == "shiva":
        plot_shiva_original_style(pairwise, output_dir)
    else:
        LOGGER.warning("No plot style defined for run_label=%s", analysis_spec.run_label)

    LOGGER.info(
        "Completed %s | rows=%d | output=%s",
        analysis_spec.run_label,
        len(model_df),
        output_dir,
    )


def build_on_label_dataset(patient_df: pd.DataFrame, paths: Paths) -> pd.DataFrame:
    """Construct grouped survival table for on-label/off-label analysis.

    LEVEL tiers: Low < 0; 0 <= Int < 1000; High > 1000.
    On-label (On): LEVEL > 1000 AND On_label = IGAZ AND APPROVED = IGAZ.
    Off-label approved (Off): LEVEL > 1000 AND On_label = HAMIS AND APPROVED = IGAZ.
    Off+Experimental (Off+Ex): LEVEL > 1000 AND On_label = HAMIS (all approval statuses).
    SC: standard-of-care chemotherapy reference cohort.
    """
    chemo = pd.read_csv(paths.chemo_csv)
    chemo["CENSOR"] = True
    chemo["GROUP"] = "SC"

    d = patient_df.copy()
    g_int = d[(d["LEVEL"] >= 0) & (d["LEVEL"] < 1000)].copy()
    g_h = d[(d["LEVEL"] > 1000)].copy()
    g_h_on = d[(d["LEVEL"] > 1000) & (d["On_label"] == "IGAZ") & (d["APPROVED"] == "IGAZ")].copy()
    g_h_off = d[(d["LEVEL"] > 1000) & (d["On_label"] == "HAMIS") & (d["APPROVED"] == "IGAZ")].copy()
    g_h_off_ex = d[(d["LEVEL"] > 1000) & (d["On_label"] == "HAMIS")].copy()

    g_int["GROUP"] = "Int"
    g_h["GROUP"] = "H"
    g_h_on["GROUP"] = "H On"
    g_h_off["GROUP"] = "H Off"
    g_h_off_ex["GROUP"] = "H Off+Ex"

    combined = pd.concat([g_int, g_h, g_h_on, g_h_off, g_h_off_ex, chemo], ignore_index=True)
    return combined


def build_ranking_dataset(source_df: pd.DataFrame) -> pd.DataFrame:
    """Assign dense per-patient rank groups for ranking analysis.

    Dense rank: 1 = highest LEVEL within patient; ranks >=7 pooled as '7-10'.
    """
    d = clean_base_patient_table(normalize_id_columns(source_df, "Model", "COMPOUND"))
    return prepare_rank_groups(d)


def build_shiva_dataset(patient_df: pd.DataFrame, paths: Paths) -> pd.DataFrame:
    """Construct grouped survival table for SHIVA-tier analysis.

    SHIVA tiers by LEVEL: Low < 0; 0 <= Int < 1000; High >= 1000.
    All: all targeted compound rows combined.
    SC: standard-of-care chemotherapy reference cohort.
    Non_assigned: rows without molecular profiling assignment.
    """
    chemo = pd.read_csv(paths.chemo_csv)
    if "Treatment" in chemo.columns and "COMPOUND" not in chemo.columns:
        chemo = chemo.rename(columns={"Treatment": "COMPOUND"})
    chemo["CENSOR"] = True
    chemo = chemo[["Model", "COMPOUND", "TimeToDouble", "CENSOR"]].copy()
    chemo["GROUP"] = "SC"

    non_assigned = pd.read_excel(paths.non_assigned_xlsx)
    non_assigned["CENSOR"] = True
    non_assigned["GROUP"] = "Non_assigned"

    d = patient_df.copy()
    g_low = d[d["LEVEL"] < 0].copy()
    g_int = d[(d["LEVEL"] >= 0) & (d["LEVEL"] < 1000)].copy()
    g_high = d[d["LEVEL"] > 1000].copy()
    g_low["GROUP"] = "Low"
    g_int["GROUP"] = "Int"
    g_high["GROUP"] = "High"

    source_all = patient_df.copy()
    source_all["GROUP"] = "All"

    combined = pd.concat([g_low, g_int, g_high, chemo, source_all, non_assigned], ignore_index=True)
    return combined


def main() -> None:
    args = parse_args()
    timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
    output_root = args.output_dir.resolve() / f"HR_PDX_run_{timestamp}"
    output_root.mkdir(parents=True, exist_ok=True)

    alternate_source = (
        args.source_alternate.resolve() if args.source_alternate is not None else None
    )
    paths = Paths(
        project_root=REPO_ROOT,
        source_primary=args.source.resolve(),
        chemo_csv=DATA_DIR / "chemo_data_2024_12_06.csv",
        non_assigned_xlsx=DATA_DIR / "non_assigned_merged.xlsx",
    )

    selected_source, source_meta = resolve_and_validate_source(
        paths.source_primary, alternate_source
    )
    patient_source = load_patient_table(selected_source)

    on_label_spec = AnalysisSpec(
        run_label="on_label",
        source_script_stem="HR_On_label_final",
        legacy_name="HR_On_label_off_label+exp_plus_chemo_combined",
        order_list=[
            "H vs. SC",
            "Int vs. SC",
            "H On vs. H Off",
            "H On vs. H Off+Ex",
            "H On vs. SC",
            "H On vs. Int",
            "H Off vs. SC",
            "H Off vs. Int",
            "H Off vs. H Off+Ex",
            "H Off+Ex vs. SC",
            "H Off+Ex vs. Int",
        ],
    )
    ranking_spec = AnalysisSpec(
        run_label="ranking",
        source_script_stem="HR_ranking_final_OI_Jun2O26_rerank",
        legacy_name="HR_primary_ranking_rerank",
        order_list=[
            "1 vs. 2",
            "1 vs. 3",
            "1 vs. 4",
            "1 vs. 5",
            "1 vs. 6",
            "1 vs. 7-10",
        ],
    )
    shiva_spec = AnalysisSpec(
        run_label="shiva",
        source_script_stem="HR_SHIVA_final",
        legacy_name="SHIVA_HR_",
        order_list=[
            "High vs. Int",
            "High vs. Low",
            "Int vs. Low",
            "Non_assigned vs. High",
            "Non_assigned vs. SC",
            "Non_assigned vs. Int",
            "Non_assigned vs. Low",
            "High vs. All",
            "Int vs. All",
            "All vs. Low",
        ],
    )

    datasets = {
        "on_label": build_on_label_dataset(patient_source, paths),
        "ranking": build_ranking_dataset(patient_source),
        "shiva": build_shiva_dataset(patient_source, paths),
    }

    analysis_specs = [on_label_spec, ranking_spec, shiva_spec]
    fixed_output_names = {
        "on_label": "HR_utility",
        "ranking": "HR_Ranking",
        "shiva": "HR_SHIVA",
    }
    analysis_folders: Dict[str, str] = {}
    for spec in analysis_specs:
        subfolder_name = fixed_output_names.get(
            spec.run_label,
            filesystem_safe_name(f"{spec.source_script_stem}__{spec.legacy_name}"),
        )
        analysis_output_dir = output_root / subfolder_name
        analysis_output_dir.mkdir(parents=True, exist_ok=True)
        analysis_folders[spec.run_label] = str(analysis_output_dir)

        run_cox_pipeline(spec, datasets[spec.run_label], analysis_output_dir)

    metadata = {
        "seed": SEED,
        "generated_at": datetime.now().isoformat(),
        "output_root": str(output_root),
        "analysis_subfolders": analysis_folders,
        **source_meta,
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    LOGGER.info("All analyses completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
