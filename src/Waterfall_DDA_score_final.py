"""Generate all-in-one DDA score waterfall plots coloured by mRECIST response.

Produces one plot per endpoint configuration (BestAvgResponse and
Response_Best_Response).

Usage
-----
python Waterfall_DDA_score_final.py [--source PATH] [--project-root DIR] [--run-stamp STAMP]

Outputs are written to <project-root>/score_distributions/ (shared with
scores_dist_plots.py and Case_level_distribution.R). Pass the same --run-stamp to
all three scripts if coordinating a single pipeline run.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from plot_typography import apply_arial_font, composite_cell_in, panel_font_sizes

apply_arial_font()

SCRIPT_NAME = "Waterfall_DDA_score_final"
RUN_FOLDER_PREFIX = "score_distributions"

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_SOURCE_SHEET = "PCT_DRUG_RESPONSE"

# Each entry produces one waterfall figure (sort column + matching mRECIST hue).
ENDPOINT_CONFIGS = [
    {
        "endpoint": "BestAvgResponse",
        "sort_column": "BestAvgResponse",
        "response_column": "Response_BestAvgResponse",
    },
    {
        "endpoint": "Response_Best_Response",
        "sort_column": "BestResponse",
        "response_column": "Response_Best_Response",
    },
]

RESPONSE_CUSTOM_PALETTE = {
    "PR": "#1E90FF",
    "SD": "#1ecbff",
    "CR": "#1e4fff",
    "PD": "#FF0000",
}

# Typography: 2x Fig. 2d base size (18 pt), +10% for composite legibility.
FIGURE_FONT_SIZE = 36
LEGEND_MARKER_SIZE = 30  # 3x prior 10 pt mRECIST colour swatches
SAVE_PAD_INCHES = 0.12

# Publication panel titles (distinct from internal endpoint column names).
DISPLAY_TITLES = {
    "Response_Best_Response": "Best response",
    "BestAvgResponse": "Best average response",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an all-in-one DDA score waterfall plot coloured by mRECIST response."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE_XLSX,
        help=f"Input workbook (default: data/{DEFAULT_SOURCE_XLSX.name})",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=_SCRIPT_DIR,
        help="Output root directory (timestamped subfolder created inside).",
    )
    parser.add_argument(
        "--run-stamp",
        type=str,
        default=None,
        help=(
            "Run folder timestamp suffix (default: current time). "
            "Use the same value across scores_dist_plots.py, "
            "Waterfall_DDA_score_final.py, and Case_level_distribution.R."
        ),
    )
    return parser.parse_args()


def build_output_dirs(project_root: Path, run_stamp: str) -> Dict[str, Path]:
    run_root = project_root / RUN_FOLDER_PREFIX
    figures_dir = run_root / "figures"
    tables_dir = run_root / "tables"
    run_root.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    return {"run_root": run_root, "figures": figures_dir, "tables": tables_dir}


def harmonize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    supplementary_renames = {
        "PCT_ID": "Model",
        "DDA score": "LEVEL",
        "TimeToDouble (days)": "TimeToDouble",
        "BestResponse (%)": "BestResponse",
        "BestAvgResponse (%)": "BestAvgResponse",
        "Response BestResponse": "Response_Best_Response",
        "Response BestAvgResponse": "Response_BestAvgResponse",
        "DDA_score_tier": "SHIVA_scores",
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
    return out


def read_source_table(source_xlsx: Path, sheet: str = DEFAULT_SOURCE_SHEET) -> pd.DataFrame:
    if source_xlsx.suffix.lower() == ".csv":
        return pd.read_csv(source_xlsx)
    return pd.read_excel(source_xlsx, sheet_name=sheet)


def prepare_data(source_xlsx: Path) -> pd.DataFrame:
    if not source_xlsx.exists():
        raise FileNotFoundError(f"Input workbook not found: {source_xlsx}")

    df = harmonize_column_names(read_source_table(source_xlsx))

    if "Treatment type" in df.columns:
        df = df[df["Treatment type"].astype(str).str.strip().str.lower() == "single"].copy()
    if "On_label" in df.columns:
        df = df[df["On_label"].astype(str).str.strip().str.lower() != "chemo"].copy()

    if "SHIVA_scores" not in df.columns:
        df["SHIVA_scores"] = pd.NA
        df.loc[df["LEVEL"] < 0, "SHIVA_scores"] = "LOW"
        df.loc[(df["LEVEL"] >= 0) & (df["LEVEL"] < 1000), "SHIVA_scores"] = "INTERMEDIATE"
        df.loc[df["LEVEL"] >= 1000, "SHIVA_scores"] = "HIGH"

    if "Response_Best_Response" not in df.columns and "ResponseCategory" in df.columns:
        df["Response_Best_Response"] = (
            df["ResponseCategory"].astype(str).str.split("-->").str[0].str.strip()
        )

    required_columns = ["LEVEL", "SHIVA_scores"]
    for config in ENDPOINT_CONFIGS:
        required_columns.extend([config["sort_column"], config["response_column"]])
    required_columns = list(dict.fromkeys(required_columns))

    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(
            "Input workbook is missing required columns: "
            f"{missing_columns}. Present columns: {list(df.columns)}"
        )

    return df


def create_waterfall_plot(
    df: pd.DataFrame,
    output_path: Path,
    *,
    endpoint: str,
    sort_column: str,
    response_column: str,
) -> Path:
    # X-axis order follows LEVEL descending: highest DDA score at index 0 (left).
    df_sorted = df.sort_values(
        by=["LEVEL", sort_column], ascending=[False, False]
    ).reset_index(drop=True)
    df_sorted["ResponseColor"] = (
        df_sorted[response_column]
        .astype(str)
        .str.strip()
        .str.upper()
        .map(RESPONSE_CUSTOM_PALETTE)
        .fillna("#808080")
    )

    plt.rcParams.update(
        {
            "font.size": FIGURE_FONT_SIZE,
            "axes.titlesize": FIGURE_FONT_SIZE,
            "axes.labelsize": FIGURE_FONT_SIZE,
            "legend.fontsize": FIGURE_FONT_SIZE,
            "xtick.labelsize": FIGURE_FONT_SIZE,
            "ytick.labelsize": FIGURE_FONT_SIZE,
        }
    )

    fig = plt.figure(figsize=(22, 12))
    ax = fig.add_subplot(111)

    ax.bar(
        x=df_sorted.index,
        height=df_sorted["LEVEL"],
        color=df_sorted["ResponseColor"],
        width=0.8,
    )

    ax.axhline(0, color="black", linestyle="-", linewidth=1)
    ax.set_yscale("symlog")
    ax.set_yticks(
        [-100000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000, 100000]
    )

    display_title = DISPLAY_TITLES.get(endpoint, endpoint)
    ax.set_xlabel("Treatment", fontsize=FIGURE_FONT_SIZE)
    ax.set_ylabel("DDA Score (log)", fontsize=FIGURE_FONT_SIZE)
    ax.set_title(display_title, fontsize=FIGURE_FONT_SIZE)
    ax.title.set_position([0.5, 1.05])
    ax.set_xticks(np.arange(0, 1300, 100))
    ax.tick_params(axis="both", labelsize=FIGURE_FONT_SIZE)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            label=label,
            markerfacecolor=color,
            markersize=LEGEND_MARKER_SIZE,
        )
        for label, color in RESPONSE_CUSTOM_PALETTE.items()
    ]
    ax.legend(
        handles=legend_handles,
        title="mRECIST",
        fontsize=FIGURE_FONT_SIZE,
        title_fontsize=FIGURE_FONT_SIZE,
        frameon=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=SAVE_PAD_INCHES,
    )
    plt.close(fig)

    return output_path


RANK_GROUP_ORDER = [str(i) for i in range(1, 7)] + ["7-10"]
DCR_STATES = {"CR", "PR", "SD"}
ORR_STATES = {"CR", "PR"}

# Fig. 3d/3e occupy a half-width cell of the 4x6 Fig. 3 composite grid. Matching
# the cell aspect means the panel fills its slot, and the font sizes below are
# derived so the text lands at the shared on-page point sizes.
FIG3_GRID = dict(nrows=4, ncols=6, rowspan=1, colspan=3)
_FIG3_CELL_W, _FIG3_CELL_H = composite_cell_in(**FIG3_GRID)
FIG3_WIDTH = 8.0
FIG3_HEIGHT = FIG3_WIDTH * _FIG3_CELL_H / _FIG3_CELL_W
FIG3_FONTS = panel_font_sizes(FIG3_WIDTH, FIG3_HEIGHT, **FIG3_GRID)

FIG3_AXIS_TEXT = FIG3_FONTS["tick"]
FIG3_AXIS_LABEL = FIG3_FONTS["label"]
FIG3_TITLE = FIG3_FONTS["title"]
FIG3_ANNOTATION = FIG3_FONTS["annotation"]
FIG3_LEGEND_TEXT = FIG3_FONTS["legend"]
FIG3_LEGEND_TITLE = FIG3_FONTS["legend"]
FIG3_RANK_GAP = 16.0 * (FIG3_WIDTH / 22)
FIG3_Y_TICKS = [-100000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000, 100000]

FIG3_MARGINS = dict(left=0.115, right=0.99, top=0.85, bottom=0.12)
# Rank groups sit close together on the right of the axis, so the DCR/ORR header
# carries bare percentages with a single row label instead of per-group prefixes.
FIG3_ORR_ROW_Y = 1.02
FIG3_DCR_ROW_Y = 1.085
FIG3_TITLE_PAD = 40

FIG3_BYRANK_TITLES = {
    "Best_Response": "BestResponse",
    "BestAvgResponse": "BestAvgResponse",
}


def _symlog_tick_labels(values):
    labels = []
    for v in values:
        if v == 0:
            labels.append("0")
            continue
        exp = int(round(np.log10(abs(v))))
        if np.isclose(abs(v), 10**exp):
            labels.append(rf"$-10^{{{exp}}}$" if v < 0 else rf"$10^{{{exp}}}$")
        else:
            labels.append(str(v))
    return labels


def _assign_rank_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Dense-rank LEVEL within each Model; pool ranks 7+ as '7-10'."""
    out = df.copy()
    out["_rank"] = (
        out.groupby("Model", group_keys=False)["LEVEL"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )
    out["RankGroup"] = out["_rank"].apply(lambda r: str(r) if r <= 6 else "7-10")
    return out


def _compute_group_rates(df: pd.DataFrame, response_column: str) -> pd.DataFrame:
    resp = df[response_column].astype(str).str.strip().str.upper()
    rows = []
    for g in RANK_GROUP_ORDER:
        sub = resp[df["RankGroup"] == g]
        n = len(sub)
        if n == 0:
            rows.append({"RankGroup": g, "n": 0, "DCR": 0.0, "ORR": 0.0})
            continue
        rows.append({
            "RankGroup": g,
            "n": n,
            "DCR": 100.0 * sub.isin(DCR_STATES).mean(),
            "ORR": 100.0 * sub.isin(ORR_STATES).mean(),
        })
    return pd.DataFrame(rows)


def create_byrank_dcr_waterfall(
    df: pd.DataFrame,
    output_path: Path,
    *,
    response_column: str,
    column_label: str,
) -> Path:
    """By-rank DDA-score waterfall with DCR/ORR annotations above top axis."""
    ranked = _assign_rank_groups(df)
    plot_df = ranked[ranked["RankGroup"].isin(RANK_GROUP_ORDER)].copy()
    resp_upper = plot_df[response_column].astype(str).str.strip().str.upper()
    plot_df["ResponseColor"] = resp_upper.map(RESPONSE_CUSTOM_PALETTE).fillna("#808080")
    plot_df["_rsort"] = plot_df["RankGroup"].map(
        {g: i for i, g in enumerate(RANK_GROUP_ORDER)}
    )
    plot_df = plot_df.sort_values(
        ["_rsort", "LEVEL"], ascending=[True, False]
    ).reset_index(drop=True)

    bar_w = 0.8
    parts, x_cursor, sep_xs = [], 0.0, []
    for gi, g in enumerate(RANK_GROUP_ORDER):
        gdf = plot_df[plot_df["RankGroup"] == g].copy()
        if gdf.empty:
            continue
        n = len(gdf)
        gdf["_x"] = x_cursor + np.arange(n, dtype=float)
        parts.append(gdf)
        x_cursor += n
        if gi < len(RANK_GROUP_ORDER) - 1:
            has_later = any(
                not plot_df[plot_df["RankGroup"] == RANK_GROUP_ORDER[j]].empty
                for j in range(gi + 1, len(RANK_GROUP_ORDER))
            )
            if has_later:
                sep_xs.append(x_cursor + FIG3_RANK_GAP / 2)
                x_cursor += FIG3_RANK_GAP

    bar_df = pd.concat(parts, ignore_index=True)

    # create_waterfall_plot() leaves its own (much larger) sizes in rcParams.
    plt.rcParams.update({
        "font.size": FIG3_AXIS_TEXT,
        "axes.titlesize": FIG3_TITLE,
        "axes.labelsize": FIG3_AXIS_LABEL,
        "legend.fontsize": FIG3_LEGEND_TEXT,
        "xtick.labelsize": FIG3_AXIS_TEXT,
        "ytick.labelsize": FIG3_AXIS_TEXT,
    })

    fig = plt.figure(figsize=(FIG3_WIDTH, FIG3_HEIGHT))
    ax = fig.add_subplot(111)

    for xs in sep_xs:
        ax.axvline(xs, color="#d8d8d8", linestyle="-", linewidth=2, zorder=0)

    ax.bar(x=bar_df["_x"], height=bar_df["LEVEL"], color=bar_df["ResponseColor"],
           width=bar_w, zorder=2)
    ax.axhline(0, color="black", linestyle="-", linewidth=1, zorder=1)
    ax.set_yscale("symlog")
    ax.set_yticks(FIG3_Y_TICKS)
    ax.set_yticklabels(_symlog_tick_labels(FIG3_Y_TICKS))

    group_centers = bar_df.groupby("RankGroup", sort=False)["_x"].mean()
    ax.set_xticks(group_centers.values)
    ax.set_xticklabels(group_centers.index, fontsize=FIG3_AXIS_TEXT)
    ax.tick_params(axis="y", labelsize=FIG3_AXIS_TEXT)

    rank_rates = _compute_group_rates(plot_df, response_column)
    for label, y, key in (("DCR (%)", FIG3_DCR_ROW_Y, "DCR"), ("ORR (%)", FIG3_ORR_ROW_Y, "ORR")):
        ax.text(
            -0.015, y, label,
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FIG3_ANNOTATION, zorder=4,
        )
        for _, row in rank_rates.iterrows():
            g = row["RankGroup"]
            if g not in group_centers.index or row["n"] == 0:
                continue
            ax.text(
                group_centers[g], y, f"{row[key]:.1f}",
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=FIG3_ANNOTATION, zorder=4,
            )

    ax.set_xlabel("Drug Rank Within Tumor", fontsize=FIG3_AXIS_LABEL)
    ax.set_ylabel("DDA Score (log)", fontsize=FIG3_AXIS_LABEL)
    ax.set_title(
        FIG3_BYRANK_TITLES.get(column_label, column_label),
        fontsize=FIG3_TITLE, pad=FIG3_TITLE_PAD,
    )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    legend_handles = [
        Line2D([0], [0], marker="s", color="w", label=lab,
               markerfacecolor=col, markersize=FIG3_LEGEND_TEXT * 0.6)
        for lab, col in RESPONSE_CUSTOM_PALETTE.items()
    ]
    # Legend inside the axes, beside the y-axis at the start of rank 1.
    ax.legend(
        handles=legend_handles, title="mRECIST",
        fontsize=FIG3_LEGEND_TEXT, title_fontsize=FIG3_LEGEND_TITLE,
        frameon=True, fancybox=False, edgecolor="none", framealpha=0.9,
        loc="lower left", bbox_to_anchor=(0.0, 0.0),
        ncol=1, columnspacing=0.8,
        handletextpad=0.3, borderpad=0.35, labelspacing=0.25,
    )

    x_left = float(bar_df["_x"].min() - bar_w / 2)
    x_right = float(bar_df["_x"].max() + bar_w / 2)
    ax.set_xlim(x_left, x_right)

    xaxis_trans = ax.get_xaxis_transform()
    edge_color = ax.spines["left"].get_edgecolor()
    edge_lw = ax.spines["left"].get_linewidth()
    for g in RANK_GROUP_ORDER:
        gdf = bar_df[bar_df["RankGroup"] == g]
        if gdf.empty:
            continue
        x_min = float(gdf["_x"].min() - bar_w / 2)
        x_max = float(gdf["_x"].max() + bar_w / 2)
        ax.plot([x_min, x_max], [0, 0], transform=xaxis_trans,
                color=edge_color, linewidth=edge_lw, clip_on=False, zorder=3)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(**FIG3_MARGINS)
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def save_run_metadata(
    run_root: Path,
    source_xlsx: Path,
    produced_files: List[Path],
    run_stamp: str,
    row_count: int,
    endpoints: List[str],
) -> Path:
    metadata = {
        "script_name": SCRIPT_NAME,
        "run_stamp": run_stamp,
        "source_xlsx": str(source_xlsx.resolve()),
        "rows_used": row_count,
        "endpoints": endpoints,
        "produced_files": [str(path) for path in produced_files],
    }
    metadata_path = run_root / f"run_metadata_{SCRIPT_NAME}.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata_path


def main() -> None:
    args = parse_args()
    source_xlsx = args.source
    project_root = args.project_root

    run_stamp = args.run_stamp or datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    date_stamp = datetime.now().strftime("%Y_%m_%d")
    output_dirs = build_output_dirs(project_root=project_root, run_stamp=run_stamp)

    patient_data = prepare_data(source_xlsx=source_xlsx)

    output_paths: List[Path] = []
    endpoint_names: List[str] = []
    for config in ENDPOINT_CONFIGS:
        endpoint = config["endpoint"]
        if endpoint == "Response_Best_Response":
            panel = "Fig2e"
        elif endpoint == "BestAvgResponse":
            panel = "Fig2f"
        else:
            panel = "Fig2"
        output_path = (
            output_dirs["figures"]
            / f"{panel}_PDX_DDA_score_{endpoint}_all_inone_waterfall_Plot_{date_stamp}.png"
        )
        create_waterfall_plot(
            patient_data,
            output_path=output_path,
            endpoint=endpoint,
            sort_column=config["sort_column"],
            response_column=config["response_column"],
        )
        output_paths.append(output_path)
        endpoint_names.append(endpoint)

    byrank_configs = [
        ("Fig3d", "Response_Best_Response", "Best_Response"),
        ("Fig3e", "Response_BestAvgResponse", "BestAvgResponse"),
    ]
    for panel, resp_col, label in byrank_configs:
        if resp_col not in patient_data.columns:
            continue
        byrank_path = (
            output_dirs["figures"]
            / f"{panel}_F2B_PDX_DDA_waterfall_{label}_{date_stamp}.png"
        )
        create_byrank_dcr_waterfall(
            patient_data, byrank_path,
            response_column=resp_col, column_label=label,
        )
        output_paths.append(byrank_path)
        endpoint_names.append(f"byrank_{label}")

    metadata_path = save_run_metadata(
        run_root=output_dirs["run_root"],
        source_xlsx=source_xlsx,
        produced_files=output_paths,
        run_stamp=run_stamp,
        row_count=len(patient_data),
        endpoints=endpoint_names,
    )

    print("Analysis completed.")
    print(f"Run folder:               {output_dirs['run_root']}")
    print(f"Run metadata:             {metadata_path}")
    for output_path in output_paths:
        print(f"Generated figure:         {output_path}")


if __name__ == "__main__":
    main()
