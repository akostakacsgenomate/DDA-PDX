"""Boxplot and violin plot of alteration counts per tumor by classification.

Produces Driver / VUS / Non-Driver / All boxplots and matching violins (solid
median) from the VARIANTS sheet. All is the per-tumor sum of Driver + VUS +
Non-Driver only (other CLASSIFICATION labels are excluded).

Usage
-----
python box_plot_final.py [--source PATH] [--sheet SHEET] [--output-dir DIR]

Outputs are written to <output-dir>/variant_counts/ (default: script directory).
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_typography import apply_arial_font

apply_arial_font()


# ---------------------------------------------------------------------------
# Repository-relative paths (no hard-coded absolute paths)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

BOX_PLOT_OUTPUT_DIRNAME = "variant_counts"
DEFAULT_SOURCE_SHEET = "VARIANTS"
DEFAULT_NAME = "Fig2b_PDX_all"

DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"

WES_LIST = ["Driver", "VUS", "Non-Driver", "All"]
# Fig. 2b: three core classes + All (= per-tumor sum of those three only).
VIOLIN_LIST = ["Driver", "VUS", "Non-Driver", "All"]
CUSTOM_PALETTE = ["#A9CDE8", "#C5E3B8", "#9BB3C4", "#FFC4A0"]
VIOLIN_PALETTE = CUSTOM_PALETTE[: len(VIOLIN_LIST)]
# Match Fig. 2d compound distribution typography, reduced 20% for panel balance.
FIGURE_FONT_SIZE = 14
# Compact layout so panel b matches pie-chart panels in the Fig. 2 composite.
FIGURE_SIZE = (4.2, 5.8)
SAVE_PAD_INCHES = 0.12
Y_MAX = 130
Y_TICKS = [0, 10, 30, 50, 70, 90, 110, 130]
VIOLIN_BREAK_Y = 100
VIOLIN_UPPER_COMPRESS = 4
VIOLIN_LOWER_TICKS = list(range(0, 100, 10))
MEAN_LINEWIDTH_PT = 4
HALF_BAR = 0.28
CLASSIFICATION_LABELS = {
    "DRIVER": "Driver",
    "VUS": "VUS",
    "NON_DRIVER": "Non-Driver",
    "All": "All",
}
# Keep only these raw labels after VUS_IN_DRIVER → VUS merge (violin + boxplot).
# All is built exclusively from this set (never GENOMIC_MARKER / BIOMARKER / etc.).
CORE_CLASSIFICATIONS = frozenset({"DRIVER", "VUS", "NON_DRIVER"})
EXCLUDED_CLASSIFICATIONS = frozenset(
    {"GENOMIC_MARKER", "BIOMARKER", "NON_CONFIRMED_DRIVER"}
)

LOGGER = logging.getLogger("box_plot_final")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Boxplot and violin plot of alteration counts per tumor by classification."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help=(
            f"Input workbook override containing a VARIANTS sheet "
            f"(default: data/{DEFAULT_SOURCE_XLSX.name})."
        ),
    )
    parser.add_argument(
        "--sheet",
        type=str,
        default=DEFAULT_SOURCE_SHEET,
        help=f"Worksheet name (default: {DEFAULT_SOURCE_SHEET!r}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help=f"Parent output directory; results in <output-dir>/{BOX_PLOT_OUTPUT_DIRNAME}/.",
    )
    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help="Date stamp for output filenames (default: today's YYYY_MM_DD).",
    )
    return parser.parse_args()


def resolve_source_path(explicit: Optional[Path] = None) -> Path:
    """Resolve the variants workbook."""
    if explicit is not None:
        path = explicit.resolve()
        if not path.exists():
            raise FileNotFoundError(f"Source file not found: {path}")
        return path

    if DEFAULT_SOURCE_XLSX.exists() and _sheet_exists(DEFAULT_SOURCE_XLSX, DEFAULT_SOURCE_SHEET):
        return DEFAULT_SOURCE_XLSX.resolve()
    if DEFAULT_SOURCE_XLSX.exists():
        return DEFAULT_SOURCE_XLSX.resolve()
    raise FileNotFoundError(f"Source file not found: {DEFAULT_SOURCE_XLSX}")


def _sheet_exists(path: Path, sheet_name: str) -> bool:
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return False
    return sheet_name in xl.sheet_names


def load_variants(source_path: Path, sheet_name: str) -> pd.DataFrame:
    if not _sheet_exists(source_path, sheet_name):
        available = pd.ExcelFile(source_path).sheet_names
        raise ValueError(
            f"Sheet {sheet_name!r} not found in {source_path}. "
            f"Available sheets: {available}"
        )
    variants = pd.read_excel(source_path, sheet_name=sheet_name)
    required = {"PCT_ID", "CLASSIFICATION"}
    missing = required - set(variants.columns)
    if missing:
        raise ValueError(f"VARIANTS sheet missing columns: {sorted(missing)}")
    return variants


def prepare_boxplot_data(variants: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return long-form plot data and per-class descriptive statistics.

    All = per-tumor count of Driver + VUS + Non-Driver alterations only.
    """
    data = variants.copy()
    data = data[~data["CLASSIFICATION"].isin(EXCLUDED_CLASSIFICATIONS)]
    # Merge VUS-in-driver alterations into the VUS group (figure categories).
    data["CLASSIFICATION"] = data["CLASSIFICATION"].replace({"VUS_IN_DRIVER": "VUS"})
    data = data[data["CLASSIFICATION"].isin(CORE_CLASSIFICATIONS)].copy()
    dropped = sorted(
        set(variants["CLASSIFICATION"].dropna().astype(str)) - CORE_CLASSIFICATIONS - {"VUS_IN_DRIVER"}
    )
    if dropped:
        LOGGER.info(
            "Excluded from Driver/VUS/Non-Driver/All counts: %s",
            ", ".join(dropped),
        )

    wes_data = data.groupby(["PCT_ID", "CLASSIFICATION"]).size().reset_index(name="COUNT")
    # All is strictly the sum of the three core classes per tumor.
    all_data = data.groupby("PCT_ID").size().reset_index(name="COUNT")
    all_data["CLASSIFICATION"] = "All"

    box_plot_data = pd.concat([wes_data, all_data], ignore_index=True)
    box_plot_data["CLASSIFICATION"] = box_plot_data["CLASSIFICATION"].replace(
        CLASSIFICATION_LABELS
    )

    descriptive_statistics = (
        box_plot_data.groupby("CLASSIFICATION")["COUNT"]
        .agg(median="median", average="mean", sd="std")
        .reindex(WES_LIST)
        .reset_index()
    )
    return box_plot_data, descriptive_statistics


def violin_plot_data(box_plot_data: pd.DataFrame) -> pd.DataFrame:
    """Driver / VUS / Non-Driver / All (All = sum of the three core classes)."""
    return box_plot_data[box_plot_data["CLASSIFICATION"].isin(VIOLIN_LIST)].copy()


def _style_count_axes(ax: plt.Axes, *, limit_y: bool = True) -> None:
    if limit_y:
        ax.set_ylim(0, Y_MAX)
        ax.set_yticks(Y_TICKS)
    else:
        ax.set_ylim(bottom=0)
    ax.tick_params(axis="y", labelsize=FIGURE_FONT_SIZE)
    ax.set_xticks(range(len(WES_LIST)))
    ax.set_xticklabels(WES_LIST, rotation=90, fontsize=FIGURE_FONT_SIZE, ha="center")
    ax.set_xlabel(None)
    ax.set_ylabel("Number of alterations per tumor", fontsize=FIGURE_FONT_SIZE)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(right=False, top=False)


def save_boxplot(box_plot_data: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    sns.boxplot(
        x="CLASSIFICATION",
        y="COUNT",
        data=box_plot_data,
        palette=CUSTOM_PALETTE,
        order=WES_LIST,
        ax=ax,
    )
    _style_count_axes(ax)
    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=SAVE_PAD_INCHES,
    )
    plt.close(fig)


def _draw_violin_contents(ax: plt.Axes, plot_data: pd.DataFrame) -> None:
    """Draw violins with solid median bars on each category."""
    palette = dict(zip(VIOLIN_LIST, VIOLIN_PALETTE))
    sns.violinplot(
        x="CLASSIFICATION",
        y="COUNT",
        data=plot_data,
        hue="CLASSIFICATION",
        palette=palette,
        order=VIOLIN_LIST,
        hue_order=VIOLIN_LIST,
        legend=False,
        inner=None,
        cut=0,
        density_norm="width",
        linewidth=1,
        linecolor="black",
        saturation=1,
        ax=ax,
    )
    median = (
        plot_data.groupby("CLASSIFICATION")["COUNT"].median().reindex(VIOLIN_LIST)
    )
    for x, cat in enumerate(VIOLIN_LIST):
        y = median.get(cat)
        if y is None or pd.isna(y):
            continue
        ax.hlines(
            float(y),
            x - HALF_BAR,
            x + HALF_BAR,
            colors="black",
            linewidth=MEAN_LINEWIDTH_PT,
            linestyles="solid",
            zorder=4,
            capstyle="butt",
        )


def _upper_yticks(break_y: float, y_top: float) -> list[float]:
    step = 100
    start = (int(break_y) // step + 1) * step
    ticks = [float(t) for t in range(start, int(y_top) + 1, step)]
    if not ticks or ticks[-1] < y_top - 1:
        ticks.append(float(y_top))
    return ticks


def _add_y_axis_break_marks(ax_top: plt.Axes, ax_bottom: plt.Axes) -> None:
    """Diagonal marks on the left y-spine to show the axis break."""
    d = 0.015
    kwargs = dict(color="black", linewidth=1, clip_on=False, solid_capstyle="butt")
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), transform=ax_bottom.transAxes, **kwargs)
    ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes, **kwargs)


def _add_close_ylabel(fig: plt.Figure, ax_top: plt.Axes, ax_bottom: plt.Axes) -> None:
    """Place the y-label just left of the tick labels, without overlapping them."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    top_pos = ax_top.get_position()
    bot_pos = ax_bottom.get_position()
    y_mid = 0.5 * (top_pos.y1 + bot_pos.y0)
    spine_px = float(ax_bottom.transAxes.transform((0.0, 0.0))[0])
    leftmost_px = spine_px
    for ax in (ax_top, ax_bottom):
        for lab in ax.get_yticklabels():
            if not lab.get_text():
                continue
            leftmost_px = min(leftmost_px, float(lab.get_window_extent(renderer=renderer).x0))
    tick_extent_pt = (spine_px - leftmost_px) * 72.0 / fig.dpi
    gap_pt = 3.0
    ax_bottom.annotate(
        "Number of alterations per tumor",
        xy=(0.0, y_mid),
        xycoords=("axes fraction", "figure fraction"),
        xytext=(-(tick_extent_pt + gap_pt + FIGURE_FONT_SIZE / 2.0), 0.0),
        textcoords="offset points",
        rotation=90,
        va="center",
        ha="center",
        fontsize=FIGURE_FONT_SIZE,
        annotation_clip=False,
        clip_on=False,
    )


def save_violinplot(box_plot_data: pd.DataFrame, output_path: Path) -> None:
    """Violin plot (Driver/VUS/Non-Driver/All) with mean lines and axis break."""
    plot_data = violin_plot_data(box_plot_data)
    if plot_data.empty:
        raise ValueError(
            "No Driver/VUS/Non-Driver/All rows available for the violin plot."
        )
    data_max = float(plot_data["COUNT"].max())
    y_top = max(VIOLIN_BREAK_Y + 50, math.ceil(data_max / 50) * 50)
    lower_span = float(VIOLIN_BREAK_Y)
    upper_span = (y_top - VIOLIN_BREAK_Y) / VIOLIN_UPPER_COMPRESS

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=FIGURE_SIZE,
        gridspec_kw={"height_ratios": [upper_span, lower_span]},
    )
    _draw_violin_contents(ax_top, plot_data)
    _draw_violin_contents(ax_bottom, plot_data)

    for ax in (ax_top, ax_bottom):
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks(range(len(VIOLIN_LIST)))
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=FIGURE_FONT_SIZE, pad=2, right=False, top=False)
        ax.tick_params(axis="x", labelsize=FIGURE_FONT_SIZE)
        ax.set_autoscaley_on(False)

    ax_top.spines["top"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.tick_params(bottom=False, labelbottom=False, top=False)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.set_xticklabels(VIOLIN_LIST, rotation=90, fontsize=FIGURE_FONT_SIZE, ha="center")

    ax_bottom.set_ylim(0, VIOLIN_BREAK_Y)
    ax_bottom.set_yticks(VIOLIN_LOWER_TICKS)
    ax_top.set_ylim(VIOLIN_BREAK_Y, y_top)
    ax_top.set_yticks(_upper_yticks(VIOLIN_BREAK_Y, y_top))

    _add_y_axis_break_marks(ax_top, ax_bottom)
    fig.subplots_adjust(hspace=0.12, left=0.20, right=0.97, top=0.98, bottom=0.16)
    _add_close_ylabel(fig, ax_top, ax_bottom)
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=SAVE_PAD_INCHES,
    )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    date_stamp = args.date or datetime.now().strftime("%Y_%m_%d")
    source_path = resolve_source_path(args.source)
    sheet_name = args.sheet.strip() if args.sheet and args.sheet.strip() else DEFAULT_SOURCE_SHEET

    output_root = args.output_dir.resolve() / BOX_PLOT_OUTPUT_DIRNAME
    output_root.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Loading variants from %s (sheet=%s)", source_path, sheet_name)
    variants = load_variants(source_path, sheet_name)
    box_plot_data, descriptive_statistics = prepare_boxplot_data(variants)

    LOGGER.info("Descriptive statistics:\n%s", descriptive_statistics.to_string(index=False))

    stats_path = output_root / f"{DEFAULT_NAME}_descriptive_statistics_classification_{date_stamp}.xlsx"
    box_path = output_root / f"{DEFAULT_NAME}_boxplot_classification_{date_stamp}.jpeg"
    violin_path = output_root / f"{DEFAULT_NAME}_violinplot_classification_{date_stamp}.jpeg"

    descriptive_statistics.to_excel(stats_path, index=False)
    save_boxplot(box_plot_data, box_path)
    save_violinplot(box_plot_data, violin_path)

    metadata: Dict[str, object] = {
        "generated_at": datetime.now().isoformat(),
        "source_path": str(source_path),
        "source_sheet": sheet_name,
        "output_root": str(output_root),
        "outputs": {
            "boxplot": str(box_path),
            "violinplot": str(violin_path),
            "descriptive_statistics": str(stats_path),
        },
        "categories": WES_LIST,
        "violin_categories": VIOLIN_LIST,
        "n_variant_rows": int(len(variants)),
        "n_boxplot_rows": int(len(box_plot_data)),
        "n_violin_rows": int(len(violin_plot_data(box_plot_data))),
        "violin_means": {
            str(k): float(v)
            for k, v in (
                violin_plot_data(box_plot_data)
                .groupby("CLASSIFICATION")["COUNT"]
                .mean()
                .reindex(VIOLIN_LIST)
                .items()
            )
        },
        "violin_medians": {
            str(k): float(v)
            for k, v in (
                violin_plot_data(box_plot_data)
                .groupby("CLASSIFICATION")["COUNT"]
                .median()
                .reindex(VIOLIN_LIST)
                .items()
            )
        },
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    LOGGER.info("box_plot_final completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
