"""Generate Fig. 1 data-processing flowchart.

Reproduces the manuscript panel describing how the PDX clinical trial (PCT)
dataset of Gao et al. was filtered into the validation cohort: MTA monotherapy
filtering, molecular-profile filtering, DDA score generation, DDA-score
filtering, and the final predictivity evaluation.

Geometry is expressed in the pixel coordinate system of the reference artwork
(863 x 1024, y increasing downward) so that box positions, table column
alignment, arrow routing and font sizes match the manuscript layout exactly.

Usage
-----
python flowchart_final.py [--output-dir DIR] [--date YYYY_MM_DD]

Outputs are written to <output-dir>/flowchart/ (default: script directory).
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

from plot_typography import apply_arial_font

apply_arial_font()


# ---------------------------------------------------------------------------
# Repository-relative paths (no hard-coded absolute paths)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

FLOWCHART_OUTPUT_DIRNAME = "flowchart"
OUTPUT_BASENAME = "Fig1_data_processing_flowchart"
DEFAULT_VARIANTS_SHEET = "VARIANTS"
DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
CHEMO_DATA_CSV = DATA_DIR / "chemo_data_2024_12_06.csv"

LOGGER = logging.getLogger("flowchart_final")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Reference artwork canvas: 1 unit == 1 px at 100 px/inch.
CANVAS_W, CANVAS_H = 863.0, 1024.0
PX_PER_INCH = 100.0
FIGURE_SIZE = (CANVAS_W / PX_PER_INCH, CANVAS_H / PX_PER_INCH)

BULLET = "\u2022"

FONT_STACK = ["Arial", "DejaVu Sans"]

BOX_LW = 1.0
RULE_LW = 1.0
ARROW_LW = 1.0
ARROW_SCALE = 11

FS_TITLE = 10.0
FS_BODY = 10.0
FS_SUBHEAD = 9.5

# Shared numeric column centres for both exclusion tables.
COL_T_EXCLUDED = 505.0
COL_T_REMAINING = 575.0
COL_U_EXCLUDED = 648.0
COL_U_REMAINING = 718.0
COL_DIVIDER_X = 612.0
GROUP_T_CENTER = 540.0
GROUP_U_CENTER = 683.0
RULE_T_SPAN = (478.0, 606.0)
RULE_U_SPAN = (620.0, 751.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the Fig. 1 data-processing flowchart."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help=f"Parent output directory; results in <output-dir>/{FLOWCHART_OUTPUT_DIRNAME}/.",
    )
    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help="Date stamp for output filenames (default: today's YYYY_MM_DD).",
    )
    return parser.parse_args()


def _sheet_exists(path: Path, sheet_name: str) -> bool:
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return False
    return sheet_name in xl.sheet_names


def resolve_variants_source_path(explicit: Optional[Path] = None) -> Path:
    if explicit is not None:
        path = explicit.resolve()
        if not path.exists():
            raise FileNotFoundError(f"Source file not found: {path}")
        return path

    if DEFAULT_SOURCE_XLSX.exists() and _sheet_exists(DEFAULT_SOURCE_XLSX, DEFAULT_VARIANTS_SHEET):
        return DEFAULT_SOURCE_XLSX.resolve()
    if DEFAULT_SOURCE_XLSX.exists():
        return DEFAULT_SOURCE_XLSX.resolve()
    raise FileNotFoundError(f"VARIANTS source workbook not found: {DEFAULT_SOURCE_XLSX}")


def _format_count(value: int) -> str:
    return f"{value:,}"


def load_flowchart_stats(variants_source: Optional[Path] = None) -> Dict[str, int]:
    """Compute variant and chemotherapy counts used in the Fig. 1 flowchart."""
    source_path = resolve_variants_source_path(variants_source)
    variants = pd.read_excel(source_path, sheet_name=DEFAULT_VARIANTS_SHEET)
    variant_type = variants["VARIANT"].astype(str).str.upper()
    is_cnv = variant_type.isin(["AMPLIFICATION", "LOSS"])
    validation_tumors = set(variants["PCT_ID"].astype(str).str.lower())

    stats = {
        "n_mutations": int((~is_cnv).sum()),
        "n_cnvs": int(is_cnv.sum()),
        "n_variants_total": int(len(variants)),
        "n_validation_tumors": int(variants["PCT_ID"].nunique()),
        "n_chemo_pct_total": 218,
        "n_chemo_with_molecular_profile": 204,
        "variants_source": str(source_path),
    }

    if CHEMO_DATA_CSV.exists():
        chemo = pd.read_csv(CHEMO_DATA_CSV)
        stats["n_chemo_pct_total"] = int(len(chemo))
        missing_profile = chemo["molprofil"].astype(str).str.contains("HI", case=False, na=False)
        in_validation = chemo["Model"].astype(str).str.lower().isin(validation_tumors)
        stats["n_chemo_with_molecular_profile"] = int((in_validation & ~missing_profile).sum())

    return stats


def _box(ax: plt.Axes, x0: float, y0: float, x1: float, y1: float) -> None:
    ax.add_patch(
        Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            linewidth=BOX_LW,
            edgecolor="black",
            facecolor="white",
            zorder=2,
        )
    )


def _dashed_box(ax: plt.Axes, x0: float, y0: float, x1: float, y1: float) -> None:
    ax.add_patch(
        Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            linewidth=BOX_LW,
            edgecolor="black",
            facecolor="white",
            linestyle="--",
            zorder=2,
        )
    )


def _text(
    ax: plt.Axes,
    x: float,
    y: float,
    label: str,
    *,
    ha: str = "center",
    size: float = FS_BODY,
    bold: bool = False,
) -> None:
    ax.text(
        x,
        y,
        label,
        ha=ha,
        va="center",
        fontsize=size,
        fontweight="bold" if bold else "normal",
        fontfamily=FONT_STACK,
        color="black",
        zorder=3,
    )


def _line(ax: plt.Axes, x0: float, y0: float, x1: float, y1: float, lw: float = RULE_LW) -> None:
    ax.plot([x0, x1], [y0, y1], color="black", lw=lw, solid_capstyle="butt", zorder=2)


def _elbow_arrow(ax: plt.Axes, points: Sequence[Tuple[float, float]]) -> None:
    """Draw an orthogonal connector; only the final segment carries the head."""
    for (x0, y0), (x1, y1) in zip(points[:-2], points[1:-1]):
        _line(ax, x0, y0, x1, y1, lw=ARROW_LW)
    (xa, ya), (xb, yb) = points[-2], points[-1]
    ax.add_patch(
        FancyArrowPatch(
            (xa, ya),
            (xb, yb),
            arrowstyle="-|>",
            mutation_scale=ARROW_SCALE,
            linewidth=ARROW_LW,
            color="black",
            shrinkA=0,
            shrinkB=0,
            joinstyle="miter",
            zorder=2,
        )
    )


def _table_header(
    ax: plt.Axes,
    text_left_x: float,
    group_y: float,
    subhead_y: float,
    rule_y: float,
) -> None:
    """Group headers, their underlines, and the excluded/remaining sub-headers."""
    _text(ax, GROUP_T_CENTER, group_y, "Treatments (n)", size=FS_TITLE, bold=True)
    _text(ax, GROUP_U_CENTER, group_y, "Tumors (n)", size=FS_TITLE, bold=True)
    _line(ax, RULE_T_SPAN[0], rule_y, RULE_T_SPAN[1], rule_y)
    _line(ax, RULE_U_SPAN[0], rule_y, RULE_U_SPAN[1], rule_y)

    _text(ax, text_left_x, subhead_y, "Exclusion reason", ha="left", size=FS_TITLE, bold=True)
    _text(ax, COL_T_EXCLUDED, subhead_y, "excluded", size=FS_SUBHEAD)
    _text(ax, COL_T_REMAINING, subhead_y, "remaining", size=FS_SUBHEAD)
    _text(ax, COL_U_EXCLUDED, subhead_y, "excluded", size=FS_SUBHEAD)
    _text(ax, COL_U_REMAINING, subhead_y, "remaining", size=FS_SUBHEAD)


def _table_row(
    ax: plt.Axes,
    reason_lines: Sequence[str],
    line_ys: Sequence[float],
    reason_x: float,
    values: Sequence[str],
    value_y: float,
    *,
    bold_values: Sequence[bool] = (False, False, False, False),
) -> None:
    for label, y in zip(reason_lines, line_ys):
        _text(ax, reason_x, y, label, ha="left")
    columns = (COL_T_EXCLUDED, COL_T_REMAINING, COL_U_EXCLUDED, COL_U_REMAINING)
    for x, value, bold in zip(columns, values, bold_values):
        _text(ax, x, value_y, value, bold=bold)


def plot_flowchart(output_dir: Path, date_stamp: str, stats: Dict[str, int]) -> Tuple[Path, Path]:
    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="white")
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0.0, CANVAS_W)
    ax.set_ylim(CANVAS_H, 0.0)
    ax.axis("off")

    # ---- PCT dataset -----------------------------------------------------
    _box(ax, 271, 45, 549, 106)
    _text(ax, 410, 60, "PCT dataset", bold=True)
    _text(ax, 410, 78, "n=4,758 treatments")
    _text(ax, 410, 95, "n=282 tumors with outcome data")

    # ---- Tumor molecular profiles ---------------------------------------
    _box(ax, 50, 153, 290, 266)
    _text(ax, 59, 172, "Tumor molecular profiles filtered", ha="left", bold=True)
    _text(ax, 59, 190, "for 1,534 cancer genes", ha="left", bold=True)
    _text(ax, 59, 213, "retained:", ha="left")
    _text(ax, 59, 232, BULLET, ha="left")
    _text(ax, 76, 232, f"{_format_count(stats['n_mutations'])} mutations", ha="left")
    _text(ax, 59, 250, BULLET, ha="left")
    _text(ax, 76, 250, f"{_format_count(stats['n_cnvs'])} CNVs", ha="left")

    # ---- Filtered for MTA monotherapies ---------------------------------
    _box(ax, 331, 194, 756, 381)
    _text(ax, 341, 212, "Filtered for MTA monotherapies", ha="left", bold=True)
    _table_header(ax, 341, group_y=242, subhead_y=267, rule_y=254)
    _line(ax, COL_DIVIDER_X, 256, COL_DIVIDER_X, 381)

    _text(ax, 349, 287, BULLET, ha="left")
    _table_row(
        ax,
        ("combination", "therapies"),
        (287, 305),
        365,
        ("1,279", "3,479", "25", "256"),
        296,
    )
    _text(ax, 349, 325, BULLET, ha="left")
    _table_row(
        ax,
        ("chemotherapies",),
        (325,),
        365,
        (_format_count(stats["n_chemo_pct_total"]), "3,261", "0", "256"),
        325,
    )
    _text(ax, 349, 344, BULLET, ha="left")
    _table_row(
        ax,
        ("untreated", "control data"),
        (344, 362),
        365,
        ("226", "3,035", "1", "255"),
        353,
    )

    # ---- Digital Drug Assignment ----------------------------------------
    _box(ax, 50, 375, 290, 460)
    _text(ax, 170, 391, "Digital Drug Assignment (DDA)", bold=True)
    _text(ax, 59, 412, "1. Molecular profiles upload", ha="left")
    _text(ax, 59, 430, "2. Data processing", ha="left")
    _text(ax, 59, 448, "3. DDA scores generated (n=70,460)", ha="left")

    # ---- Filtered for DDA scores ----------------------------------------
    _box(ax, 322, 456, 756, 678)
    _text(ax, 332, 475, "Filtered for DDA scores", ha="left", bold=True)
    _table_header(ax, 332, group_y=507, subhead_y=529, rule_y=519)
    _line(ax, COL_DIVIDER_X, 521, COL_DIVIDER_X, 678)

    _text(ax, 337, 554, BULLET, ha="left")
    _table_row(
        ax,
        ("gastric (no", "molecular data)"),
        (554, 572),
        348,
        ("698", "2,337", "64", "191"),
        563,
    )
    _text(ax, 337, 588, BULLET, ha="left")
    _table_row(
        ax,
        ("others without", "molecular profile"),
        (588, 606),
        348,
        ("139", "2,198", "13", "178"),
        597,
    )
    _text(ax, 337, 625, BULLET, ha="left")
    _table_row(
        ax,
        ("treatments not", "associated to the", "molecular profile"),
        (625, 645, 664),
        348,
        ("1,047", "1,151", "0", "178"),
        645,
        bold_values=(False, True, False, True),
    )

    # ---- Validation set --------------------------------------------------
    validation_x0, validation_y0, validation_x1, validation_y1 = 281.0, 765.0, 525.0, 877.0
    _box(ax, validation_x0, validation_y0, validation_x1, validation_y1)
    _text(ax, 403, 791, "Validation set", bold=True)
    _text(ax, 403, 813, "1151 monotherapy (with DDA scores)")
    _text(ax, 403, 832, "+ treatment outcome data points")
    _text(ax, 403, 851, "(178 tumors, 23 MTAs)")

    # ---- Additional benchmark (chemotherapy reference) -------------------
    benchmark_x0, benchmark_x1 = 538.0, 710.0
    _dashed_box(ax, benchmark_x0, validation_y0, benchmark_x1, validation_y1)
    _text(
        ax,
        (benchmark_x0 + benchmark_x1) / 2.0,
        (validation_y0 + validation_y1) / 2.0,
        (
            "Additional\n"
            "benchmark\n"
            f"{_format_count(stats['n_chemo_with_molecular_profile'])} chemotherapy\n"
            "outcome data points"
        ),
        size=FS_SUBHEAD,
    )

    # ---- Statistical evaluation -----------------------------------------
    _box(ax, 281, 931, 525, 1004)
    _text(ax, 290, 948, "Statistical evaluation of predictivity", ha="left", bold=True)
    _text(ax, 290, 971, BULLET, ha="left")
    _text(ax, 318, 971, "survival", ha="left")
    _text(ax, 290, 990, BULLET, ha="left")
    _text(ax, 318, 990, "response", ha="left")
    _text(ax, 397, 980, "vs. DDA score", ha="left")

    # ---- Connectors ------------------------------------------------------
    # PCT -> molecular profiles (left elbow off the box side)
    _elbow_arrow(ax, [(271, 78), (170, 78), (170, 153)])
    # PCT -> MTA monotherapies (elbow off the box right edge)
    _elbow_arrow(ax, [(549, 78), (600, 78), (600, 150), (550, 150), (550, 194)])
    # molecular profiles -> DDA
    _elbow_arrow(ax, [(170, 266), (170, 375)])
    # DDA -> merge lane feeding the DDA-score filter
    _elbow_arrow(ax, [(290, 413), (533, 413)])
    # MTA monotherapies -> DDA-score filter through the same lane
    _elbow_arrow(ax, [(533, 381), (533, 456)])
    # DDA-score filter -> validation set (stepped)
    _elbow_arrow(ax, [(605, 678), (605, 723), (403, 723), (403, 765)])
    # validation set -> statistical evaluation
    _elbow_arrow(ax, [(403, 877), (403, 931)])

    png_path = output_dir / f"{OUTPUT_BASENAME}_{date_stamp}.png"
    jpeg_path = output_dir / f"{OUTPUT_BASENAME}_{date_stamp}.jpeg"
    fig.savefig(png_path, dpi=300, facecolor="white")
    fig.savefig(jpeg_path, dpi=300, facecolor="white")
    plt.close(fig)
    return png_path, jpeg_path


def main() -> None:
    args = parse_args()
    date_stamp = args.date or datetime.now().strftime("%Y_%m_%d")

    output_root = args.output_dir.resolve() / FLOWCHART_OUTPUT_DIRNAME
    output_root.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Computing flowchart statistics from supplementary data")
    stats = load_flowchart_stats()
    LOGGER.info(
        "Variant counts: %s mutations, %s CNVs (%s tumors); chemo: %s total PCT, %s w/ molecular profile",
        _format_count(stats["n_mutations"]),
        _format_count(stats["n_cnvs"]),
        stats["n_validation_tumors"],
        stats["n_chemo_pct_total"],
        stats["n_chemo_with_molecular_profile"],
    )

    LOGGER.info("Rendering Fig. 1 flowchart with font stack %s", FONT_STACK)
    png_path, jpeg_path = plot_flowchart(output_root, date_stamp, stats)

    metadata: Dict[str, object] = {
        "generated_at": datetime.now().isoformat(),
        "figure_panel": "Fig. 1",
        "description": "Data processing flowchart (PCT filtering to validation cohort)",
        "canvas_px": [CANVAS_W, CANVAS_H],
        "font_stack": FONT_STACK,
        "output_root": str(output_root),
        "flowchart_stats": stats,
        "outputs": {
            "figure_png": str(png_path),
            "figure_jpeg": str(jpeg_path),
        },
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    LOGGER.info("flowchart_final completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
