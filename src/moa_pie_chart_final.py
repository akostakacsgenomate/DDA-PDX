"""Generate Fig. 2c MoA / target pie-of-pie chart.

Primary mechanism-of-action pie with the RTK slice exploded and linked to an
RTK-target sub-pie, count labels on wedges, and legends below.

Usage
-----
python src/moa_pie_chart_final.py

Output is written to src/moa_distribution/.
"""
from __future__ import annotations
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import ConnectionPatch, Patch
from matplotlib.transforms import ScaledTranslation
from plot_typography import apply_arial_font
apply_arial_font()
# ---------------------------------------------------------------------------
# Repository-relative paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
OUTPUT_DIR = _SCRIPT_DIR / "moa_distribution"
OUTPUT_BASENAME = "Fig2c_MoA_piechart"
DATE_STAMP = datetime.now().strftime("%Y_%m_%d")
# Primary pie: signaling mechanism / MoA class (unique MTA targets)
MOA_LABELS: List[str] = [
    "RTK",
    "MAPK",
    "PI3K",
    "CDK",
    "JAK",
    "HR",
    "SMO",
    "MDM2",
    "HSP90",
]
MOA_SIZES: List[int] = [10, 4, 3, 1, 1, 1, 1, 1, 1]
# Soft pastels aligned with pie_chart_final.py (Fig. 2a), one step darker
MOA_COLORS: List[str] = [
    "#D7E8A9",  # RTK (lung light green)
    "#E5A5B9",  # MAPK (breast pink)
    "#9AB8DC",  # PI3K (light colon blue)
    "#8EC4C4",  # CDK (soft teal; not used elsewhere in main pie)
    "#B5CF98",  # JAK (soft green)
    "#E9BB4C",  # HR (unknown gold)
    "#E387A7",  # SMO (custom pink)
    "#EBCE78",  # MDM2 (light yellow)
    "#62839A",  # HSP90 (skin blue-gray)
]
# Secondary pie: RTK inhibitors by target (sums to RTK = 10)
RTK_LABELS: List[str] = [
    "EGFR",
    "FGFR",
    "ERBB2",
    "ERBB3",
    "NTRK",
    "ALK",
    "IGFR",
    "MET",
]
RTK_SIZES: List[int] = [2, 2, 1, 1, 1, 1, 1, 1]
# Sub-pie legend: 2 columns x 4 rows (left-to-right, top-to-bottom)
RTK_LEGEND_ORDER: List[str] = [
    "EGFR",
    "ERBB2",
    "NTRK",
    "IGFR",
    "FGFR",
    "ERBB3",
    "ALK",
    "MET",
]
# RTK sub-pie colors - soft pastels aligned with pie_chart_final.py, one step darker
RTK_COLORS: List[str] = [
    "#CEBFEA",  # EGFR (pancreas lavender)
    "#E5A5B9",  # FGFR (breast pink)
    "#A8D4C8",  # ERBB2 (soft mint; not used in main pie)
    "#62839A",  # ERBB3 (skin blue-gray)
    "#3870C1",  # NTRK (colon blue)
    "#998CE7",  # ALK (purple; swapped with MET)
    "#E387A7",  # IGFR (custom palette pink)
    "#E9BB4C",  # MET (gold; swapped with ALK)
]
FIGURE_SIZE = (12.0, 5.8)
# Typography: +20% vs prior draft, aligned close to Fig. 2a (23 pt wedge / 21 pt legend).
PRIMARY_LABEL_FONT_SIZE = 22
SUB_LABEL_FONT_SIZE = 20
LEGEND_FONT_SIZE = 20
RTK_LEGEND_Y_OFFSET_POINTS = 2.0  # nudge sub-pie key 2 pt downward
RTK_EXPLODE_POINTS = 5.0  # gap between RTK wedge and main pie (typographic points)
AX1_WIDTH_FRAC = 0.48
AX1_XLIM_SPAN = 2.50  # ax1 xlim: -1.20 to +1.30
PRIMARY_RADIUS = 1.0
SECONDARY_RADIUS = 0.92

def _rtk_explode_fraction() -> float:
    """Convert RTK_EXPLODE_POINTS to matplotlib pie explode units."""
    ax_width_in = FIGURE_SIZE[0] * AX1_WIDTH_FRAC
    data_units_per_inch = AX1_XLIM_SPAN / ax_width_in
    offset_data = (RTK_EXPLODE_POINTS / 72.0) * data_units_per_inch
    return offset_data / PRIMARY_RADIUS

def make_count_autopct(sizes: Sequence[int]):
    """Label wedges with absolute counts (black numerals on pastel wedges)."""
    state = {"i": 0}
    def _fmt(_pct: float) -> str:
        idx = state["i"]
        state["i"] = idx + 1
        return str(int(sizes[idx]))
    return _fmt

def ensure_output_directory() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR

def write_count_tables(output_dir: Path, date_stamp: str) -> Path:
    moa = pd.DataFrame({"MoA_class": MOA_LABELS, "n_targets": MOA_SIZES})
    rtk = pd.DataFrame({"RTK_target": RTK_LABELS, "n_targets": RTK_SIZES})
    path = output_dir / f"{OUTPUT_BASENAME}_counts_{date_stamp}.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        moa.to_excel(writer, sheet_name="MoA_primary", index=False)
        rtk.to_excel(writer, sheet_name="RTK_secondary", index=False)
    return path

def _rtk_start_angle_facing_right() -> float:
    """Place the RTK wedge so its bisector points toward the secondary pie (+x)."""
    rtk_span = 360.0 * MOA_SIZES[0] / float(sum(MOA_SIZES))
    # First wedge starts at startangle and grows CCW; bisector at start + span/2.
    # Want bisector near 0 deg (east / toward secondary pie).
    return -0.5 * rtk_span

def _wedge_edge_xy(
    center: Tuple[float, float],
    radius: float,
    theta_deg: float,
) -> Tuple[float, float]:
    theta = np.deg2rad(theta_deg)
    return center[0] + radius * np.cos(theta), center[1] + radius * np.sin(theta)

def plot_moa_pies(output_dir: Path, date_stamp: str) -> Path:
    """Draw an Excel-style pie-of-pie matching the manuscript Fig. 2c layout."""
    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="white")
    # Main pie left; RTK sub-pie upper-right - close but not touching
    ax1 = fig.add_axes([0.02, 0.12, 0.48, 0.82])
    ax2 = fig.add_axes([0.52, 0.40, 0.34, 0.50])
    explode = [_rtk_explode_fraction()] + [0.0] * (len(MOA_SIZES) - 1)
    startangle = _rtk_start_angle_facing_right()
    wedges1, _, autotexts1 = ax1.pie(
        MOA_SIZES,
        colors=MOA_COLORS,
        explode=explode,
        startangle=startangle,
        counterclock=True,
        autopct=make_count_autopct(MOA_SIZES),
        pctdistance=0.55,
        wedgeprops={"linewidth": 1.5, "edgecolor": "white"},
        textprops={"fontsize": PRIMARY_LABEL_FONT_SIZE, "fontweight": "bold", "color": "black"},
    )
    for t in autotexts1:
        t.set_color("black")
        t.set_fontsize(PRIMARY_LABEL_FONT_SIZE)
        t.set_fontweight("bold")
    _, _, autotexts2 = ax2.pie(
        RTK_SIZES,
        colors=RTK_COLORS,
        startangle=90,
        counterclock=True,
        autopct=make_count_autopct(RTK_SIZES),
        pctdistance=0.55,
        radius=SECONDARY_RADIUS,
        wedgeprops={"linewidth": 1.5, "edgecolor": "white", "zorder": 2},
        textprops={"fontsize": SUB_LABEL_FONT_SIZE, "fontweight": "bold", "color": "black"},
    )
    for t in autotexts2:
        t.set_color("black")
        t.set_fontsize(SUB_LABEL_FONT_SIZE)
        t.set_fontweight("bold")
        t.set_zorder(3)
    ax1.set_xlim(-1.20, 1.30)
    ax1.set_ylim(-1.10, 1.15)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax2.set_xlim(-1.05, 1.05)
    ax2.set_ylim(-1.05, 1.05)
    ax2.set_aspect("equal")
    ax2.axis("off")
    # Connecting lines: RTK radial edges -> sub-pie rim at 12 and 6 o'clock
    rtk = wedges1[0]
    theta1, theta2 = float(rtk.theta1), float(rtk.theta2)
    bisector = 0.5 * (theta1 + theta2)
    rtk_explode = _rtk_explode_fraction()
    rtk_center = (
        rtk_explode * PRIMARY_RADIUS * np.cos(np.deg2rad(bisector)),
        rtk_explode * PRIMARY_RADIUS * np.sin(np.deg2rad(bisector)),
    )
    edge_a = _wedge_edge_xy(rtk_center, PRIMARY_RADIUS, theta1)
    edge_b = _wedge_edge_xy(rtk_center, PRIMARY_RADIUS, theta2)
    if edge_a[1] >= edge_b[1]:
        p_12, p_6 = edge_a, edge_b
    else:
        p_12, p_6 = edge_b, edge_a
    q_12 = (0.0, float(SECONDARY_RADIUS))   # 12 o'clock
    q_6 = (0.0, float(-SECONDARY_RADIUS))  # 6 o'clock
    for p, q in ((p_12, q_12), (p_6, q_6)):
        con = ConnectionPatch(
            xyA=p,
            xyB=q,
            coordsA="data",
            coordsB="data",
            axesA=ax1,
            axesB=ax2,
            color="black",
            lw=1.1,
            zorder=5,
            clip_on=False,
        )
        ax2.add_artist(con)
    # Legends attached under each pie (minimal gap)
    handles_moa = [
        Patch(facecolor=c, edgecolor="white", label=lab)
        for c, lab in zip(MOA_COLORS, MOA_LABELS)
    ]
    ax1.legend(
        handles=handles_moa,
        loc="upper center",
        bbox_to_anchor=(0.48, 0.00),
        ncol=3,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        handlelength=1.0,
        handleheight=1.0,
        columnspacing=1.0,
        labelspacing=0.25,
        borderaxespad=0.0,
    )
    rtk_color_map = dict(zip(RTK_LABELS, RTK_COLORS))
    handles_rtk = [
        Patch(facecolor=rtk_color_map[lab], edgecolor="white", label=lab)
        for lab in RTK_LEGEND_ORDER
    ]
    # RTK color key: 2 cols x 4 rows, directly under the sub-pie
    rtk_legend_transform = ax2.transAxes + ScaledTranslation(
        0.0, -RTK_LEGEND_Y_OFFSET_POINTS / 72.0, fig.dpi_scale_trans
    )
    ax2.legend(
        handles=handles_rtk,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.04),
        bbox_transform=rtk_legend_transform,
        ncol=2,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        handlelength=1.0,
        handleheight=1.0,
        columnspacing=1.2,
        labelspacing=0.25,
        borderaxespad=0.0,
    )
    output_path = output_dir / f"{OUTPUT_BASENAME}_{date_stamp}.png"
    fig.savefig(
        output_path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.08
    )
    jpeg_path = output_dir / f"{OUTPUT_BASENAME}_{date_stamp}.jpeg"
    fig.savefig(
        jpeg_path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.08
    )
    plt.close(fig)
    return output_path

def main() -> None:
    if sum(RTK_SIZES) != MOA_SIZES[MOA_LABELS.index("RTK")]:
        raise ValueError(
            f"RTK sub-pie sum ({sum(RTK_SIZES)}) must equal RTK primary count "
            f"({MOA_SIZES[MOA_LABELS.index('RTK')]})"
        )
    date_stamp = DATE_STAMP
    output_dir = ensure_output_directory()
    fig_path = plot_moa_pies(output_dir, date_stamp)
    counts_path = write_count_tables(output_dir, date_stamp)
    metadata: Dict[str, object] = {
        "generated_at": datetime.now().isoformat(),
        "figure_panel": "Fig. 2c",
        "layout": "pie-of-pie (RTK exploded + connected RTK sub-pie)",
        "output_dir": str(output_dir.resolve()),
        "outputs": {
            "figure_png": str(fig_path),
            "figure_jpeg": str(fig_path.with_suffix(".jpeg")),
            "counts_xlsx": str(counts_path),
        },
        "moa_primary": dict(zip(MOA_LABELS, MOA_SIZES)),
        "rtk_secondary": dict(zip(RTK_LABELS, RTK_SIZES)),
    }
    (output_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(f"Saved: {fig_path}")
    print(f"Saved: {fig_path.with_suffix('.jpeg')}")
    print(f"Saved: {counts_path}")
    print(f"Output folder: {output_dir.resolve()}")

if __name__ == "__main__":
    main()
