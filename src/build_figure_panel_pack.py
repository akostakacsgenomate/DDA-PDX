"""Assemble numbered/lettered figure panels into per-figure folders.

Collects the current figure panel images produced by the analysis scripts and
packages them for manuscript submission.

Usage
-----
python build_figure_panel_pack.py [--output-dir DIR]

Outputs are written to <output-dir>/figure_panels/ (default: script directory):

  FIGURE_PANEL_INVENTORY.xlsx  - comprehensive panel list
  FigN/FigN_composite_A4.pdf   - vector composite (A4 layout)
  FigN/FigN_composite_A4.png   - lossless raster composite (300 dpi)
  FigN/FigN_composite_A4.jpeg  - high-quality raster composite (300 dpi, q=95)
  FigN/FigN{letter}_*.jpeg     - individual panel copies
  FigN/FigN{letter}_*.pdf      - individual panel copies

An absolute path pointer is written to build_figure_panel_pack_output_folder.txt.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import ScaledTranslation
from PIL import Image, ImageDraw, ImageFont

from plot_typography import (
    A4_H_IN,
    A4_W_IN,
    GRID_BOTTOM,
    GRID_HSPACE,
    GRID_LEFT,
    GRID_RIGHT,
    GRID_TOP,
    GRID_WSPACE,
    apply_arial_font,
)

apply_arial_font()


# ---------------------------------------------------------------------------
# Repository-relative paths (no hard-coded absolute paths)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent

PANEL_PACK_OUTPUT_DIRNAME = "figure_panels"
INVENTORY_FILENAME = "FIGURE_PANEL_INVENTORY.xlsx"

# A4 size in inches (shared with plot_typography so panel scripts can size fonts)
A4_W, A4_H = A4_W_IN, A4_H_IN
COMPOSITE_PDF_DPI = 250
COMPOSITE_RASTER_DPI = 300

LOGGER = logging.getLogger("build_figure_panel_pack")
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble numbered/lettered figure panels into per-figure folders."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help=f"Parent output directory; results in <output-dir>/{PANEL_PACK_OUTPUT_DIRNAME}/.",
    )
    return parser.parse_args()


@dataclass(frozen=True)
class PanelSpec:
    figure: int
    letter: str
    description: str
    source: Path
    layout_row: int
    layout_col: int
    colspan: int = 1
    rowspan: int = 1

    @property
    def panel_id(self) -> str:
        if self.letter:
            return f"Fig{self.figure}{self.letter}"
        return f"Fig{self.figure}"


def _first_existing(*candidates: Path) -> Optional[Path]:
    for p in candidates:
        if p is not None and p.exists() and p.is_file():
            return p
    return None


def _latest_matching(directory: Path, pattern: str) -> Optional[Path]:
    if not directory.is_dir():
        return None
    matches = sorted(directory.glob(pattern), reverse=True)
    return matches[0] if matches else None


def _scores_figures_dir(src: Path) -> Optional[Path]:
    """Return figures/ from the score_distributions folder."""
    figures_dir = src / "score_distributions" / "figures"
    return figures_dir if figures_dir.is_dir() else None


def _find_in_scores_dirs(src: Path, pattern: str) -> Optional[Path]:
    """Search score_distributions/figures/ for a file."""
    figures_dir = src / "score_distributions" / "figures"
    return _latest_matching(figures_dir, pattern)


def _km_ranking_figures_dir(src: Path) -> Optional[Path]:
    """Return the Kaplan-Meier ranking output folder."""
    km_root = src / "kaplan_meier" / "ranking"
    return km_root if km_root.is_dir() else None


def _dcr_orr_figures_dir(src: Path) -> Optional[Path]:
    """Return figures/ from the dcr_orr_db output folder."""
    figures_dir = src / "dcr_orr_db" / "figures"
    return figures_dir if figures_dir.is_dir() else None


def _latest_clinical_utility_run(src: Path) -> Optional[Path]:
    """Return the clinical_utility_ranking output folder."""
    cur_root = src / "clinical_utility_ranking"
    return cur_root if cur_root.is_dir() else None


def _hr_ranking_forest_plot(src: Path) -> Optional[Path]:
    """Return forest_plot.jpg from the hazard-ratio analysis."""
    candidate = src / "hazard_ratios" / "forest_plot.jpg"
    return candidate if candidate.is_file() else None


def _find_clinical_utility_waterfall_plot(src: Path, *patterns: str) -> Optional[Path]:
    """Find a clinical-utility plot matching any pattern."""
    cur_root = src / "clinical_utility_ranking"
    if not cur_root.is_dir():
        return None
    for pattern in patterns:
        for match in sorted(cur_root.glob(pattern), reverse=True):
            if "tilted" not in match.name.lower():
                return match
    return None


def build_panel_specs() -> List[PanelSpec]:
    src = _SCRIPT_DIR
    fig3 = src / "figure_panels" / "Fig3"
    fig4 = src / "figure_panels" / "Fig4"
    km = _km_ranking_figures_dir(src) or Path()
    dcr = _dcr_orr_figures_dir(src) or Path()
    cu_run = _latest_clinical_utility_run(src)
    pub = cu_run if cu_run is not None else Path()
    wf = cu_run if cu_run is not None else Path()
    dm = src / "data_modeling" / "figures"
    moa = src / "moa_distribution"
    box = src / "variant_counts"
    pie = src / "tumor_site"
    flow = src / "flowchart"

    panels: List[PanelSpec] = []

    # Fig 1 - a: graphical abstract (top); b: data-processing flowchart (below)
    # Composite is sized to half A4 height.
    fig_panels = src / "figure_panels"
    data_dir = src.parent / "data"
    fig1a = _first_existing(
        data_dir / "Fig1a_PDX_Graphical_abstract.jpg",
        fig_panels / "Fig1" / "Fig1a_PDX_Graphical_abstract.jpg",
        fig_panels / "Fig1a_PDX_Graphical_abstract.jpg",
        fig_panels / "1A_PDX Graphical abstract.jpg",
    )
    fig1b_candidates = sorted(
        flow.glob("Fig1_data_processing_flowchart_*.jpeg"), reverse=True
    )
    fig1b = _first_existing(
        *fig1b_candidates,
        _latest_matching(fig_panels / "Fig1", "Fig1_data_processing_flowchart_*.jpeg"),
    )
    if fig1a is not None:
        panels.append(
            PanelSpec(
                1,
                "a",
                "PDX graphical abstract (DDA scoring + PDX validation overview)",
                fig1a,
                0,
                0,
                colspan=1,
                rowspan=1,
            )
        )
    if fig1b is not None:
        panels.append(
            PanelSpec(
                1,
                "b",
                "Data-processing flowchart (PCT filtering to analysis cohort)",
                fig1b,
                1,
                0,
                colspan=1,
                rowspan=1,
            )
        )

    # Fig 2 - a b c / d / e f (sources: latest script outputs under src/)
    fig2_map = [
        (
            "a",
            "Primary tumor site distribution (n=178)",
            _first_existing(
                _latest_matching(pie, "Fig2a_distribution_Primary_Tumor_Site_*.jpeg"),
            ),
            0,
            0,
        ),
        (
            "b",
            "Alteration counts per tumor by Driver / VUS / Non-Driver / All (violin, median)",
            _first_existing(
                _latest_matching(box, "Fig2b_PDX_all_violinplot_classification_*.jpeg"),
                _latest_matching(box, "Fig2b_PDX_all_boxplot_classification_*.jpeg"),
            ),
            0,
            1,
        ),
        (
            "c",
            "MoA / target pie with RTK sub-pie",
            _first_existing(
                _latest_matching(moa, "Fig2c_MoA_piechart_*.jpeg"),
                _latest_matching(moa, "Fig2c_MoA_piechart_*.png"),
            ),
            0,
            2,
        ),
        (
            "d",
            "Individual DDA drug score distribution (compound plot)",
            _first_existing(
                _find_in_scores_dirs(src, "Fig2d_PDX_compound_plot_*.jpg"),
                _find_in_scores_dirs(src, "PDX_compound_plot_*.jpg"),
            ),
            1,
            0,
        ),
        (
            "e",
            "DDA-score waterfall colored by BestResponse",
            _first_existing(
                _find_in_scores_dirs(
                    src,
                    "Fig2e_PDX_DDA_score_Response_Best_Response_all_inone_waterfall_Plot_*.png",
                ),
            ),
            2,
            0,
        ),
        (
            "f",
            "DDA-score waterfall colored by BestAvgResponse",
            _first_existing(
                _find_in_scores_dirs(
                    src,
                    "Fig2f_PDX_DDA_score_BestAvgResponse_all_inone_waterfall_Plot_*.png",
                ),
            ),
            2,
            1,
        ),
    ]
    for letter, desc, path, row, col in fig2_map:
        if path is None:
            continue
        colspan = 3 if letter == "d" else 1
        # e and f share row with 2 cols - use colspan 1 each but grid is 2 for bottom
        if letter in {"e", "f"}:
            # remap to 2-column bottom row: e=0, f=1 in a 2-col sense via colspan later
            pass
        panels.append(
            PanelSpec(2, letter, desc, path, row, col, colspan=colspan, rowspan=1)
        )

    # Fig 3 - a b c / d e / f g / h i
    fig3_map = [
        (
            "a",
            "Best response by DDA score ranks (BestResponse) + Cochran-Armitage",
            _first_existing(
                dcr / "Fig3a_F456_mRECIST_by_rank_ORR_DCR_Response_Best_Response.png",
                fig3 / "Fig3a_F456_mRECIST_by_rank_ORR_DCR_Response_Best_Response.png",
            ),
            0,
            0,
        ),
        (
            "b",
            "Best response by DDA score ranks (BestAvgResponse) + Cochran-Armitage",
            _first_existing(
                dcr / "Fig3b_F456_mRECIST_by_rank_ORR_DCR_Response_BestAvgResponse.png",
                fig3 / "Fig3b_F456_mRECIST_by_rank_ORR_DCR_Response_BestAvgResponse.png",
            ),
            0,
            1,
        ),
        (
            "c",
            "Durable benefit by rank + Cochran-Armitage",
            _first_existing(
                dcr / "Fig3c_F467_Durable_Benefit_by_rank_CA_trend.png",
                fig3 / "Fig3c_F467_Durable_Benefit_by_rank_CA_trend.png",
            ),
            0,
            2,
        ),
        (
            "d",
            "By-rank DDA-score waterfalls (BestResponse)",
            _first_existing(
                _find_in_scores_dirs(src, "Fig3d_F2B_PDX_DDA_waterfall_Best_Response_*.png"),
                fig3 / "Fig3d_F2B_PDX_DDA_waterfall_Best_Response_2026_06_02_.png",
            ),
            1,
            0,
        ),
        (
            "e",
            "By-rank DDA-score waterfalls (BestAvgResponse)",
            _first_existing(
                _find_in_scores_dirs(src, "Fig3e_F2B_PDX_DDA_waterfall_BestAvgResponse_*.png"),
                fig3 / "Fig3e_F2_PDX_DDA_waterfall_BestAvgResponse_2026_06_02_.png",
            ),
            1,
            1,
        ),
        (
            "f",
            "PFS Kaplan-Meier by DDA rank groups",
            _first_existing(
                _latest_matching(km, "all_Kmplot.jpg"),
                fig3 / "Fig3f_all_Kmplot.jpg",
            ),
            2,
            0,
        ),
        (
            "g",
            "Pairwise HR forest plot (rank strata)",
            _first_existing(
                _hr_ranking_forest_plot(src),
                fig3 / "Fig3g_forest_plot.jpg",
            ),
            2,
            1,
        ),
        (
            "h",
            "SCAM vs GAM avg PFS vs DDA rank",
            _first_existing(
                dm / "Fig3h_SCAM_vs_GAM_all_ranks.png",
                fig3 / "Fig3h_SCAM_vs_GAM_all_ranks.png",
            ),
            3,
            0,
        ),
        (
            "i",
            "LMM effect sizes vs top rank",
            _first_existing(
                dm / "Fig3i_F8_Effect_forest_enhanced.png",
                fig3 / "Fig3i_F8_Effect_forest_enhanced.png",
            ),
            3,
            1,
        ),
    ]
    for letter, desc, path, row, col in fig3_map:
        if path is None:
            continue
        panels.append(PanelSpec(3, letter, desc, path, row, col, colspan=1))

    # Fig 4 - a b / c d / e / f g
    fig4_map = [
        (
            "a",
            "KM Top-rank vs SC",
            _first_existing(
                _latest_matching(km, "Fig4a_*Kmplot_1_vs_SC*.jpg"),
                _latest_matching(km, "all_Kmplot_1_vs_SC*.jpg"),
            ),
            0,
            0,
        ),
        (
            "b",
            "KM Bottom ranks 7-10 vs SC",
            _first_existing(
                _latest_matching(km, "Fig4b_*Kmplot_7-10_vs_SC*.jpg"),
                _latest_matching(km, "all_Kmplot_7-10_vs_SC*.jpg"),
            ),
            0,
            1,
        ),
        (
            "c",
            "On-label / Off-label / Experimental composition by DDA rank",
            _first_existing(pub / "Fig4c_pub_rank_utility_composition.png"),
            1,
            0,
        ),
        (
            "d",
            "KM Top (no on-label) vs all SC (unstratified)",
            _first_existing(pub / "Fig4d_pub_km_no_onlabel_rank1_vs_SC.png"),
            1,
            1,
        ),
        (
            "e",
            "Paired PFS boxplot: rank vs matched SC",
            _first_existing(
                _find_clinical_utility_waterfall_plot(
                    src,
                    "Fig4e_*paired_boxplot*.jpg",
                    "clinical_utility_rank_mPFS_vs_SC_paired_boxplot_*.jpg",
                ),
                _latest_matching(fig4, "Fig4e_*paired_boxplot*.jpg"),
                _latest_matching(fig4, "Fig4e_clinical_utility_rank_mPFS_vs_SC_paired_boxplot_*.jpg"),
            ),
            2,
            0,
        ),
        (
            "f",
            "Rank1 off/exp vs on-label lower-rank pool PFS boxplot",
            _first_existing(
                _find_clinical_utility_waterfall_plot(
                    src,
                    "Fig4f_*combined_pool_mpfs_boxplot*.jpg",
                    "clinical_utility_q1_combined_pool_mpfs_boxplot_*.jpg",
                ),
                _latest_matching(fig4, "Fig4f_*combined_pool_mpfs_boxplot*.jpg"),
                _latest_matching(fig4, "Fig4f_clinical_utility_q1_combined_pool_mpfs_boxplot_*.jpg"),
            ),
            3,
            0,
        ),
        (
            "g",
            "KM Non assigned vs Bottom ranks",
            _first_existing(
                _latest_matching(km, "Fig4g_*Kmplot_Non_assigned_vs_7-10*.jpg"),
                _latest_matching(km, "all_Kmplot_Non_assigned_vs_7-10*.jpg"),
            ),
            3,
            1,
        ),
    ]
    for letter, desc, path, row, col in fig4_map:
        if path is None:
            continue
        colspan = 2 if letter == "e" else 1
        panels.append(PanelSpec(4, letter, desc, path, row, col, colspan=colspan))

    return panels


def write_inventory(panels: Sequence[PanelSpec], output_root: Path) -> Path:
    rows = []
    for p in panels:
        rows.append(
            {
                "figure": p.figure,
                "panel": p.letter or "(single)",
                "panel_id": p.panel_id,
                "description": p.description,
                "source_path": str(p.source),
                "source_exists": p.source.exists(),
                "layout_row": p.layout_row,
                "layout_col": p.layout_col,
                "colspan": p.colspan,
                "folder": f"{PANEL_PACK_OUTPUT_DIRNAME}/Fig{p.figure}/",
            }
        )
    df = pd.DataFrame(rows)
    inventory_path = output_root / INVENTORY_FILENAME
    with pd.ExcelWriter(inventory_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="All_panels", index=False)
        for fig in sorted(df["figure"].unique()):
            sub = df[df["figure"] == fig]
            sub.to_excel(writer, sheet_name=f"Fig{fig}", index=False)
        summary = (
            df.groupby("figure")
            .agg(n_panels=("panel_id", "count"), missing=("source_exists", lambda s: int((~s).sum())))
            .reset_index()
        )
        summary.to_excel(writer, sheet_name="Summary", index=False)
    return inventory_path


def _open_rgb(path: Path) -> Image.Image:
    im = Image.open(path)
    if im.mode in ("RGBA", "P"):
        bg = Image.new("RGB", im.size, (255, 255, 255))
        rgba = im.convert("RGBA")
        bg.paste(rgba, mask=rgba.split()[-1])
        return bg
    return im.convert("RGB")


def _save_fig1_three_quarter_a4(
    im: Image.Image,
    *,
    letter: str,
    jpeg_path: Path,
    pdf_path: Path,
    png_path: Path,
    alias_stem: Optional[Path] = None,
) -> None:
    """Export one Fig. 1 panel as a separate exact 3/4-A4 page (210 × 222.75 mm).

    Raster pages are composed at 300 dpi to the exact pixel size; PDF page box is
    set to the same physical size (no bbox_inches='tight' shrinkage).
    """
    page_w_mm, page_h_mm = 210.0, 297.0 * 0.75
    page_w_in, page_h_in = page_w_mm / 25.4, page_h_mm / 25.4
    dpi = COMPOSITE_RASTER_DPI
    canvas_w = int(round(page_w_in * dpi))
    canvas_h = int(round(page_h_in * dpi))

    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    # Keep aspect; scale to the largest fit inside a small margin.
    margin = int(round(0.02 * min(canvas_w, canvas_h)))
    max_w, max_h = canvas_w - 2 * margin, canvas_h - 2 * margin
    src = im.convert("RGB")
    scale = min(max_w / src.width, max_h / src.height)
    new_w = max(1, int(round(src.width * scale)))
    new_h = max(1, int(round(src.height * scale)))
    fitted = src.resize((new_w, new_h), Image.Resampling.LANCZOS)
    x0 = (canvas_w - new_w) // 2
    y0 = (canvas_h - new_h) // 2
    canvas.paste(fitted, (x0, y0))

    # Panel letter in the top-left corner (manuscript style).
    draw = ImageDraw.Draw(canvas)
    tag = f"Fig. 1{letter}" if letter else "Fig. 1"
    font = None
    for candidate in (
        "C:/Windows/Fonts/arialbd.ttf",
        "C:/Windows/Fonts/Arial.ttf",
        "arialbd.ttf",
        "Arial.ttf",
    ):
        try:
            font = ImageFont.truetype(candidate, size=max(28, int(round(dpi * 12 / 72))))
            break
        except OSError:
            continue
    if font is None:
        font = ImageFont.load_default()
    pad = max(6, int(round(dpi * 0.04)))
    bbox = draw.textbbox((0, 0), tag, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    box = [pad - 4, pad - 4, pad + tw + 8, pad + th + 8]
    draw.rectangle(box, fill=(255, 255, 255))
    draw.text((pad, pad), tag, fill=(0, 0, 0), font=font)

    canvas.save(png_path, format="PNG", optimize=True)
    canvas.save(jpeg_path, format="JPEG", quality=95, optimize=True)

    # Vector-friendly single-page PDF with an exact 3/4-A4 mediabox.
    fig = plt.figure(figsize=(page_w_in, page_h_in), facecolor="white")
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.imshow(canvas, interpolation="nearest")
    ax.axis("off")
    fig.savefig(
        pdf_path,
        format="pdf",
        dpi=dpi,
        facecolor="white",
        bbox_inches=None,
        pad_inches=0,
    )
    plt.close(fig)

    if alias_stem is not None:
        for src_p, ext in ((pdf_path, ".pdf"), (png_path, ".png"), (jpeg_path, ".jpeg")):
            alias_stem.with_suffix(ext).write_bytes(src_p.read_bytes())


def save_individual_copies(panel: PanelSpec, dest_dir: Path) -> Tuple[Path, Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{panel.panel_id}_{panel.source.stem}"
    # avoid double Fig prefix in stem if source already starts with panel id
    if panel.source.stem.startswith(panel.panel_id):
        stem = panel.source.stem
    elif panel.letter and panel.source.stem.startswith(f"Fig{panel.figure}_"):
        stem = f"{panel.panel_id}_{panel.source.stem[len(f'Fig{panel.figure}_'):]}"

    jpeg_path = dest_dir / f"{stem}.jpeg"
    pdf_path = dest_dir / f"{stem}.pdf"
    png_path = dest_dir / f"{stem}.png"

    im = _open_rgb(panel.source)

    # Fig. 1A / 1B: separate exact 3/4-A4 pages (not full A4, not the half-A4 composite).
    if panel.figure == 1 and panel.letter in {"a", "b"}:
        alias = dest_dir / f"Fig1{panel.letter}_3quarter_A4"
        _save_fig1_three_quarter_a4(
            im,
            letter=panel.letter,
            jpeg_path=jpeg_path,
            pdf_path=pdf_path,
            png_path=png_path,
            alias_stem=alias,
        )
        return jpeg_path, pdf_path

    im.save(jpeg_path, format="JPEG", quality=95, optimize=True)

    # Single-page PDF of the panel image (fit to A4 width with margin)
    fig = plt.figure(figsize=(A4_W, A4_H))
    ax = fig.add_axes([0.04, 0.04, 0.92, 0.92])
    ax.imshow(im)
    ax.axis("off")
    nice = f"Fig. {panel.figure}{panel.letter}" if panel.letter else f"Fig. {panel.figure}"
    ax.text(
        0.0,
        1.01,
        nice,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        fontfamily="Arial",
        va="bottom",
        ha="left",
    )
    fig.savefig(pdf_path, format="pdf", dpi=200, facecolor="white")
    plt.close(fig)
    return jpeg_path, pdf_path


def _grid_for_figure(fig_num: int) -> Tuple[int, int, Dict[str, Tuple[int, int, int, int]]]:
    """Return nrows, ncols, and panel->(r0,c0,rs,cs) spans for GridSpec."""
    if fig_num == 1:
        # Half-A4 stack: graphical abstract (a) above flowchart (b).
        return 2, 1, {"a": (0, 0, 1, 1), "b": (1, 0, 1, 1)}
    if fig_num == 2:
        # row0: a b c; row1: d spanning 3; row2: e f spanning via 3-col with e at 0-1, f at 2? 
        # Better: 3 cols; d colspan 3; e colspan 1.5 -> use 2 cols on bottom by e at (2,0) cs=1 and f at (2,1) but need 2-col row
        # Use 6-col grid for flexibility
        # a:0-2, b:2-4, c:4-6; d:0-6; e:0-3; f:3-6
        return (
            3,
            6,
            {
                "a": (0, 0, 1, 2),
                "b": (0, 2, 1, 2),
                "c": (0, 4, 1, 2),
                "d": (1, 0, 1, 6),
                "e": (2, 0, 1, 3),
                "f": (2, 3, 1, 3),
            },
        )
    if fig_num == 3:
        # a b c; d e; f g; h i  using 6-col
        return (
            4,
            6,
            {
                "a": (0, 0, 1, 2),
                "b": (0, 2, 1, 2),
                "c": (0, 4, 1, 2),
                "d": (1, 0, 1, 3),
                "e": (1, 3, 1, 3),
                "f": (2, 0, 1, 3),
                "g": (2, 3, 1, 3),
                "h": (3, 0, 1, 3),
                "i": (3, 3, 1, 3),
            },
        )
    if fig_num == 4:
        return (
            4,
            4,
            {
                "a": (0, 0, 1, 2),
                "b": (0, 2, 1, 2),
                "c": (1, 0, 1, 2),
                "d": (1, 2, 1, 2),
                "e": (2, 0, 1, 4),
                "f": (3, 0, 1, 2),
                "g": (3, 2, 1, 2),
            },
        )
    return 1, 1, {}


def _detect_yaxis_x(im: Image.Image) -> int:
    """Find the x-pixel position of the y-axis line in a plot image.

    Scans left-to-right looking for the first column with a continuous run of
    dark pixels spanning >= 12 % of the image height.  If the very first hit is
    at the image edge (likely a plot frame), it is skipped and the next one is
    returned.
    """
    arr = np.array(im.convert("L"))
    h, w = arr.shape
    search_w = int(w * 0.40)
    min_run = int(h * 0.12)
    skip_edge = max(3, int(w * 0.005))

    first_hit: Optional[int] = None
    for x in range(search_w):
        col = arr[:, x] < 160
        if not col.any():
            continue
        padded = np.concatenate([[False], col, [False]])
        d = np.diff(padded.astype(np.int8))
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0]
        if len(starts) > 0 and int((ends - starts).max()) >= min_run:
            if first_hit is None:
                first_hit = x
                if x > skip_edge:
                    return x
            else:
                return x
    return first_hit if first_hit is not None else 0


def _pad_image_left(im: Image.Image, pad_px: int) -> Image.Image:
    """Return a copy of *im* with *pad_px* white pixels prepended on the left."""
    if pad_px <= 0:
        return im
    w, h = im.size
    canvas = Image.new("RGB", (w + pad_px, h), (255, 255, 255))
    canvas.paste(im, (pad_px, 0))
    return canvas


def _equalise_sizes(letters: List[str], images: Dict[str, Image.Image]) -> None:
    """Pad images (in-place) to have identical dimensions (max w x max h)."""
    sizes = [(l, images[l].size) for l in letters if l in images]
    if len(sizes) < 2:
        return
    max_w = max(s[0] for _, s in sizes)
    max_h = max(s[1] for _, s in sizes)
    for letter, (w, h) in sizes:
        if w < max_w or h < max_h:
            canvas = Image.new("RGB", (max_w, max_h), (255, 255, 255))
            canvas.paste(images[letter], (0, 0))
            images[letter] = canvas


def _align_column_yaxes(
    letters: List[str],
    images: Dict[str, Image.Image],
) -> None:
    """Pad images on the left (in-place) so y-axes within a column group align."""
    info: Dict[str, Tuple[int, int]] = {}
    for l in letters:
        im = images.get(l)
        if im is None:
            continue
        info[l] = (_detect_yaxis_x(im), im.size[0])

    if len(info) < 2:
        return

    fracs = {l: yax / w for l, (yax, w) in info.items()}
    target = max(fracs.values())

    for l, (yax, w) in info.items():
        if target - fracs[l] < 0.005:
            continue
        pad = int((target * w - yax) / (1 - target))
        if pad > 2:
            images[l] = _pad_image_left(images[l], pad)
            LOGGER.info("  %s: +%d px left (y-axis align)", l, pad)


def _align_fig4_panels(images: Dict[str, Image.Image], wspace: float) -> None:
    """Align all Fig 4 panels: KM plots on identical canvas, y-axes aligned per column.

    Solves three constraints simultaneously:
    - All KM plots (a, b, d, g) end up with identical pixel dimensions
      so they display at the same visual size in the composite.
    - Y-axes within each column group sit at the same fraction,
      producing a single vertical line per column in the composite.
    - The full-width panel (4e) aligns with the left column group
      without needing separate padding.
    """
    km_letters = ["a", "b", "d", "g"]
    right_col_kms = ["b", "d", "g"]
    hw = 2 + (2 - 1) * wspace          # half-width span in grid units
    fw = 4 + (4 - 1) * wspace          # full-width span in grid units

    # 1) Equalise KM pixel sizes (pad bottom-right to max w x max h)
    _equalise_sizes(km_letters, images)
    base_w, base_h = images[km_letters[0]].size

    # 2) Detect y-axis positions on the equalised KM images
    yax: Dict[str, int] = {}
    for l in km_letters:
        if l in images:
            yax[l] = _detect_yaxis_x(images[l])
    LOGGER.info("  KM y-axis positions (equalised): %s", yax)

    # 3) Determine left-column target fraction, accounting for 4e alignment.
    #    Converting 4e's fraction into the half-column equivalent tells us
    #    the minimum left-column fraction that keeps 4e aligned without
    #    needing any (impossible) right-side padding.
    left_km_frac = yax.get("a", 0) / base_w
    e_equiv_frac = 0.0
    im_e = images.get("e")
    if im_e is not None:
        e_yax = _detect_yaxis_x(im_e)
        e_frac = e_yax / im_e.size[0]
        e_equiv_frac = e_frac * fw / hw
        LOGGER.info("  4e y-axis frac=%.4f -> half-col equiv=%.4f", e_frac, e_equiv_frac)

    left_target_frac = max(left_km_frac, e_equiv_frac)

    # 4) Right-column target = max y-axis pixel among right-column KMs
    right_target = max((yax[l] for l in right_col_kms if l in yax), default=0)

    # 5) Compute shifts: left-column KMs shift to match left_target_frac on base_w
    left_target_px = int(round(left_target_frac * base_w))
    shifts: Dict[str, int] = {}
    shifts["a"] = left_target_px - yax.get("a", 0)
    for l in right_col_kms:
        if l in yax:
            shifts[l] = right_target - yax[l]

    # 6) Common canvas width = base + largest shift (keeps all content visible)
    max_shift = max(shifts.values(), default=0)
    canvas_w = base_w + max_shift

    # 7) Place every KM on the same-size canvas
    for l in km_letters:
        if l not in images:
            continue
        shift = shifts.get(l, 0)
        new = Image.new("RGB", (canvas_w, base_h), (255, 255, 255))
        new.paste(images[l], (shift, 0))
        images[l] = new
        LOGGER.info("  KM %s: shift=%d px, canvas %dx%d", l, shift, canvas_w, base_h)

    # 8) Compute the final left-column fraction on the canvas
    left_frac = (yax.get("a", 0) + shifts.get("a", 0)) / canvas_w
    LOGGER.info("  left-col target frac: %.4f", left_frac)

    # 9) Align non-KM left-column panels (c, f) to match left_frac
    for l in ["c", "f"]:
        im = images.get(l)
        if im is None:
            continue
        cur_yax = _detect_yaxis_x(im)
        cur_frac = cur_yax / im.size[0]
        if left_frac - cur_frac < 0.003:
            continue
        w = im.size[0]
        pad = int(round((left_frac * w - cur_yax) / (1 - left_frac)))
        if pad > 2:
            images[l] = _pad_image_left(images[l], pad)
            LOGGER.info("  %s: +%d px left (align to left frac %.4f)", l, pad, left_frac)

    # 10) Align 4e (full-width row) with the left column
    if im_e is not None:
        target_e_frac = left_frac * hw / fw
        e_yax = _detect_yaxis_x(images["e"])
        e_cur_frac = e_yax / images["e"].size[0]
        delta = target_e_frac - e_cur_frac
        if abs(delta) > 0.003:
            w = images["e"].size[0]
            if delta > 0:
                pad = int(round((target_e_frac * w - e_yax) / (1 - target_e_frac)))
                if pad > 2:
                    images["e"] = _pad_image_left(images["e"], pad)
                    LOGGER.info("  e: +%d px left (align frac %.4f)", pad, target_e_frac)
            else:
                pad_right = int(round(e_yax / target_e_frac - w))
                if pad_right > 2:
                    ew, eh = images["e"].size
                    canvas = Image.new("RGB", (ew + pad_right, eh), (255, 255, 255))
                    canvas.paste(images["e"], (0, 0))
                    images["e"] = canvas
                    LOGGER.info("  e: +%d px right (align frac %.4f)", pad_right, target_e_frac)
        else:
            LOGGER.info("  e: already aligned (frac %.4f, target %.4f)", e_cur_frac, target_e_frac)


def _compute_height_ratios(
    nrows: int,
    ncols: int,
    spans: Dict[str, Tuple[int, int, int, int]],
    images: Dict[str, Image.Image],
    wspace: float,
) -> List[float]:
    """Row height ratios proportional to the tallest image in each row.

    For each row, finds the panel whose image needs the most vertical space
    (given its column span and aspect ratio) and uses that as the row height.
    This prevents oversized rows that would leave images floating mid-cell.
    """
    ratios: List[float] = []
    for r in range(nrows):
        max_h = 0.0
        for letter, (r0, c0, rs, cs) in spans.items():
            if r0 != r or rs != 1:
                continue
            im = images.get(letter)
            if im is None:
                continue
            w, h = im.size
            eff_w = cs + (cs - 1) * wspace
            max_h = max(max_h, (h / w) * eff_w)
        ratios.append(max_h if max_h > 0 else 1.0)
    return ratios


def build_composite_a4(
    fig_num: int, panels: Sequence[PanelSpec], dest_pdf: Path
) -> Tuple[Path, Path, Path]:
    nrows, ncols, spans = _grid_for_figure(fig_num)
    by_letter = {p.letter: p for p in panels}

    # Pre-load panel images so we can size rows to their content.
    images: Dict[str, Image.Image] = {}
    for letter in spans:
        panel = by_letter.get(letter)
        if panel is None and letter == "" and panels:
            panel = panels[0]
        if panel is not None and panel.source.exists():
            images[letter] = _open_rgb(panel.source)

    # --- Figure-specific alignment ---
    if fig_num == 3:
        _align_column_yaxes(["d", "f", "h"], images)
    if fig_num == 4:
        _align_fig4_panels(images, GRID_WSPACE)

    height_ratios = _compute_height_ratios(nrows, ncols, spans, images, GRID_WSPACE)
    if fig_num == 1:
        # Equal share of the half-A4 stack for graphical abstract + flowchart.
        height_ratios = [1.0, 1.0]

    # Fig. 1 is a two-panel stack designed for the top half of an A4 page.
    fig_h = (A4_H * 0.5) if fig_num == 1 else A4_H
    fig = plt.figure(figsize=(A4_W, fig_h), facecolor="white")

    # Slightly tighter vertical gap for the half-page Fig. 1 stack.
    hspace = 0.06 if fig_num == 1 else GRID_HSPACE
    top = 0.96 if fig_num == 1 else GRID_TOP
    bottom = 0.04 if fig_num == 1 else GRID_BOTTOM

    gs = fig.add_gridspec(
        nrows,
        ncols,
        left=GRID_LEFT,
        right=GRID_RIGHT,
        top=top,
        bottom=bottom,
        wspace=GRID_WSPACE,
        hspace=hspace,
        height_ratios=height_ratios,
    )

    for letter, (r0, c0, rs, cs) in spans.items():
        im = images.get(letter)
        if im is None:
            continue
        ax = fig.add_subplot(gs[r0 : r0 + rs, c0 : c0 + cs])
        ax.imshow(im)
        ax.axis("off")
        ax.set_anchor("NW")

        # Fine-tune panel horizontal placement (positive = left, negative = right)
        fig3_shift_pt = {"e": -4.0}
        fig4_shift_pt = {"c": 1.0, "e": 5.0, "f": -1.0, "g": 12.0}

        shift_pt = 0.0
        if fig_num == 3 and letter in fig3_shift_pt:
            shift_pt = fig3_shift_pt[letter]
        elif fig_num == 4 and letter in fig4_shift_pt:
            shift_pt = fig4_shift_pt[letter]

        if shift_pt:
            pos = ax.get_position()
            shift_frac = shift_pt / 72.0 / A4_W
            ax.set_position([pos.x0 - shift_frac, pos.y0, pos.width, pos.height])

        panel = by_letter.get(letter)
        if panel is None and letter == "" and panels:
            panel = panels[0]
        tag = panel.letter.upper() if panel and panel.letter else ""
        if tag:
            lift_pt = 2.0 if panel.letter in {"d", "e"} else 0.0
            text_transform = ax.transAxes + ScaledTranslation(
                0.0, lift_pt / 72.0, fig.dpi_scale_trans
            )
            ax.text(
                0.01,
                0.99,
                tag,
                transform=text_transform,
                fontsize=9,
                fontweight="bold",
                fontfamily="Arial",
                va="top",
                ha="left",
                color="black",
                clip_on=False,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85),
            )

    dest_pdf.parent.mkdir(parents=True, exist_ok=True)
    dest_png = dest_pdf.with_suffix(".png")
    dest_jpeg = dest_pdf.with_suffix(".jpeg")
    fig.savefig(dest_pdf, format="pdf", dpi=COMPOSITE_PDF_DPI, facecolor="white")
    fig.savefig(
        dest_png,
        format="png",
        dpi=COMPOSITE_RASTER_DPI,
        facecolor="white",
        pil_kwargs={"compress_level": 3},
    )
    _open_rgb(dest_png).save(dest_jpeg, format="JPEG", quality=95, optimize=True)
    plt.close(fig)
    return dest_pdf, dest_png, dest_jpeg


def write_pack_readme(output_root: Path) -> Path:
    readme = output_root / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Figure panel pack",
                "",
                "Per-figure folders with:",
                "",
                "- `FigN_composite_A4.pdf` - vector composite on one A4 page",
                "- `FigN_composite_A4.png` - lossless raster composite (300 dpi)",
                "- `FigN_composite_A4.jpeg` - high-quality raster composite (300 dpi, q=95)",
                "- `FigN{letter}_*.jpeg` / `.pdf` - individual panel exports",
                "- `EFig2/` - extended Fig 2 (leave-one-out SCAM vs GAM panels a–d)",
                "",
                f"Master inventory: [`{INVENTORY_FILENAME}`]({INVENTORY_FILENAME})",
                "",
                "Regenerate:",
                "",
                "```bash",
                "python src/build_figure_panel_pack.py",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return readme


def _latest_lou_scam_gam_dir(src: Path) -> Optional[Path]:
    """Leave-one-out SCAM vs GAM collect folder."""
    cand = src / "data_modeling" / "leave_one_out" / "SCAM_vs_GAM_figures"
    return cand if cand.is_dir() else None


def _crop_plot_whitespace(im: Image.Image, *, pad: int = 8, threshold: int = 250) -> Image.Image:
    """Trim near-white margins so composite rows can sit close together."""
    arr = np.asarray(im.convert("L"))
    mask = arr < threshold
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    left = max(0, int(xs.min()) - pad)
    top = max(0, int(ys.min()) - pad)
    right = min(im.size[0], int(xs.max()) + 1 + pad)
    bottom = min(im.size[1], int(ys.max()) + 1 + pad)
    return im.crop((left, top, right, bottom))


def _collapse_white_bands(
    im: Image.Image,
    *,
    max_band: int = 24,
    threshold: int = 250,
    empty_frac: float = 0.985,
    top_frac: float = 0.28,
    bottom_frac: float = 0.18,
    collapse_top: bool = False,
) -> Image.Image:
    """Shrink tall near-empty bands near the bottom (and optionally top).

    Top collapse is off by default so ggplot titles keep clear air above the
    plot panel and do not sit on the grid / n= labels.
    """
    arr = np.asarray(im.convert("L"))
    empty = (arr >= threshold).mean(axis=1) >= empty_frac
    keep = np.ones(arr.shape[0], dtype=bool)
    n = arr.shape[0]
    top_lim = int(round(n * top_frac))
    bot_lim = int(round(n * (1.0 - bottom_frac)))
    i = 0
    while i < n:
        if not empty[i]:
            i += 1
            continue
        j = i + 1
        while j < n and empty[j]:
            j += 1
        run = j - i
        in_top = j <= top_lim
        in_bot = i >= bot_lim
        if in_top and not collapse_top:
            i = j
            continue
        if (in_top or in_bot) and run > max_band:
            keep_start = i + (run - max_band) // 2
            keep_end = keep_start + max_band
            keep[i:keep_start] = False
            keep[keep_end:j] = False
        i = j
    if keep.all():
        return im
    rgb = np.asarray(im)
    return Image.fromarray(rgb[keep], mode=im.mode)


def _pad_top_for_letter(im: Image.Image, pad_px: int) -> Image.Image:
    """White strip above the plot so panel letters sit clear of axes/titles."""
    if pad_px <= 0:
        return im
    out = Image.new("RGB", (im.size[0], im.size[1] + pad_px), (255, 255, 255))
    out.paste(im, (0, pad_px))
    return out


def _title_ink_end(gray: np.ndarray) -> Optional[int]:
    """Bottom row of the centered ggplot title (first contiguous top ink run)."""
    h, w = gray.shape
    center = gray[: int(h * 0.14), int(w * 0.28) : int(w * 0.72)]
    ink = (center < 200).mean(axis=1) > 0.02
    if not ink.any():
        return None
    start = int(np.argmax(ink))
    end = start
    while end + 1 < ink.size and ink[end + 1]:
        end += 1
    return end


def _clear_orphaned_nlabel_tops(im: Image.Image) -> Image.Image:
    """Remove clipped ``)`` / dash fragments left above high n= labels.

    Rotated ``(n = …)`` labels on tall points can extend into the title band.
    ``_insert_space_below_title`` then splits those tips from the rest of the
    string, leaving floating ``)`` and ``=`` strokes under the title. Only
    clear short ink in a shallow band just under the title so real n= glyphs
    further down are left intact.
    """
    arr = np.array(im, copy=True)
    gray = np.asarray(im.convert("L"))
    h, w = gray.shape
    title_end = _title_ink_end(gray)
    if title_end is None:
        return im
    # Shallow band only: tips that sit in the title/margin split zone.
    y0 = title_end + 1
    y1 = min(h, title_end + 55)
    x_zones = (
        (int(w * 0.10), int(w * 0.32)),
        (int(w * 0.32), int(w * 0.68)),
    )
    for x0, x1 in x_zones:
        zone = gray[y0:y1, x0:x1] < 110
        if not zone.any():
            continue
        # White out all dark ink in this shallow band (orphans / title crumbs).
        # Real n= bodies sit lower, below the panel top.
        ys, xs = np.where(zone)
        arr[y0 + ys, x0 + xs, :] = 255
    return Image.fromarray(arr, mode=im.mode)


def _insert_space_below_title(im: Image.Image, extra_px: int = 160) -> Image.Image:
    """Push plot content down so the title does not sit on n= labels / grid."""
    if extra_px <= 0:
        return im
    im = _clear_orphaned_nlabel_tops(im)
    arr = np.asarray(im)
    gray = np.asarray(im.convert("L"))
    h, w = gray.shape
    end = _title_ink_end(gray)
    if end is None:
        return im
    # Keep split tight under the title so we do not bisect residual n= ink.
    split = min(h - 1, end + 4)
    # Full-width wipe of the few pad rows under the title (catches center dashes).
    arr = np.array(arr, copy=True)
    arr[end + 1 : split + 1, :, :] = 255
    white = np.full((extra_px, w, arr.shape[2]), 255, dtype=arr.dtype)
    out = np.concatenate([arr[: split + 1], white, arr[split + 1 :]], axis=0)
    return Image.fromarray(out, mode=im.mode)


def _arial_font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    names = (
        ("arialbd.ttf", "Arial Bold.ttf", "DejaVuSans-Bold.ttf")
        if bold
        else ("arial.ttf", "Arial.ttf", "DejaVuSans.ttf")
    )
    win_fonts = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    candidates = []
    for name in names:
        candidates.append(name)
        candidates.append(str(win_fonts / name))
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _relabel_plot_title(
    im: Image.Image,
    new_title: str,
    *,
    font_size: Optional[int] = None,
) -> Image.Image:
    """Replace the centered ggplot title (top of LOU SCAM/GAM panels).

    ``font_size`` defaults to the native ggplot title size measured on these
    3900×3000 exports (~121 px Arial Bold), so overlays match All / CDK4/6 /
    panPI3K.
    """
    out = im.copy()
    w, h = out.size
    # Title band measured on 3900x3000 LOU exports (~y 125–211).
    y0 = max(0, int(round(h * 0.037)))
    y1 = min(h, int(round(h * 0.075)))
    if font_size is None:
        font_size = max(28, int(round((y1 - y0) * 1.15)))  # ~121 on stock exports
    font = _arial_font(font_size, bold=True)
    # Full-width clear so leftover title ink / clipped n= tips cannot remain
    # under the new title (visible as floating dashes or ')').
    draw = ImageDraw.Draw(out)
    draw.rectangle([0, y0, w, y1], fill=(255, 255, 255))
    bbox = draw.textbbox((0, 0), new_title, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    tx = (w - tw) / 2.0
    ty = y0 + ((y1 - y0) - th) / 2.0 - bbox[1]
    draw.text((tx, ty), new_title, fill=(0, 0, 0), font=font)
    return out


# Native ggplot title size on LOU 3900×3000 SCAM/GAM exports (Arial Bold px).
_EFIG2_TITLE_FONT_PX = 121
# White header above each panel for a–d letters. Sized so a ~15 pt letter on the
# A4-width composite fits entirely above the ggplot title (not on axes/plot).
_EFIG2_LETTER_HEADER_PX = 240
_EFIG2_PANEL_TITLES = {
    "a": "All",
    "b": "CDK4/6",
    "c": "MAP2K1/2",
    "d": "panPI3K",
}


def build_efig2_pack(output_root: Path) -> Optional[Dict[str, object]]:
    """Assemble Extended Fig 2: LOU SCAM vs GAM (All, CDK4/6, MAP2K1/2, panPI3K)."""
    lou = _latest_lou_scam_gam_dir(_SCRIPT_DIR)
    if lou is None:
        LOGGER.warning("EFig2 skipped: leave-one-out SCAM/GAM folder not found.")
        return None

    sources = {
        "a": lou / "EFig2a_SCAM_vs_GAM_leave_out_baseline_ALL.png",
        "b": lou / "EFig2b_SCAM_vs_GAM_leave_out_CDK4_6.png",
        "c": lou / "EFig2c_SCAM_vs_GAM_leave_out_MAP2K1_2.png",
        "d": lou / "EFig2d_SCAM_vs_GAM_leave_out_panPI3K.png",
    }
    missing = [f"EFig2{k}" for k, p in sources.items() if not p.is_file()]
    if missing:
        LOGGER.warning("EFig2 skipped; missing panels: %s", ", ".join(missing))
        return None

    dest_dir = output_root / "EFig2"
    dest_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("EFig2 -> %s", dest_dir)

    images: Dict[str, Image.Image] = {}
    panel_outputs: Dict[str, Dict[str, str]] = {}
    for letter, src in sources.items():
        im = _open_rgb(src)
        # Uniform title size; keep title–plot gap; reserve header for letters.
        im = _relabel_plot_title(
            im, _EFIG2_PANEL_TITLES[letter], font_size=_EFIG2_TITLE_FONT_PX
        )
        im = _collapse_white_bands(_crop_plot_whitespace(im), collapse_top=False)
        im = _insert_space_below_title(im, extra_px=200)
        im = _pad_top_for_letter(im, _EFIG2_LETTER_HEADER_PX)
        images[letter] = im

        stem = f"EFig2{letter}"
        png_path = dest_dir / f"{stem}.png"
        jpeg_path = dest_dir / f"{stem}.jpeg"
        pdf_path = dest_dir / f"{stem}.pdf"
        im.save(png_path, format="PNG", compress_level=3)
        im.save(jpeg_path, format="JPEG", quality=95, optimize=True)

        fig = plt.figure(figsize=(A4_W, A4_H), facecolor="white")
        ax = fig.add_axes([0.04, 0.04, 0.92, 0.92])
        ax.imshow(im)
        ax.axis("off")
        ax.text(
            0.0,
            1.01,
            f"EFig. 2{letter}",
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            fontfamily="Arial",
            va="bottom",
            ha="left",
        )
        fig.savefig(pdf_path, format="pdf", dpi=200, facecolor="white")
        plt.close(fig)
        LOGGER.info("  %s: %s | %s | %s", stem, png_path.name, jpeg_path.name, pdf_path.name)
        panel_outputs[stem] = {
            "source": str(src),
            "png": str(png_path),
            "jpeg": str(jpeg_path),
            "pdf": str(pdf_path),
        }

    # Compact 2x2 collage: modest column gap, 100 px between rows.
    order = (("a", "b"), ("c", "d"))
    max_w = max(im.size[0] for im in images.values())
    max_h = max(im.size[1] for im in images.values())
    col_gap = 24
    row_gap = 100
    outer = 6
    canvas_w = outer * 2 + max_w * 2 + col_gap
    canvas_h = outer * 2 + max_h * 2 + row_gap
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    for r, row_letters in enumerate(order):
        for c, letter in enumerate(row_letters):
            im = images[letter]
            if im.size != (max_w, max_h):
                # Scale to a shared cell size so residual crop differences
                # do not open an artificial gap between rows.
                im = im.resize((max_w, max_h), Image.Resampling.LANCZOS)
            x = outer + c * (max_w + col_gap)
            y = outer + r * (max_h + row_gap)
            canvas.paste(im, (x, y))

    composite_pdf = dest_dir / "EFig2_composite_A4.pdf"
    composite_png = dest_dir / "EFig2_composite_A4.png"
    composite_jpeg = dest_dir / "EFig2_composite_A4.jpeg"
    fig_w = A4_W
    fig_h = A4_W * canvas_h / canvas_w
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(canvas)
    ax.axis("off")
    # Panel letters sit in the reserved white header above each plot (not on axes).
    panel_letter_pt = 15.0
    for r, row_letters in enumerate(order):
        for c, letter in enumerate(row_letters):
            x0 = outer + c * (max_w + col_gap)
            y0 = outer + r * (max_h + row_gap)
            # Anchor near the top of the header strip; header height is chosen so
            # the full glyph stays above the panel title.
            lx = (x0 + 0.012 * max_w) / canvas_w
            ly = 1.0 - (y0 + 10) / canvas_h
            ax.text(
                lx,
                ly,
                letter.upper(),
                transform=ax.transAxes,
                fontsize=panel_letter_pt,
                fontweight="bold",
                fontfamily="Arial",
                va="top",
                ha="left",
                color="black",
                clip_on=False,
                bbox=dict(
                    boxstyle="square,pad=0.12",
                    facecolor="white",
                    edgecolor="none",
                    alpha=1.0,
                ),
            )
    fig.savefig(composite_pdf, format="pdf", dpi=COMPOSITE_PDF_DPI, facecolor="white")
    fig.savefig(
        composite_png,
        format="png",
        dpi=COMPOSITE_RASTER_DPI,
        facecolor="white",
        pil_kwargs={"compress_level": 3},
    )
    _open_rgb(composite_png).save(composite_jpeg, format="JPEG", quality=95, optimize=True)
    plt.close(fig)

    # Mirror composite at figure_panels root (same pattern as Fig1–Fig4).
    root_jpeg = output_root / "EFig2_composite_A4.jpeg"
    root_jpeg.write_bytes(composite_jpeg.read_bytes())
    LOGGER.info(
        "  composite: %s | %s | %s",
        composite_pdf.name,
        composite_png.name,
        composite_jpeg.name,
    )
    return {
        "folder": str(dest_dir),
        "composite_pdf": str(composite_pdf),
        "composite_png": str(composite_png),
        "composite_jpeg": str(composite_jpeg),
        "root_composite_jpeg": str(root_jpeg),
        "panels": panel_outputs,
    }


def main() -> None:
    args = parse_args()
    output_root = args.output_dir.resolve() / PANEL_PACK_OUTPUT_DIRNAME
    output_root.mkdir(parents=True, exist_ok=True)

    panels = build_panel_specs()
    missing = [p.panel_id for p in panels if not p.source.exists()]
    if missing:
        LOGGER.warning("Panels with missing source images: %s", ", ".join(missing))

    inventory_path = write_inventory(panels, output_root)
    LOGGER.info("Wrote inventory: %s", inventory_path)

    by_fig: Dict[int, List[PanelSpec]] = {}
    for p in panels:
        by_fig.setdefault(p.figure, []).append(p)

    outputs: Dict[str, object] = {"inventory": str(inventory_path), "figures": {}}
    for fig_num, fig_panels in sorted(by_fig.items()):
        dest_dir = output_root / f"Fig{fig_num}"
        dest_dir.mkdir(parents=True, exist_ok=True)
        LOGGER.info("Fig %d -> %s", fig_num, dest_dir)

        panel_outputs: Dict[str, Dict[str, str]] = {}
        for p in sorted(fig_panels, key=lambda x: x.letter or " "):
            jpeg_p, pdf_p = save_individual_copies(p, dest_dir)
            LOGGER.info("  %s: %s | %s", p.panel_id, jpeg_p.name, pdf_p.name)
            panel_outputs[p.panel_id] = {
                "source": str(p.source),
                "jpeg": str(jpeg_p),
                "pdf": str(pdf_p),
            }

        composite_pdf = dest_dir / f"Fig{fig_num}_composite_A4.pdf"
        composite_pdf, composite_png, composite_jpeg = build_composite_a4(
            fig_num, fig_panels, composite_pdf
        )
        LOGGER.info(
            "  composite: %s | %s | %s",
            composite_pdf.name,
            composite_png.name,
            composite_jpeg.name,
        )

        root_jpeg = output_root / f"Fig{fig_num}_composite_A4.jpeg"
        root_jpeg.write_bytes(composite_jpeg.read_bytes())
        LOGGER.info("  root jpeg: %s", root_jpeg.name)

        outputs["figures"][f"Fig{fig_num}"] = {
            "folder": str(dest_dir),
            "composite_pdf": str(composite_pdf),
            "composite_png": str(composite_png),
            "composite_jpeg": str(composite_jpeg),
            "root_composite_jpeg": str(root_jpeg),
            "panels": panel_outputs,
        }

    efig2 = build_efig2_pack(output_root)
    if efig2 is not None:
        outputs["figures"]["EFig2"] = efig2

    readme = write_pack_readme(output_root)
    outputs["readme"] = str(readme)

    metadata: Dict[str, object] = {
        "generated_at": datetime.now().isoformat(),
        "source_root": str(_SCRIPT_DIR),
        "output_root": str(output_root),
        "outputs": outputs,
        "n_figures": len(by_fig),
        "n_panels": len(panels),
        "missing_panels": missing,
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    LOGGER.info("build_figure_panel_pack completed. Output: %s", output_root)


if __name__ == "__main__":
    main()
