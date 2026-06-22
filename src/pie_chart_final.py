"""Generate the PDX primary tumour site distribution pie chart.

Counts unique PDX models by primary tumour site after excluding
chemotherapy and tubulin-inhibitor reference treatments. The figure
displays model counts and percentages for each tumour site category.

Usage
-----
Run from the repository root:

python src/pie_chart_final.py

Output is written to a script-local output directory (created on first run).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch


# ---------------------------------------------------------------------------
# Repository-relative paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

SOURCE_FILE = DATA_DIR / "pdx_curve_metrics_single_treatments_dda_scores.xlsx"
OUTPUT_DIR = _SCRIPT_DIR / "Primary_Tumor_Site_distribution"

DATE_STAMP = datetime.now().strftime("%Y_%m_%d")
OUTPUT_BASENAME = "distribution_Primary_Tumor_Site"

REQUIRED_COLUMNS = [
    "Model",
    "Primary Tumor Site",
    "Treatment target",
]

# Display order and colours for tumour site categories (manuscript figure).
SITE_ORDER: List[str] = ["lung", "skin", "unknown", "pancreas", "colon", "breast"]

SITE_COLORS: Dict[str, str] = {
    "breast": "#F692B6",
    "colon": "#0F52BA",
    "pancreas": "#D1B9FD",
    "lung": "#E2FC9A",
    "skin": "#476F8A",
    "unknown": "#FDBA1A",
}

# Wedge labels rendered in black for contrast on light-coloured segments.
BLACK_LABEL_SITES = {"unknown", "lung", "breast", "pancreas"}

FIGURE_SIZE = (13, 8)
LABEL_FONT_SIZE = 19
LEGEND_FONT_SIZE = 17
PIE_START_ANGLE = -13
WEDGE_EXPLODE = 0.02


def harmonize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Align column names with the shared PDX analysis schema."""
    out = df.copy()
    if "Treatment target" not in out.columns and "TARGET" in out.columns:
        out = out.rename(columns={"TARGET": "Treatment target"})
    return out


def load_and_preprocess_data(source_xlsx: Path) -> pd.DataFrame:
    """Load treatment data and retain non-reference targeted therapy rows."""
    df = harmonize_column_names(pd.read_excel(source_xlsx))

    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            f"Input workbook is missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    df = df[
        ~df["Treatment target"].str.contains(
            r"\bchemotherapy\b", case=False, na=False
        )
    ].copy()
    df = df[
        ~df["Treatment target"].str.contains(r"\bTubulin\b", case=False, na=False)
    ].copy()

    df["Primary Tumor Site"] = df["Primary Tumor Site"].replace({"UNKNOWN": "unknown"})
    return df


def compute_site_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """Count unique models per primary tumour site and compute percentages."""
    grouped = (
        df.groupby("Primary Tumor Site")["Model"]
        .nunique()
        .reset_index()
        .rename(columns={"Model": "n_models"})
    )
    total_models = grouped["n_models"].sum()
    grouped["percentage"] = (grouped["n_models"] / total_models) * 100
    return grouped


def order_sites(
    labels: Sequence[str], sizes: Sequence[int]
) -> Tuple[List[str], List[int]]:
    """Sort site labels and counts according to the manuscript display order."""
    label_list = list(labels)
    size_list = list(sizes)
    sorted_indices = [SITE_ORDER.index(label) for label in label_list]
    return (
        [label_list[i] for i in sorted_indices],
        [size_list[i] for i in sorted_indices],
    )


def make_autopct_formatter(
    sizes: Sequence[int], labels: Sequence[str]
) -> Callable[[float], str]:
    """Return a pie-chart label formatter (percentage and model count per wedge)."""
    wedge_index = {"value": 0}

    def autopct(pct: float) -> str:
        idx = wedge_index["value"]
        wedge_index["value"] += 1
        absolute = int(round(pct / 100.0 * sum(sizes)))
        if labels[idx].lower() == "unknown":
            return ""
        return f"{pct:.1f}%\n(n={absolute})"

    return autopct


def ensure_output_directory() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


def plot_primary_tumor_site_distribution(
    grouped_data: pd.DataFrame, output_dir: Path
) -> Path:
    """Render and save the primary tumour site distribution pie chart."""
    labels, sizes = order_sites(
        grouped_data["Primary Tumor Site"].tolist(),
        grouped_data["n_models"].tolist(),
    )
    colors = [SITE_COLORS.get(site, "#808080") for site in labels]
    explode = [WEDGE_EXPLODE] * len(labels)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    _, _, autotexts = ax.pie(
        sizes,
        autopct=make_autopct_formatter(sizes, labels),
        startangle=PIE_START_ANGLE,
        pctdistance=0.7,
        explode=explode,
        colors=colors,
        labeldistance=0.8,
        textprops={"fontsize": LABEL_FONT_SIZE, "weight": "bold"},
    )

    # Annotate the single-model unknown category outside the wedge label.
    unknown_row = grouped_data[grouped_data["Primary Tumor Site"] == "unknown"]
    if not unknown_row.empty:
        unknown_pct = unknown_row["percentage"].iloc[0]
        unknown_n = int(unknown_row["n_models"].iloc[0])
        ax.text(
            1.2,
            0.8,
            f"{unknown_pct:.1f}%\n(n={unknown_n})",
            fontsize=LABEL_FONT_SIZE,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.5"),
            transform=ax.transAxes,
        )

    for label, autotext in zip(labels, autotexts):
        autotext.set_color("black" if label in BLACK_LABEL_SITES else "white")

    handles = [Patch(color=color) for color in colors[::-1]]
    ax.legend(
        handles=handles,
        labels=labels[::-1],
        bbox_to_anchor=(0.90, 0.89),
        fontsize=LEGEND_FONT_SIZE,
        frameon=False,
    )

    plt.tight_layout()
    output_path = output_dir / f"{OUTPUT_BASENAME}_{DATE_STAMP}.jpeg"
    plt.savefig(output_path, dpi=300)
    plt.close()
    return output_path


def main() -> None:
    if not SOURCE_FILE.exists():
        raise FileNotFoundError(
            f"Source workbook not found: {SOURCE_FILE}\n"
            f"Expected location: {DATA_DIR / SOURCE_FILE.name}"
        )

    output_dir = ensure_output_directory()
    cleaned_data = load_and_preprocess_data(SOURCE_FILE)
    site_counts = compute_site_distribution(cleaned_data)
    output_path = plot_primary_tumor_site_distribution(site_counts, output_dir)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
