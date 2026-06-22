"""Generate DDA score distribution and waterfall plots for PDX treatment cohorts.

Produces four publication figures:
  1. Compound-level strip plot: DDA score by compound, coloured by tier.
  2. Patient-level strip plot: DDA score by patient model, colour- and marker-coded
     by compound.
  3. Waterfall plot (Figure S2A): BestResponse stratified by SHIVA tier.
  4. Waterfall plot (Figure S2B): BestAvgResponse stratified by SHIVA tier.

Distribution tier thresholds (compound and patient plots):
  Low:          LEVEL < 0
  Intermediate: 0 <= LEVEL <= 1000
  High:         LEVEL > 1000

Waterfall SHIVA tier thresholds:
  LOW:          LEVEL < 0
  INTERMEDIATE: 0 <= LEVEL < 1000
  HIGH:         LEVEL >= 1000

Chemotherapy and Tubulin-inhibitor rows are excluded from distribution plots.
Duplicate (Model, COMPOUND, LEVEL) combinations are removed (first kept).
Waterfall plots use the full workbook without those filters.

Usage
-----
python scores_dist_plots.py [--source PATH] [--project-root DIR] [--seed INT]

Random seed (--seed, default 42) controls jitter in strip plots.
Outputs are written to <project-root>/scores_dist_plots_<timestamp>/.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


SCRIPT_NAME = "scores_dist_plots"

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

DEFAULT_SOURCE_XLSX = DATA_DIR / "pdx_curve_metrics_single_treatments_dda_scores.xlsx"

DCR_ORR_FONT = 18
SYMLOG_XTICKS_COMPOUND = [
    -100000,
    -10000,
    -1000,
    -100,
    -10,
    -1,
    1,
    10,
    100,
    1000,
    10000,
    100000,
]
SYMLOG_XTICKS_PATIENT = [-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000]

MARKER_STYLES = ["o", "s", "D", "^", "v", "<", ">", "p", "h", "H", "*"]
HEX_COLORS = [
    "#332288",
    "#88CCEE",
    "#44AA99",
    "#117733",
    "#999933",
    "#DDCC77",
    "#CC6677",
    "#882255",
    "#AA4499",
    "#661100",
    "#6699CC",
    "#EECC88",
    "#EE8866",
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#0077BB",
    "#EE7733",
    "#009988",
]
MARKER_NAMES = {
    "o": "circle",
    "s": "square",
    "D": "diamond",
    "^": "triangle_up",
    "v": "triangle_down",
    "<": "triangle_left",
    ">": "triangle_right",
    "p": "pentagon",
    "h": "hexagon1",
    "H": "hexagon2",
    "*": "star",
}

DISTRIBUTION_REQUIRED_COLUMNS = [
    "Model",
    "COMPOUND",
    "LEVEL",
    "TimeToDouble",
    "Treatment target",
    "ResponseCategory",
]

SHIVA_ORDER = ["HIGH", "INTERMEDIATE", "LOW"]
SHIVA_PALETTE = {
    "HIGH": "#50C878",
    "INTERMEDIATE": "#ffbc15",
    "LOW": "#f71c6c",
}


@dataclass(frozen=True)
class EndpointConfig:
    """Configuration for one waterfall endpoint panel set."""

    figure_label: str
    name: str
    upper_cutoff: float
    lower_cutoff: float
    y_min: float = -100
    y_max: float = 100


ENDPOINTS: tuple[EndpointConfig, ...] = (
    EndpointConfig(
        figure_label="Figure S2A",
        name="BestResponse",
        upper_cutoff=35,
        lower_cutoff=-50,
    ),
    EndpointConfig(
        figure_label="Figure S2B",
        name="BestAvgResponse",
        upper_cutoff=30,
        lower_cutoff=-20,
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate DDA score distribution and waterfall plots for PDX cohorts."
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
    parser.add_argument("--seed", type=int, default=42, help="Random seed for jitter.")
    return parser.parse_args()


def symlog_tick_labels(values: Iterable[int]) -> List[str]:
    labels: List[str] = []
    for value in values:
        if value == 0:
            labels.append("0")
            continue
        exponent = int(round(np.log10(abs(value))))
        if np.isclose(abs(value), 10**exponent):
            labels.append(rf"$-10^{{{exponent}}}$" if value < 0 else rf"$10^{{{exponent}}}$")
        else:
            labels.append(str(value))
    return labels


def build_output_dirs(project_root: Path, run_stamp: str) -> Dict[str, Path]:
    run_root = project_root / f"{SCRIPT_NAME}_{run_stamp}"
    figures_dir = run_root / "figures"
    tables_dir = run_root / "tables"
    run_root.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    return {"run_root": run_root, "figures": figures_dir, "tables": tables_dir}


def harmonize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "Treatment target" not in out.columns and "TARGET" in out.columns:
        out = out.rename(columns={"TARGET": "Treatment target"})
    if "Treatment" in out.columns and "COMPOUND" not in out.columns:
        out = out.rename(columns={"Treatment": "COMPOUND"})
    return out


def load_distribution_data(source_xlsx: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess data for distribution strip plots.

    Excludes chemotherapy and Tubulin-inhibitor rows and removes duplicate
    (Model, COMPOUND, LEVEL) combinations (first occurrence kept).

    Returns (cleaned_data, excluded_duplicates).
    """
    if not source_xlsx.exists():
        raise FileNotFoundError(f"Input workbook not found: {source_xlsx}")

    data = harmonize_column_names(pd.read_excel(source_xlsx))
    missing_columns = [column for column in DISTRIBUTION_REQUIRED_COLUMNS if column not in data.columns]
    if missing_columns:
        raise ValueError(
            "Input workbook is missing required columns: "
            f"{missing_columns}. Present columns: {list(data.columns)}"
        )

    data["Model"] = data["Model"].astype(str).str.replace("-", "", regex=False).str.lower()
    data["COMPOUND"] = data["COMPOUND"].astype(str).str.lower()

    data = data[data["LEVEL"] != "#HIÁNYZIK"].copy()
    data["LEVEL"] = pd.to_numeric(data["LEVEL"], errors="coerce")
    data["TimeToDouble"] = pd.to_numeric(data["TimeToDouble"], errors="coerce")
    data = data.dropna(subset=["LEVEL"]).copy()
    data["CENSOR"] = True

    data = data[
        ~data["Treatment target"].str.contains(r"\bchemotherapy\b", case=False, na=False)
    ].copy()
    data = data[~data["Treatment target"].str.contains(r"\bTubulin\b", case=False, na=False)].copy()

    dedup_subset = ["Model", "COMPOUND", "LEVEL"]
    duplicate_mask = data.duplicated(subset=dedup_subset, keep="first")
    excluded_duplicates = data[duplicate_mask].copy()
    cleaned_data = data[~duplicate_mask].copy()
    return cleaned_data, excluded_duplicates


def level_category(value: float) -> str:
    if value > 1000:
        return "High"
    if value >= 0:
        return "Intermediate"
    return "Low"


def configure_distribution_style(seed: int) -> None:
    np.random.seed(seed)
    sns.set_theme(style="white", context="talk")
    plt.rcParams.update(
        {
            "font.size": DCR_ORR_FONT,
            "axes.titlesize": DCR_ORR_FONT,
            "axes.labelsize": DCR_ORR_FONT,
            "legend.fontsize": DCR_ORR_FONT,
            "figure.dpi": 150,
            "savefig.dpi": 600,
        }
    )


def plot_compound_distribution(cleaned_data: pd.DataFrame, output_dir: Path, date_stamp: str) -> Path:
    """Strip plot of DDA score by compound, ordered by intra-compound score range."""
    df = cleaned_data.copy()
    df["Tiers"] = df["LEVEL"].apply(level_category)

    compound_order = (
        df.groupby("COMPOUND")["LEVEL"].agg(lambda x: x.max() - x.min()).sort_values(ascending=False).index
    )
    tier_palette = {"High": "#4fd260", "Intermediate": "#ffbc15", "Low": "#ff687a"}

    plt.figure(figsize=(16, 9))
    sns.stripplot(
        data=df,
        y="COMPOUND",
        x="LEVEL",
        hue="Tiers",
        palette=tier_palette,
        dodge=False,
        jitter=True,
        order=compound_order,
        size=8,
        alpha=0.80,
    )
    plt.legend(
        title="Tiers",
        loc="lower right",
        bbox_to_anchor=(1.26, 0.80),
        borderaxespad=0,
        frameon=True,
        fontsize=DCR_ORR_FONT,
        title_fontsize=DCR_ORR_FONT,
    )
    plt.xscale("symlog")
    plt.xticks(
        ticks=SYMLOG_XTICKS_COMPOUND,
        labels=symlog_tick_labels(SYMLOG_XTICKS_COMPOUND),
        rotation=90,
        fontsize=DCR_ORR_FONT,
    )
    plt.axvline(x=1000, color="black", linestyle=":", linewidth=2)
    plt.axvline(x=0, color="black", linestyle="-", linewidth=2)
    plt.xlabel("DDA score", fontsize=DCR_ORR_FONT)
    plt.ylabel(None)
    plt.yticks(fontsize=DCR_ORR_FONT)
    plt.tight_layout()

    output_path = output_dir / f"PDX_compound_plot_{date_stamp}.jpg"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()
    return output_path


def make_compound_style_map(compounds: Iterable[str]) -> Dict[str, Dict[str, str]]:
    sorted_compounds = sorted(set(compounds), key=str)
    num_combinations = len(MARKER_STYLES) * len(HEX_COLORS)
    if len(sorted_compounds) > num_combinations:
        raise ValueError("Insufficient marker-colour combinations. Extend MARKER_STYLES or HEX_COLORS.")

    styles: Dict[str, Dict[str, str]] = {}
    for index, compound in enumerate(sorted_compounds):
        styles[compound] = {
            "marker": MARKER_STYLES[index % len(MARKER_STYLES)],
            "color": HEX_COLORS[(index // len(MARKER_STYLES)) % len(HEX_COLORS)],
        }
    return styles


def plot_patient_distribution(
    cleaned_data: pd.DataFrame,
    figure_dir: Path,
    table_dir: Path,
    date_stamp: str,
) -> Tuple[Path, Path]:
    """Strip plot of DDA score by patient model with compound-specific markers."""
    df = cleaned_data.copy()
    compound_styles = make_compound_style_map(df["COMPOUND"].tolist())

    map_df = pd.DataFrame(
        [
            {
                "COMPOUND": compound,
                "marker": style["marker"],
                "marker_name": MARKER_NAMES[style["marker"]],
                "color_hex": style["color"],
            }
            for compound, style in sorted(compound_styles.items(), key=lambda item: str(item[0]))
        ]
    )
    map_path = table_dir / "COMPOUND_symbol_color_map.csv"
    map_df.to_csv(map_path, index=False)

    model_order = (
        df.groupby("Model")["LEVEL"].agg(lambda x: x.max() - x.min()).sort_values(ascending=False).index
    )
    model_counts = df.groupby("Model")["LEVEL"].count()

    marker_size = 13
    plt.figure(figsize=(13.0 * 1.5, 28 * 1.5))
    for compound, style in compound_styles.items():
        sns.stripplot(
            data=df[df["COMPOUND"] == compound],
            y="Model",
            x="LEVEL",
            marker=style["marker"],
            dodge=True,
            jitter=True,
            order=model_order,
            color=style["color"],
            size=marker_size,
            alpha=0.70,
        )

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker=style["marker"],
            color="w",
            label=compound,
            markerfacecolor=style["color"],
            markersize=marker_size,
        )
        for compound, style in compound_styles.items()
    ]
    plt.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.05, 1.0),
        title="COMPOUND",
        fontsize=DCR_ORR_FONT,
        title_fontsize=DCR_ORR_FONT,
    )

    ytick_labels = [f"{model} (n={model_counts[model]})" for model in model_order]
    plt.yticks(ticks=range(len(model_order)), labels=ytick_labels)
    plt.xlabel("DDA score", labelpad=15, fontsize=DCR_ORR_FONT)
    plt.xscale("symlog")
    plt.xticks(
        ticks=SYMLOG_XTICKS_PATIENT,
        labels=symlog_tick_labels(SYMLOG_XTICKS_PATIENT),
        rotation=90,
        fontsize=DCR_ORR_FONT,
    )
    plt.axvline(x=1000, color="black", linestyle=":", linewidth=2)
    plt.axvline(x=0, color="black", linestyle="-", linewidth=2)

    axis = plt.gca()
    axis.xaxis.set_ticks_position("top")
    axis.xaxis.set_label_position("top")
    plt.yticks(fontsize=DCR_ORR_FONT)
    plt.tight_layout()

    output_path = figure_dir / f"PDX_score_dist_compound_plot_by_case_{date_stamp}.jpg"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()
    return output_path, map_path


def ensure_required_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {', '.join(missing)}")


def assign_shiva_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Assign SHIVA tier labels from DDA LEVEL score for waterfall plots."""
    categorized = df.copy()

    categorized["SHIVA_scores"] = pd.NA
    categorized.loc[categorized["LEVEL"] < 0, "SHIVA_scores"] = "LOW"
    categorized.loc[
        (categorized["LEVEL"] >= 0) & (categorized["LEVEL"] < 1000), "SHIVA_scores"
    ] = "INTERMEDIATE"
    categorized.loc[categorized["LEVEL"] >= 1000, "SHIVA_scores"] = "HIGH"

    categorized = categorized.dropna(subset=["SHIVA_scores"]).copy()
    return categorized


def add_segmented_line(
    ax: plt.Axes,
    y_position: float,
    segments: list[tuple[float, float]],
    labels: list[str],
    font_size: int,
    line_width: float,
) -> None:
    """Draw a segmented arrow annotation below a waterfall panel."""
    for i, (start, end) in enumerate(segments):
        arrow = mpatches.FancyArrowPatch(
            (start, y_position),
            (end, y_position),
            arrowstyle="<|-|>",
            mutation_scale=10,
            color="black",
            linewidth=line_width,
        )
        ax.add_patch(arrow)

        midpoint = (start + end) / 2
        ax.text(
            x=midpoint,
            y=y_position + 5,
            s=labels[i],
            ha="center",
            fontsize=font_size,
            color="black",
        )


def make_waterfall_plot(df: pd.DataFrame, cfg: EndpointConfig, output_dir: Path) -> Path:
    """Create and save a multi-panel waterfall plot for one endpoint."""
    ensure_required_columns(df, [cfg.name, "SHIVA_scores"])

    df_sorted = df.sort_values(by=cfg.name, ascending=False).reset_index(drop=True)

    plt.rcParams.update({"font.size": 31})
    font_size = 33
    line_width = 2

    fig = plt.figure(figsize=(23, 18))
    axes: list[plt.Axes] = []

    for panel_index, score in enumerate(SHIVA_ORDER):
        group_df = df_sorted[df_sorted["SHIVA_scores"] == score].reset_index(drop=True)
        if group_df.empty:
            continue

        n_total = len(group_df)
        n_upper = int((group_df[cfg.name] > cfg.upper_cutoff).sum())
        n_middle = int(
            ((group_df[cfg.name] > cfg.lower_cutoff) & (group_df[cfg.name] <= cfg.upper_cutoff)).sum()
        )
        n_lower = int((group_df[cfg.name] <= cfg.lower_cutoff).sum())

        if n_total == 0:
            pct_upper = pct_middle = pct_lower = 0
        else:
            pct_upper = round((n_upper / n_total) * 100)
            pct_middle = round((n_middle / n_total) * 100)
            pct_lower = round((n_lower / n_total) * 100)

        ax = fig.add_subplot(len(SHIVA_ORDER), 1, panel_index + 1)
        axes.append(ax)

        bar_width = 0.25
        bar_spacing = 3
        x_positions = [idx * (bar_width + bar_spacing) for idx in group_df.index]

        sns.barplot(
            x=x_positions,
            y=group_df[cfg.name],
            hue=group_df["SHIVA_scores"],
            palette=SHIVA_PALETTE,
            ax=ax,
            dodge=False,
            width=0.8,
        )

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        ax.grid(axis="y", linestyle="--", alpha=0.7)
        ax.set_title(f"{score} (n={n_total})", fontsize=font_size)
        ax.set_ylabel(cfg.name, fontsize=font_size)
        ax.axhline(0, color="black", linestyle="-", linewidth=line_width)
        ax.axhline(cfg.upper_cutoff, color="black", linestyle="--", linewidth=line_width)
        ax.axhline(cfg.lower_cutoff, color="black", linestyle="--", linewidth=line_width)

        ax.set_yticks([-100, -75, -50, -20, 0, 30, 50, 75, 100])
        ax.set_ylim(cfg.y_min, cfg.y_max)
        ax.set_xticks([])
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        segments = [
            (0, n_upper),
            (n_upper, n_upper + n_middle),
            (n_upper + n_middle, n_total),
        ]
        labels = [
            f"{n_upper} ({pct_upper}%)",
            f"{n_middle} ({pct_middle}%)",
            f"{n_lower} ({pct_lower}%)",
        ]
        add_segmented_line(
            ax=ax,
            y_position=-95,
            segments=segments,
            labels=labels,
            font_size=font_size,
            line_width=line_width,
        )

    if axes:
        axes[-1].set_xlabel("Treatment", labelpad=20, fontsize=font_size)

    fig.text(0.01, 0.995, cfg.figure_label, ha="left", va="top", fontsize=24, fontweight="bold")
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)

    run_date = datetime.now().strftime("%Y_%m_%d")
    output_file = output_dir / f"PDX_waterfall_{cfg.figure_label.replace(' ', '_')}_{run_date}_{cfg.name}.png"
    fig.savefig(output_file, dpi=300)
    plt.close(fig)
    return output_file


def save_reproducibility_files(
    run_root: Path,
    table_dir: Path,
    source_xlsx: Path,
    cleaned_data: pd.DataFrame,
    excluded_duplicates: pd.DataFrame,
    produced_files: List[Path],
    seed: int,
    run_stamp: str,
) -> Tuple[Path, Path]:
    duplicates_path = table_dir / "excluded_duplicates_model_compound_level.csv"
    excluded_duplicates.to_csv(duplicates_path, index=False)

    metadata = {
        "script_name": SCRIPT_NAME,
        "run_stamp": run_stamp,
        "source_xlsx": str(source_xlsx.resolve()),
        "random_seed": seed,
        "rows_after_cleaning": int(len(cleaned_data)),
        "excluded_duplicate_rows": int(len(excluded_duplicates)),
        "unique_models": int(cleaned_data["Model"].nunique()),
        "unique_compounds": int(cleaned_data["COMPOUND"].nunique()),
        "produced_files": [str(path) for path in produced_files + [duplicates_path]],
    }
    metadata_path = run_root / "run_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return duplicates_path, metadata_path


def main() -> None:
    args = parse_args()
    source_xlsx = args.source
    project_root = args.project_root
    seed = args.seed

    run_stamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    date_stamp = datetime.now().strftime("%Y_%m_%d")
    output_dirs = build_output_dirs(project_root=project_root, run_stamp=run_stamp)

    configure_distribution_style(seed=seed)
    cleaned_data, excluded_duplicates = load_distribution_data(source_xlsx=source_xlsx)

    produced_files: List[Path] = []
    produced_files.append(
        plot_compound_distribution(
            cleaned_data=cleaned_data,
            output_dir=output_dirs["figures"],
            date_stamp=date_stamp,
        )
    )
    patient_plot_path, map_path = plot_patient_distribution(
        cleaned_data=cleaned_data,
        figure_dir=output_dirs["figures"],
        table_dir=output_dirs["tables"],
        date_stamp=date_stamp,
    )
    produced_files.extend([patient_plot_path, map_path])

    raw_df = pd.read_excel(source_xlsx)
    ensure_required_columns(raw_df, ["LEVEL", "BestResponse", "BestAvgResponse"])
    categorized_df = assign_shiva_groups(raw_df)

    for endpoint_cfg in ENDPOINTS:
        produced_files.append(
            make_waterfall_plot(categorized_df, endpoint_cfg, output_dirs["figures"])
        )

    duplicates_path, metadata_path = save_reproducibility_files(
        run_root=output_dirs["run_root"],
        table_dir=output_dirs["tables"],
        source_xlsx=source_xlsx,
        cleaned_data=cleaned_data,
        excluded_duplicates=excluded_duplicates,
        produced_files=produced_files,
        seed=seed,
        run_stamp=run_stamp,
    )

    print("Analysis completed.")
    print(f"Run folder:               {output_dirs['run_root']}")
    print(f"Run metadata:             {metadata_path}")
    print(f"Excluded duplicates table:{duplicates_path}")
    for file_path in produced_files:
        print(f"Generated figure/table:   {file_path}")


if __name__ == "__main__":
    main()
