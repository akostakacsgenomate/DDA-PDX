"""Generate DDA score distribution plots for PDX treatment cohorts.

Produces a compound-level strip plot: DDA score by compound, coloured by score sign.

Compound plot colours points by DDA score sign and magnitude:
  negative scores shade from dark red (most negative) to light red (near zero);
  positive scores shade from light blue (near zero) to dark blue (most positive).

Chemotherapy and Tubulin-inhibitor rows are excluded from distribution plots.
Duplicate (Model, COMPOUND, LEVEL) combinations are removed (first kept).

Usage
-----
python scores_dist_plots.py [--source PATH] [--project-root DIR] [--seed INT] [--run-stamp STAMP]

Random seed (--seed, default 42) controls jitter in strip plots.
Outputs are written to <project-root>/score_distributions/ (shared with
Waterfall_DDA_score_final.py and Case_level_distribution.R). Pass the same
--run-stamp to all three scripts if coordinating a single pipeline run.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plot_typography import apply_arial_font

apply_arial_font()


SCRIPT_NAME = "scores_dist_plots"
RUN_FOLDER_PREFIX = "score_distributions"

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "data"

DEFAULT_SOURCE_XLSX = DATA_DIR / "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_SOURCE_SHEET = "PCT_DRUG_RESPONSE"

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

DISTRIBUTION_REQUIRED_COLUMNS = [
    "Model",
    "COMPOUND",
    "LEVEL",
    "TimeToDouble",
    "Treatment target",
    "ResponseCategory",
]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate DDA score distribution plots for PDX cohorts."
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
    if "Treatment" in out.columns and "COMPOUND" not in out.columns:
        out = out.rename(columns={"Treatment": "COMPOUND"})
    if "ResponseCategory" not in out.columns and "Response_Best_Response" in out.columns:
        out["ResponseCategory"] = out["Response_Best_Response"]
    return out


def read_source_table(source_xlsx: Path, sheet: str = DEFAULT_SOURCE_SHEET) -> pd.DataFrame:
    if source_xlsx.suffix.lower() == ".csv":
        return pd.read_csv(source_xlsx)
    return pd.read_excel(source_xlsx, sheet_name=sheet)


def load_distribution_data(source_xlsx: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess data for distribution strip plots.

    Excludes chemotherapy and Tubulin-inhibitor rows and removes duplicate
    (Model, COMPOUND, LEVEL) combinations (first occurrence kept).

    Returns (cleaned_data, excluded_duplicates).
    """
    if not source_xlsx.exists():
        raise FileNotFoundError(f"Input workbook not found: {source_xlsx}")

    data = harmonize_column_names(read_source_table(source_xlsx))
    missing_columns = [column for column in DISTRIBUTION_REQUIRED_COLUMNS if column not in data.columns]
    if missing_columns:
        raise ValueError(
            "Input workbook is missing required columns: "
            f"{missing_columns}. Present columns: {list(data.columns)}"
        )

    if "Treatment type" in data.columns:
        data = data[data["Treatment type"].astype(str).str.strip().str.lower() == "single"].copy()
    if "On_label" in data.columns:
        data = data[data["On_label"].astype(str).str.strip().str.lower() != "chemo"].copy()

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


def _symlog_transform(values: np.ndarray, linthresh: float = 1.0) -> np.ndarray:
    """Match matplotlib symlog spacing for colour normalization."""
    values = np.asarray(values, dtype=float)
    out = np.zeros_like(values)
    nonzero = values != 0
    sign = np.sign(values[nonzero])
    magnitude = np.abs(values[nonzero])
    out[nonzero] = sign * linthresh * np.log10(magnitude / linthresh + 1.0)
    return out


def dda_level_colors(levels: np.ndarray, linthresh: float = 1.0) -> np.ndarray:
    """Map DDA scores to red (negative) or blue (positive) shading by magnitude."""
    levels = np.asarray(levels, dtype=float)
    transformed = _symlog_transform(levels, linthresh=linthresh)
    colors = np.empty((len(levels), 4))

    neg_mask = levels < 0
    if neg_mask.any():
        neg_transformed = transformed[neg_mask]
        neg_min = neg_transformed.min()
        span = -neg_min
        shade = (neg_transformed - neg_min) / span if span > 0 else np.ones_like(neg_transformed)
        colors[neg_mask] = plt.cm.Reds(0.85 - 0.50 * shade)

    pos_mask = levels >= 0
    if pos_mask.any():
        pos_transformed = transformed[pos_mask]
        pos_max = pos_transformed.max()
        shade = pos_transformed / pos_max if pos_max > 0 else np.zeros_like(pos_transformed)
        colors[pos_mask] = plt.cm.Blues(0.35 + 0.55 * shade)

    return colors


def plot_compound_distribution(cleaned_data: pd.DataFrame, output_dir: Path, date_stamp: str) -> Path:
    """Strip plot of DDA score by compound, ordered by intra-compound score range."""
    df = cleaned_data.copy()

    compound_order = (
        df.groupby("COMPOUND")["LEVEL"].agg(lambda x: x.max() - x.min()).sort_values(ascending=False).index
    )
    compound_order_list = list(compound_order)[::-1]
    compound_to_y = {compound: index for index, compound in enumerate(compound_order_list)}

    point_size = 8
    fig, ax = plt.subplots(figsize=(16, 9))
    y_base = df["COMPOUND"].map(compound_to_y).to_numpy()
    y_jitter = y_base + np.random.uniform(-0.25, 0.25, size=len(df))
    point_colors = dda_level_colors(df["LEVEL"].to_numpy())

    ax.scatter(
        df["LEVEL"],
        y_jitter,
        c=point_colors,
        s=point_size**2,
        alpha=0.80,
        linewidths=0,
        zorder=3,
    )
    ax.set_yticks(range(len(compound_order_list)))
    ax.set_yticklabels(compound_order_list, fontsize=DCR_ORR_FONT)
    ax.set_xscale("symlog")
    ax.set_xticks(SYMLOG_XTICKS_COMPOUND)
    ax.set_xticklabels(
        symlog_tick_labels(SYMLOG_XTICKS_COMPOUND),
        rotation=90,
        fontsize=DCR_ORR_FONT,
    )
    ax.axvline(x=0, color="black", linestyle="-", linewidth=2)
    ax.set_xlabel("DDA score", fontsize=DCR_ORR_FONT)
    ax.set_ylabel(None)
    fig.tight_layout()

    output_path = output_dir / f"Fig2d_PDX_compound_plot_{date_stamp}.jpg"
    fig.savefig(
        output_path,
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.12,
    )
    plt.close(fig)
    return output_path


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
    metadata_path = run_root / f"run_metadata_{SCRIPT_NAME}.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return duplicates_path, metadata_path


def main() -> None:
    args = parse_args()
    source_xlsx = args.source
    project_root = args.project_root
    seed = args.seed

    run_stamp = args.run_stamp or datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
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
