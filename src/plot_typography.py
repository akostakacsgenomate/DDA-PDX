"""Shared Arial font setup and composite-page font scaling.

Two responsibilities:

1. ``apply_arial_font`` makes Arial the preferred sans-serif without overriding
   ``font.family``, so existing per-script rcParams and layouts are preserved.
2. ``panel_font_sizes`` converts "how large should this text look on the printed
   A4 composite" into "what point size must the source figure use". Panels are
   authored at very different figure sizes and then shrunk by different amounts
   when placed into the composite grid, so identical point sizes in two source
   scripts do not produce identical text on the page.

Also provides shared on-plot *P*-value formatters. Prefer matplotlib mathtext
(``$...$``) over Unicode superscripts so Arial missing-glyph boxes do not appear.

The grid constants below must match the GridSpec used by
``build_figure_panel_pack.build_composite_a4``.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import math

A4_W_IN, A4_H_IN = 8.27, 11.69

GRID_LEFT, GRID_RIGHT = 0.04, 0.98
GRID_TOP, GRID_BOTTOM = 0.95, 0.03
GRID_WSPACE, GRID_HSPACE = 0.08, 0.12

# Target text height on the final A4 sheet, in points. Nature's floor is 5 pt.
ON_PAGE_PT: Dict[str, float] = {
    "annotation": 5.6,
    "legend": 6.5,
    "tick": 7.0,
    "label": 7.0,
    "title": 7.5,
}


def apply_arial_font() -> None:
    """Make Arial the first-choice sans-serif and ensure editable PDF text."""
    plt.rcParams.update({
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def composite_cell_in(
    nrows: int, ncols: int, rowspan: int = 1, colspan: int = 1
) -> Tuple[float, float]:
    """Size in inches of one composite grid cell (spans included)."""
    grid_w = (GRID_RIGHT - GRID_LEFT) * A4_W_IN
    grid_h = (GRID_TOP - GRID_BOTTOM) * A4_H_IN
    col_w = grid_w / (ncols + (ncols - 1) * GRID_WSPACE)
    row_h = grid_h / (nrows + (nrows - 1) * GRID_HSPACE)
    return (
        colspan * col_w + (colspan - 1) * GRID_WSPACE * col_w,
        rowspan * row_h + (rowspan - 1) * GRID_HSPACE * row_h,
    )


def composite_scale(
    fig_w_in: float,
    fig_h_in: float,
    nrows: int,
    ncols: int,
    rowspan: int = 1,
    colspan: int = 1,
) -> float:
    """Linear shrink factor applied to a source figure by the composite."""
    cell_w, cell_h = composite_cell_in(nrows, ncols, rowspan, colspan)
    rendered_w = min(cell_w, cell_h * fig_w_in / fig_h_in)
    return rendered_w / fig_w_in


def panel_font_sizes(
    fig_w_in: float,
    fig_h_in: float,
    nrows: int,
    ncols: int,
    rowspan: int = 1,
    colspan: int = 1,
) -> Dict[str, float]:
    """Source-figure point sizes that render at ``ON_PAGE_PT`` on the composite."""
    scale = composite_scale(fig_w_in, fig_h_in, nrows, ncols, rowspan, colspan)
    return {role: pt / scale for role, pt in ON_PAGE_PT.items()}


def _p_mantissa_exp(p_value: Any, digits: int = 1) -> Optional[Tuple[str, int]]:
    """Return (mantissa, exponent) for one-decimal scientific notation, or None."""
    try:
        p = float(p_value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(p):
        return None
    text = f"{p:.{int(digits)}e}"
    mantissa, exp = text.split("e")
    return mantissa, int(exp)


def format_p_sci_mathtext(p_value: Any, digits: int = 1) -> str:
    """Value-only mathtext, e.g. ``$5.6 \\\\times 10^{-2}$``."""
    parts = _p_mantissa_exp(p_value, digits=digits)
    if parts is None:
        return "n/a"
    mantissa, exp = parts
    return rf"${mantissa} \times 10^{{{exp}}}$"


def format_p_equals_mathtext(p_value: Any, digits: int = 1) -> str:
    """Italic *P* equality in mathtext, e.g. ``$\\mathit{P} = 5.6 \\\\times 10^{-2}$``."""
    parts = _p_mantissa_exp(p_value, digits=digits)
    if parts is None:
        return r"$\mathit{P}$ = n/a"
    mantissa, exp = parts
    return rf"$\mathit{{P}} = {mantissa} \times 10^{{{exp}}}$"


def format_p_sci_plain(p_value: Any, digits: int = 1) -> str:
    """ASCII fallback for tables/logs (not for Arial on-plot superscripts)."""
    parts = _p_mantissa_exp(p_value, digits=digits)
    if parts is None:
        return "n/a"
    mantissa, exp = parts
    return f"{mantissa}e{exp:+d}"
