## Shared Arial font setup and composite-page font scaling.
##
## Source this file and call setup_arial_font() once after loading ggplot2.
## Does NOT replace theme_minimal or any other theme - only sets the default
## text family so existing plot layouts are preserved.
##
## panel_font_sizes() mirrors plot_typography.py: it converts "how large should
## this text look on the printed A4 composite" into "what point size must the
## source figure use", because panels are authored at different figure sizes and
## then shrunk by different amounts when placed into the composite grid.

setup_arial_font <- function() {
  ggplot2::theme_update(text = ggplot2::element_text(family = "sans"))
}

A4_W_IN <- 8.27
A4_H_IN <- 11.69

GRID_LEFT <- 0.04
GRID_RIGHT <- 0.98
GRID_TOP <- 0.95
GRID_BOTTOM <- 0.03
GRID_WSPACE <- 0.08
GRID_HSPACE <- 0.12

## Target text height on the final A4 sheet, in points. Nature's floor is 5 pt.
ON_PAGE_PT <- list(
  annotation = 5.6,
  legend = 6.5,
  tick = 7.0,
  label = 7.0,
  title = 7.5
)

composite_cell_in <- function(nrows, ncols, rowspan = 1, colspan = 1) {
  grid_w <- (GRID_RIGHT - GRID_LEFT) * A4_W_IN
  grid_h <- (GRID_TOP - GRID_BOTTOM) * A4_H_IN
  col_w <- grid_w / (ncols + (ncols - 1) * GRID_WSPACE)
  row_h <- grid_h / (nrows + (nrows - 1) * GRID_HSPACE)
  list(
    width = colspan * col_w + (colspan - 1) * GRID_WSPACE * col_w,
    height = rowspan * row_h + (rowspan - 1) * GRID_HSPACE * row_h
  )
}

composite_scale <- function(fig_w_in, fig_h_in, nrows, ncols, rowspan = 1, colspan = 1) {
  cell <- composite_cell_in(nrows, ncols, rowspan, colspan)
  rendered_w <- min(cell$width, cell$height * fig_w_in / fig_h_in)
  rendered_w / fig_w_in
}

panel_font_sizes <- function(fig_w_in, fig_h_in, nrows, ncols, rowspan = 1, colspan = 1) {
  scale <- composite_scale(fig_w_in, fig_h_in, nrows, ncols, rowspan, colspan)
  lapply(ON_PAGE_PT, function(pt) pt / scale)
}

## ggplot2 geom_text()/annotate() sizes are in mm, theme text sizes are in pt.
pt_to_mm <- function(pt) pt / 2.845276

## ---------------------------------------------------------------------------
## On-plot P-value formatting (plotmath; one decimal; italic P)
## Example rendered form: italic(P) = 5.6 × 10^-2
## ---------------------------------------------------------------------------
.p_mantissa_exp <- function(p, digits = 1L) {
  if (length(p) == 0L || is.na(p) || !is.finite(as.numeric(p))) {
    return(NULL)
  }
  text <- formatC(as.numeric(p), format = "e", digits = as.integer(digits)[1])
  parts <- strsplit(text, "e", fixed = TRUE)[[1]]
  list(mantissa = parts[[1]], exp = as.integer(parts[[2]]))
}

## Value-only plotmath fragment: 5.6 %*% 10^{-2}
fmt_pval_sci_plotmath <- function(p, digits = 1L) {
  parts <- .p_mantissa_exp(p, digits = digits)
  if (is.null(parts)) return("'n/a'")
  paste0("'", parts$mantissa, "' %*% 10^{", parts$exp, "}")
}

## Full equality: paste(italic(P), ' = ', 5.6 %*% 10^{-2})
fmt_pval_equals_plotmath <- function(p, digits = 1L, op = "=") {
  parts <- .p_mantissa_exp(p, digits = digits)
  if (is.null(parts)) {
    return(paste0("paste(italic(P), ' ", op, " n/a')"))
  }
  paste0(
    "paste(italic(P), ' ", op, " ', ",
    "'", parts$mantissa, "' %*% 10^{", parts$exp, "}",
    ")"
  )
}
