# DCR_ORR_DB.R
# Fig. 3a–c: DCR, ORR, and durable-benefit panels.
#   - Fig3a_F456_mRECIST_by_rank_ORR_DCR_Response_Best_Response
#   - Fig3b_F456_mRECIST_by_rank_ORR_DCR_Response_BestAvgResponse
#   - Fig3c_F467_Durable_Benefit_by_rank_CA_trend
#
# SCAM vs GAM (Fig. 3h) and the LMM effect forest (Fig. 3i) are written by
# src/data_modelling.R into data_modeling/figures/.
#
# Uses plot_data_cache.rds from a completed data_modelling.R run.
# Output:
#   <output-dir>/dcr_orr_db/figures/
#   <output-dir>/dcr_orr_db/run_metadata.json
#
# Usage (from repository root):
#   Rscript src/DCR_ORR_DB.R
#   Rscript src/DCR_ORR_DB.R --output-dir src --cache path/to/plot_data_cache.rds

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
})

## ---------------------------------------------------------------------------
## Repository-relative paths + CLI (aligned with KM_PDX.py / other src scripts)
## ---------------------------------------------------------------------------
.resolve_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE)))
  }
  ofile <- NULL
  for (i in seq_along(sys.frames())) {
    frame_ofile <- sys.frames()[[i]]$ofile
    if (!is.null(frame_ofile) && nzchar(frame_ofile)) {
      ofile <- frame_ofile
      break
    }
  }
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

.parse_cli_args <- function(args) {
  out <- list(output_dir = NULL, cache = NULL, ca_helper = NULL)
  i <- 1L
  while (i <= length(args)) {
    if (identical(args[[i]], "--output-dir") && i < length(args)) {
      out$output_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(args[[i]], "--cache") && i < length(args)) {
      out$cache <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(args[[i]], "--ca-helper") && i < length(args)) {
      out$ca_helper <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

.resolve_latest_normal_rank_cache <- function(data_modeling_root) {
  if (!dir.exists(data_modeling_root)) {
    return(NULL)
  }
  candidates <- list.files(
    data_modeling_root,
    pattern = "plot_data_cache\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(candidates) == 0L) {
    return(NULL)
  }
  # Prefer normal-rank modelling runs (exclude leave-one-out).
  preferred <- candidates[
    !grepl("leave_one_out", candidates, ignore.case = TRUE)
  ]
  pool <- if (length(preferred) > 0L) preferred else candidates
  info <- file.info(pool)
  pool[order(info$mtime, decreasing = TRUE)][[1]]
}

file_meta <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  meta <- list(path = path_norm, exists = file.exists(path))
  if (isTRUE(meta$exists)) {
    info <- file.info(path)
    meta$size_bytes <- as.character(as.integer(info$size))
    meta$mtime <- format(info$mtime, "%Y-%m-%dT%H:%M:%S")
    meta$md5 <- tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  }
  meta
}

write_json_simple <- function(path, named_list) {
  escape_json <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", as.character(x), perl = TRUE)
    x <- gsub("\"", "\\\\\"", x, perl = TRUE)
    x <- gsub("\n", "\\\\n", x, perl = TRUE)
    x <- gsub("\r", "", x, perl = TRUE)
    x
  }
  lines <- character(0)
  lines <- c(lines, "{")
  nms <- names(named_list)
  for (i in seq_along(named_list)) {
    key <- nms[[i]]
    val <- named_list[[i]]
    comma <- if (i < length(named_list)) "," else ""
    if (is.null(val)) {
      lines <- c(lines, sprintf('  "%s": null%s', key, comma))
    } else if (is.logical(val) && length(val) == 1L) {
      lines <- c(lines, sprintf('  "%s": %s%s', key, if (isTRUE(val)) "true" else "false", comma))
    } else if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
      lines <- c(lines, sprintf('  "%s": %s%s', key, as.character(val), comma))
    } else if (is.list(val) && !is.null(names(val))) {
      inner <- vapply(names(val), function(k) {
        sprintf('"%s": "%s"', k, escape_json(val[[k]]))
      }, character(1))
      lines <- c(lines, sprintf('  "%s": {%s}%s', key, paste(inner, collapse = ", "), comma))
    } else if (length(val) > 1L) {
      items <- paste0('"', escape_json(val), '"', collapse = ", ")
      lines <- c(lines, sprintf('  "%s": [%s]%s', key, items, comma))
    } else {
      lines <- c(lines, sprintf('  "%s": "%s"%s', key, escape_json(val), comma))
    }
  }
  lines <- c(lines, "}")
  writeLines(lines, path, useBytes = TRUE)
}

DCR_ORR_DB_OUTPUT_DIRNAME <- "dcr_orr_db"
SCRIPT_DIR <- .resolve_script_dir()
source(file.path(SCRIPT_DIR, "plot_typography.R"))
setup_arial_font()
REPO_ROOT <- if (basename(SCRIPT_DIR) == "src") dirname(SCRIPT_DIR) else SCRIPT_DIR
DATA_MODELING_OUTPUT_ROOT <- file.path(SCRIPT_DIR, "data_modeling")

cli_args <- .parse_cli_args(commandArgs(trailingOnly = TRUE))
output_parent <- if (!is.null(cli_args$output_dir)) {
  normalizePath(cli_args$output_dir, winslash = "/", mustWork = FALSE)
} else {
  SCRIPT_DIR
}
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)

run_date <- format(Sys.Date(), "%Y_%m_%d")
output_root <- file.path(output_parent, DCR_ORR_DB_OUTPUT_DIRNAME)
run_dir <- output_root
fig_dir <- file.path(run_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

cache_path <- if (!is.null(cli_args$cache)) {
  normalizePath(cli_args$cache, winslash = "/", mustWork = FALSE)
} else {
  .resolve_latest_normal_rank_cache(DATA_MODELING_OUTPUT_ROOT)
}
ca_helper_path <- if (!is.null(cli_args$ca_helper)) {
  normalizePath(cli_args$ca_helper, winslash = "/", mustWork = FALSE)
} else {
  file.path(SCRIPT_DIR, "cochran_armitage_trend.R")
}

## Logging (same spirit as KM_PDX.py / data_modelling.R)
timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
step_log_path <- file.path(run_dir, "step_log.txt")
analysis_log_path <- file.path(run_dir, "analysis_log.txt")
heartbeat_path <- file.path(run_dir, "last_activity.txt")
writeLines(paste0("=== DCR_ORR_DB run started: ", timestamp_now(), " ==="), step_log_path)
writeLines(paste0("=== DCR_ORR_DB run started: ", timestamp_now(), " ==="), analysis_log_path)
writeLines(paste0(timestamp_now(), " | Run started."), heartbeat_path)

log_step <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  cat(line, "\n")
  write(line, step_log_path, append = TRUE)
  write(paste0(timestamp_now(), " | ", message), heartbeat_path)
}

log_note <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  cat(line, "\n")
  write(line, analysis_log_path, append = TRUE)
  write(line, step_log_path, append = TRUE)
  write(paste0(timestamp_now(), " | ", message), heartbeat_path)
}

## ---------------------------------------------------------------------------
## Plot sizing / fonts
## ---------------------------------------------------------------------------
FIG_WIDTH <- 13
FIG_HEIGHT <- 10
FONT_SCALE_ALL <- 1.40
fsz <- function(x) x * FONT_SCALE_ALL

## F99 panels (Fig. 3a/3b/3c) occupy a third-width cell of the 4x6 Fig. 3
## composite grid. Matching the cell aspect makes the panel fill its slot, and
## the point sizes below are derived so text lands at the shared on-page sizes.
FIG_F99_GRID <- list(nrows = 4, ncols = 6, rowspan = 1, colspan = 2)
FIG_F99_CELL <- composite_cell_in(4, 6, 1, 2)
FIG_F99_WIDTH <- 8
FIG_F99_HEIGHT <- FIG_F99_WIDTH * FIG_F99_CELL$height / FIG_F99_CELL$width
F99_FONT <- panel_font_sizes(FIG_F99_WIDTH, FIG_F99_HEIGHT, 4, 6, 1, 2)

log_note(paste0("Output root -> ", normalizePath(output_root, winslash = "/", mustWork = FALSE)))
log_note(paste0("Run folder  -> ", normalizePath(run_dir, winslash = "/", mustWork = FALSE)))
log_note(paste0("Figures     -> ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE)))

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
space_plot_separators <- function(s) {
  if (length(s) == 0L) return(s)
  out <- as.character(s)
  out <- gsub("([A-Za-z\\)\\%\\|])(=)([^=<>])", "\\1 = \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(<)(?!=)", "\\1 < \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(>)(?!=)", "\\1 > \\3", out, perl = TRUE)
  gsub("\\s{2,}", " ", out, perl = TRUE)
}

quote_plotmath_text <- function(s) {
  paste0("'", gsub("'", "''", space_plot_separators(as.character(s)), fixed = TRUE), "'")
}

sprintf_ca_trend_label <- function(z, p, prefix = NULL) {
  left <- sprintf("Cochran-Armitage: Z = %s, ", format(round(z, 3), nsmall = 3))
  if (!is.null(prefix) && nzchar(as.character(prefix)[1])) {
    left <- paste0(as.character(prefix)[1], ": ", left)
  }
  paste0(
    "paste(",
    quote_plotmath_text(left),
    ", ",
    fmt_pval_equals_plotmath(p, digits = 1L),
    ")"
  )
}

sig_stars_from_p <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 5e-4) return("***")
  if (p < 5e-3) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

build_adjacent_prop_annotations <- function(
  summary_df,
  rank_col,
  success_col,
  total_col,
  trend_df = NULL,
  trend_rank_col = NULL,
  trend_value_col = NULL,
  layer_offset = 4,
  span_step = 2.7,
  y_cap = 106
) {
  if (!all(c(rank_col, success_col, total_col) %in% names(summary_df))) {
    return(data.frame())
  }
  ord <- summary_df %>%
    dplyr::filter(!is.na(.data[[rank_col]]), !is.na(.data[[success_col]]), !is.na(.data[[total_col]])) %>%
    dplyr::arrange(.data[[rank_col]])
  if (nrow(ord) < 2L) return(data.frame())

  comp_pairs <- tibble::tibble(x = rep(1L, nrow(ord) - 1L), xend = 2L:nrow(ord)) %>%
    dplyr::mutate(x = pmin(x, xend), xend = pmax(x, xend)) %>%
    dplyr::distinct(x, xend) %>%
    dplyr::mutate(span = xend - x) %>%
    dplyr::arrange(span, x, xend)

  comps <- lapply(seq_len(nrow(comp_pairs)), function(i) {
    i1 <- comp_pairs$x[i]
    i2 <- comp_pairs$xend[i]
    s1 <- as.numeric(ord[[success_col]][i1])
    n1 <- as.numeric(ord[[total_col]][i1])
    s2 <- as.numeric(ord[[success_col]][i2])
    n2 <- as.numeric(ord[[total_col]][i2])
    p_val <- NA_real_
    if (is.finite(s1) && is.finite(n1) && is.finite(s2) && is.finite(n2) &&
        n1 > 0 && n2 > 0 && s1 >= 0 && s2 >= 0 && s1 <= n1 && s2 <= n2) {
      p_val <- tryCatch(
        stats::fisher.test(matrix(c(s1, n1 - s1, s2, n2 - s2), nrow = 2, byrow = TRUE))$p.value,
        error = function(e) NA_real_
      )
    }
    tibble::tibble(x = i1, xend = i2, span = i2 - i1, p_value = p_val, label = sig_stars_from_p(p_val))
  })
  comp_df <- dplyr::bind_rows(comps)
  if (nrow(comp_df) == 0L) return(comp_df)
  idx_p <- which(is.finite(comp_df$p_value))
  if (length(idx_p) > 0L) {
    comp_df$p_fdr <- NA_real_
    comp_df$p_fdr[idx_p] <- stats::p.adjust(comp_df$p_value[idx_p], method = "BH")
    comp_df$label <- vapply(comp_df$p_fdr, sig_stars_from_p, character(1))
  }

  use_trend <- !is.null(trend_df) &&
    !is.null(trend_rank_col) &&
    !is.null(trend_value_col) &&
    all(c(trend_rank_col, trend_value_col) %in% names(trend_df))
  if (use_trend) {
    tr <- trend_df %>%
      dplyr::select(
        rank_val = dplyr::all_of(trend_rank_col),
        trend_val = dplyr::all_of(trend_value_col)
      ) %>%
      dplyr::filter(!is.na(rank_val), !is.na(trend_val))
    lookup <- stats::setNames(as.numeric(tr$trend_val), as.character(tr$rank_val))
    comp_df <- comp_df %>%
      dplyr::mutate(
        y_anchor = pmax(
          as.numeric(lookup[as.character(ord[[rank_col]][x])]),
          as.numeric(lookup[as.character(ord[[rank_col]][xend])]),
          na.rm = TRUE
        )
      )
    n_br <- nrow(comp_df)
    stack_base <- max(comp_df$y_anchor, na.rm = TRUE) + layer_offset
    comp_df$y <- pmin(y_cap, stack_base + (seq_len(n_br) - 1L) * span_step)
  } else {
    prop_pct <- 100 * ord[[success_col]] / pmax(ord[[total_col]], 1)
    base_y <- max(55, min(88, max(prop_pct, na.rm = TRUE) + 6))
    if (nrow(comp_df) == 1L) {
      comp_df$y <- base_y
    } else {
      step_y <- max(3, min(8, (y_cap - base_y) / (nrow(comp_df) - 1L)))
      comp_df$y <- pmin(y_cap, base_y + (seq_len(nrow(comp_df)) - 1L) * step_y)
    }
  }
  comp_df %>% dplyr::mutate(y_low = pmax(0, y - 1.8), x_mid = (x + xend) / 2)
}

add_pairwise_signif_brackets <- function(
  plot_obj,
  ann_df,
  text_size = pt_to_mm(F99_FONT$annotation),
  text_y_gap = 0.8,
  line_width = 0.35,
  line_color = "black",
  text_color = "black"
) {
  if (is.null(ann_df) || nrow(ann_df) == 0L) return(plot_obj)
  label_size <- pmax(0.55, text_size)
  for (i in seq_len(nrow(ann_df))) {
    plot_obj <- plot_obj +
      ggplot2::annotate("segment", x = ann_df$x[i], xend = ann_df$xend[i], y = ann_df$y[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate("segment", x = ann_df$x[i], xend = ann_df$x[i], y = ann_df$y_low[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate("segment", x = ann_df$xend[i], xend = ann_df$xend[i], y = ann_df$y_low[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate(
        "text",
        x = ann_df$x_mid[i],
        y = ann_df$y[i] + text_y_gap,
        label = ann_df$label[i],
        size = label_size,
        fontface = "bold",
        colour = text_color,
        vjust = 0
      )
  }
  plot_obj
}

save_fig <- function(plot_obj, path, width, height, dpi = 300) {
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
  log_note(paste0("Saved: ", basename(path)))
  base_no_ext <- sub("\\.[^.]*$", "", basename(path))
  svg_path <- file.path(dirname(path), paste0(base_no_ext, ".svg"))
  ggplot2::ggsave(svg_path, plot_obj, width = width, height = height, dpi = 320, device = "svg")
  log_note(paste0("Saved: ", basename(svg_path)))
  invisible(c(path, svg_path))
}

## ---------------------------------------------------------------------------
## Load cache + CA helper
## ---------------------------------------------------------------------------
if (is.null(cache_path) || !file.exists(cache_path)) {
  stop(
    "plot_data_cache.rds not found.\n",
    "Run data_modelling.R first, or pass --cache path/to/plot_data_cache.rds.\n",
    "Looked under: ", DATA_MODELING_OUTPUT_ROOT,
    call. = FALSE
  )
}
if (!file.exists(ca_helper_path)) {
  stop(
    "cochran_armitage_trend.R not found: ", ca_helper_path,
    "\nExpected beside this script in src/, or pass --ca-helper PATH.",
    call. = FALSE
  )
}

log_step(paste0("Loading CA helper: ", ca_helper_path))
source(ca_helper_path, local = FALSE)
ca_registry_reset()

log_step(paste0("Loading plot cache: ", cache_path))
plot_cache <- readRDS(cache_path)
required_fields <- c("df_all_ranks", "rank_levels", "max_rank")
missing_fields <- setdiff(required_fields, names(plot_cache))
if (length(missing_fields) > 0L) {
  stop("Plot cache missing: ", paste(missing_fields, collapse = ", "), call. = FALSE)
}
log_note(paste0("Plot cache OK (fields: ", paste(names(plot_cache), collapse = ", "), ")"))

df_all_ranks <- as.data.frame(plot_cache$df_all_ranks)
rank_levels <- as.integer(unlist(plot_cache$rank_levels))
if ("rank_labels" %in% names(plot_cache)) {
  rank_labels <- as.character(unlist(plot_cache$rank_labels))
} else {
  rank_labels <- as.character(rank_levels)
}
max_rank <- as.integer(unlist(plot_cache$max_rank))[1]

theme_pub <- ggplot2::theme_minimal(base_size = fsz(14)) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = fsz(16), color = "black"),
    plot.subtitle = ggplot2::element_text(size = fsz(13), color = "black"),
    axis.title = ggplot2::element_text(size = fsz(14), color = "black"),
    axis.text = ggplot2::element_text(size = fsz(12), color = "black"),
    axis.ticks = ggplot2::element_line(color = "black"),
    axis.ticks.length = grid::unit(0.15, "cm"),
    panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", color = "black")
  )

theme_pub_nogrid <- theme_pub +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "grey85")
  )

## Overrides for the F99 panels (Fig. 3a/3b/3c). The legend sits below the plot
## so the Cochran-Armitage captions can use the full panel width.
## Title and CA results text are shifted upward to clear the plot / brackets.
F99_LIFT_PT <- 5
## Extra headroom inside the panel so CA lines sit above the signif. brackets.
F99_Y_EXPAND_TOP <- 0.22
theme_f99 <- ggplot2::theme(
  axis.ticks.x = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(
    size = F99_FONT$title,
    margin = ggplot2::margin(t = F99_LIFT_PT, b = F99_LIFT_PT),
    face = "bold",
    hjust = 0.5
  ),
  plot.subtitle = ggplot2::element_text(
    size = F99_FONT$annotation,
    margin = ggplot2::margin(t = 0, b = 6 + F99_LIFT_PT),
    hjust = 0.5
  ),
  plot.margin = ggplot2::margin(t = 10 + F99_LIFT_PT, r = 5, b = 5, l = 5),
  axis.title = ggplot2::element_text(size = F99_FONT$label, color = "black"),
  axis.title.x = ggplot2::element_text(size = F99_FONT$label, margin = ggplot2::margin(t = 8)),
  axis.text = ggplot2::element_text(size = F99_FONT$tick, color = "black"),
  legend.position = "bottom",
  legend.box = "vertical",
  legend.direction = "horizontal",
  legend.box.just = "left",
  legend.key.height = grid::unit(0.5, "cm"),
  legend.box.spacing = grid::unit(0.15, "cm"),
  legend.margin = ggplot2::margin(t = 0, b = 0),
  legend.spacing.y = grid::unit(0.05, "cm"),
  legend.spacing.x = grid::unit(0.3, "cm"),
  legend.text = ggplot2::element_text(size = F99_FONT$legend, color = "black"),
  legend.title = ggplot2::element_text(size = F99_FONT$legend, color = "black")
)

mrecist_colors <- c(
  "CR" = "#1565c0",
  "PR" = "#bbdefb",
  "SD" = "#d4d0f5",
  "PD" = "#d9d9d9"
)


## Shared: summarize mRECIST / ORR / DCR by response column
## ---------------------------------------------------------------------------
summarize_mrecist_panel <- function(df, response_col, panel_id_suffix, fill_title, plot_title, plot_subtitle = NULL) {
  if (!response_col %in% colnames(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = paste0(panel_id_suffix, ": ", response_col, " missing")))
  }

  mrecist_by_rank <- df %>%
    dplyr::mutate(resp = trimws(toupper(as.character(.data[[response_col]])))) %>%
    dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
    dplyr::group_by(rank_pooled, resp) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(rank_pooled) %>%
    dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(resp = factor(resp, levels = c("PD", "SD", "PR", "CR")))

  dcr_by_rank <- df %>%
    dplyr::mutate(resp = trimws(toupper(as.character(.data[[response_col]])))) %>%
    dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
    dplyr::mutate(is_dcr = resp %in% c("CR", "PR", "SD")) %>%
    dplyr::group_by(rank_pooled) %>%
    dplyr::summarise(n_dcr = sum(is_dcr, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      dcr_pct = 100 * n_dcr / n,
      se_pct = dplyr::if_else(n > 0, 100 * sqrt((dcr_pct / 100) * (1 - dcr_pct / 100) / n), 0)
    )

  orr_by_rank <- df %>%
    dplyr::mutate(resp = trimws(toupper(as.character(.data[[response_col]])))) %>%
    dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
    dplyr::group_by(rank_pooled) %>%
    dplyr::summarise(
      n_orr = sum(resp %in% c("CR", "PR"), na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      orr_pct = dplyr::if_else(n > 0, 100 * n_orr / n, 0),
      se_pct = dplyr::if_else(n > 0, 100 * sqrt((orr_pct / 100) * (1 - orr_pct / 100) / n), 0)
    )

  if (nrow(mrecist_by_rank) == 0L || nrow(dcr_by_rank) == 0L || nrow(orr_by_rank) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = paste0(panel_id_suffix, ": insufficient data")))
  }

  dcr_ordered <- dcr_by_rank %>% dplyr::arrange(rank_pooled)
  orr_ordered <- orr_by_rank %>% dplyr::arrange(rank_pooled)
  trend_df_dcr <- data.frame(
    rank_pooled = rank_levels,
    trend_fit_dcr = stats::predict(stats::lm(dcr_pct ~ rank_pooled, data = dcr_ordered), newdata = data.frame(rank_pooled = rank_levels))
  )
  trend_df_orr <- data.frame(
    rank_pooled = rank_levels,
    trend_fit_orr = stats::predict(stats::lm(orr_pct ~ rank_pooled, data = orr_ordered), newdata = data.frame(rank_pooled = rank_levels))
  )

  ca_ann_dcr <- ca_trend_annotate(
    scores = dcr_ordered$rank_pooled, n = dcr_ordered$n, success = dcr_ordered$n_dcr,
    metric = "DCR", panel_id = panel_id_suffix, analysis = "normal_rank", stratification = "rank_pooled"
  )
  ca_ann_orr <- ca_trend_annotate(
    scores = orr_ordered$rank_pooled, n = orr_ordered$n, success = orr_ordered$n_orr,
    metric = "ORR", panel_id = panel_id_suffix, analysis = "normal_rank", stratification = "rank_pooled"
  )

  p <- ggplot2::ggplot(
    mrecist_by_rank,
    ggplot2::aes(
      x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
      y = pct,
      fill = resp
    )
  ) +
    ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
    ggplot2::geom_line(
      data = trend_df_dcr,
      ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_dcr, group = 1, color = "CA (DCR%)"),
      linewidth = 1, linetype = "solid", inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = trend_df_orr,
      ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_orr, group = 1, color = "CA (ORR%)"),
      linewidth = 1, linetype = "solid", inherit.aes = FALSE
    ) +
    ## vjust reduced by ~5 pt / annotation size so CA lines sit above the bars.
    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = ca_ann_orr$label, parse = TRUE,
      hjust = 1.02, vjust = 2.10 - (F99_LIFT_PT / max(F99_FONT$annotation, 1)),
      size = pt_to_mm(F99_FONT$annotation), colour = "black"
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = ca_ann_dcr$label, parse = TRUE,
      hjust = 1.02, vjust = 1.00 - (F99_LIFT_PT / max(F99_FONT$annotation, 1)),
      size = pt_to_mm(F99_FONT$annotation), colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = mrecist_colors, name = fill_title) +
    ggplot2::scale_color_manual(
      values = c("CA (DCR%)" = "#008000", "CA (ORR%)" = "#FFA500"),
      breaks = c("CA (DCR%)", "CA (ORR%)"),
      name = NULL
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = ggplot2::expansion(mult = c(0, F99_Y_EXPAND_TOP))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Drug Rank Within Tumor", y = "Best Response (%)", title = plot_title, subtitle = plot_subtitle) +
    theme_pub_nogrid +
    theme_f99
  if (is.null(plot_title) && (is.null(plot_subtitle) || !nzchar(as.character(plot_subtitle)[1]))) {
    p <- p + ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )
  }

  ann_orr <- build_adjacent_prop_annotations(
    summary_df = orr_ordered, rank_col = "rank_pooled", success_col = "n_orr", total_col = "n",
    trend_df = trend_df_orr, trend_rank_col = "rank_pooled", trend_value_col = "trend_fit_orr",
    layer_offset = 3.0, span_step = 5.0
  )
  ann_dcr <- build_adjacent_prop_annotations(
    summary_df = dcr_ordered, rank_col = "rank_pooled", success_col = "n_dcr", total_col = "n",
    trend_df = trend_df_dcr, trend_rank_col = "rank_pooled", trend_value_col = "trend_fit_dcr",
    layer_offset = 7.0, span_step = 5.0
  )
  p <- add_pairwise_signif_brackets(p, ann_orr)
  p <- add_pairwise_signif_brackets(p, ann_dcr)
  p
}

## ---------------------------------------------------------------------------
## 1) Fig3a Best_Response
## ---------------------------------------------------------------------------
p456 <- summarize_mrecist_panel(
  df = df_all_ranks,
  response_col = "RESPONSE_BEST_RESPONSE",
  panel_id_suffix = "F456_mRECIST_Best_Response",
  fill_title = "mRECIST\n(BestResponse)",
  plot_title = NULL,
  plot_subtitle = NULL
)
log_step("1/3 Writing Fig3a Best_Response panel")
save_fig(p456, file.path(fig_dir, "Fig3a_F456_mRECIST_by_rank_ORR_DCR_Response_Best_Response.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)

## ---------------------------------------------------------------------------
## 2) Fig3b BestAvgResponse
## ---------------------------------------------------------------------------
p456_avg <- summarize_mrecist_panel(
  df = df_all_ranks,
  response_col = "RESPONSE_BESTAVGRESPONSE",
  panel_id_suffix = "F456_mRECIST_BestAvgResponse",
  fill_title = "mRECIST\n(BestAvgResponse)",
  plot_title = NULL,
  plot_subtitle = NULL
)
log_step("2/3 Writing Fig3b BestAvgResponse panel")
save_fig(p456_avg, file.path(fig_dir, "Fig3b_F456_mRECIST_by_rank_ORR_DCR_Response_BestAvgResponse.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)

## ---------------------------------------------------------------------------
## 3) Fig3c Durable Benefit (TIME_TO_DOUBLE > 42)
## ---------------------------------------------------------------------------
db_colors <- c("DB" = "#81c784", "NonDB" = "#d9d9d9")
db_df_rank <- df_all_ranks %>%
  dplyr::filter(is.finite(TIME_TO_DOUBLE), !is.na(rank_pooled)) %>%
  dplyr::mutate(is_db = TIME_TO_DOUBLE > 42)

if (nrow(db_df_rank) > 0L) {
  db_by_rank <- db_df_rank %>%
    dplyr::group_by(rank_pooled) %>%
    dplyr::summarise(n_db = sum(is_db, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(db_pct = dplyr::if_else(n > 0, 100 * n_db / n, 0)) %>%
    dplyr::arrange(rank_pooled)

  db_comp_rank <- db_df_rank %>%
    dplyr::mutate(db_status = dplyr::if_else(is_db, "DB", "NonDB")) %>%
    dplyr::group_by(rank_pooled, db_status) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::right_join(
      expand.grid(rank_pooled = rank_levels, db_status = c("NonDB", "DB"), stringsAsFactors = FALSE),
      by = c("rank_pooled", "db_status")
    ) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
    dplyr::group_by(rank_pooled) %>%
    dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(db_status = factor(db_status, levels = c("NonDB", "DB")))

  trend_df_db_rank <- data.frame(
    rank_pooled = rank_levels,
    trend_fit_db = stats::predict(stats::lm(db_pct ~ rank_pooled, data = db_by_rank), newdata = data.frame(rank_pooled = rank_levels))
  )

  ca_ann_f467 <- ca_trend_annotate(
    scores = db_by_rank$rank_pooled, n = db_by_rank$n, success = db_by_rank$n_db,
    metric = "DB", panel_id = "F467_Durable_Benefit_by_rank", analysis = "normal_rank", stratification = "rank_pooled"
  )

  p467 <- ggplot2::ggplot(
    db_comp_rank,
    ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = db_status)
  ) +
    ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
    ggplot2::geom_line(
      data = trend_df_db_rank,
      ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_db, group = 1, color = "CA (DB%)"),
      linewidth = 1, linetype = "solid", inherit.aes = FALSE
    ) +
    ## Align CA (DB%) annotation with the top CA line used on Fig3a/3b (DCR).
    ggplot2::annotate(
      "text", x = Inf, y = Inf, label = ca_ann_f467$label, parse = TRUE,
      hjust = 1.02, vjust = 1.00 - (F99_LIFT_PT / max(F99_FONT$annotation, 1)),
      size = pt_to_mm(F99_FONT$annotation), colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = db_colors, name = "Durable Benefit\n(>42 days)") +
    ggplot2::scale_color_manual(values = c("CA (DB%)" = "black"), breaks = c("CA (DB%)"), name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = ggplot2::expansion(mult = c(0, F99_Y_EXPAND_TOP))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Drug Rank Within Tumor", y = "Durable Benefit (%)", title = NULL, subtitle = NULL) +
    theme_pub_nogrid +
    theme_f99 +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )

  ann_f467_db <- build_adjacent_prop_annotations(
    summary_df = db_by_rank, rank_col = "rank_pooled", success_col = "n_db", total_col = "n",
    trend_df = trend_df_db_rank, trend_rank_col = "rank_pooled", trend_value_col = "trend_fit_db",
    layer_offset = 3.0, span_step = 5.0
  )
  p467 <- add_pairwise_signif_brackets(p467, ann_f467_db)
} else {
  p467 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F467: Durable Benefit by ranks not available")
}

log_step("3/3 Writing Fig3c Durable Benefit panel")
save_fig(p467, file.path(fig_dir, "Fig3c_F467_Durable_Benefit_by_rank_CA_trend.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)

## Export CA comparison tables into the dated run folder
log_step("Exporting Cochran-Armitage registry")
ca_registry_export(run_dir)

## Reproducibility metadata + output pointer (KM_PDX.py pattern)
produced_figures <- sort(list.files(fig_dir, full.names = FALSE))
metadata <- list(
  script_name = "DCR_ORR_DB.R",
  seed = 42,
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  run_date = run_date,
  output_root = normalizePath(output_root, winslash = "/", mustWork = FALSE),
  run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
  figures_dir = normalizePath(fig_dir, winslash = "/", mustWork = FALSE),
  cache_meta = file_meta(cache_path),
  ca_helper = normalizePath(ca_helper_path, winslash = "/", mustWork = FALSE),
  produced_figures = produced_figures
)
write_json_simple(file.path(output_root, "run_metadata.json"), metadata)
write_json_simple(file.path(run_dir, "run_metadata.json"), metadata)
write_json_simple(file.path(run_dir, "analysis_info.json"), list(
  run_key = "dcr_orr_db_figures",
  date = paste0(run_date, "_"),
  cache_path = normalizePath(cache_path, winslash = "/", mustWork = FALSE),
  cache_meta = file_meta(cache_path)
))

log_note(paste0("run_metadata.json written under ", output_root))
log_step(paste0("DCR_ORR_DB completed. Output: ", normalizePath(output_root, winslash = "/", mustWork = FALSE)))
log_note(paste0("Done. Figures written to: ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE)))
