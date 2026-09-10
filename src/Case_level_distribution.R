## =============================================================================
## Case_level_distribution.R
## Case-level DDA score distribution panel: score by case in rank order,
## with a companion rank panel.
##
## Input (same as scores_dist_plots.py):
##   ../data/Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx (sheet PCT_DRUG_RESPONSE)
##
## Output (shared run folder with scores_dist_plots.py and Waterfall_DDA_score_final.py):
##   <project-root>/score_distributions/figures/
##
## Usage:
##   Rscript src/Case_level_distribution.R
##   Rscript src/Case_level_distribution.R --source path/to/workbook.xlsx
##   Rscript src/Case_level_distribution.R --project-root src --run-stamp 2026_07_06__14_45_43
##
## Pass the same --run-stamp to all three scripts to collect every figure in one folder.
## =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

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
  out <- list(project_root = NULL, run_stamp = NULL, source = NULL)
  i <- 1L
  while (i <= length(args)) {
    if (identical(args[[i]], "--project-root") && i < length(args)) {
      out$project_root <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(args[[i]], "--run-stamp") && i < length(args)) {
      out$run_stamp <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(args[[i]], "--source") && i < length(args)) {
      out$source <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

.resolve_input_path <- function(explicit_source = NULL, data_dir, default_xlsx) {
  if (!is.null(explicit_source)) {
    path <- normalizePath(explicit_source, winslash = "/", mustWork = FALSE)
    if (!file.exists(path)) {
      stop("Input file not found: ", path, call. = FALSE)
    }
    return(path)
  }
  path <- file.path(data_dir, default_xlsx)
  if (!file.exists(path)) {
    stop(
      "Input workbook not found: ", path,
      "\nPlace ", default_xlsx, " in data/ and run from the repository root.",
      call. = FALSE
    )
  }
  path
}

RUN_FOLDER_PREFIX <- "score_distributions"
SCRIPT_DIR <- .resolve_script_dir()
REPO_ROOT <- if (basename(SCRIPT_DIR) == "src") dirname(SCRIPT_DIR) else SCRIPT_DIR
DATA_DIR <- file.path(REPO_ROOT, "data")
DEFAULT_INPUT_XLSX <- "Takacs_et_al_Supplementary_File_1_2026_07_06.xlsx"
DEFAULT_INPUT_SHEET <- "PCT_DRUG_RESPONSE"

cli_args <- .parse_cli_args(commandArgs(trailingOnly = TRUE))
input_path <- .resolve_input_path(cli_args$source, DATA_DIR, DEFAULT_INPUT_XLSX)
message("Using input workbook: ", input_path)
project_root <- if (!is.null(cli_args$project_root)) {
  normalizePath(cli_args$project_root, winslash = "/", mustWork = TRUE)
} else {
  SCRIPT_DIR
}
run_stamp <- if (!is.null(cli_args$run_stamp)) {
  cli_args$run_stamp
} else {
  format(Sys.time(), "%Y_%m_%d__%H_%M_%S")
}
date_stamp <- format(Sys.Date(), "%Y_%m_%d")
run_root <- file.path(project_root, RUN_FOLDER_PREFIX)
figures_dir <- file.path(run_root, "figures")
tables_dir <- file.path(run_root, "tables")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggnewscale)
})
source(file.path(SCRIPT_DIR, "plot_typography.R"))
setup_arial_font()
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required.", call. = FALSE)
}
library(patchwork)

# F22B panel saver; output path aligned with scores_dist_plots.py.
save_f22b_panel <- function(
    model_levels,
    score_dist_df,
    model_order_span_sd,
    row_step_sd,
    has_site_col_f9a,
    model_site_tbl_f9a,
    build_score_dist_plot,
    p_rank_for_panel,
    site_levels_panel,
    site_levels_f9a,
    site_color_map_panel,
    left_margin_pt_sd,
    fig_w_in,
    fig_h_in,
    fig_w_in_sd,
    fig_h_in_sd,
    figures_dir,
    date_stamp
) {
  if (missing(model_levels) || is.null(model_levels)) {
    message("F22B panel skipped: F8 rank order (model_levels) not available.")
    return(invisible(NULL))
  }

  models_in_score_dist <- unique(as.character(score_dist_df$Model))
  model_order_sd_rank <- as.character(model_levels)
  model_order_sd_rank <- model_order_sd_rank[model_order_sd_rank %in% models_in_score_dist]
  model_order_sd_rank <- c(
    model_order_sd_rank,
    setdiff(model_order_span_sd, model_order_sd_rank)
  )
  if (length(model_order_sd_rank) == 0L) {
    stop("F22B panel skipped: no overlapping models between F9 and F8 rank order.", call. = FALSE)
  }

  n_models_sd_rank <- length(model_order_sd_rank)
  score_dist_df_rank <- score_dist_df |>
    dplyr::mutate(
      Model = factor(as.character(.data$Model), levels = rev(model_order_sd_rank)),
      y_id = as.integer(factor(as.character(.data$Model), levels = rev(model_order_sd_rank))) * row_step_sd
    )
  axis_model_labels_sd_rank <- model_order_sd_rank
  model_site_tbl_f9a_rank <- if (has_site_col_f9a && nrow(model_site_tbl_f9a) > 0L) {
    model_site_tbl_f9a |>
      dplyr::mutate(
        Model = factor(as.character(.data$Model), levels = rev(model_order_sd_rank)),
        y_id = as.integer(factor(as.character(.data$Model), levels = rev(model_order_sd_rank))) * row_step_sd
      ) |>
      dplyr::filter(!is.na(.data$Model))
  } else {
    model_site_tbl_f9a
  }

  p_score_dist_a_rank <- if (has_site_col_f9a && nrow(model_site_tbl_f9a_rank) > 0L) {
    build_score_dist_plot(
      site_tbl = model_site_tbl_f9a_rank,
      show_legend = TRUE,
      df_in = score_dist_df_rank,
      axis_labels_in = axis_model_labels_sd_rank,
      n_models_in = n_models_sd_rank,
      show_vline_neg_1000 = FALSE,
      reverse_x = TRUE
    )
  } else {
    build_score_dist_plot(
      site_tbl = NULL,
      show_legend = TRUE,
      df_in = score_dist_df_rank,
      axis_labels_in = axis_model_labels_sd_rank,
      n_models_in = n_models_sd_rank,
      show_vline_neg_1000 = FALSE,
      reverse_x = TRUE
    )
  }

  p_f9_for_panel_b <- p_score_dist_a_rank +
    ggplot2::guides(shape = "none", colour = "none", fill = "none") +
    ggplot2::theme(plot.margin = ggplot2::margin(16, 8, 16, left_margin_pt_sd))
  if (length(site_levels_panel) > 0L && length(site_levels_f9a) > 0L) {
    p_f9_for_panel_b <- p_f9_for_panel_b +
      ggplot2::scale_fill_manual(
        values = site_color_map_panel,
        breaks = site_levels_panel,
        limits = site_levels_panel,
        drop = FALSE,
        name = "Primary Tumor Site",
        guide = "none"
      )
  }

  p_panel_f8_f9_b <- p_f9_for_panel_b + p_rank_for_panel +
    patchwork::plot_layout(ncol = 2, widths = c(1.0, 1.0), guides = "collect") &
    ggplot2::theme(legend.position = "right")

  f22b_panel_path <- file.path(figures_dir, sprintf("F22B_Score_ranking_panel_%s.png", date_stamp))
  ggplot2::ggsave(
    filename = f22b_panel_path,
    plot = p_panel_f8_f9_b,
    width = min(36, fig_w_in + fig_w_in_sd + 3.2),
    height = max(fig_h_in, fig_h_in_sd),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  message("Saved combined panel F22B (F9 in rank order + F8A): ", f22b_panel_path)
  invisible(f22b_panel_path)
}

# ---- Colour helpers (shared Primary Tumor Site palette) ----
CHM_VIVID_PALETTE <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#FFD700", "#F58231",
  "#911EB4", "#00BFFF", "#FF00FF", "#ADFF2F", "#FF1493",
  "#008B8B", "#8B00FF", "#9A6324", "#FF8C00", "#800000",
  "#00FA9A", "#808000", "#1E90FF", "#000075", "#696969",
  "#DC143C", "#00CED1", "#32CD32", "#FF4500", "#8A2BE2"
)
CHM_CLASS_DISTINCT_PALETTE <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
  "#AA3377", "#332288", "#009988", "#E69F00", "#9467BD"
)
chm_vivid_map <- function(levels_chr, rotate = 0L) {
  n <- length(levels_chr)
  if (n == 0L) return(character(0))
  pal_n <- length(CHM_VIVID_PALETTE)
  idx <- ((seq_len(n) - 1L + as.integer(rotate)) %% pal_n) + 1L
  stats::setNames(CHM_VIVID_PALETTE[idx], levels_chr)
}
chm_site_map <- function(levels_chr) {
  n <- length(levels_chr)
  if (n == 0L) return(character(0))
  site_cols <- chm_vivid_map(levels_chr, rotate = 0L)
  if ("Unknown" %in% names(site_cols)) site_cols[["Unknown"]] <- "#BDBDBD"
  if ("Multiple" %in% names(site_cols)) site_cols[["Multiple"]] <- "#4D4D4D"
  if ("Breast" %in% names(site_cols)) site_cols[["Breast"]] <- "#2E7D32"
  site_cols
}

compound_shape_pool_f8 <- c(16, 17, 15, 18, 8, 3, 4, 7, 0, 2, 5, 6)
compound_color_pool_f8 <- unique(c(CHM_VIVID_PALETTE, CHM_CLASS_DISTINCT_PALETTE))
build_compound_style_maps <- function(compound_levels_lc) {
  compound_levels_lc <- sort(unique(compound_levels_lc))
  n_comp <- length(compound_levels_lc)
  if (n_comp == 0L) {
    return(list(shape_map = character(0), color_map = character(0)))
  }
  list(
    shape_map = stats::setNames(
      compound_shape_pool_f8[((seq_len(n_comp) - 1L) %% length(compound_shape_pool_f8)) + 1L],
      compound_levels_lc
    ),
    color_map = stats::setNames(
      compound_color_pool_f8[((seq_len(n_comp) - 1L) %% length(compound_color_pool_f8)) + 1L],
      compound_levels_lc
    )
  )
}

ranking_viz_colors <- c(
  "1" = "#FFD700", "2" = "#C0C0C0", "3" = "#CD7F32", "4" = "#C2185B",
  "5" = "#64B5F6", "6" = "#D8F3DC", "7-10" = "#E8D5FF"
)

# ---- Read and normalise input ----
input_sheets <- readxl::excel_sheets(input_path)
compounds_sheet_name <- if (DEFAULT_INPUT_SHEET %in% input_sheets) {
  DEFAULT_INPUT_SHEET
} else if ("Compounds" %in% input_sheets) {
  "Compounds"
} else {
  input_sheets[[1]]
}
message("Reading: ", input_path, " (sheet = '", compounds_sheet_name, "')")
raw_data <- as.data.frame(readxl::read_excel(input_path, sheet = compounds_sheet_name))
raw_names_norm <- toupper(trimws(gsub("\\s+", "_", names(raw_data))))
keep_cols <- !duplicated(raw_names_norm) & !grepl("^UNNAMED", raw_names_norm)
raw_data <- raw_data[, keep_cols, drop = FALSE]
names(raw_data) <- raw_names_norm[keep_cols]

if ("PCT_ID" %in% names(raw_data) && !"MODEL" %in% names(raw_data)) {
  raw_data <- raw_data |> dplyr::rename(MODEL = PCT_ID)
}
if ("DDA_SCORE" %in% names(raw_data) && !"LEVEL" %in% names(raw_data)) {
  raw_data <- raw_data |> dplyr::rename(LEVEL = DDA_SCORE)
}
if ("TARGET" %in% names(raw_data) && !"TREATMENT_TARGET" %in% names(raw_data)) {
  raw_data <- raw_data |> dplyr::rename(TREATMENT_TARGET = TARGET)
}

col_builtin_cands <- c(
  "COMPOUND_BUILT-IN_GROUP_TYPE",
  "COMPOUND_BUILTIN_GROUP_TYPE",
  "COMPOUND_BUILT_IN_GROUP_TYPE"
)
hit_builtin <- col_builtin_cands[col_builtin_cands %in% names(raw_data)]
if (length(hit_builtin)) {
  col_builtin_group <- hit_builtin[[1]]
  v <- tolower(as.character(raw_data[[col_builtin_group]]))
  raw_data <- raw_data[is.na(v) | !grepl("hemothe", v, ignore.case = TRUE), , drop = FALSE]
}
if ("ON_LABEL" %in% names(raw_data)) {
  raw_data <- raw_data[
    tolower(trimws(as.character(raw_data$ON_LABEL))) != "chemo",
    ,
    drop = FALSE
  ]
}
if ("TREATMENT_TYPE" %in% names(raw_data)) {
  raw_data <- raw_data[
    tolower(trimws(as.character(raw_data$TREATMENT_TYPE))) == "single",
    ,
    drop = FALSE
  ]
}
if ("TREATMENT_TARGET" %in% names(raw_data)) {
  tt <- tolower(as.character(raw_data$TREATMENT_TARGET))
  raw_data <- raw_data[
    !grepl("chemotherapy|tubulin", tt, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
}

if ("REPORT_ID" %in% names(raw_data) && !"MODEL" %in% names(raw_data)) {
  raw_data <- raw_data |> dplyr::rename(MODEL = REPORT_ID)
}
if (!"MODEL" %in% names(raw_data)) {
  stop("Input must contain MODEL or REPORT_ID.", call. = FALSE)
}
raw_data <- raw_data |> dplyr::rename(Model = MODEL)

resolve_col <- function(candidates, df_cols) {
  hit <- candidates[candidates %in% df_cols]
  if (length(hit) == 0) NA_character_ else hit[1]
}
col_compound_dedup <- resolve_col(c("COMPOUND", "DRUG", "COMPOUND_NAME"), names(raw_data))
if (is.na(col_compound_dedup)) {
  stop("COMPOUND (or DRUG / COMPOUND_NAME) column required.", call. = FALSE)
}

raw_data <- raw_data |>
  dplyr::mutate(
    .source_row_id = dplyr::row_number(),
    Model = trimws(as.character(.data$Model)),
    COMPOUND = trimws(as.character(.data[[col_compound_dedup]])),
    LEVEL = suppressWarnings(as.numeric(.data$LEVEL))
  )
exact_dup_keep_ids <- raw_data |>
  dplyr::group_by(.data$Model, .data$COMPOUND, .data$LEVEL) |>
  dplyr::slice_min(.data$.source_row_id, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::pull(.data$.source_row_id)
raw_data_after_exact <- raw_data |> dplyr::filter(.data$.source_row_id %in% exact_dup_keep_ids)
model_compound_keep_ids <- raw_data_after_exact |>
  dplyr::group_by(.data$Model, .data$COMPOUND) |>
  dplyr::slice_max(.data$LEVEL, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::pull(.data$.source_row_id)
raw_data <- raw_data_after_exact |>
  dplyr::filter(.data$.source_row_id %in% model_compound_keep_ids) |>
  dplyr::select(-.data$.source_row_id)

# ---- F8: ranking panel (models ordered by n compounds descending) ----
raw_data_f8 <- raw_data
if ("PRIMARY_TUMOR_SITE" %in% names(raw_data_f8)) {
  raw_data_f8 <- raw_data_f8 |>
    dplyr::mutate(
      PRIMARY_TUMOR_SITE = trimws(as.character(.data$PRIMARY_TUMOR_SITE)),
      PRIMARY_TUMOR_SITE = dplyr::if_else(
        is.na(.data$PRIMARY_TUMOR_SITE) | .data$PRIMARY_TUMOR_SITE == "",
        NA_character_,
        .data$PRIMARY_TUMOR_SITE
      )
    ) |>
    dplyr::filter(!is.na(.data$PRIMARY_TUMOR_SITE), toupper(.data$PRIMARY_TUMOR_SITE) != "UNKNOWN")
}

rank_viz_base <- raw_data_f8 |>
  dplyr::mutate(
    Model = gsub("-", "", tolower(trimws(as.character(.data$Model)))),
    LEVEL = suppressWarnings(as.numeric(.data$LEVEL))
  ) |>
  dplyr::filter(!is.na(.data$LEVEL)) |>
  dplyr::group_by(Model) |>
  dplyr::filter(dplyr::n() >= 2L) |>
  dplyr::mutate(rank_bin = dplyr::row_number(dplyr::desc(LEVEL))) |>
  dplyr::ungroup()
if (nrow(rank_viz_base) == 0L) {
  stop("No models with >= 2 compounds and valid LEVEL.", call. = FALSE)
}

model_order_tbl <- rank_viz_base |>
  dplyr::count(Model, name = "n_rank_rows") |>
  dplyr::arrange(dplyr::desc(.data$n_rank_rows), .data$Model)
n_models_rank_viz <- nrow(model_order_tbl)
model_levels <- model_order_tbl |> dplyr::pull(Model) |> as.character()
model_label_map_f8 <- stats::setNames(model_levels, model_levels)

rank_viz_x_step <- 0.05
bracket <- rank_viz_x_step * 0.5
rank_col_x <- numeric(7L)
rank_col_x[1L:6L] <- 1 + (seq_len(6L) - 1L) * rank_viz_x_step
x_right_5_6 <- rank_col_x[6L] + bracket
tail_x_spread <- rank_viz_x_step
max_tail_any <- rank_viz_base |>
  dplyr::group_by(Model) |>
  dplyr::summarise(n_tail = sum(rank_bin >= 7L), .groups = "drop") |>
  dplyr::pull(n_tail) |>
  max(0L)
tail_point_inset <- max(0.018, rank_viz_x_step * 0.65)
tail_start_x <- x_right_5_6 + tail_point_inset
tail_band_width <- if (max_tail_any <= 1L) 0 else (max_tail_any - 1L) * tail_x_spread
rank_col_x[7L] <- tail_start_x + tail_band_width / 2
tail_box_xmin <- x_right_5_6
tail_box_xmax <- tail_start_x + tail_band_width + bracket + 0.06
row_step_f8 <- 1.19
rank_viz_boxes <- data.frame(
  xmin = c(rank_col_x[1L] - bracket, rank_col_x[3L] - bracket, rank_col_x[5L] - bracket, tail_box_xmin),
  xmax = c(rank_col_x[2L] + bracket, rank_col_x[4L] + bracket, rank_col_x[6L] + bracket, tail_box_xmax),
  ymin = 0.5 * row_step_f8,
  ymax = (n_models_rank_viz + 0.5) * row_step_f8
)
rank_col_x_num <- as.numeric(rank_col_x)
rank_col_mid <- (rank_col_x_num[-length(rank_col_x_num)] + rank_col_x_num[-1L]) / 2
x_left_1_2 <- rank_col_x[1L] - bracket
rank_column_bands_f8 <- tibble::tibble(
  rank_key = factor(c("1", "2", "3", "4", "5", "6", "7-10"), levels = c("1", "2", "3", "4", "5", "6", "7-10")),
  xmin = c(x_left_1_2, rank_col_mid[1L:5L], tail_box_xmin),
  xmax = c(rank_col_mid[1L:5L], x_right_5_6, tail_box_xmax),
  ymin = 0.5 * row_step_f8,
  ymax = (n_models_rank_viz + 0.5) * row_step_f8
)
rank_col_separators_f8 <- c(rank_col_mid[1L:5L], x_right_5_6)
rank_viz_pts <- rank_viz_base |>
  dplyr::arrange(Model, rank_bin) |>
  dplyr::group_by(Model) |>
  dplyr::mutate(
    rank_key = dplyr::if_else(rank_bin <= 6L, as.character(rank_bin), "7-10"),
    tail_idx = dplyr::if_else(rank_bin >= 7L, cumsum(rank_bin >= 7L), NA_integer_),
    x_plot = dplyr::if_else(
      rank_bin <= 6L,
      rank_col_x[rank_bin],
      tail_start_x + (tail_idx - 1L) * tail_x_spread
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    rank_key = factor(.data$rank_key, levels = c("1", "2", "3", "4", "5", "6", "7-10")),
    Model = factor(.data$Model, levels = model_levels),
    COMPOUND_LC = tolower(trimws(as.character(.data$COMPOUND))),
    y_id = as.integer(factor(as.character(.data$Model), levels = rev(model_levels))) * row_step_f8
  )

compound_breaks_f8 <- sort(unique(rank_viz_pts$COMPOUND_LC))
f8_style_maps <- build_compound_style_maps(compound_breaks_f8)
compound_shape_map_f8 <- f8_style_maps$shape_map
compound_color_map_f8 <- f8_style_maps$color_map
top_model_f8 <- model_levels[1L]
top_bar_y_f8 <- (n_models_rank_viz + 1.02) * row_step_f8
top_bar_pts_f8 <- rank_viz_pts |>
  dplyr::filter(as.character(.data$Model) == top_model_f8) |>
  dplyr::mutate(y_id = top_bar_y_f8)
x_left_f8 <- rank_col_x[1L] - bracket - 0.012
x_right_f8 <- tail_box_xmax + 0.05
site_strip_x_f8a <- x_left_f8 - 0.032
site_strip_xmin_f8a <- site_strip_x_f8a - 0.012
site_strip_xmax_f8a <- site_strip_x_f8a + 0.012
has_site_col_f8a <- "PRIMARY_TUMOR_SITE" %in% names(rank_viz_base)
model_site_tbl_f8a <- if (has_site_col_f8a) {
  rank_viz_base |>
    dplyr::mutate(
      PRIMARY_TUMOR_SITE = trimws(as.character(.data$PRIMARY_TUMOR_SITE)),
      PRIMARY_TUMOR_SITE = dplyr::if_else(
        is.na(.data$PRIMARY_TUMOR_SITE) | .data$PRIMARY_TUMOR_SITE == "",
        "Unknown",
        .data$PRIMARY_TUMOR_SITE
      )
    ) |>
    dplyr::group_by(Model) |>
    dplyr::summarise(
      site = if (dplyr::n_distinct(.data$PRIMARY_TUMOR_SITE) == 1L) {
        dplyr::first(.data$PRIMARY_TUMOR_SITE)
      } else {
        "Multiple"
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      y_id = as.integer(factor(.data$Model, levels = rev(model_levels))) * row_step_f8
    )
} else {
  tibble::tibble()
}
site_levels_f8a <- sort(unique(model_site_tbl_f8a$site))
site_color_map_f8a <- chm_site_map(site_levels_f8a)
y_lab_pt <- 2 * max(4.8, min(7.2, 520 / n_models_rank_viz))
left_margin_pt <- 10 + min(280, max(nchar(model_levels), na.rm = TRUE) * 4.4)
fig_h_in <- max(8, min(400, n_models_rank_viz * 0.13))
fig_w_in <- max(11, min(28, 10.5 + min(10, max(nchar(model_levels), na.rm = TRUE) * 0.06)))

build_rank_plot <- function(rank_pts_df, axis_model_levels, top_bar_df, site_tbl = NULL) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rank_column_bands_f8,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$rank_key),
      inherit.aes = FALSE, colour = NA, alpha = 0.18, show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = rank_col_separators_f8, colour = "black", linewidth = 0.25) +
    ggplot2::geom_rect(
      data = rank_viz_boxes,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.35
    ) +
    ggplot2::geom_point(
      data = rank_pts_df,
      ggplot2::aes(x = .data$x_plot, y = .data$y_id, shape = .data$COMPOUND_LC, colour = .data$COMPOUND_LC),
      stroke = 0.5, size = 3.0, alpha = 0.95, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = ranking_viz_colors, breaks = c("1", "2", "3", "4", "5", "6", "7-10"), guide = "none") +
    ggplot2::scale_shape_manual(values = compound_shape_map_f8, breaks = compound_breaks_f8, name = "COMPOUND") +
    ggplot2::scale_colour_manual(values = compound_color_map_f8, breaks = compound_breaks_f8, name = "COMPOUND", guide = "none") +
    ggplot2::scale_x_continuous(
      position = "top",
      breaks = rank_col_x,
      labels = c("1", "2", "3", "4", "5", "6", "7 - 10"),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0.012, 0.07))
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_models_rank_viz) * row_step_f8,
      labels = rev(axis_model_levels),
      expand = ggplot2::expansion(add = c(0.35, 0.65))
    ) +
    ggplot2::coord_cartesian(
      xlim = c(x_left_f8, x_right_f8),
      ylim = c(0.5 * row_step_f8, (n_models_rank_viz + 0.80) * row_step_f8),
      clip = "off"
    ) +
    ggplot2::labs(x = "Rank group", y = NULL, title = NULL) +
    ggplot2::theme_bw(base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.text.x = ggplot2::element_text(colour = "black", size = 15, vjust = 0.65, margin = ggplot2::margin(b = 0)),
      axis.text.y = ggplot2::element_text(colour = "black", size = y_lab_pt, margin = ggplot2::margin(r = 4, l = 0)),
      axis.title = ggplot2::element_text(colour = "black", size = 24),
      plot.title = ggplot2::element_text(hjust = 0.5, colour = "black", size = 18),
      legend.position = c(0.985, 1.0),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(16, 110, 16, left_margin_pt)
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        ncol = 1,
        override.aes = list(
          shape = unname(compound_shape_map_f8[compound_breaks_f8]),
          colour = unname(compound_color_map_f8[compound_breaks_f8]),
          stroke = 0.5, size = 3.0, alpha = 0.95
        )
      )
    )
  if (!is.null(site_tbl) && nrow(site_tbl) > 0L) {
    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = site_strip_xmin_f8a, xmax = site_strip_xmax_f8a,
          ymin = 0.5 * row_step_f8, ymax = (n_models_rank_viz + 0.5) * row_step_f8
        ),
        inherit.aes = FALSE, fill = "#F6F6F6", colour = "black", linewidth = 0.3
      ) +
      ggplot2::geom_point(
        data = site_tbl,
        ggplot2::aes(x = site_strip_x_f8a, y = .data$y_id, fill = .data$site),
        shape = 22, colour = "black", stroke = 0.25, size = 2.9, alpha = 1
      ) +
      ggplot2::scale_fill_manual(values = site_color_map_f8a, breaks = site_levels_f8a, name = "Primary Tumor Site") +
      ggplot2::coord_cartesian(
        xlim = c(site_strip_xmin_f8a - 0.01, x_right_f8),
        ylim = c(0.5 * row_step_f8, (n_models_rank_viz + 0.80) * row_step_f8),
        clip = "off"
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          ncol = 1,
          override.aes = list(shape = 22, size = 3.4, colour = "black", alpha = 1)
        )
      )
  }
  p
}

p_rank_viz <- if (has_site_col_f8a && nrow(model_site_tbl_f8a) > 0L) {
  build_rank_plot(rank_viz_pts, unname(model_label_map_f8[model_levels]), top_bar_pts_f8, site_tbl = model_site_tbl_f8a)
} else {
  build_rank_plot(rank_viz_pts, unname(model_label_map_f8[model_levels]), top_bar_pts_f8, site_tbl = NULL)
}

# ---- F9: DDA score distribution by case ----
score_dist_df <- raw_data |>
  dplyr::mutate(
    LEVEL_RAW = suppressWarnings(as.numeric(LEVEL)),
    COMPOUND_LC = tolower(trimws(as.character(COMPOUND))),
    Model = gsub("-", "", tolower(trimws(as.character(Model))))
  ) |>
  dplyr::filter(
    !is.na(LEVEL_RAW),
    !is.na(Model), trimws(Model) != "",
    !is.na(COMPOUND_LC), COMPOUND_LC != ""
  )
if ("PRIMARY_TUMOR_SITE" %in% names(score_dist_df)) {
  score_dist_df <- score_dist_df |>
    dplyr::mutate(
      PRIMARY_TUMOR_SITE = trimws(as.character(.data$PRIMARY_TUMOR_SITE)),
      PRIMARY_TUMOR_SITE = dplyr::if_else(
        is.na(.data$PRIMARY_TUMOR_SITE) | .data$PRIMARY_TUMOR_SITE == "",
        NA_character_,
        .data$PRIMARY_TUMOR_SITE
      )
    ) |>
    dplyr::filter(!is.na(.data$PRIMARY_TUMOR_SITE), toupper(.data$PRIMARY_TUMOR_SITE) != "UNKNOWN")
}
if (nrow(score_dist_df) == 0L) {
  stop("No valid Model/COMPOUND/LEVEL rows for score distribution.", call. = FALSE)
}

model_span_sd <- score_dist_df |>
  dplyr::group_by(Model) |>
  dplyr::summarise(level_span = max(LEVEL_RAW, na.rm = TRUE) - min(LEVEL_RAW, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$level_span), .data$Model)
model_order_span_sd <- model_span_sd |> dplyr::pull(Model) |> as.character()
model_order_sd <- model_order_span_sd[model_order_span_sd %in% unique(score_dist_df$Model)]
if (length(model_order_sd) > 0L) {
  score_dist_df <- score_dist_df |> dplyr::filter(as.character(.data$Model) %in% model_order_sd)
}
model_counts_sd <- score_dist_df |> dplyr::count(Model, name = "n_compounds")

compound_pool_panel <- unique(c(as.character(rank_viz_pts$COMPOUND_LC), as.character(score_dist_df$COMPOUND_LC)))
score_style_maps <- build_compound_style_maps(compound_pool_panel)
compound_shape_map_sd <- score_style_maps$shape_map
compound_color_map_sd <- score_style_maps$color_map
compound_breaks_sd <- names(compound_shape_map_sd)
row_step_sd <- row_step_f8
model_label_tbl <- score_dist_df |>
  dplyr::distinct(Model) |>
  dplyr::mutate(
    Model = factor(as.character(.data$Model), levels = rev(model_order_sd)),
    Model_lab = as.character(.data$Model)
  ) |>
  dplyr::mutate(model_chr = as.character(.data$Model)) |>
  dplyr::arrange(match(.data$model_chr, model_order_sd))
model_label_map <- stats::setNames(model_label_tbl$Model_lab, model_label_tbl$model_chr)
score_dist_df <- score_dist_df |>
  dplyr::left_join(model_counts_sd, by = "Model") |>
  dplyr::mutate(
    Model = factor(as.character(.data$Model), levels = rev(model_order_sd)),
    Model_lab = as.character(.data$Model),
    LEVEL_PLOT = -.data$LEVEL_RAW,
    y_id = as.integer(factor(as.character(.data$Model), levels = rev(model_order_sd))) * row_step_sd
  )
axis_model_labels_sd <- unname(model_label_map[model_order_sd])
x_breaks_sd_display <- c(10000, 1000, 100, 10, 1, -1, -10, -100, -1000, -10000)
pow10_label <- function(x) {
  vapply(x, function(v) {
    if (is.na(v) || v == 0) {
      "0"
    } else if (abs(v) == 1) {
      if (v < 0) "-10^{0}" else "10^{0}"
    } else {
      sign_txt <- if (v < 0) "-" else ""
      exp_val <- as.integer(round(log10(abs(v))))
      sprintf("%s10^{%d}", sign_txt, exp_val)
    }
  }, character(1))
}
n_models_sd <- length(unique(as.character(score_dist_df$Model)))
y_lab_pt_sd <- y_lab_pt
left_margin_pt_sd <- left_margin_pt
fig_h_in_sd <- fig_h_in * 1.35
fig_w_in_sd <- fig_w_in
has_site_col_f9a <- "PRIMARY_TUMOR_SITE" %in% names(score_dist_df)
model_site_tbl_f9a <- if (has_site_col_f9a) {
  score_dist_df |>
    dplyr::transmute(Model = as.character(.data$Model), site = trimws(as.character(.data$PRIMARY_TUMOR_SITE))) |>
    dplyr::mutate(site = dplyr::if_else(is.na(.data$site) | .data$site == "", "Unknown", .data$site)) |>
    dplyr::count(Model, site, name = "n_site") |>
    dplyr::group_by(Model) |>
    dplyr::arrange(dplyr::desc(.data$n_site), .data$site, .by_group = TRUE) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Model = factor(.data$Model, levels = rev(model_order_sd)),
      y_id = as.integer(factor(as.character(.data$Model), levels = rev(model_order_sd))) * row_step_sd
    ) |>
    dplyr::filter(!is.na(.data$Model))
} else {
  tibble::tibble()
}
site_levels_f9a <- sort(unique(model_site_tbl_f9a$site))
site_color_map_f9a <- chm_site_map(site_levels_f9a)

x_rng_sd <- range(c(score_dist_df$LEVEL_PLOT, -x_breaks_sd_display), na.rm = TRUE)
x_trans_sd <- scales::pseudo_log_trans(base = 10, sigma = 1)
x_rng_sd_t <- x_trans_sd$transform(x_rng_sd)
x_span_sd_t <- max(1e-6, diff(x_rng_sd_t))
site_strip_x_t_sd <- x_rng_sd_t[1] - 0.045 * x_span_sd_t
site_strip_half_w_t_sd <- 0.012 * x_span_sd_t
site_strip_xmin_sd <- x_trans_sd$inverse(site_strip_x_t_sd - site_strip_half_w_t_sd)
site_strip_xmax_sd <- x_trans_sd$inverse(site_strip_x_t_sd + site_strip_half_w_t_sd)
site_strip_x_sd <- x_trans_sd$inverse(site_strip_x_t_sd)
x_left_sd <- site_strip_xmin_sd
x_right_sd <- x_trans_sd$inverse(x_rng_sd_t[2] + 0.02 * x_span_sd_t)

build_score_dist_plot <- function(
    site_tbl = NULL,
    show_legend = TRUE,
    df_in = score_dist_df,
    axis_labels_in = axis_model_labels_sd,
    n_models_in = n_models_sd,
    show_vline_neg_1000 = TRUE,
    vline_neg_1000_linetype = "dotted",
    reverse_x = FALSE
) {
  if (isTRUE(reverse_x)) {
    trans_x_local <- scales::trans_new(
      name = "rev_pseudo_log",
      transform = function(x) -scales::pseudo_log_trans(base = 10, sigma = 1)$transform(x),
      inverse = function(x) scales::pseudo_log_trans(base = 10, sigma = 1)$inverse(-x),
      domain = c(-Inf, Inf)
    )
    site_strip_x_t_local <- x_rng_sd_t[2] + 0.045 * x_span_sd_t
    site_strip_xmin_local <- x_trans_sd$inverse(site_strip_x_t_local - site_strip_half_w_t_sd)
    site_strip_xmax_local <- x_trans_sd$inverse(site_strip_x_t_local + site_strip_half_w_t_sd)
    site_strip_x_local <- x_trans_sd$inverse(site_strip_x_t_local)
    x_left_local <- x_trans_sd$inverse(x_rng_sd_t[1] - 0.02 * x_span_sd_t)
    x_right_local <- site_strip_xmax_local
    x_scale_expand_local <- ggplot2::expansion(mult = c(0, 0))
  } else {
    trans_x_local <- scales::pseudo_log_trans(base = 10, sigma = 1)
    site_strip_xmin_local <- site_strip_xmin_sd
    site_strip_xmax_local <- site_strip_xmax_sd
    site_strip_x_local <- site_strip_x_sd
    x_left_local <- x_left_sd
    x_right_local <- x_right_sd
    x_scale_expand_local <- ggplot2::expansion(mult = c(0.05, 0.05))
  }
  p <- ggplot2::ggplot(
    df_in,
    ggplot2::aes(x = .data$LEVEL_PLOT, y = .data$y_id, shape = .data$COMPOUND_LC, colour = .data$COMPOUND_LC)
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 123),
      size = 3.0, alpha = 0.95, stroke = 0.5
    ) +
    ggplot2::scale_shape_manual(values = compound_shape_map_sd, breaks = compound_breaks_sd, name = "COMPOUND") +
    ggplot2::scale_colour_manual(values = compound_color_map_sd, breaks = compound_breaks_sd, name = "COMPOUND") +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_models_in) * row_step_sd,
      labels = rev(axis_labels_in),
      expand = ggplot2::expansion(add = c(0.35, 0.65))
    ) +
    ggplot2::scale_x_continuous(
      trans = trans_x_local,
      breaks = -x_breaks_sd_display,
      labels = function(x) parse(text = pow10_label(-x)),
      position = "top",
      expand = x_scale_expand_local
    ) +
    (if (isTRUE(show_vline_neg_1000)) {
      ggplot2::geom_vline(xintercept = -1000, colour = "black", linetype = vline_neg_1000_linetype, linewidth = 0.5)
    } else {
      NULL
    }) +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linetype = "solid", linewidth = 0.5) +
    ggplot2::coord_cartesian(
      xlim = c(x_left_local, x_right_local),
      ylim = c(0.5 * row_step_sd, (n_models_in + 0.80) * row_step_sd),
      clip = "off"
    ) +
    ggplot2::labs(x = "DDA score", y = "Model", title = NULL) +
    ggplot2::theme_bw(base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.text.x = ggplot2::element_text(angle = 0, colour = "black", size = 15, vjust = 0.5, hjust = 0.5, margin = ggplot2::margin(b = 0)),
      axis.text.y = ggplot2::element_text(colour = "black", size = y_lab_pt_sd, hjust = 1, margin = ggplot2::margin(r = 0, l = 0, t = 0, b = 0)),
      axis.title = ggplot2::element_text(colour = "black", size = 24),
      plot.title = ggplot2::element_text(hjust = 0.5, colour = "black", size = 18),
      legend.position = c(0.985, 1.0),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(16, 110, 16, left_margin_pt_sd)
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        ncol = 1,
        override.aes = list(
          shape = unname(compound_shape_map_sd[compound_breaks_sd]),
          colour = unname(compound_color_map_sd[compound_breaks_sd]),
          alpha = 1, size = 3.0, stroke = 0.5
        )
      ),
      colour = "none"
    )
  if (!is.null(site_tbl) && nrow(site_tbl) > 0L) {
    p <- p +
      ggnewscale::new_scale_fill() +
      ggplot2::annotate(
        "rect",
        xmin = site_strip_xmin_local, xmax = site_strip_xmax_local,
        ymin = 0.5 * row_step_sd, ymax = (n_models_in + 0.5) * row_step_sd,
        fill = "#F6F6F6", colour = "black", linewidth = 0.3
      ) +
      ggplot2::geom_point(
        data = site_tbl,
        ggplot2::aes(x = site_strip_x_local, y = .data$y_id, fill = .data$site),
        shape = 22, colour = "black", stroke = 0.25, size = 2.9, alpha = 1
      ) +
      ggplot2::scale_fill_manual(values = site_color_map_f9a, breaks = site_levels_f9a, name = "Primary Tumor Site") +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          ncol = 1,
          override.aes = list(shape = 22, size = 3.4, colour = "black", alpha = 1)
        )
      )
  }
  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

site_levels_panel <- sort(unique(c(site_levels_f8a, site_levels_f9a)))
site_color_map_panel <- chm_site_map(site_levels_panel)
p_rank_for_panel <- p_rank_viz +
  ggplot2::scale_shape_manual(
    values = compound_shape_map_sd, breaks = compound_breaks_sd,
    limits = compound_breaks_sd, drop = FALSE, name = "COMPOUND"
  ) +
  ggplot2::scale_colour_manual(
    values = compound_color_map_sd, breaks = compound_breaks_sd,
    limits = compound_breaks_sd, drop = FALSE, name = "COMPOUND", guide = "none"
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend(
      ncol = 1,
      override.aes = list(
        shape = unname(compound_shape_map_sd[compound_breaks_sd]),
        colour = unname(compound_color_map_sd[compound_breaks_sd]),
        stroke = 0.5, size = 3.0, alpha = 0.95
      )
    )
  ) +
  ggplot2::theme(plot.margin = ggplot2::margin(16, 8, 16, left_margin_pt))
if (length(site_levels_panel) > 0L && length(site_levels_f8a) > 0L) {
  p_rank_for_panel <- p_rank_for_panel +
    ggplot2::scale_fill_manual(
      values = site_color_map_panel, breaks = site_levels_panel,
      limits = site_levels_panel, drop = FALSE, name = "Primary Tumor Site"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = 1,
        override.aes = list(shape = 22, size = 3.4, colour = "black", alpha = 1)
      )
    )
}

# ---- F22B: case-level distribution panel ----
output_path <- save_f22b_panel(
  model_levels = model_levels,
  score_dist_df = score_dist_df,
  model_order_span_sd = model_order_span_sd,
  row_step_sd = row_step_sd,
  has_site_col_f9a = has_site_col_f9a,
  model_site_tbl_f9a = model_site_tbl_f9a,
  build_score_dist_plot = build_score_dist_plot,
  p_rank_for_panel = p_rank_for_panel,
  site_levels_panel = site_levels_panel,
  site_levels_f9a = site_levels_f9a,
  site_color_map_panel = site_color_map_panel,
  left_margin_pt_sd = left_margin_pt_sd,
  fig_w_in = fig_w_in,
  fig_h_in = fig_h_in,
  fig_w_in_sd = fig_w_in_sd,
  fig_h_in_sd = fig_h_in_sd,
  figures_dir = figures_dir,
  date_stamp = date_stamp
)
message("Run folder: ", run_root)
