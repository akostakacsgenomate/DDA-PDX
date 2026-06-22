## =============================================================================
## data_modelling.R
## Reproducible, standalone analysis script for peer review.
##
## Linear mixed models, SCAM/GAM regression models, and Cochran-Armitage trend tests.
##
## Publication outputs: SCAM residual diagnostics, effect forest, SCAM vs GAM,
## mRECIST response panels, durable benefit (Cochran-Armitage), and leave-one-out
## SCAM vs GAM panels per treatment target (section 3).
##
## Usage (from repository root):
##   Rscript src/data_modelling.R
##   PDX_PLOT_ONLY=1 Rscript src/data_modelling.R   # Step 9 figures only (requires plot_data_cache.rds)
## =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

## Repository-relative paths (standalone; run from repo root: Rscript src/data_modelling.R)
.resolve_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SCRIPT_DIR <- .resolve_script_dir()
REPO_ROOT <- if (basename(SCRIPT_DIR) == "src") dirname(SCRIPT_DIR) else SCRIPT_DIR
DATA_DIR <- file.path(REPO_ROOT, "data")
DEFAULT_INPUT_XLSX <- "pdx_curve_metrics_single_treatments_dda_scores.xlsx"
DATA_MODELING_OUTPUT_ROOT <- file.path(SCRIPT_DIR, "data_modeling_outputs")
if (!file.exists(file.path(DATA_DIR, DEFAULT_INPUT_XLSX))) {
  stop(paste0(
    "Input file not found: ", file.path(DATA_DIR, DEFAULT_INPUT_XLSX),
    "\nPlace the workbook in data/ and run from the repository root."
  ), call. = FALSE)
}


## SECTION 1: SHIVA-tier data modeling analysis


input_path <- file.path(DATA_DIR, DEFAULT_INPUT_XLSX)
## All analysis below uses all available rows with non-missing TimeToDouble values.
## Publication mixed-model contrasts: each SHIVA tier vs High (reference) on tier-classified PDX rows only.
## Inference: lmerTest (Satterthwaite p-values from summary); forest-plot CIs via profile likelihood (Wald fallback).
output_root <- DATA_MODELING_OUTPUT_ROOT
date <- format(Sys.Date(), "%Y_%m_%d")
name <- "shiva_tier_data_modeling"

set.seed(42)
PLOT_ONLY_MODE <- FALSE
if (identical(Sys.getenv("PDX_PLOT_ONLY"), "1")) PLOT_ONLY_MODE <- TRUE
SAVE_PLOT_DATA_CACHE <- TRUE

# Name of the response column (after renaming: spaces -> '_' and uppercasing).
# Default for earlier runs was 'TIMETODOUBLE'. Currently set to 'DAY_LAST'.
# Example alternative (for later use): 'TimeToDouble' -> internal name 'TIMETODOUBLE'
# RESPONSE_COL_NAME <- "TIMETODOUBLE"
RESPONSE_COL_NAME <- "TIMETODOUBLE"
# Label used in plot titles/axes for the response; kept simple so it always matches the chosen column.
RESPONSE_COL_LABEL <- RESPONSE_COL_NAME

# SHIVA ranks: three groups by LEVEL (all comparisons use High as reference).
# Low: LEVEL < 0; Intermediate: 0 <= LEVEL < 1000; High: LEVEL >= 1000.
SHIVA_GROUP_LOW  <- "Low"
SHIVA_GROUP_INT  <- "Intermediate"
SHIVA_GROUP_HIGH <- "High"
SHIVA_RANK_LEVELS <- c(1L, 2L, 3L)   # 1 = High, 2 = Intermediate, 3 = Low
SHIVA_RANK_LABELS <- c(SHIVA_GROUP_HIGH, SHIVA_GROUP_INT, SHIVA_GROUP_LOW)

# Performance metrics: ground truth from Response_BestAvgResponse.
#   truth = 1 (True)  when response in GROUND_TRUTH_TRUE (e.g. ORR: CR, PR)
#   truth = 0 (False) when response in GROUND_TRUTH_FALSE (e.g. ORR: PD, SD)
#   Prediction (SHIVA): High = predicted true (1, good); Low = predicted false (0).
GROUND_TRUTH_TRUE <- c("CR", "PR", "SD")   # ORR responders
GROUND_TRUTH_FALSE <- c("PD")   # non-responders (for ORR)
# Response type for performance metric titles, plot names, and ground_response in Excel.
RESPONSE_TYPE <- "DCR"

# Figure layout controls (publication style)
FIG_WIDTH <- 13
FIG_HEIGHT <- 10

# Per-plot figure sizes (edit here)
FIG_SCAM_MONO_WIDTH <- FIG_WIDTH
FIG_SCAM_MONO_HEIGHT <- FIG_HEIGHT
FIG_F22_JPG_WIDTH <- FIG_WIDTH
FIG_F22_JPG_HEIGHT <- FIG_HEIGHT
FIG_EFFECTS_LEGACY_WIDTH <- FIG_WIDTH
FIG_EFFECTS_LEGACY_HEIGHT <- FIG_HEIGHT
FIG_F2A_WIDTH <- FIG_WIDTH
FIG_F2A_HEIGHT <- FIG_HEIGHT
FIG_F24_WIDTH <- FIG_WIDTH
FIG_F24_HEIGHT <- FIG_HEIGHT
FIG_F2_WIDTH <- FIG_WIDTH
FIG_F2_HEIGHT <- FIG_HEIGHT
FIG_F3_WIDTH <- FIG_WIDTH
FIG_F3_HEIGHT <- FIG_HEIGHT
FIG_F25_WIDTH <- FIG_WIDTH
FIG_F25_HEIGHT <- FIG_HEIGHT
FIG_F4B_WIDTH <- FIG_WIDTH
FIG_F4B_HEIGHT <- FIG_HEIGHT
FIG_F6_WIDTH <- FIG_WIDTH
FIG_F6_HEIGHT <- FIG_HEIGHT
FIG_F22_WIDTH <- FIG_WIDTH
FIG_F22_HEIGHT <- FIG_HEIGHT
FIG_F22B_WIDTH <- FIG_WIDTH
FIG_F22B_HEIGHT <- FIG_HEIGHT
FIG_F22B_MARGIN <- ggplot2::margin(t = 40, r = 18, b = 12, l = 12, unit = "pt")
F22B_NLABEL_Y_PAD <- 8
FIG_F22C_WIDTH <- FIG_WIDTH
FIG_F22C_HEIGHT <- FIG_HEIGHT
FIG_F8_WIDTH <- FIG_WIDTH
FIG_F8_HEIGHT <- FIG_HEIGHT
FIG_F88_WIDTH <- 8
FIG_F88_HEIGHT <- 5
FIG_F99_WIDTH <- 8
FIG_F99_HEIGHT <- 5
FIG_F99B_WIDTH <- 8
FIG_F99B_HEIGHT <- 5
FIG_F99C_WIDTH <- 8
FIG_F99C_HEIGHT <- 5
FIG_F34_WIDTH <- 10
FIG_F34_HEIGHT <- 6
## F100: two heatmaps (DCR % | ORR %) by COMPOUND x SHIVA group; width fixed, height scales with n compounds
FIG_F100_WIDTH <- 14
FIG_F100_HEIGHT_BASE <- 7
FIG_F100_HEIGHT_PER_COMPOUND <- 0.32
FIG_F77_WIDTH <- 24
FIG_F77_HEIGHT <- 7
FIG_F23_WIDTH <- FIG_WIDTH
FIG_F23_HEIGHT <- FIG_HEIGHT
PANEL_HORIZONTAL_GAP <- 0.15  # F23: small gap between A-B and C-D
PANEL_VERTICAL_GAP <- 0.0    # F23: gap between top row (A-B) and bottom row (C-D)
FONT_SCALE_ALL <- 1.40

fsz <- function(x) x * FONT_SCALE_ALL

## Step timeout: if a step runs longer than this (seconds), abort (stuck detection). 0 = no limit.
STEP_TIMEOUT_SEC <- 0L
if (STEP_TIMEOUT_SEC > 0) {
  message("Step timeout enabled: ", STEP_TIMEOUT_SEC, "s per step (stuck detection).")
}

## Create date-named folder; save everything there (create if not present)
analysis_dir <- file.path(output_root, paste(date, name, sep = "_"))
fig_dir <- file.path(analysis_dir, "figures")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
plot_data_cache_path <- file.path(analysis_dir, "plot_data_cache.rds")


timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

step_log_path <- file.path(analysis_dir, "step_log.txt")
analysis_log_path <- file.path(analysis_dir, "analysis_log.txt")
heartbeat_path <- file.path(analysis_dir, "last_activity.txt")
writeLines(paste0("=== Analysis run started: ", timestamp_now(), " ==="), step_log_path)
writeLines(paste0("=== Analysis run started: ", timestamp_now(), " ==="), analysis_log_path)
writeLines(paste0(timestamp_now(), " | Run started."), heartbeat_path)

.step_counter <- 0L

log_step <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  cat(line, "\n")
  write(line, step_log_path, append = TRUE)
}

log_note <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  write(line, analysis_log_path, append = TRUE)
}

validate_plot_cache <- function(cache_obj) {
  required_fields <- c(
    "df", "df_all_ranks", "scam_out", "wiggle_out", "effects_out",
    "rank_levels", "max_rank"
  )
  missing_fields <- setdiff(required_fields, names(cache_obj))
  if (length(missing_fields) > 0) {
    stop(paste0("Plot cache missing required fields: ", paste(missing_fields, collapse = ", ")))
  }
  if (!is.data.frame(cache_obj$df) || !is.data.frame(cache_obj$df_all_ranks)) {
    stop("Plot cache validation failed: df/df_all_ranks are not data frames.")
  }
  invisible(TRUE)
}

## Ensure plot-facing text uses spaced separators, e.g. "p = 0.05", "n = 10", "EDF = 1.2".
space_plot_separators <- function(s) {
  if (length(s) == 0L) return(s)
  out <- as.character(s)
  out <- gsub("([A-Za-z\\)\\%\\|])(=)([^=<>])", "\\1 = \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(<)(?!=)", "\\1 < \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(>)(?!=)", "\\1 > \\3", out, perl = TRUE)
  gsub("\\s{2,}", " ", out, perl = TRUE)
}

run_step <- function(step_name, fn) {
  .step_counter <<- .step_counter + 1L
  n <- .step_counter
  start_time <- Sys.time()
  ts_start <- timestamp_now()

  start_msg <- paste0("========== Step ", n, ": ", step_name, " - START at ", ts_start, " ==========")
  cat(start_msg, "\n")
  write(start_msg, step_log_path, append = TRUE)
  write(paste0(ts_start, " | Step ", n, " START: ", step_name), heartbeat_path)

  if (STEP_TIMEOUT_SEC > 0) {
    setTimeLimit(elapsed = STEP_TIMEOUT_SEC, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  }
  out <- NULL
  err <- NULL
  tryCatch(
    {
      out <- fn()
    },
    error = function(e) {
      err <<- conditionMessage(e)
    }
  )
  if (!is.null(err)) {
    fail_msg <- paste0("========== Step ", n, ": ", step_name, " - FAIL at ", timestamp_now(), " ==========\n", err)
    cat(fail_msg, "\n")
    write(fail_msg, step_log_path, append = TRUE)
    write(paste0(timestamp_now(), " | Step ", n, " FAILED: ", step_name), heartbeat_path)
    stop(err)
  }

  ts_end <- timestamp_now()
  dur_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  ok_msg <- paste0("========== Step ", n, ": ", step_name, " - OK at ", ts_end, " (duration ", dur_sec, "s) ==========")
  cat(ok_msg, "\n")
  write(ok_msg, step_log_path, append = TRUE)
  write(paste0(ts_end, " | Step ", n, " OK (", dur_sec, "s): ", step_name), heartbeat_path)
  log_step(paste0("Heartbeat: ", ts_end, " - Step ", n, " completed"))
  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), "p = NA", paste0("p = ", format.pval(p, digits = 3, eps = 1e-4)))
}

fmt_p_sci1 <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", formatC(p, format = "e", digits = 0))
}

fmt_p_sci1_cap <- function(p, min_p = 1e-100) {
  if (is.na(p)) return("p = NA")
  if (p < min_p) return(paste0("p < ", format(min_p, scientific = TRUE)))
  paste0("p = ", formatC(p, format = "e", digits = 0))
}

fmt_pval_sci2 <- function(p) {
  if (is.na(p)) return("NA")
  formatC(p, format = "e", digits = 2)
}

fmt_pval_sci2_cap <- function(p, min_p = 1e-30) {
  if (is.na(p)) return("NA")
  if (p < min_p) return("<1e-30")
  formatC(p, format = "e", digits = 2)
}

## F22/F22b corner annotations: shared SCAM/GAM smooth p formatting (matches normal-ranks script).
## For p < cap_below (underflow / numeric zero is also caught because 0 < cap_below), shows p < floor.
## sci_digits = 2L => two decimals after the leading digit in scientific notation, e.g. "p = 4.19e-21".
fmt_sci_e <- function(x, sci_digits = 2L) {
  formatC(as.numeric(x), format = "e", digits = sci_digits)
}

fmt_edf_plot <- function(edf) {
  format(round(as.numeric(edf), 3), nsmall = 3)
}

fmt_smooth_p_plot_label <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  if (is.na(p) || !is.finite(as.numeric(p))) return("p = NA")
  p <- as.numeric(p)
  if (p < cap_below) {
    return(paste0("p < ", fmt_sci_e(cap_below, sci_digits)))
  }
  paste0("p = ", fmt_sci_e(p, sci_digits))
}

fmt_fdr <- function(p) {
  if (is.na(p)) return("FDR = NA")
  if (p < 0.001) return("FDR < 0.001")
  if (p < 0.01) return("FDR < 0.01")
  if (p < 0.05) return("FDR < 0.05")
  "FDR >= 0.05"
}

## On-figure labels only: render "p" as plotmath italic(p) (ggplot parse = TRUE).
num_str_to_plotmath_rhs <- function(s, mantissa_digits = 2L) {
  s <- trimws(as.character(s))
  if (!nzchar(s) || grepl("^[Nn][Aa]$", s)) return("NA")
  if (grepl("^[<>]", s)) return(s)
  m <- regmatches(s, regexec("^([0-9.+\\-]+)[eE]([-+]?\\d+)$", s, perl = TRUE))[[1]]
  if (length(m) == 3L) {
    mant <- sprintf(paste0("%.", mantissa_digits, "f"), as.numeric(m[2]))
    return(paste0("'", mant, "' %*% 10^{", as.integer(m[3]), "}"))
  }
  s
}

quote_plotmath_text <- function(s) {
  paste0("'", gsub("'", "''", space_plot_separators(as.character(s)), fixed = TRUE), "'")
}

## DCR/ORR Cochran-Armitage corner text: italic plotmath p; value = fmt_pval_sci2 (2-decimal e-notation).
sprintf_ca_trend_label <- function(z, p, prefix = NULL) {
  left <- sprintf("Cochran-Armitage: Z = %s, ", format(round(z, 3), nsmall = 3))
  if (!is.null(prefix) && nzchar(as.character(prefix)[1])) {
    left <- paste0(as.character(prefix)[1], ": ", left)
  }
  paste0(
    "paste(",
    quote_plotmath_text(left),
    ", italic(p), ' = ",
    gsub("'", "''", fmt_pval_sci2(p), fixed = TRUE),
    "')"
  )
}


## Cochran-Armitage trend test helpers (DCR / ORR / Durable Benefit panels)

if (!exists(".ca_registry", inherits = FALSE)) {
  .ca_registry <- new.env(parent = emptyenv())
  assign("rows", list(), envir = .ca_registry)
}

ca_registry_reset <- function() {
  assign("rows", list(), envir = .ca_registry)
}

ca_registry_add <- function(row) {
  if (!is.list(row)) stop("ca_registry_add expects a list.")
  rows <- get("rows", envir = .ca_registry)
  assign("rows", c(rows, list(row)), envir = .ca_registry)
}

ca_registry_get_rows <- function() {
  get("rows", envir = .ca_registry)
}

ca_registry_as_tibble <- function() {
  rows <- ca_registry_get_rows()
  if (length(rows) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(rows)
}

#' Manual Cochran-Armitage Z (normal approximation), matching legacy script.
ca_trend_manual <- function(success, n, scores) {
  success <- as.numeric(success)
  n <- as.numeric(n)
  scores <- as.numeric(scores)
  if (length(success) != length(n) || length(n) != length(scores)) {
    stop("success, n, and scores must have the same length.")
  }
  ok <- !is.na(success) & !is.na(n) & !is.na(scores) & n > 0
  success <- success[ok]
  n <- n[ok]
  scores <- scores[ok]
  if (length(success) < 2L) {
    return(list(
      z = 0, p = 1, chi2 = 0, p_overall = NA_real_, s_bar = NA_real_,
      method = "manual_Cochran-Armitage (insufficient groups)"
    ))
  }
  r_tot <- sum(success)
  n_tot <- sum(n)
  p_overall <- r_tot / n_tot
  s_bar <- sum(n * scores) / n_tot
  num <- sum((scores - s_bar) * (success - n * p_overall))
  den <- sqrt(p_overall * (1 - p_overall) * sum(n * (scores - s_bar)^2))
  z <- if (is.finite(den) && den > 0) num / den else 0
  p <- 2 * (1 - stats::pnorm(abs(z)))
  list(
    z = z,
    p = p,
    chi2 = z^2,
    p_overall = p_overall,
    s_bar = s_bar,
    method = "manual_Cochran-Armitage (normal approximation, two-sided)"
  )
}

#' Standard R implementation: prop.trend.test; Z = sign(trend) * sqrt(chi-square).
ca_trend_standard <- function(success, n, scores) {
  success <- as.numeric(success)
  n <- as.numeric(n)
  scores <- as.numeric(scores)
  if (length(success) != length(n) || length(n) != length(scores)) {
    stop("success, n, and scores must have the same length.")
  }
  ok <- !is.na(success) & !is.na(n) & !is.na(scores) & n > 0
  success <- success[ok]
  n <- n[ok]
  scores <- scores[ok]
  if (length(success) < 2L) {
    return(list(
      z = 0, p = 1, chi2 = 0,
      method = "stats::prop.trend.test (insufficient groups)",
      alternative = NA_character_
    ))
  }
  p_overall <- sum(success) / sum(n)
  s_bar <- sum(n * scores) / sum(n)
  manual_num <- sum((scores - s_bar) * (success - n * p_overall))
  pt <- tryCatch(
    stats::prop.trend.test(success, n, scores),
    error = function(e) NULL
  )
  if (is.null(pt)) {
    man <- ca_trend_manual(success, n, scores)
    return(list(
      z = man$z,
      p = man$p,
      chi2 = man$chi2,
      method = "stats::prop.trend.test failed; manual Cochran-Armitage fallback",
      alternative = "two.sided",
      data_name = NA_character_
    ))
  }
  z <- if (manual_num == 0) 0 else sign(manual_num) * sqrt(as.numeric(pt$statistic))
  list(
    z = z,
    p = pt$p.value,
    chi2 = as.numeric(pt$statistic),
    method = pt$method,
    alternative = pt$alternative,
    data_name = pt$data.name
  )
}

#' Run manual and standard CA; return comparison row for registry / export.
ca_trend_compare <- function(success, n, scores,
                             metric = NA_character_,
                             panel_id = NA_character_,
                             analysis = NA_character_,
                             stratification = NA_character_) {
  man <- ca_trend_manual(success, n, scores)
  std <- ca_trend_standard(success, n, scores)
  z_diff <- std$z - man$z
  p_diff <- std$p - man$p
  tibble::tibble(
    analysis = analysis,
    stratification = stratification,
    panel_id = panel_id,
    metric = metric,
    n_groups = length(success),
    scores = paste(scores, collapse = ","),
    n_per_group = paste(n, collapse = ","),
    success_per_group = paste(success, collapse = ","),
    p_overall = man$p_overall,
    s_bar_scores = man$s_bar,
    z_manual = man$z,
    p_manual = man$p,
    chi2_manual = man$chi2,
    z_standard = std$z,
    p_standard = std$p,
    chi2_standard = std$chi2,
    z_diff = z_diff,
    p_diff = p_diff,
    z_match = isTRUE(all.equal(man$z, std$z, tolerance = 1e-6)),
    p_match = isTRUE(all.equal(man$p, std$p, tolerance = 1e-8)),
    method_manual = man$method,
    method_standard = std$method,
    alternative_standard = std$alternative
  )
}

#' Build plot annotation using standard package; register comparison metadata.
ca_trend_annotate <- function(scores, n, success,
                             metric = NULL,
                             panel_id = NA_character_,
                             analysis = NA_character_,
                             stratification = NA_character_) {
  cmp <- ca_trend_compare(
    success = success, n = n, scores = scores,
    metric = ifelse(is.null(metric), NA_character_, as.character(metric)[1]),
    panel_id = panel_id,
    analysis = analysis,
    stratification = stratification
  )
  ca_registry_add(as.list(cmp[1, ]))
  if (!exists("sprintf_ca_trend_label", mode = "function")) {
    stop("sprintf_ca_trend_label must be defined before calling ca_trend_annotate (define in parent script).")
  }
  list(
    z = cmp$z_standard[[1]],
    p = cmp$p_standard[[1]],
    label = sprintf_ca_trend_label(cmp$z_standard[[1]], cmp$p_standard[[1]], metric),
    comparison = cmp
  )
}

ca_registry_export <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cmp <- ca_registry_as_tibble()
  if (nrow(cmp) == 0L) {
    warning("Cochran-Armitage registry is empty; no CSV exported.")
    return(invisible(cmp))
  }
  inputs <- cmp %>%
    dplyr::select(
      analysis, stratification, panel_id, metric,
      n_groups, scores, n_per_group, success_per_group, p_overall, s_bar_scores
    )
  out_cmp <- cmp %>%
    dplyr::select(
      dplyr::any_of(c(
        "analysis", "stratification", "panel_id", "metric",
        "z_manual", "p_manual", "chi2_manual",
        "z_standard", "p_standard", "chi2_standard",
        "z_diff", "p_diff", "z_match", "p_match",
        "method_manual", "method_standard", "alternative_standard"
      ))
    )
  utils::write.csv(inputs, file.path(out_dir, "cochran_armitage_trend_inputs.csv"), row.names = FALSE)
  utils::write.csv(out_cmp, file.path(out_dir, "cochran_armitage_trend_comparison.csv"), row.names = FALSE)
  saveRDS(
    list(
      comparison = cmp,
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      standard_package = "stats::prop.trend.test",
      manual_note = "Legacy in-script normal approximation (Z from numerator/denominator formula)"
    ),
    file.path(out_dir, "cochran_armitage_trend_metadata.rds")
  )
  invisible(cmp)
}

plain_plot_label_line_to_plotmath <- function(s) {
  s <- space_plot_separators(as.character(s)[1])
  if (!grepl("[Pp]\\s*[=<>]", s)) return(quote_plotmath_text(s))
  chunks <- character(0)
  rest <- s
  while (nzchar(rest) && grepl("[Pp]\\s*[=<>]", rest)) {
    loc <- regexpr("[Pp]\\s*[=<>]", rest, perl = TRUE)[1]
    if (loc > 1L) {
      chunks <- c(chunks, quote_plotmath_text(substr(rest, 1L, loc - 1L)))
    }
    rest <- substr(rest, loc, nchar(rest))
    m <- regexpr("^[Pp]\\s*(=|<)\\s*([^,\\)\\n]+)", rest, perl = TRUE)
    if (m[1L] == -1L) {
      chunks <- c(chunks, quote_plotmath_text(rest))
      break
    }
    ml <- attr(m, "match.length")
    token <- substr(rest, 1L, ml)
    if (grepl("^[Pp]\\s*<", token, perl = TRUE)) {
      rhs <- sub("^[Pp]\\s*<\\s*", "", token, ignore.case = TRUE)
      chunks <- c(chunks, paste0("italic(p) < ", num_str_to_plotmath_rhs(rhs)))
    } else {
      rhs <- sub("^[Pp]\\s*=\\s*", "", token, ignore.case = TRUE)
      chunks <- c(chunks, paste0("italic(p) == ", num_str_to_plotmath_rhs(rhs)))
    }
    rest <- substr(rest, ml + 1L, nchar(rest))
  }
  if (nzchar(rest)) chunks <- c(chunks, quote_plotmath_text(rest))
  if (length(chunks) == 1L) return(chunks[[1L]])
  paste0("paste(", paste(chunks, collapse = ", "), ")")
}

plain_plot_label_to_plotmath <- function(s) {
  lines <- strsplit(as.character(s)[1], "\n", fixed = TRUE)[[1]]
  exprs <- vapply(lines, plain_plot_label_line_to_plotmath, character(1))
  if (length(exprs) == 1L) return(exprs[[1L]])
  paste0("paste(", paste(exprs, collapse = ", '\n', "), ")")
}

label_uses_italic_p <- function(s) {
  grepl("italic\\(p\\)", plain_plot_label_to_plotmath(s), fixed = TRUE)
}

as_ggplot_label <- function(s) {
  s <- space_plot_separators(s)
  pm <- plain_plot_label_to_plotmath(s)
  if (!label_uses_italic_p(s)) return(s)
  parse(text = pm)[[1]]
}

fmt_p_plot <- function(p) {
  plain_plot_label_to_plotmath(fmt_p(p))
}

fmt_smooth_p_plot_label_plot <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  plain_plot_label_to_plotmath(fmt_smooth_p_plot_label(p, cap_below = cap_below, sci_digits = sci_digits))
}

## F22/F22b on-figure: italic(p) with literal 2-digit e-notation (not %*% 10^{} plotmath).
fmt_smooth_p_f22_plotmath <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  plain <- fmt_smooth_p_plot_label(p, cap_below = cap_below, sci_digits = sci_digits)
  if (grepl("^[Pp]\\s*<", plain)) {
    rhs <- trimws(sub("^[Pp]\\s*<\\s*", "", plain))
    return(paste0("paste(italic(p), ' < ", gsub("'", "''", rhs, fixed = TRUE), "')"))
  }
  if (grepl("^[Pp]\\s*=", plain)) {
    rhs <- trimws(sub("^[Pp]\\s*=\\s*", "", plain))
    return(paste0("paste(italic(p), ' = ", gsub("'", "''", rhs, fixed = TRUE), "')"))
  }
  quote_plotmath_text(plain)
}

sig_stars_from_p <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 5e-4) return("***")
  if (p < 5e-3) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

## F22/F22b: GAM on top, SCAM below — two annotate layers (plotmath newline is unreliable).
## top_right (SHIVA): GAM top at y = 50; SCAM one line below (F22_SHIVA_CAPTION_LINE_GAP).
F22_SHIVA_CAPTION_Y <- 50
F22_SHIVA_CAPTION_LINE_GAP <- 2.05
## F22b oxfordblue only: captions higher and 2 pt smaller than default F22b scaled size.
F22B_SHIVA_CAPTION_Y <- 56
F22B_SHIVA_CAPTION_LINE_GAP <- F22_SHIVA_CAPTION_LINE_GAP
F22B_SHIVA_CAPTION_SIZE_MULT <- 1.5 * 1.2
F22B_SHIVA_CAPTION_SIZE <- fsz(3.8 * F22B_SHIVA_CAPTION_SIZE_MULT - 2)

f22_caption_y_mid <- function(corner = c("top_right", "f22b_ranks"), y = NULL) {
  corner <- match.arg(corner)
  if (!is.null(y)) return(as.numeric(y)[1])
  if (corner == "f22b_ranks") 10 else F22_SHIVA_CAPTION_Y
}

f22_one_line_caption_layer <- function(model, edf, p_val, y, size = fsz(3.8), vjust = 0.5) {
  left <- paste0(sprintf("%-4s", model), " EDF = ", fmt_edf_plot(edf), ", ")
  lab <- paste0("paste(", quote_plotmath_text(left), ", ", fmt_smooth_p_f22_plotmath(p_val), ")")
  ggplot2::annotate(
    "text", x = Inf, y = y, label = lab,
    parse = TRUE, hjust = 1.02, vjust = vjust, size = size, color = "black"
  )
}

f22_scam_gam_caption_gam_layer <- function(
  gam_edf, gam_p,
  corner = c("top_right", "f22b_ranks"),
  size = fsz(3.8),
  y = NULL,
  line_gap = 3.5
) {
  y_mid <- f22_caption_y_mid(corner, y)
  gap <- as.numeric(line_gap)[1]
  shiva_stack <- !is.null(y) && is.finite(gap) && gap > 0 &&
    (isTRUE(all.equal(y, F22_SHIVA_CAPTION_Y)) || isTRUE(all.equal(y, F22B_SHIVA_CAPTION_Y)))
  if (shiva_stack) {
    f22_one_line_caption_layer("GAM", gam_edf, gam_p, y_mid, size, vjust = 1)
  } else if (!is.finite(gap) || gap <= 0) {
    f22_one_line_caption_layer("GAM", gam_edf, gam_p, y_mid, size, vjust = 0.5)
  } else {
    f22_one_line_caption_layer("GAM", gam_edf, gam_p, y_mid + gap / 2, size, vjust = 0.5)
  }
}

f22_scam_gam_caption_scam_layer <- function(
  scam_edf, scam_p,
  corner = c("top_right", "f22b_ranks"),
  size = fsz(3.8),
  y = NULL,
  line_gap = 3.5
) {
  y_mid <- f22_caption_y_mid(corner, y)
  gap <- as.numeric(line_gap)[1]
  shiva_stack <- !is.null(y) && is.finite(gap) && gap > 0 &&
    (isTRUE(all.equal(y, F22_SHIVA_CAPTION_Y)) || isTRUE(all.equal(y, F22B_SHIVA_CAPTION_Y)))
  if (shiva_stack) {
    f22_one_line_caption_layer("SCAM", scam_edf, scam_p, y_mid - gap, size, vjust = 1)
  } else if (!is.finite(gap) || gap <= 0) {
    f22_one_line_caption_layer("SCAM", scam_edf, scam_p, y_mid, size, vjust = 0.5)
  } else {
    f22_one_line_caption_layer("SCAM", scam_edf, scam_p, y_mid - gap / 2, size, vjust = 0.5)
  }
}

log_scam_gam_pvalues <- function(scam_out, wiggle_out, analysis_dir, pipeline_tag) {
  cap_scam <- fmt_smooth_p_plot_label(scam_out$p_value)
  cap_gam <- fmt_smooth_p_plot_label(wiggle_out$p_value)
  lines <- c(
    paste0("pipeline=", pipeline_tag),
    paste0("SCAM: p=", scam_out$p_value, " | edf=", scam_out$edf),
    paste0("GAM:  p=", wiggle_out$p_value, " | edf=", wiggle_out$edf),
    paste0("F22 caption labels: SCAM=", cap_scam, " GAM=", cap_gam)
  )
  out_path <- file.path(analysis_dir, "scam_gam_pvalues.txt")
  writeLines(lines, out_path)
  row1 <- data.frame(
    pipeline = pipeline_tag,
    scam_p = scam_out$p_value,
    scam_edf = scam_out$edf,
    gam_p = wiggle_out$p_value,
    gam_edf = wiggle_out$edf,
    scam_f22_label = cap_scam,
    gam_f22_label = cap_gam,
    stringsAsFactors = FALSE
  )
  write.csv(row1, file.path(analysis_dir, "scam_gam_pvalues.csv"), row.names = FALSE)
  for (ln in lines) log_note(paste0("SCAM/GAM p-values: ", ln))
}

## All-vs-all pairwise Fisher-test brackets (mirrors build_adjacent_prop_annotations style
## from PDX_crossover_normal_ranks_*.R, but generates every i<j pair instead of first-vs-rest).
## y positions stack above trend line (if provided) using layer_offset + span_step * k; same
## visual style as the F99D / F456 / F432 brackets in the normal-ranks file.
build_allvsall_prop_annotations <- function(
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
  n_grp <- nrow(ord)
  if (n_grp < 2L) return(data.frame())

  ## Brackets: every pair i < j (all vs all).
  comp_pairs <- expand.grid(x = seq_len(n_grp), xend = seq_len(n_grp)) %>%
    tibble::as_tibble() %>%
    dplyr::filter(x < xend) %>%
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
    tibble::tibble(
      x = i1,
      xend = i2,
      span = i2 - i1,
      p_value = p_val,
      label = sig_stars_from_p(p_val)
    )
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
  comp_df %>%
    dplyr::mutate(
      y_low = pmax(0, y - 1.8),
      x_mid = (x + xend) / 2
    )
}

add_pairwise_signif_brackets <- function(
  plot_obj,
  ann_df,
  text_size = fsz(3),
  text_y_gap = 0.0,
  line_width = 0.55,
  bracket_x_pad = 0.14,
  line_color = "black",
  text_color = "black"
) {
  if (is.null(ann_df) || nrow(ann_df) == 0L) return(plot_obj)
  label_pt <- pmax(1.2, text_size - 0.4)
  for (i in seq_len(nrow(ann_df))) {
    x_left <- ann_df$x[i] - bracket_x_pad
    x_right <- ann_df$xend[i] + bracket_x_pad
    plot_obj <- plot_obj +
      ggplot2::annotate("segment", x = x_left, xend = x_right, y = ann_df$y[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate("segment", x = x_left, xend = x_left, y = ann_df$y_low[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate("segment", x = x_right, xend = x_right, y = ann_df$y_low[i], yend = ann_df$y[i], linewidth = line_width, colour = line_color) +
      ggplot2::annotate(
        "text",
        x = ann_df$x_mid[i],
        y = ann_df$y[i] + text_y_gap,
        label = ann_df$label[i],
        size = label_pt,
        fontface = "bold",
        colour = text_color,
        vjust = 0
      )
  }
  plot_obj
}

n_workers <- parallel::detectCores(logical = TRUE)
if (is.na(n_workers) || n_workers < 1) n_workers <- 1L

run_step("Step 1: Load packages and parallel config", function() {
  required_pkgs <- c(
    "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "tibble",
    "lme4", "lmerTest", "scam", "mgcv", "writexl", "parallel",
    "patchwork", "gridExtra", "RColorBrewer", "survival", "survminer"
  )
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  log_note(paste0("R version: ", R.version.string))
  tryCatch(log_note(paste0("sessionInfo: ", paste(capture.output(sessionInfo()), collapse = "\n"))), error = function(e) log_note("sessionInfo capture skipped"))

  for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")) {
    tryCatch(do.call(Sys.setenv, setNames(list(as.character(n_workers)), v)), error = function(e) NULL)
  }
  log_note(paste0("Parallel config: using ", n_workers, " workers (all logical cores). BLAS/OpenMP threads set to ", n_workers, "."))
})

if (!isTRUE(PLOT_ONLY_MODE)) {
df <- run_step("Step 2: Read and standardize data", function() {
  df_raw <- readxl::read_xlsx(input_path)
  df_local <- df_raw %>%
    rename_with(~ stringr::str_replace_all(., "\\s+", "_")) %>%
    rename_with(~ toupper(.))

  required_candidates <- list(
    REPORT_ID = c("REPORT_ID", "PATIENT", "PATIENT_ID", "SAMPLE_ID", "MODEL", "MODEL_ID"),
    COMPOUND = c("COMPOUND", "DRUG", "COMPOUND_NAME"),
    SCORE = c("LEVEL", "SCORE", "PREDICTIVE_SCORE", "PRED_SCORE"),
    TIME_TO_DOUBLE = c(RESPONSE_COL_NAME, "TIMETODOUBLE", "TIME_TO_DOUBLE", "TIME_TO_DOUBLING", "TTD")
  )

  resolve_col <- function(candidates, df_cols) {
    hit <- candidates[candidates %in% df_cols]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  col_map <- sapply(required_candidates, resolve_col, df_cols = colnames(df_local))
  if (any(is.na(col_map))) {
    stop(paste("Missing required columns:", paste(names(col_map)[is.na(col_map)], collapse = ", ")))
  }

  has_response <- "RESPONSE_BESTAVGRESPONSE" %in% colnames(df_local)
  has_response_best <- "RESPONSE_BEST_RESPONSE" %in% colnames(df_local)
  df_local <- df_local %>%
    transmute(
      REPORT_ID = as.character(.data[[col_map["REPORT_ID"]]]),
      COMPOUND = as.character(.data[[col_map["COMPOUND"]]]),
      SCORE = as.numeric(.data[[col_map["SCORE"]]]),
      TIME_TO_DOUBLE = as.numeric(.data[[col_map["TIME_TO_DOUBLE"]]]),
      RESPONSE_BESTAVGRESPONSE = if (has_response) as.character(.data[["RESPONSE_BESTAVGRESPONSE"]]) else NA_character_,
      RESPONSE_BEST_RESPONSE = if (has_response_best) as.character(.data[["RESPONSE_BEST_RESPONSE"]]) else NA_character_
    ) %>%
    filter(!is.na(REPORT_ID), !is.na(COMPOUND), !is.na(SCORE), !is.na(TIME_TO_DOUBLE))
  log_note(paste0("Using all non-missing TimeToDouble rows: ", nrow(df_local), " rows, ", dplyr::n_distinct(df_local$REPORT_ID), " patients."))
  df_local
})
log_note(paste0("Post Step 2: ", nrow(df), " rows, ", dplyr::n_distinct(df$REPORT_ID), " patients."))


df_all_ranks <- run_step("Step 3: SHIVA ranks based on LEVEL (Low / Intermediate / High)", function() {
  ## SHIVA ranks: assign each row to Low (LEVEL < 0), Intermediate (0 <= LEVEL < 1000), or High (LEVEL >= 1000).
  ## rank_pooled: 1 = High (reference), 2 = Intermediate, 3 = Low. rank_scaled: 1 = best, 0 = worst.
  all_ranks <- df %>%
    mutate(
      LEVEL = SCORE,
      SHIVA_GROUP = dplyr::case_when(
        LEVEL < 0                       ~ SHIVA_GROUP_LOW,
        LEVEL >= 0 & LEVEL < 1000       ~ SHIVA_GROUP_INT,
        LEVEL >= 1000                   ~ SHIVA_GROUP_HIGH,
        TRUE                            ~ NA_character_
      )
    ) %>%
    filter(!is.na(SHIVA_GROUP)) %>%
    mutate(
      rank_pooled = dplyr::case_when(
        SHIVA_GROUP == SHIVA_GROUP_HIGH ~ 1L,
        SHIVA_GROUP == SHIVA_GROUP_INT  ~ 2L,
        SHIVA_GROUP == SHIVA_GROUP_LOW  ~ 3L,
        TRUE ~ NA_integer_
      ),
      rank_bin = rank_pooled,
      rank_scaled = dplyr::case_when(
        rank_pooled == 1L ~ 1,
        rank_pooled == 2L ~ 0.5,
        rank_pooled == 3L ~ 0,
        TRUE ~ NA_real_
      )
    )

  if (!"rank_pooled" %in% colnames(all_ranks)) {
    stop("Step 3 failed to create rank_pooled (SHIVA ranks).")
  }
  if (!"rank_scaled" %in% colnames(all_ranks)) {
    all_ranks <- all_ranks %>%
      mutate(rank_scaled = dplyr::case_when(rank_pooled == 1L ~ 1, rank_pooled == 2L ~ 0.5, rank_pooled == 3L ~ 0, TRUE ~ NA_real_))
  }

  coverage <- all_ranks %>%
    group_by(REPORT_ID) %>%
    summarise(n_ranks = n(), .groups = "drop")
  write.csv(coverage, file.path(analysis_dir, "rank_coverage_by_patient.csv"), row.names = FALSE)

  write.csv(all_ranks, file.path(analysis_dir, "data_all_ranks_ranks.csv"), row.names = FALSE)

  if (nrow(all_ranks) == 0) {
    stop("No data after computing SHIVA ranks. Cannot proceed.")
  }
  log_note(paste0("SHIVA ranks: ", nrow(all_ranks), " rows, ", dplyr::n_distinct(all_ranks$REPORT_ID), " patients (High = reference)."))
  all_ranks
})

## SHIVA ranks: rank_pooled 1 = High (reference), 2 = Intermediate, 3 = Low; rank_scaled 1 = best, 0 = worst.
rank_levels <- SHIVA_RANK_LEVELS
rank_labels <- SHIVA_RANK_LABELS
max_rank <- length(rank_levels)
n_unique_ranks <- max_rank
log_note(paste0("Ranks scope: ", paste(rank_labels, collapse = ", "), " (reference = High)."))

## SHIVA group colors (High, Intermediate, Low) — same hex as KM_SHIVA_final.py / KM plots
rank_colors_shiva <- c("1" = "#50C878", "2" = "#ffbc15", "3" = "#f71c6c")
names(rank_colors_shiva) <- as.character(rank_levels)


scam_out <- run_step("Step 4: SCAM monotonicity test on ranks", function() {
  # rank_scaled: 1 = High (best), 0.5 = Intermediate, 0 = Low (worst).
  scam_k <- 4L
  model_df <- df_all_ranks %>%
    mutate(
      REPORT_ID = factor(REPORT_ID),
      rank_ordered = factor(rank_pooled, levels = rank_levels, labels = rank_labels, ordered = TRUE)
    ) %>%
    as.data.frame()

  # Guard against missing/invalid rank fields in incoming data.
  if (!"rank_pooled" %in% colnames(model_df) || all(is.na(model_df$rank_pooled))) {
    stop(paste0("Step 4 input is missing rank_pooled. Columns: ", paste(colnames(model_df), collapse = ", ")))
  }
  if (!"rank_scaled" %in% colnames(model_df) || all(is.na(model_df$rank_scaled))) {
    model_df <- model_df %>%
      dplyr::mutate(rank_scaled = dplyr::case_when(rank_pooled == 1L ~ 1, rank_pooled == 2L ~ 0.5, rank_pooled == 3L ~ 0, TRUE ~ NA_real_))
    log_note("Step 4: rank_scaled was missing/invalid and was reconstructed from rank_pooled.")
  }

  model_df <- model_df %>% dplyr::filter(is.finite(rank_scaled), is.finite(TIME_TO_DOUBLE))
  build_step4_fallback <- function(reason_msg) {
    obs_rank_single <- model_df %>%
      group_by(.data$rank_pooled) %>%
      summarise(
        mean_ttd = mean(TIME_TO_DOUBLE, na.rm = TRUE),
        se_ttd = dplyr::if_else(n() > 1, sd(TIME_TO_DOUBLE, na.rm = TRUE) / sqrt(n()), 0),
        n = n(),
        .groups = "drop"
      ) %>%
      rename(rank_int = rank_pooled)
    if (nrow(obs_rank_single) == 0L) {
      stop("Step 4 failed: no finite TIME_TO_DOUBLE rows after filtering.")
    }
    overall_mean <- mean(model_df$TIME_TO_DOUBLE, na.rm = TRUE)
    rank_scaled_full <- c(1, 0.5, 0)[rank_levels]
    pred_grid_ranks_single <- tibble::tibble(
      rank_int = rank_levels,
      rank_scaled = rank_scaled_full,
      pred_ttd = rep(overall_mean, length(rank_levels)),
      se = 0,
      ci_low = rep(overall_mean, length(rank_levels)),
      ci_high = rep(overall_mean, length(rank_levels))
    )
    pred_grid_single <- tibble::tibble(
      rank_scaled = seq(min(rank_scaled_full), max(rank_scaled_full), length.out = 50),
      pred_ttd = rep(overall_mean, 50),
      se = 0,
      ci_low = rep(overall_mean, 50),
      ci_high = rep(overall_mean, 50)
    )
    log_note(reason_msg)
    list(
      model = NULL,
      model_df = model_df,
      s_table = tibble::tibble(edf = NA_real_, p.value = NA_real_, term = "s(rank_scaled)"),
      edf = NA_real_,
      p_value = NA_real_,
      pred_grid = pred_grid_single,
      pred_grid_ranks = pred_grid_ranks_single,
      obs_rank = obs_rank_single
    )
  }

  n_unique_rank_scaled <- dplyr::n_distinct(model_df$rank_scaled)
  if (n_unique_rank_scaled < 2L) {
    return(build_step4_fallback("Step 4: only one SHIVA level observed; SCAM smooth skipped and flat fallback used."))
  }

  # For 3 discrete SHIVA groups, scam::scam with bs="mpi" is unstable at k=3
  # (can throw spline range errors). k>=4 is empirically stable here.
  scam_k <- max(4L, min(8L, n_unique_rank_scaled + 1L))
  model <- tryCatch(
    scam::scam(
      TIME_TO_DOUBLE ~ s(rank_scaled, bs = "mpi", k = scam_k) + s(REPORT_ID, bs = "re"),
      data = model_df
    ),
    error = function(e) e
  )
  if (inherits(model, "error")) {
    return(build_step4_fallback(paste0("Step 4: SCAM fit failed; fallback used. Reason: ", conditionMessage(model))))
  }
  sm <- summary(model)
  capture.output(sm, file = file.path(analysis_dir, "model_scam_all_ranks_summary.txt"))

  st_raw <- sm$s.table
  if (is.null(st_raw) || (is.matrix(st_raw) && (nrow(st_raw) == 0L || ncol(st_raw) == 0L))) {
    edf <- NA_real_
    p_val <- NA_real_
    s_tab <- tibble::tibble(edf = edf, p.value = p_val, term = "s(rank_scaled)")
  } else {
    s_tab <- as.data.frame(st_raw)
    s_tab$term <- rownames(s_tab)
    rownames(s_tab) <- NULL
    nms <- names(s_tab)
    if (!is.null(nms) && length(nms) == ncol(s_tab)) {
      names(s_tab) <- make.names(nms)
    }
    edf <- if ("edf" %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab$edf[1] else NA_real_
    pcol <- grep("^p", names(s_tab), ignore.case = TRUE, value = TRUE)[1]
    p_val <- if (!is.na(pcol) && pcol %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab[[pcol]][1] else NA_real_
  }

  rank_scaled_obs_min <- min(model_df$rank_scaled, na.rm = TRUE)
  rank_scaled_obs_max <- max(model_df$rank_scaled, na.rm = TRUE)
  pred_seq <- seq(rank_scaled_obs_min, rank_scaled_obs_max, length.out = 50)
  pred_grid <- tibble::tibble(
    rank_scaled = pred_seq,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred <- tryCatch(
    predict(model, newdata = pred_grid, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)"),
    error = function(e) e
  )
  if (inherits(pred, "error")) {
    return(build_step4_fallback(paste0("Step 4: SCAM predict(grid) failed; fallback used. Reason: ", conditionMessage(pred))))
  }
  pred_grid <- pred_grid %>%
    mutate(
      pred_ttd = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  rank_scaled_mid <- c(1, 0.5, 0)[rank_levels]
  rank_scaled_mid <- pmin(rank_scaled_obs_max, pmax(rank_scaled_obs_min, rank_scaled_mid))
  pred_grid_ranks <- tibble::tibble(
    rank_scaled = rank_scaled_mid,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred_dec <- tryCatch(
    predict(model, newdata = pred_grid_ranks, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)"),
    error = function(e) e
  )
  if (inherits(pred_dec, "error")) {
    return(build_step4_fallback(paste0("Step 4: SCAM predict(ranks) failed; fallback used. Reason: ", conditionMessage(pred_dec))))
  }
  pred_grid_ranks <- pred_grid_ranks %>%
    mutate(
      rank_int = rank_levels,
      pred_ttd = as.numeric(pred_dec$fit),
      se = as.numeric(pred_dec$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  obs_rank <- model_df %>%
    group_by(.data$rank_pooled) %>%
    summarise(
      mean_ttd = mean(TIME_TO_DOUBLE, na.rm = TRUE),
      se_ttd = sd(TIME_TO_DOUBLE, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    rename(rank_int = rank_pooled)

  fig_w <- FIG_SCAM_MONO_WIDTH
  fig_h <- FIG_SCAM_MONO_HEIGHT
  plot_mono <- ggplot2::ggplot(obs_rank, ggplot2::aes(factor(rank_int, levels = rank_levels, labels = rank_labels), mean_ttd)) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd, ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd),
      data = obs_rank, width = 0.15
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd, label = paste0("(n = ", n, ")")),
      data = obs_rank, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(2.8)
    ) +
    ggplot2::geom_line(data = pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd), color = "#1f78b4", linewidth = 1, group = 1) +
    ggplot2::geom_ribbon(data = pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high),
                         inherit.aes = FALSE, alpha = 0.2, fill = "#1f78b4", group = 1) +
    ggplot2::theme_minimal(base_size = fsz(12)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 30, r = 5, b = 5, l = 5)) +
    ggplot2::labs(
      title = as_ggplot_label(paste0("SCAM monotonicity (SHIVA groups, reference = High): EDF = ", round(edf, 3), ", ", fmt_p(p_val))),
      x = "SHIVA group",
      y = paste("Avg", RESPONSE_COL_LABEL)
    )

  list(
    model = model,
    model_df = model_df,
    s_table = s_tab,
    edf = edf,
    p_value = p_val,
    pred_grid = pred_grid,
    pred_grid_ranks = pred_grid_ranks,
    obs_rank = obs_rank
  )
})

wiggle_out <- run_step("Step 4a: Unconstrained GAM visualization (wiggly reference)", function() {
  model_df <- df_all_ranks %>%
    mutate(REPORT_ID = factor(REPORT_ID)) %>%
    as.data.frame()

  model_df <- model_df %>% dplyr::filter(is.finite(rank_scaled), is.finite(TIME_TO_DOUBLE))
  build_step4a_fallback <- function(reason_msg) {
    pred_grid_gam <- scam_out$pred_grid_ranks %>%
      dplyr::select(rank_scaled, rank_int, pred_ttd, se, ci_low, ci_high)
    log_note(reason_msg)
    list(
      model = NULL,
      edf = NA_real_,
      p_value = NA_real_,
      pred_grid = pred_grid_gam,
      result_label = reason_msg
    )
  }

  n_unique_rank_scaled <- dplyr::n_distinct(model_df$rank_scaled)
  if (n_unique_rank_scaled < 2L) {
    return(build_step4a_fallback("Inconclusive: only one SHIVA level observed; GAM smooth skipped and SCAM fallback grid reused."))
  }

  # With discrete SHIVA levels, k must not exceed unique x-levels.
  gam_k <- max(2L, min(10L, n_unique_rank_scaled))
  gam_model <- tryCatch(
    mgcv::gam(
      TIME_TO_DOUBLE ~ s(rank_scaled, bs = "tp", k = gam_k) + s(REPORT_ID, bs = "re"),
      data = model_df,
      method = "REML"
    ),
    error = function(e) e
  )
  if (inherits(gam_model, "error")) {
    return(build_step4a_fallback(paste0("Inconclusive: GAM fit failed; SCAM fallback grid reused. Reason: ", conditionMessage(gam_model))))
  }
  gam_sm <- summary(gam_model)
  capture.output(gam_sm, file = file.path(analysis_dir, "model_gam_unconstrained_all_ranks_summary.txt"))

  gam_st <- gam_sm$s.table
  gam_edf <- NA_real_
  gam_p <- NA_real_
  if (!is.null(gam_st) && nrow(gam_st) >= 1L && ncol(gam_st) >= 1L) {
    gam_edf <- if ("edf" %in% colnames(gam_st)) gam_st[1, "edf"] else NA_real_
    pcol <- grep("^p", colnames(gam_st), ignore.case = TRUE, value = TRUE)[1]
    gam_p <- if (!is.na(pcol) && pcol %in% colnames(gam_st)) gam_st[1, pcol] else NA_real_
  }

  gam_result_label <- if (is.na(gam_p) || is.na(gam_edf)) {
    "Inconclusive: GAM summary statistics unavailable."
  } else if (gam_p < 0.05 && gam_edf > 1.1) {
    "Significant nonlinearity: unconstrained GAM indicates a wiggly rank-response pattern."
  } else if (gam_p < 0.05 && gam_edf <= 1.1) {
    "Significant but near-linear trend: unconstrained GAM supports rank association with limited wiggle."
  } else {
    "No strong evidence of nonlinear rank effect in the unconstrained GAM."
  }

  rank_scaled_obs_min <- min(model_df$rank_scaled, na.rm = TRUE)
  rank_scaled_obs_max <- max(model_df$rank_scaled, na.rm = TRUE)
  rank_scaled_mid <- if (max_rank > 1L) 1 - (rank_levels - 1) / (max_rank - 1) else rep(1, length(rank_levels))
  rank_scaled_mid <- pmin(rank_scaled_obs_max, pmax(rank_scaled_obs_min, rank_scaled_mid))
  pred_grid_gam <- tibble::tibble(
    rank_scaled = rank_scaled_mid,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred_gam <- tryCatch(
    predict(gam_model, newdata = pred_grid_gam, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)"),
    error = function(e) e
  )
  if (inherits(pred_gam, "error")) {
    return(build_step4a_fallback(paste0("Inconclusive: GAM predict failed; SCAM fallback grid reused. Reason: ", conditionMessage(pred_gam))))
  }
  pred_grid_gam <- pred_grid_gam %>%
    mutate(
      rank_int = rank_levels,
      pred_ttd = as.numeric(pred_gam$fit),
      se = as.numeric(pred_gam$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  fig_w <- FIG_F22_JPG_WIDTH
  fig_h <- FIG_F22_JPG_HEIGHT
  tryCatch({
    obs_rank_4a <- scam_out$obs_rank
    y_min_4a <- min(obs_rank_4a$mean_ttd - 1.96 * obs_rank_4a$se_ttd, scam_out$pred_grid_ranks$ci_low, pred_grid_gam$ci_low, na.rm = TRUE)
    y_max_4a <- max(obs_rank_4a$mean_ttd + 1.96 * obs_rank_4a$se_ttd, scam_out$pred_grid_ranks$ci_high, pred_grid_gam$ci_high, na.rm = TRUE)
    y_pad_4a <- 0.16 * (y_max_4a - y_min_4a + 1e-8)
    y_lab_off_4a <- 0.04 * (y_max_4a - y_min_4a + 1e-8)
    p_wiggle_base <- ggplot2::ggplot(obs_rank_4a, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
      ggplot2::geom_point(size = 1.8, color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd), width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd + y_lab_off_4a, label = paste0("(n = ", n, ")")),
        data = obs_rank_4a, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(3.7), color = "black"
      ) +
      ggplot2::geom_ribbon(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
        inherit.aes = FALSE, alpha = 0.12) +
      ggplot2::geom_line(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"), linewidth = 1.05) +
      ggplot2::geom_ribbon(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
        inherit.aes = FALSE, alpha = 0.08) +
      ggplot2::geom_line(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"),
        linewidth = 1.05, linetype = "22") +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = "#2166ac", "GAM" = "#b2182b")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = "#2166ac", "GAM" = "#b2182b"), guide = "none") +
      ggplot2::theme_minimal(base_size = fsz(15)) +
      ggplot2::theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.text = ggplot2::element_text(color = "black", size = fsz(12)),
        plot.title = ggplot2::element_text(color = "black", face = "bold", size = fsz(16)),
        plot.subtitle = ggplot2::element_text(color = "black", size = fsz(13)),
        axis.title = ggplot2::element_text(color = "black", size = fsz(14)),
        axis.text = ggplot2::element_text(color = "black", size = fsz(12)),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.ticks.length = grid::unit(0.2, "cm"),
        plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 10)),
        plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 12)),
        plot.margin = ggplot2::margin(t = 30, r = 8, b = 8, l = 8)
      ) +
      ggplot2::labs(
        title = NULL,
        subtitle = NULL,
        x = "DDA score tiers",
        y = "Avg PFS (days)"
      ) +
      ggplot2::coord_cartesian(ylim = c(y_min_4a, y_max_4a + y_pad_4a), clip = "off")
    p_wiggle <- p_wiggle_base +
      f22_scam_gam_caption_gam_layer(gam_edf, gam_p, corner = "top_right", y = F22_SHIVA_CAPTION_Y, line_gap = F22_SHIVA_CAPTION_LINE_GAP, size = fsz(3.8)) +
      f22_scam_gam_caption_scam_layer(scam_out$edf, scam_out$p_value, corner = "top_right", y = F22_SHIVA_CAPTION_Y, line_gap = F22_SHIVA_CAPTION_LINE_GAP, size = fsz(3.8))

    p_wiggle_ox <- p_wiggle_base +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4"), guide = "none") +
      f22_scam_gam_caption_gam_layer(
        gam_edf, gam_p, corner = "top_right",
        y = F22B_SHIVA_CAPTION_Y, line_gap = F22B_SHIVA_CAPTION_LINE_GAP, size = F22B_SHIVA_CAPTION_SIZE
      ) +
      f22_scam_gam_caption_scam_layer(
        scam_out$edf, scam_out$p_value, corner = "top_right",
        y = F22B_SHIVA_CAPTION_Y, line_gap = F22B_SHIVA_CAPTION_LINE_GAP, size = F22B_SHIVA_CAPTION_SIZE
      )
  }, error = function(e) {
    log_note(paste0("Step 4a: Skipped SCAM vs GAM plot (plot/ggsave error): ", conditionMessage(e)))
  })
  write.csv(pred_grid_gam, file.path(analysis_dir, "gam_unconstrained_pred_by_rank.csv"), row.names = FALSE)
  writeLines(
    c(
      paste0("Timestamp: ", timestamp_now()),
      paste0("GAM EDF: ", round(gam_edf, 4)),
      paste0("GAM p-value: ", format.pval(gam_p, digits = 4, eps = 1e-4))
    ),
    file.path(analysis_dir, "gam_unconstrained_interpretation.txt")
  )

  list(
    model = gam_model,
    edf = gam_edf,
    p_value = gam_p,
    pred_grid = pred_grid_gam
  )
})

log_scam_gam_pvalues(scam_out, wiggle_out, analysis_dir, "SHIVA_tier_data_modeling")

## Step 4b inactivated for a while (do not delete - will be needed shortly).
# scam_sensitivity <- run_step("Step 4b: SCAM k-sensitivity (multiple k in parallel)", function() {
#   model_df <- df_all_ranks %>%
#     mutate(REPORT_ID = factor(REPORT_ID)) %>%
#     as.data.frame()
#
#   fit_scam_k <- function(k_val) {
#     m <- scam::scam(
#       TIME_TO_DOUBLE ~ s(rank_scaled, bs = "mpi", k = k_val) + s(REPORT_ID, bs = "re"),
#       data = model_df
#     )
#     st <- summary(m)$s.table
#     edf_k <- NA_real_
#     pval <- NA_real_
#     if (!is.null(st) && nrow(st) >= 1L && ncol(st) >= 1L) {
#       edf_k <- if ("edf" %in% colnames(st)) st[1, "edf"] else NA_real_
#       pcol <- grep("^p", colnames(st), ignore.case = TRUE, value = TRUE)[1]
#       pval <- if (!is.na(pcol) && pcol %in% colnames(st)) st[1, pcol] else NA_real_
#     }
#     tibble::tibble(k = k_val, edf = edf_k, p_value = pval)
#   }
#
#   ## Use only k=7 for speed; lower k = more conservative (rigid) smooth = more trustworthy
#   k_vals <- c(min(7L, n_unique_ranks))
#   n_k <- min(n_workers, length(k_vals))
#   if (n_k > 1) {
#     cl <- parallel::makeCluster(n_k)
#     on.exit(parallel::stopCluster(cl), add = TRUE)
#     parallel::clusterExport(cl, "model_df", envir = environment())
#     parallel::clusterEvalQ(cl, library(scam))
#     sensitivity_list <- parallel::parLapply(cl, k_vals, fit_scam_k)
#   } else {
#     sensitivity_list <- lapply(k_vals, fit_scam_k)
#   }
#   sens_df <- dplyr::bind_rows(sensitivity_list)
#   write.csv(sens_df, file.path(analysis_dir, "scam_k_sensitivity.csv"), row.names = FALSE)
#   log_note(paste0("SCAM k-sensitivity: ", paste(sprintf("k=%d EDF=%.3f p=%s", sens_df$k, sens_df$edf, format.pval(sens_df$p_value, digits = 2)), collapse = "; ")))
#   sens_df
# })
scam_sensitivity <- tibble::tibble(k = NA_integer_, edf = NA_real_, p_value = NA_real_)

## Mixed model (publication): TIME_TO_DOUBLE ~ rank_f + (1|REPORT_ID); High tier = reference.
## P-values: lmerTest summary (Satterthwaite). Forest CIs: confint(profile), Wald if profiling fails.
effects_out <- run_step("Mixed-model effect sizes", function() {
  df_eff <- df_all_ranks %>% mutate(rank_f = factor(rank_pooled, levels = rank_levels))
  model <- lmerTest::lmer(TIME_TO_DOUBLE ~ rank_f + (1 | REPORT_ID), data = df_eff)
  sm <- summary(model)
  capture.output(sm, file = file.path(analysis_dir, "model_lmer_all_ranks_summary.txt"))

  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>%
    filter(grepl("^rank_f", term)) %>%
    mutate(
      rank_int = suppressWarnings(as.integer(gsub("^rank_f", "", term))),
      estimate = Estimate,
      std_error = `Std. Error`,
      p_value = `Pr(>|t|)`
    ) %>%
    select(rank_int, estimate, std_error, p_value)

  ci <- tryCatch(
    suppressMessages(confint(model, method = "profile")),
    error = function(e) {
      log_note(paste0("Profile CIs failed (", conditionMessage(e), "); falling back to Wald."))
      suppressMessages(confint(model, method = "Wald"))
    }
  )
  ci_df <- as.data.frame(ci)
  ci_df$term <- rownames(ci_df)
  rownames(ci_df) <- NULL
  cn <- colnames(ci_df)
  cilo <- cn[grep("2\\.5|lower", cn, ignore.case = TRUE)[1]]
  cihi <- cn[grep("97\\.5|upper", cn, ignore.case = TRUE)[1]]
  if (is.na(cilo)) cilo <- cn[1]
  if (is.na(cihi)) cihi <- cn[2]
  ci_df <- ci_df %>%
    filter(grepl("^rank_f", term)) %>%
    mutate(rank_int = as.integer(gsub("^rank_f", "", term)), ci_low = .data[[cilo]], ci_high = .data[[cihi]]) %>%
    select(rank_int, ci_low, ci_high)

  effects <- coef_df %>%
    left_join(ci_df, by = "rank_int") %>%
    bind_rows(tibble::tibble(rank_int = 1L, estimate = 0, std_error = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_)) %>%
    arrange(rank_int) %>%
    mutate(
      p_label = ifelse(is.na(p_value), "ref", fmt_p(p_value)),
      p_label_plot = ifelse(is.na(p_value), "ref", fmt_p_plot(p_value))
    )
  idx_eff <- !is.na(effects$p_value)
  effects$p_adj_fdr <- NA_real_
  effects$p_adj_bonf <- NA_real_
  if (any(idx_eff)) {
    effects$p_adj_fdr[idx_eff] <- p.adjust(effects$p_value[idx_eff], method = "BH")
    effects$p_adj_bonf[idx_eff] <- p.adjust(effects$p_value[idx_eff], method = "bonferroni")
  }

  write.csv(effects, file.path(analysis_dir, "effect_sizes_all_ranks_vs_top_rank.csv"), row.names = FALSE)

  x_max <- max(effects$ci_high, na.rm = TRUE)
  x_pad <- 0.06 * (max(effects$ci_high, na.rm = TRUE) - min(effects$ci_low, na.rm = TRUE) + 1e-8)
  eff_height <- max(6, 4 + max_rank * 0.2)
  label_size_eff <- max(2, min(3.5, 3.5 - max_rank * 0.02))
  p_fx <- ggplot2::ggplot(effects, ggplot2::aes(x = estimate, y = factor(rank_int, levels = rev(rank_levels), labels = rev(rank_labels)))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0.2, orientation = "y", color = "#1f78b4") +
    ggplot2::geom_point(ggplot2::aes(color = rank_int == 1L), size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(x = x_max + x_pad, label = p_label_plot), parse = TRUE, hjust = 0, size = label_size_eff) +
    ggplot2::scale_color_manual(values = c("FALSE" = "#1f78b4", "TRUE" = "black"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Effect sizes by Tiers (mixed model)",
      x = "Estimated TimeToDouble difference vs High (reference)",
      y = "SHIVA group"
    ) +
    ggplot2::coord_cartesian(xlim = c(min(effects$ci_low, na.rm = TRUE), x_max + 3 * x_pad))

  effects
})

if (isTRUE(SAVE_PLOT_DATA_CACHE)) {
  run_step("Step 8a: Save cached data bundle for plotting-only reruns", function() {
    plot_cache <- list(
      df = as.data.frame(df),
      df_all_ranks = as.data.frame(df_all_ranks),
      scam_out = scam_out,
      wiggle_out = wiggle_out,
      effects_out = as.data.frame(effects_out),
      rank_levels = rank_levels,
      rank_labels = rank_labels,
      rank_colors_shiva = rank_colors_shiva,
      max_rank = max_rank
    )
    saveRDS(plot_cache, plot_data_cache_path)
    log_note(paste0("Saved plot data cache: ", plot_data_cache_path))
  })
}
} else {
  run_step("Step 2-8a bypass: Load cached data bundle for plotting", function() {
    if (!file.exists(plot_data_cache_path)) {
      stop(paste0("PLOT_ONLY_MODE=TRUE but cache not found: ", plot_data_cache_path,
                  ". Run once with PLOT_ONLY_MODE=FALSE and SAVE_PLOT_DATA_CACHE=TRUE."))
    }
    plot_cache <- readRDS(plot_data_cache_path)
    validate_plot_cache(plot_cache)
    assign("df", as.data.frame(plot_cache$df), envir = .GlobalEnv)
    assign("df_all_ranks", as.data.frame(plot_cache$df_all_ranks), envir = .GlobalEnv)
    assign("scam_out", plot_cache$scam_out, envir = .GlobalEnv)
    assign("wiggle_out", plot_cache$wiggle_out, envir = .GlobalEnv)
    assign("effects_out", as.data.frame(plot_cache$effects_out), envir = .GlobalEnv)
    assign("rank_levels", as.integer(plot_cache$rank_levels), envir = .GlobalEnv)
    if ("rank_labels" %in% names(plot_cache)) {
      assign("rank_labels", as.character(plot_cache$rank_labels), envir = .GlobalEnv)
    } else {
      assign("rank_labels", as.character(plot_cache$rank_levels), envir = .GlobalEnv)
    }
    assign("max_rank", as.integer(plot_cache$max_rank), envir = .GlobalEnv)
    assign("n_unique_ranks", length(as.integer(plot_cache$rank_levels)), envir = .GlobalEnv)
    if ("rank_colors_shiva" %in% names(plot_cache)) {
      assign("rank_colors_shiva", plot_cache$rank_colors_shiva, envir = .GlobalEnv)
    } else {
      assign("rank_colors_shiva", setNames(c("#50C878", "#ffbc15", "#f71c6c"), as.character(plot_cache$rank_levels)), envir = .GlobalEnv)
    }
    log_note(paste0("Loaded plot data cache: ", plot_data_cache_path))
  })
}

run_step("Step 9: Publication-quality figures (peer review & PI perspective)", function() {
  ca_registry_reset()

  save_fig <- function(plot_obj, path, width, height, dpi = 300) {
    tryCatch(
      {
        # Save PNG (or whatever extension is in 'path')
        ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
        log_note(paste0("Saved: ", basename(path)))

        # Also save an SVG version in the same folder (figures)
        base_no_ext <- sub("\\.[^.]*$", "", basename(path))
        svg_path <- file.path(dirname(path), paste0(base_no_ext, ".svg"))
        ggplot2::ggsave(svg_path, plot_obj, width = width, height = height, dpi = 320, device = "svg")
        log_note(paste0("Saved: ", basename(svg_path)))
      },
      error = function(e) {
        log_note(paste0("SKIP figure ", basename(path), " (error): ", conditionMessage(e)))
      }
    )
  }

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

  ## mRECIST / Durable Benefit: horizontal y-axis aid lines only (no vertical grid).
  theme_pub_nogrid <- theme_pub +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "grey85")
    )

  # Keep F23 panel text at current approved sizes.
  theme_f23_lock <- ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = 13, color = "black"),
    axis.title = ggplot2::element_text(size = 14, color = "black"),
    axis.text = ggplot2::element_text(size = 12, color = "black"),
    legend.text = ggplot2::element_text(size = 12, color = "black"),
    legend.title = ggplot2::element_text(size = 12, color = "black")
  )


  ## SCAM residual diagnostics
  res_scam <- residuals(scam_out$model)
  pred_scam <- fitted(scam_out$model)
  res_df <- tibble::tibble(fitted = pred_scam, residual = res_scam)
  qq_df <- tibble::tibble(residual = res_scam)
  p6a <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = residual)) +
    ggplot2::stat_qq(alpha = 0.5) +
    ggplot2::stat_qq_line(color = "#b2182b", linewidth = 0.8) +
    ggplot2::labs(title = "A. Normal Q-Q plot", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_pub
  p6b <- ggplot2::ggplot(res_df, ggplot2::aes(fitted, residual)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_smooth(se = TRUE, color = "#b2182b", alpha = 0.2) +
    ggplot2::labs(title = "B. Residuals vs fitted", x = paste("Fitted", RESPONSE_COL_LABEL), y = "Residual") +
    theme_pub
  p6 <- p6a + p6b + patchwork::plot_layout(ncol = 2)

  ## SCAM residuals vs fitted (points colored by COMPOUND | GROUP, side by side)
  md_f6b <- scam_out$model_df
  if (is.null(md_f6b) || !is.data.frame(md_f6b)) {
    md_f6b <- df_all_ranks %>%
      dplyr::mutate(
        REPORT_ID = factor(REPORT_ID),
        rank_ordered = factor(rank_pooled, levels = rank_levels, labels = rank_labels, ordered = TRUE)
      ) %>%
      as.data.frame()
    if (!"rank_scaled" %in% colnames(md_f6b) || all(is.na(md_f6b$rank_scaled))) {
      md_f6b <- md_f6b %>%
        dplyr::mutate(rank_scaled = dplyr::case_when(rank_pooled == 1L ~ 1, rank_pooled == 2L ~ 0.5, rank_pooled == 3L ~ 0, TRUE ~ NA_real_))
    }
    md_f6b <- md_f6b %>% dplyr::filter(is.finite(rank_scaled), is.finite(TIME_TO_DOUBLE))
    log_note("SCAM residuals panel: scam_out$model_df missing; reconstructed from df_all_ranks (verify row order matches SCAM fit).")
  }
  res_f6b <- residuals(scam_out$model)
  fit_f6b <- fitted(scam_out$model)
  n_res_f6b <- length(res_f6b)
  if (nrow(md_f6b) != n_res_f6b) {
    log_note(paste0("SCAM residuals panel: nrow(model_df) (", nrow(md_f6b), ") != length(residuals) (", n_res_f6b, "); subsetting to complete formula rows."))
    md_f6b <- md_f6b[
      stats::complete.cases(md_f6b[["TIME_TO_DOUBLE"]], md_f6b[["rank_scaled"]], md_f6b[["REPORT_ID"]]),
      ,
      drop = FALSE
    ]
  }
  if (nrow(md_f6b) != n_res_f6b) {
    log_note(paste0("SCAM residuals panel skipped: cannot align model_df (n=", nrow(md_f6b), ") with SCAM residuals (n=", n_res_f6b, ")."))
    p6b_residual_panels <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "SCAM residuals panel unavailable (row alignment)", size = fsz(5)) +
      ggplot2::theme_void()
  } else {
    use_explicit_group_col <- "GROUP" %in% colnames(md_f6b)
    grp_f6b <- if (use_explicit_group_col) md_f6b[["GROUP"]] else md_f6b[["SHIVA_GROUP"]]
    res_df_f6b <- tibble::tibble(
      fitted = fit_f6b,
      residual = res_f6b,
      COMPOUND = md_f6b[["COMPOUND"]],
      GROUP = grp_f6b
    ) %>%
      dplyr::mutate(COMPOUND = factor(COMPOUND))
    if (use_explicit_group_col) {
      res_df_f6b <- res_df_f6b %>% dplyr::mutate(GROUP = factor(GROUP))
    } else {
      res_df_f6b <- res_df_f6b %>% dplyr::mutate(GROUP = factor(GROUP, levels = rank_labels))
    }
    n_comp <- nlevels(res_df_f6b$COMPOUND)
    leg_ncol_comp <- max(1L, min(4L, as.integer(ceiling(sqrt(n_comp)))))
    p6b_by_compound <- ggplot2::ggplot(res_df_f6b, ggplot2::aes(x = fitted, y = residual)) +
      ggplot2::geom_point(ggplot2::aes(color = COMPOUND), alpha = 0.65, size = 1.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_smooth(se = TRUE, color = "#b2182b", alpha = 0.2, linewidth = 0.8) +
      ggplot2::labs(
        title = "A. Residuals vs fitted (by COMPOUND)",
        x = paste("Fitted", RESPONSE_COL_LABEL),
        y = "Residual",
        color = "COMPOUND"
      ) +
      theme_pub +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = leg_ncol_comp, byrow = TRUE, override.aes = list(size = 2.5, alpha = 1))) +
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold"))
    p6b_by_group <- ggplot2::ggplot(res_df_f6b, ggplot2::aes(x = fitted, y = residual)) +
      ggplot2::geom_point(ggplot2::aes(color = GROUP), alpha = 0.65, size = 1.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_smooth(se = TRUE, color = "#b2182b", alpha = 0.2, linewidth = 0.8) +
      ggplot2::labs(
        title = "B. Residuals vs fitted (by GROUP)",
        x = paste("Fitted", RESPONSE_COL_LABEL),
        y = "Residual",
        color = "GROUP"
      ) +
      theme_pub +
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold"))
    km_group_colors_f6b <- c(
      "High" = "#50C878",
      "Intermediate" = "#ffbc15",
      "Int" = "#ffbc15",
      "Low" = "#f71c6c",
      "1" = "#50C878",
      "2" = "#ffbc15",
      "3" = "#f71c6c"
    )
    group_levels_f6b <- levels(res_df_f6b$GROUP)
    group_colors_f6b <- km_group_colors_f6b[intersect(names(km_group_colors_f6b), group_levels_f6b)]
    if (length(group_colors_f6b) > 0L) {
      p6b_by_group <- p6b_by_group + ggplot2::scale_color_manual(values = group_colors_f6b, drop = FALSE, na.value = "gray50")
    } else {
      p6b_by_group <- p6b_by_group + ggplot2::scale_color_discrete(na.value = "gray50")
    }
    p6b_residual_panels <- p6b_by_compound + p6b_by_group + patchwork::plot_layout(ncol = 2)
  }

  ## F22: SCAM vs GAM with CI (standalone, blue-red and oxfordblue-purple)
  obs_rank_f22 <- scam_out$obs_rank
  y_min_f22 <- min(obs_rank_f22$mean_ttd - 1.96 * obs_rank_f22$se_ttd, scam_out$pred_grid_ranks$ci_low, wiggle_out$pred_grid$ci_low, na.rm = TRUE)
  y_max_f22 <- max(obs_rank_f22$mean_ttd + 1.96 * obs_rank_f22$se_ttd, scam_out$pred_grid_ranks$ci_high, wiggle_out$pred_grid$ci_high, na.rm = TRUE)
  y_pad_f22 <- 0.16 * (y_max_f22 - y_min_f22 + 1e-8)
  y_lab_off_f22 <- 0.04 * (y_max_f22 - y_min_f22 + 1e-8)
  make_f22_plot <- function(
    scam_color, gam_color, text_mult = 1, line_width_mult = 1, body_text_extra = 1,
    caption_corner = "top_right", ylim_lo = NULL, ylim_hi = NULL,
    caption_y = NULL, caption_gap = NULL, caption_size = NULL
  ) {
    cap_y <- if (!is.null(caption_y)) {
      caption_y
    } else if (identical(caption_corner, "top_right")) {
      F22_SHIVA_CAPTION_Y
    } else {
      NULL
    }
    cap_gap <- if (!is.null(caption_gap)) {
      caption_gap
    } else if (identical(caption_corner, "top_right")) {
      F22_SHIVA_CAPTION_LINE_GAP * text_mult * body_text_extra
    } else {
      3.5 * text_mult * body_text_extra
    }
    cap_size <- if (!is.null(caption_size)) caption_size else fsz(3.8 * text_mult * body_text_extra)
    y_lo <- if (is.null(ylim_lo)) y_min_f22 else ylim_lo
    y_hi <- if (is.null(ylim_hi)) y_max_f22 + y_pad_f22 else ylim_hi
    ggplot2::ggplot(obs_rank_f22, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
      ggplot2::geom_point(size = 1.8, color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd), width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd + y_lab_off_f22, label = paste0("(n = ", n, ")")),
        data = obs_rank_f22, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(3.7 * text_mult * body_text_extra), color = "black"
      ) +
      ggplot2::geom_ribbon(
        data = scam_out$pred_grid_ranks,
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
        inherit.aes = FALSE, alpha = 0.12
      ) +
      ggplot2::geom_line(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"), linewidth = 1.05 * line_width_mult) +
      ggplot2::geom_ribbon(
        data = wiggle_out$pred_grid,
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
        inherit.aes = FALSE, alpha = 0.08
      ) +
      ggplot2::geom_line(data = wiggle_out$pred_grid, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"), linewidth = 1.05 * line_width_mult, linetype = "22") +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = scam_color, "GAM" = gam_color)) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = scam_color, "GAM" = gam_color), guide = "none") +
      f22_scam_gam_caption_gam_layer(
        wiggle_out$edf, wiggle_out$p_value,
        corner = caption_corner, y = cap_y, size = cap_size,
        line_gap = cap_gap
      ) +
      f22_scam_gam_caption_scam_layer(
        scam_out$edf, scam_out$p_value,
        corner = caption_corner, y = cap_y, size = cap_size,
        line_gap = cap_gap
      ) +
      ggplot2::labs(
        title = NULL,
        subtitle = NULL,
        x = "DDA score tiers",
        y = "Avg PFS (days)"
      ) +
      theme_pub +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = fsz(16 * text_mult), color = "black"),
        plot.subtitle = ggplot2::element_text(size = fsz(13 * text_mult), color = "black"),
        axis.title = ggplot2::element_text(size = fsz(14 * text_mult * body_text_extra), color = "black"),
        axis.text = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.text = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.title = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      ggplot2::coord_cartesian(ylim = c(y_lo, y_hi), clip = "off", expand = FALSE)
  }
  f22b_text_mult <- 1.5 * 1.2
  f22b_n_top <- max(
    obs_rank_f22$mean_ttd + 1.96 * obs_rank_f22$se_ttd + y_lab_off_f22 * f22b_text_mult,
    na.rm = TRUE
  )
  f22b_y_lo <- max(18, y_min_f22 - 1)
  f22b_y_hi <- max(
    y_max_f22 + y_pad_f22,
    F22B_SHIVA_CAPTION_Y + F22B_SHIVA_CAPTION_LINE_GAP * f22b_text_mult + 6,
    f22b_n_top + F22B_NLABEL_Y_PAD
  )
  p22_blue_red <- make_f22_plot("#2166ac", "#b2182b")

  p22_ox_purple <- make_f22_plot(
    "#002147", "#C9B6E4", text_mult = 1.5, line_width_mult = 2, body_text_extra = 1.2,
    caption_corner = "top_right", ylim_lo = f22b_y_lo, ylim_hi = f22b_y_hi,
    caption_y = F22B_SHIVA_CAPTION_Y,
    caption_gap = F22B_SHIVA_CAPTION_LINE_GAP * f22b_text_mult,
    caption_size = F22B_SHIVA_CAPTION_SIZE
  ) +
    ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4"), breaks = c("SCAM", "GAM")) +
    ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4"), guide = "none") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed"), linewidth = c(1.2, 1.4)), reverse = FALSE)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 3)),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.width = grid::unit(1.3, "cm"),
      plot.margin = FIG_F22B_MARGIN
    )
  p22_ox_f22b <- patchwork::wrap_plots(p22_ox_purple, ncol = 1) +
    patchwork::plot_layout(ncol = 1, widths = 1, heights = 1)

  ## Effect forest — FDR-only labels (reference = High)
  effects_plot <- effects_out %>%
    mutate(
      label_f8 = ifelse(is.na(p_adj_fdr), "ref", paste0("FDR = ", formatC(p_adj_fdr, format = "e", digits = 2)))
    )
  x_eff_min <- min(effects_plot$ci_low, effects_plot$estimate, na.rm = TRUE)
  x_eff_max <- max(effects_plot$ci_high, effects_plot$estimate, na.rm = TRUE)
  x_span <- max(x_eff_max - x_eff_min, 0.1)
  x_pad <- 0.25 * x_span
  x_label_pos <- x_eff_max + x_pad
  x_plot_min <- x_eff_min - 0.05 * x_span
  x_plot_max <- x_eff_max + x_pad * 3
  label_size_f8 <- fsz(4 * 1.5 * 1.2)
  p8 <- ggplot2::ggplot(effects_plot, ggplot2::aes(x = estimate, y = factor(rank_int, levels = rev(rank_levels), labels = rev(rank_labels)))) +
    ggplot2::geom_vline(xintercept = 0, color = "gray35", linewidth = 0.9) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0.2, linewidth = 1.44, color = "black", orientation = "y") +
    ggplot2::geom_point(size = 7.21, color = "black") +
    ggplot2::geom_text(ggplot2::aes(x = x_label_pos, label = label_f8), hjust = 0, size = label_size_f8, color = "black") +
    ggplot2::labs(
      title = NULL,
      x = "Effect size (avg PFS) (mixed model)",
      y = "DDA score tiers"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_plot_min, x_plot_max), clip = "off") +
    theme_pub
  fig8_h <- 6
  p8 <- p8 + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = fsz(16 * 1.5), color = "black"),
    axis.title = ggplot2::element_text(size = fsz(14 * 1.5 * 1.2), color = "black"),
    axis.text = ggplot2::element_text(size = fsz(12 * 1.5 * 1.2), color = "black"),
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.8),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = grid::unit(0.24, "cm"),
    legend.position = "none"
  )

  if ("RESPONSE_BESTAVGRESPONSE" %in% colnames(df_all_ranks)) {
    ## Imported SHIVA-group mRECIST families (mirroring normal-rank panels):
    ## F99D/F99E/F99L/F432/F432B from RESPONSE_BESTAVGRESPONSE
    ## F456/F462/F465/F466 from RESPONSE_BEST_RESPONSE.
    make_mrecist_family_shiva <- function(response_col, legend_title, title_plot) {
      d_resp <- df_all_ranks %>%
        dplyr::mutate(resp = trimws(toupper(as.character(.data[[response_col]])))) %>%
        dplyr::filter(resp %in% c("CR", "PR", "SD", "PD"))
      if (nrow(d_resp) == 0L) return(NULL)

      mrecist_by_rank_local <- d_resp %>%
        dplyr::group_by(rank_pooled, resp) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(rank_pooled) %>%
        dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(resp = factor(resp, levels = c("PD", "SD", "PR", "CR")))

      dcr_local <- d_resp %>%
        dplyr::mutate(is_dcr = resp %in% c("CR", "PR", "SD")) %>%
        dplyr::group_by(rank_pooled) %>%
        dplyr::summarise(n_dcr = sum(is_dcr, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(dcr_pct = 100 * n_dcr / n)

      orr_local <- d_resp %>%
        dplyr::group_by(rank_pooled) %>%
        dplyr::summarise(n_orr = sum(resp %in% c("CR", "PR"), na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(orr_pct = dplyr::if_else(n > 0, 100 * n_orr / n, 0))

      if (nrow(mrecist_by_rank_local) == 0L || nrow(dcr_local) == 0L || nrow(orr_local) == 0L) return(NULL)

      dcr_local <- dcr_local %>% dplyr::arrange(rank_pooled)
      orr_local <- orr_local %>% dplyr::arrange(rank_pooled)

      # Cochran-Armitage trend-test labels for DCR/ORR annotations (matching normal-rank style).
      make_ca_label <- function(df_counts, success_col, metric_prefix, panel_id = NA_character_) {
        ord <- df_counts %>% dplyr::arrange(rank_pooled)
        ann <- ca_trend_annotate(
          scores = ord$rank_pooled, n = ord$n, success = ord[[success_col]],
          metric = metric_prefix, panel_id = panel_id,
          analysis = "SHIVA_tier", stratification = "rank_pooled"
        )
        ann$label
      }
      trend_lab_dcr_local <- make_ca_label(dcr_local, "n_dcr", "DCR", "mRECIST_family_DCR")
      trend_lab_orr_local <- make_ca_label(orr_local, "n_orr", "ORR", "mRECIST_family_ORR")

      mrecist_colors_std <- c("CR" = "#1565c0", "PR" = "#bbdefb", "SD" = "#d4d0f5", "PD" = "#d9d9d9")
      mrecist_colors_blue <- c("CR" = "#1565c0", "PR" = "#bbdefb", "SD" = "#d4d0f5", "PD" = "#d9d9d9")

      p_main <- ggplot2::ggplot(
        mrecist_by_rank_local,
        ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = resp)
      ) +
        ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
        ggplot2::geom_line(
          data = dcr_local,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = dcr_pct, group = 1, color = "CA (DCR%)"),
          linewidth = 1, inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = orr_local,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = orr_pct, group = 1, color = "CA (ORR%)"),
          linewidth = 1, inherit.aes = FALSE
        ) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_orr_local, parse = TRUE, hjust = 1.05, vjust = 2.10, size = fsz(3), colour = "black") +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_dcr_local, parse = TRUE, hjust = 1.05, vjust = 1.00, size = fsz(3), colour = "black") +
        ggplot2::scale_fill_manual(values = mrecist_colors_std, name = legend_title) +
        ggplot2::scale_color_manual(
          values = c("CA (DCR%)" = "#008000", "CA (ORR%)" = "#FFA500"),
          breaks = c("CA (DCR%)", "CA (ORR%)"),
          name = NULL
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(x = "Tiers", y = "Best Response (%)", title = title_plot) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.box = "vertical",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      ## main_trend: same as p_main but DCR/ORR lines are linear trend lines (lm fits) instead
      ## of observed-value polylines, and all-vs-all pairwise Fisher-test brackets are overlaid
      ## (style matches F99D / F456 / F432 in PDX_crossover_normal_ranks_KM_HR_19APR2O26.R).
      lm_fit_dcr_loc <- stats::lm(dcr_pct ~ rank_pooled, data = dcr_local)
      trend_df_dcr_loc <- data.frame(
        rank_pooled = rank_levels,
        trend_fit_dcr = stats::predict(lm_fit_dcr_loc, newdata = data.frame(rank_pooled = rank_levels))
      )
      lm_fit_orr_loc <- stats::lm(orr_pct ~ rank_pooled, data = orr_local)
      trend_df_orr_loc <- data.frame(
        rank_pooled = rank_levels,
        trend_fit_orr = stats::predict(lm_fit_orr_loc, newdata = data.frame(rank_pooled = rank_levels))
      )

      p_main_trend <- ggplot2::ggplot(
        mrecist_by_rank_local,
        ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = resp)
      ) +
        ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
        ggplot2::geom_line(
          data = trend_df_dcr_loc,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_dcr, group = 1, color = "CA (DCR%)"),
          linewidth = 1, linetype = "solid", inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = trend_df_orr_loc,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_orr, group = 1, color = "CA (ORR%)"),
          linewidth = 1, linetype = "solid", inherit.aes = FALSE
        ) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_orr_local, parse = TRUE, hjust = 1.05, vjust = 2.10, size = fsz(3), colour = "black") +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_dcr_local, parse = TRUE, hjust = 1.05, vjust = 1.00, size = fsz(3), colour = "black") +
        ggplot2::scale_fill_manual(values = mrecist_colors_std, name = legend_title) +
        ggplot2::scale_color_manual(
          values = c("CA (DCR%)" = "#008000", "CA (ORR%)" = "#FFA500"),
          breaks = c("CA (DCR%)", "CA (ORR%)"),
          name = NULL
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1), color = ggplot2::guide_legend(order = 2)) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(x = "Tiers", y = "Best Response (%)", title = title_plot) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.box = "vertical",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      ann_orr_loc <- build_allvsall_prop_annotations(
        summary_df = orr_local,
        rank_col = "rank_pooled",
        success_col = "n_orr",
        total_col = "n",
        trend_df = trend_df_orr_loc,
        trend_rank_col = "rank_pooled",
        trend_value_col = "trend_fit_orr",
        layer_offset = 2.2,
        span_step = 4.0
      )
      ann_dcr_loc <- build_allvsall_prop_annotations(
        summary_df = dcr_local,
        rank_col = "rank_pooled",
        success_col = "n_dcr",
        total_col = "n",
        trend_df = trend_df_dcr_loc,
        trend_rank_col = "rank_pooled",
        trend_value_col = "trend_fit_dcr",
        layer_offset = 6.8,
        span_step = 4.3
      )
      p_main_trend <- add_pairwise_signif_brackets(p_main_trend, ann_orr_loc)
      p_main_trend <- add_pairwise_signif_brackets(p_main_trend, ann_dcr_loc)

      p_bars <- ggplot2::ggplot(
        mrecist_by_rank_local,
        ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = resp)
      ) +
        ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
        ggplot2::scale_fill_manual(values = mrecist_colors_std, name = legend_title) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(x = "Tiers", y = "Best Response (%)", title = title_plot) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      p_bars_blue <- ggplot2::ggplot(
        mrecist_by_rank_local,
        ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = resp)
      ) +
        ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
        ggplot2::scale_fill_manual(values = mrecist_colors_blue, name = legend_title) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(x = "Tiers", y = "Best Response (%)", title = title_plot) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      list(main = p_main, main_trend = p_main_trend, bars = p_bars, bars_blue = p_bars_blue)
    }

    fam_avg <- make_mrecist_family_shiva(
      "RESPONSE_BESTAVGRESPONSE",
      "mRECIST\n(BestAvgResponse)",
      "Best Average Response by\nTiers (BestAvgResponse)"
    )
    if (!is.null(fam_avg)) {
      p99d <- fam_avg$main
      p99e <- fam_avg$main
      p432 <- fam_avg$main_trend +
        ggplot2::labs(
          title = "Best Response by DDA Score Tiers",
          subtitle = "(BestAvgResponse data)"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 6), hjust = 0.5)
        )
      p432b <- fam_avg$bars
    } else {
      p99e <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F99E: pooled SHIVA-group mRECIST not available")
      p432 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "pooled SHIVA-group mRECIST not available")
      p432b <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "pooled SHIVA-group mRECIST bars not available")
    }

    if ("RESPONSE_BEST_RESPONSE" %in% colnames(df_all_ranks)) {
      fam_best <- make_mrecist_family_shiva(
        "RESPONSE_BEST_RESPONSE",
        "mRECIST\n(BestResponse)",
        "Best Response by\nTiers (BestResponse)"
      )
      if (!is.null(fam_best)) {
        ## F456: linear trend lines (lm fit) + all-vs-all pairwise brackets (3 groups, 3 brackets).
        p456 <- fam_best$main_trend +
          ggplot2::labs(
            title = "Best Response by DDA Score Tiers",
            subtitle = "(BestResponse data)"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 6), hjust = 0.5)
          )
        p462 <- fam_best$bars
        p465 <- fam_best$main
        p466 <- fam_best$bars
      } else {
        p456 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "RESPONSE_BEST_RESPONSE insufficient data")
        p462 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F462: RESPONSE_BEST_RESPONSE insufficient data")
        p465 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F465: RESPONSE_BEST_RESPONSE insufficient data")
        p466 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F466: RESPONSE_BEST_RESPONSE insufficient data")
      }
    } else {
      p456 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "RESPONSE_BEST_RESPONSE column not found")
      p462 <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "F462: RESPONSE_BEST_RESPONSE column not found")
    }
  } else {
    p88 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F88: Response_BestAvgResponse not available")
    p_dcr_orr_heat <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F100: Response_BestAvgResponse not available (DCR/ORR heatmaps)")
    p99 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99: Response_BestAvgResponse not available")
    p99b <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99B: Response_BestAvgResponse not available")
    p99l <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99L: Response_BestAvgResponse not available")
    p99x <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99X: Response_BestAvgResponse not available")
    p34 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F34: Response_BestAvgResponse not available")
    p99c <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99C: Response_BestAvgResponse not available")
    p99d <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99D: Response_BestAvgResponse not available")
    p456 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Response_Best_Response not available")
    p462 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F462: Response_Best_Response not available")
  }

  ## Durable Benefit (DB) by SHIVA group
  ## DB definition: TIME_TO_DOUBLE > 42
  db_df <- df_all_ranks %>%
    dplyr::filter(is.finite(TIME_TO_DOUBLE), !is.na(rank_pooled)) %>%
    dplyr::mutate(is_db = TIME_TO_DOUBLE > 42)

  if (nrow(db_df) > 0L) {
    db_by_rank <- db_df %>%
      dplyr::group_by(rank_pooled) %>%
      dplyr::summarise(
        n_db = sum(is_db, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        db_pct = dplyr::if_else(n > 0, 100 * n_db / n, 0),
        se_pct = dplyr::if_else(n > 0, 100 * sqrt((db_pct / 100) * (1 - db_pct / 100) / n), 0)
      ) %>%
      dplyr::arrange(rank_pooled)

    db_components <- db_df %>%
      dplyr::mutate(db_status = dplyr::if_else(is_db, "DB", "Non_DB")) %>%
      dplyr::group_by(rank_pooled, db_status) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::right_join(
        expand.grid(rank_pooled = rank_levels, db_status = c("Non_DB", "DB"), stringsAsFactors = FALSE),
        by = c("rank_pooled", "db_status")
      ) %>%
      dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
      dplyr::group_by(rank_pooled) %>%
      dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(db_status = factor(db_status, levels = c("Non_DB", "DB")))

    lm_fit_db <- stats::lm(db_pct ~ rank_pooled, data = db_by_rank)
    trend_df_db <- data.frame(
      rank_pooled = rank_levels,
      trend_fit_db = stats::predict(lm_fit_db, newdata = data.frame(rank_pooled = rank_levels))
    )

    ca_ann_f467 <- ca_trend_annotate(
      scores = db_by_rank$rank_pooled, n = db_by_rank$n, success = db_by_rank$n_db,
      metric = "DB", panel_id = "Durable_Benefit_by_SHIVA_tier", analysis = "SHIVA_tier", stratification = "rank_pooled"
    )
    trend_lab_db <- ca_ann_f467$label
    ann_db <- build_allvsall_prop_annotations(
      summary_df = db_by_rank,
      rank_col = "rank_pooled",
      success_col = "n_db",
      total_col = "n",
      trend_df = trend_df_db,
      trend_rank_col = "rank_pooled",
      trend_value_col = "trend_fit_db",
      layer_offset = 4,
      span_step = 3.9,
      y_cap = 106
    )

    db_colors <- c("DB" = "#81c784", "Non_DB" = "#d9d9d9")
    p467 <- ggplot2::ggplot(
      db_components,
      ggplot2::aes(
        x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
        y = pct,
        fill = db_status
      )
    ) +
      ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
      ggplot2::geom_line(
        data = trend_df_db,
        ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_db, group = 1, color = "CA (DB%)"),
        linewidth = 1,
        inherit.aes = FALSE
      ) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_db, parse = TRUE, hjust = 1.05, vjust = 1.5, size = fsz(3), colour = "black") +
      ggplot2::scale_fill_manual(
        values = db_colors,
        name = "Durable Benefit\n(>42 days)",
        breaks = c("DB", "Non_DB"),
        labels = c("DB", "NonDB")
      ) +
      ggplot2::scale_color_manual(values = c("CA (DB%)" = "black"), breaks = c("CA (DB%)"), name = NULL) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(order = 1),
        color = ggplot2::guide_legend(order = 2)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        x = "Tiers",
        y = "Durable Benefit (%)",
        title = "Durable Benefit by DDA Score Tiers"
      ) +
      theme_pub_nogrid +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = fsz(18), margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.height = grid::unit(0.55, "cm"),
        legend.spacing.y = grid::unit(0.15, "cm"),
        legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
      )
    p467 <- add_pairwise_signif_brackets(p467, ann_db)
  } else {
    p467 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Durable Benefit by SHIVA group not available")
  }

  ca_cmp_export <- ca_registry_export(analysis_dir)
  log_note(paste0(
    "Cochran-Armitage registry exported (n=", nrow(ca_cmp_export),
    "; Z mismatches=", sum(!ca_cmp_export$z_match, na.rm = TRUE),
    "; p mismatches=", sum(!ca_cmp_export$p_match, na.rm = TRUE), ")."
  ))

  save_fig(p6, file.path(fig_dir, "SCAM_residual_diagnostics.png"), FIG_F6_WIDTH, FIG_F6_HEIGHT)
  save_fig(p6b_residual_panels, file.path(fig_dir, "SCAM_residuals_vs_fitted_by_COMPOUND_and_GROUP.png"), FIG_F6_WIDTH, FIG_F6_HEIGHT)
  save_fig(p22_ox_f22b, file.path(fig_dir, "SCAM_vs_GAM_all_ranks.png"), FIG_F22B_WIDTH, FIG_F22B_HEIGHT)
  save_fig(p8, file.path(fig_dir, "Effect_forest_enhanced.png"), FIG_F8_WIDTH, FIG_F8_HEIGHT)
  save_fig(p432, file.path(fig_dir, "mRECIST_pooled_SHIVA_group_ORR_DCR.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)
  save_fig(p432b, file.path(fig_dir, "mRECIST_pooled_SHIVA_group_bars_only.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)
  save_fig(p456, file.path(fig_dir, "mRECIST_SHIVA_group_ORR_DCR_Response_Best_Response.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)
  save_fig(p467, file.path(fig_dir, "Durable_Benefit_by_SHIVA_group_CA_trend.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)

  log_note("Step 9 complete: publication figures saved.")
  invisible(NULL)
})

log_step("All done: SHIVA-tier data modeling analysis completed.")



## SECTION 2: Normal-rank data modeling analysis


input_path <- file.path(DATA_DIR, DEFAULT_INPUT_XLSX)
## All analysis below uses all available rows with non-missing TimeToDouble values.
## Publication mixed-model contrasts: each normal rank vs rank 1 (reference) on ranked PDX rows only.
## Inference: lmerTest (Satterthwaite p-values from summary); forest-plot CIs via profile likelihood (Wald fallback).
output_root <- DATA_MODELING_OUTPUT_ROOT
date <- format(Sys.Date(), "%Y_%m_%d")
name <- "normal_rank_data_modeling"

set.seed(42)
PLOT_ONLY_MODE <- FALSE
if (identical(Sys.getenv("PDX_PLOT_ONLY"), "1")) PLOT_ONLY_MODE <- TRUE
SAVE_PLOT_DATA_CACHE <- TRUE

# Name of the response column (after renaming: spaces -> '_' and uppercasing).
# Default for earlier runs was 'TIMETODOUBLE'. Currently set to 'DAY_LAST'.
# Example alternative (for later use): 'TimeToDouble' -> internal name 'TIMETODOUBLE'
# RESPONSE_COL_NAME <- "TIMETODOUBLE"
RESPONSE_COL_NAME <- "TIMETODOUBLE"
# Label used in plot titles/axes for the response; kept simple so it always matches the chosen column.
RESPONSE_COL_LABEL <- RESPONSE_COL_NAME

# Rank pooling: ranks >= RANK_POOL_FROM are pooled into one category labeled RANK_POOL_LABEL (e.g. "7-10").
RANK_POOL_FROM <- 7L
RANK_POOL_LABEL <- "7-10"

# Performance metrics (DCR-based): ground truth from Response_BestAvgResponse; prediction from LEVEL vs threshold.
# Ground truth = 1 (True): response in GROUND_TRUTH_TRUE; = 0 (False): response in GROUND_TRUTH_FALSE.
GROUND_TRUTH_TRUE <- c("CR", "PR", "SD")   # DCR responders c("CR", "PR", "SD") 
GROUND_TRUTH_FALSE <- c("PD")              # non-responders
# Response type for performance metric titles and plot"." names (e.g. "DCR").
RESPONSE_TYPE <- "DCR"

# Figure layout controls (publication style)
FIG_WIDTH <- 13
FIG_HEIGHT <- 10

# Per-plot figure sizes (edit here)
FIG_SCAM_MONO_WIDTH <- FIG_WIDTH
FIG_SCAM_MONO_HEIGHT <- FIG_HEIGHT
FIG_F22_JPG_WIDTH <- FIG_WIDTH
FIG_F22_JPG_HEIGHT <- FIG_HEIGHT
FIG_EFFECTS_LEGACY_WIDTH <- FIG_WIDTH
FIG_EFFECTS_LEGACY_HEIGHT <- FIG_HEIGHT
FIG_F2A_WIDTH <- FIG_WIDTH
FIG_F2A_HEIGHT <- FIG_HEIGHT
FIG_F24_WIDTH <- FIG_WIDTH
FIG_F24_HEIGHT <- FIG_HEIGHT
FIG_F2_WIDTH <- FIG_WIDTH
FIG_F2_HEIGHT <- FIG_HEIGHT
FIG_F3_WIDTH <- FIG_WIDTH
FIG_F3_HEIGHT <- FIG_HEIGHT
FIG_F25_WIDTH <- FIG_WIDTH
FIG_F25_HEIGHT <- FIG_HEIGHT
FIG_F4B_WIDTH <- FIG_WIDTH
FIG_F4B_HEIGHT <- FIG_HEIGHT
FIG_F6_WIDTH <- FIG_WIDTH
FIG_F6_HEIGHT <- FIG_HEIGHT
FIG_F22_WIDTH <- FIG_WIDTH
FIG_F22_HEIGHT <- FIG_HEIGHT
FIG_F22B_WIDTH <- FIG_WIDTH
FIG_F22B_HEIGHT <- FIG_HEIGHT
FIG_F22B_MARGIN <- ggplot2::margin(t = 40, r = 18, b = 12, l = 12, unit = "pt")
F22B_NLABEL_Y_PAD <- 8
FIG_F22C_WIDTH <- FIG_WIDTH
FIG_F22C_HEIGHT <- FIG_HEIGHT
FIG_F8_WIDTH <- FIG_WIDTH
FIG_F8_HEIGHT <- FIG_HEIGHT
FIG_F88_WIDTH <- 8
FIG_F88_HEIGHT <- 5
FIG_F99_WIDTH <- 8
FIG_F99_HEIGHT <- 5
FIG_F99B_WIDTH <- 8
FIG_F99B_HEIGHT <- 5
FIG_F99C_WIDTH <- 8
FIG_F99C_HEIGHT <- 5
FIG_F34_WIDTH <- 10
FIG_F34_HEIGHT <- 6
FIG_F77_WIDTH <- 24
FIG_F77_HEIGHT <- 7
FIG_F23_WIDTH <- FIG_WIDTH
FIG_F23_HEIGHT <- FIG_HEIGHT
PANEL_HORIZONTAL_GAP <- 0.15  # F23: small gap between A-B and C-D
PANEL_VERTICAL_GAP <- 0.0    # F23: gap between top row (A-B) and bottom row (C-D)
FONT_SCALE_ALL <- 1.40

fsz <- function(x) x * FONT_SCALE_ALL

## Step timeout: if a step runs longer than this (seconds), abort (stuck detection). 0 = no limit.
STEP_TIMEOUT_SEC <- 0L
if (STEP_TIMEOUT_SEC > 0) {
  message("Step timeout enabled: ", STEP_TIMEOUT_SEC, "s per step (stuck detection).")
}

## Create date-named folder; save everything there (create if not present)
analysis_dir <- file.path(output_root, paste(date, name, sep = "_"))
fig_dir <- file.path(analysis_dir, "figures")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
plot_data_cache_path <- file.path(analysis_dir, "plot_data_cache.rds")

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

step_log_path <- file.path(analysis_dir, "step_log.txt")
analysis_log_path <- file.path(analysis_dir, "analysis_log.txt")
heartbeat_path <- file.path(analysis_dir, "last_activity.txt")
writeLines(paste0("=== Analysis run started: ", timestamp_now(), " ==="), step_log_path)
writeLines(paste0("=== Analysis run started: ", timestamp_now(), " ==="), analysis_log_path)
writeLines(paste0(timestamp_now(), " | Run started."), heartbeat_path)

.step_counter <- 0L

log_step <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  cat(line, "\n")
  write(line, step_log_path, append = TRUE)
}

log_note <- function(message) {
  line <- paste0("[", timestamp_now(), "] ", message)
  write(line, analysis_log_path, append = TRUE)
}

validate_plot_cache <- function(cache_obj) {
  required_fields <- c(
    "df", "df_all_ranks", "scam_out", "wiggle_out", "effects_out",
    "rank_levels", "max_rank"
  )
  missing_fields <- setdiff(required_fields, names(cache_obj))
  if (length(missing_fields) > 0) {
    stop(paste0("Plot cache missing required fields: ", paste(missing_fields, collapse = ", ")))
  }
  if (!is.data.frame(cache_obj$df) || !is.data.frame(cache_obj$df_all_ranks)) {
    stop("Plot cache validation failed: df/df_all_ranks are not data frames.")
  }
  invisible(TRUE)
}

## Ensure plot-facing text uses spaced separators, e.g. "p = 0.05", "n = 10", "EDF = 1.2".
space_plot_separators <- function(s) {
  if (length(s) == 0L) return(s)
  out <- as.character(s)
  out <- gsub("([A-Za-z\\)\\%\\|])(=)([^=<>])", "\\1 = \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(<)(?!=)", "\\1 < \\3", out, perl = TRUE)
  out <- gsub("([A-Za-z0-9\\)])(>)(?!=)", "\\1 > \\3", out, perl = TRUE)
  gsub("\\s{2,}", " ", out, perl = TRUE)
}

run_step <- function(step_name, fn) {
  .step_counter <<- .step_counter + 1L
  n <- .step_counter
  start_time <- Sys.time()
  ts_start <- timestamp_now()

  start_msg <- paste0("========== Step ", n, ": ", step_name, " - START at ", ts_start, " ==========")
  cat(start_msg, "\n")
  write(start_msg, step_log_path, append = TRUE)
  write(paste0(ts_start, " | Step ", n, " START: ", step_name), heartbeat_path)

  if (STEP_TIMEOUT_SEC > 0) {
    setTimeLimit(elapsed = STEP_TIMEOUT_SEC, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  }
  out <- NULL
  err <- NULL
  tryCatch(
    {
      out <- fn()
    },
    error = function(e) {
      err <<- conditionMessage(e)
    }
  )
  if (!is.null(err)) {
    fail_msg <- paste0("========== Step ", n, ": ", step_name, " - FAIL at ", timestamp_now(), " ==========\n", err)
    cat(fail_msg, "\n")
    write(fail_msg, step_log_path, append = TRUE)
    write(paste0(timestamp_now(), " | Step ", n, " FAILED: ", step_name), heartbeat_path)
    stop(err)
  }

  ts_end <- timestamp_now()
  dur_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  ok_msg <- paste0("========== Step ", n, ": ", step_name, " - OK at ", ts_end, " (duration ", dur_sec, "s) ==========")
  cat(ok_msg, "\n")
  write(ok_msg, step_log_path, append = TRUE)
  write(paste0(ts_end, " | Step ", n, " OK (", dur_sec, "s): ", step_name), heartbeat_path)
  log_step(paste0("Heartbeat: ", ts_end, " - Step ", n, " completed"))
  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), "p = NA", paste0("p = ", format.pval(p, digits = 3, eps = 1e-4)))
}

fmt_p_sci2 <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", formatC(p, format = "e", digits = 2))
}

fmt_p_sci3 <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", formatC(p, format = "e", digits = 3))
}

fmt_p_sci1 <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", formatC(p, format = "e", digits = 0))
}

fmt_p_sci1_cap <- function(p, min_p = 1e-100) {
  if (is.na(p)) return("p = NA")
  if (p < min_p) return(paste0("p < ", format(min_p, scientific = TRUE)))
  paste0("p = ", formatC(p, format = "e", digits = 0))
}

fmt_pval_sci2 <- function(p) {
  if (is.na(p)) return("NA")
  formatC(p, format = "e", digits = 2)
}

fmt_pval_sci2_cap <- function(p, min_p = 1e-30) {
  if (is.na(p)) return("NA")
  if (p < min_p) return("<1e-30")
  formatC(p, format = "e", digits = 2)
}

## F22/F22b corner annotations: single place for SCAM/GAM smooth p display.
## For p < cap_below (underflow / numeric zero is also caught because 0 < cap_below), shows p < floor.
## sci_digits = 2L => two decimals after the leading digit in scientific notation, e.g. "p = 4.19e-21".
fmt_sci_e <- function(x, sci_digits = 2L) {
  formatC(as.numeric(x), format = "e", digits = sci_digits)
}

fmt_edf_plot <- function(edf) {
  format(round(as.numeric(edf), 3), nsmall = 3)
}

fmt_smooth_p_plot_label <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  if (is.na(p) || !is.finite(as.numeric(p))) return("p = NA")
  p <- as.numeric(p)
  if (p < cap_below) {
    return(paste0("p < ", fmt_sci_e(cap_below, sci_digits)))
  }
  paste0("p = ", fmt_sci_e(p, sci_digits))
}

## On-figure labels only: render "p" as plotmath italic(p) (ggplot parse = TRUE).
num_str_to_plotmath_rhs <- function(s, mantissa_digits = 2L) {
  s <- trimws(as.character(s))
  if (!nzchar(s) || grepl("^[Nn][Aa]$", s)) return("NA")
  if (grepl("^[<>]", s)) return(s)
  m <- regmatches(s, regexec("^([0-9.+\\-]+)[eE]([-+]?\\d+)$", s, perl = TRUE))[[1]]
  if (length(m) == 3L) {
    mant <- sprintf(paste0("%.", mantissa_digits, "f"), as.numeric(m[2]))
    return(paste0("'", mant, "' %*% 10^{", as.integer(m[3]), "}"))
  }
  s
}

quote_plotmath_text <- function(s) {
  paste0("'", gsub("'", "''", space_plot_separators(as.character(s)), fixed = TRUE), "'")
}

## DCR/ORR Cochran-Armitage corner text: italic plotmath p; value = fmt_pval_sci2 (2-decimal e-notation).
sprintf_ca_trend_label <- function(z, p, prefix = NULL) {
  left <- sprintf("Cochran-Armitage: Z = %s, ", format(round(z, 3), nsmall = 3))
  if (!is.null(prefix) && nzchar(as.character(prefix)[1])) {
    left <- paste0(as.character(prefix)[1], ": ", left)
  }
  paste0(
    "paste(",
    quote_plotmath_text(left),
    ", italic(p), ' = ",
    gsub("'", "''", fmt_pval_sci2(p), fixed = TRUE),
    "')"
  )
}

plain_plot_label_line_to_plotmath <- function(s) {
  s <- space_plot_separators(as.character(s)[1])
  if (!grepl("[Pp]\\s*[=<>]", s)) return(quote_plotmath_text(s))
  chunks <- character(0)
  rest <- s
  while (nzchar(rest) && grepl("[Pp]\\s*[=<>]", rest)) {
    loc <- regexpr("[Pp]\\s*[=<>]", rest, perl = TRUE)[1]
    if (loc > 1L) {
      chunks <- c(chunks, quote_plotmath_text(substr(rest, 1L, loc - 1L)))
    }
    rest <- substr(rest, loc, nchar(rest))
    m <- regexpr("^[Pp]\\s*(=|<)\\s*([^,\\)\\n]+)", rest, perl = TRUE)
    if (m[1L] == -1L) {
      chunks <- c(chunks, quote_plotmath_text(rest))
      break
    }
    ml <- attr(m, "match.length")
    token <- substr(rest, 1L, ml)
    if (grepl("^[Pp]\\s*<", token, perl = TRUE)) {
      rhs <- sub("^[Pp]\\s*<\\s*", "", token, ignore.case = TRUE)
      chunks <- c(chunks, paste0("italic(p) < ", num_str_to_plotmath_rhs(rhs)))
    } else {
      rhs <- sub("^[Pp]\\s*=\\s*", "", token, ignore.case = TRUE)
      chunks <- c(chunks, paste0("italic(p) == ", num_str_to_plotmath_rhs(rhs)))
    }
    rest <- substr(rest, ml + 1L, nchar(rest))
  }
  if (nzchar(rest)) chunks <- c(chunks, quote_plotmath_text(rest))
  if (length(chunks) == 1L) return(chunks[[1L]])
  paste0("paste(", paste(chunks, collapse = ", "), ")")
}

plain_plot_label_to_plotmath <- function(s) {
  lines <- strsplit(as.character(s)[1], "\n", fixed = TRUE)[[1]]
  exprs <- vapply(lines, plain_plot_label_line_to_plotmath, character(1))
  if (length(exprs) == 1L) return(exprs[[1L]])
  paste0("paste(", paste(exprs, collapse = ", '\n', "), ")")
}

label_uses_italic_p <- function(s) {
  grepl("italic\\(p\\)", plain_plot_label_to_plotmath(s), fixed = TRUE)
}

as_ggplot_label <- function(s) {
  s <- space_plot_separators(s)
  pm <- plain_plot_label_to_plotmath(s)
  if (!label_uses_italic_p(s)) return(s)
  parse(text = pm)[[1]]
}

fmt_p_plot <- function(p) {
  plain_plot_label_to_plotmath(fmt_p(p))
}

fmt_smooth_p_plot_label_plot <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  plain_plot_label_to_plotmath(fmt_smooth_p_plot_label(p, cap_below = cap_below, sci_digits = sci_digits))
}

## F22/F22b on-figure: italic(p) with literal 2-digit e-notation (not %*% 10^{} plotmath).
fmt_smooth_p_f22_plotmath <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  plain <- fmt_smooth_p_plot_label(p, cap_below = cap_below, sci_digits = sci_digits)
  if (grepl("^[Pp]\\s*<", plain)) {
    rhs <- trimws(sub("^[Pp]\\s*<\\s*", "", plain))
    return(paste0("paste(italic(p), ' < ", gsub("'", "''", rhs, fixed = TRUE), "')"))
  }
  if (grepl("^[Pp]\\s*=", plain)) {
    rhs <- trimws(sub("^[Pp]\\s*=\\s*", "", plain))
    return(paste0("paste(italic(p), ' = ", gsub("'", "''", rhs, fixed = TRUE), "')"))
  }
  quote_plotmath_text(plain)
}

sig_stars_from_p <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 5e-4) return("***")
  if (p < 5e-3) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

## F22/F22b: GAM on top, SCAM below — two annotate layers (plotmath newline is unreliable).
F22B_RANKS_CAPTION_Y <- 10
F22B_RANKS_CAPTION_LINE_GAP <- 3.5

f22_caption_y_mid <- function(corner = c("top_right", "f22b_ranks"), y = NULL) {
  corner <- match.arg(corner)
  if (!is.null(y)) return(as.numeric(y)[1])
  if (corner == "f22b_ranks") F22B_RANKS_CAPTION_Y else 10
}

f22_one_line_caption_layer <- function(model, edf, p_val, y, size = fsz(3.8), vjust = 0.5) {
  left <- paste0(sprintf("%-4s", model), " EDF = ", fmt_edf_plot(edf), ", ")
  lab <- paste0("paste(", quote_plotmath_text(left), ", ", fmt_smooth_p_f22_plotmath(p_val), ")")
  ggplot2::annotate(
    "text", x = Inf, y = y, label = lab,
    parse = TRUE, hjust = 1.02, vjust = vjust, size = size, color = "black"
  )
}

f22_scam_gam_caption_gam_layer <- function(
  gam_edf, gam_p,
  corner = c("top_right", "f22b_ranks"),
  size = fsz(3.8),
  y = NULL,
  line_gap = 3.5
) {
  y_mid <- f22_caption_y_mid(corner, y)
  gap <- as.numeric(line_gap)[1]
  f22_one_line_caption_layer("GAM", gam_edf, gam_p, y_mid + gap / 2, size)
}

f22_scam_gam_caption_scam_layer <- function(
  scam_edf, scam_p,
  corner = c("top_right", "f22b_ranks"),
  size = fsz(3.8),
  y = NULL,
  line_gap = 3.5
) {
  y_mid <- f22_caption_y_mid(corner, y)
  gap <- as.numeric(line_gap)[1]
  f22_one_line_caption_layer("SCAM", scam_edf, scam_p, y_mid - gap / 2, size)
}

log_scam_gam_pvalues <- function(scam_out, wiggle_out, analysis_dir, pipeline_tag) {
  cap_scam <- fmt_smooth_p_plot_label(scam_out$p_value)
  cap_gam <- fmt_smooth_p_plot_label(wiggle_out$p_value)
  lines <- c(
    paste0("pipeline=", pipeline_tag),
    paste0("SCAM: p=", scam_out$p_value, " | edf=", scam_out$edf),
    paste0("GAM:  p=", wiggle_out$p_value, " | edf=", wiggle_out$edf),
    paste0("F22 caption labels: SCAM=", cap_scam, " GAM=", cap_gam)
  )
  out_path <- file.path(analysis_dir, "scam_gam_pvalues.txt")
  writeLines(lines, out_path)
  row1 <- data.frame(
    pipeline = pipeline_tag,
    scam_p = scam_out$p_value,
    scam_edf = scam_out$edf,
    gam_p = wiggle_out$p_value,
    gam_edf = wiggle_out$edf,
    scam_f22_label = cap_scam,
    gam_f22_label = cap_gam,
    stringsAsFactors = FALSE
  )
  write.csv(row1, file.path(analysis_dir, "scam_gam_pvalues.csv"), row.names = FALSE)
  for (ln in lines) log_note(paste0("SCAM/GAM p-values: ", ln))
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

  ## Brackets: first rank vs each other rank only (Fisher); no adjacent-only pairs (e.g. 2 vs 3).
  comp_pairs <- tibble::tibble(
    x = rep(1L, nrow(ord) - 1L),
    xend = 2L:nrow(ord)
  )
  comp_pairs <- comp_pairs %>%
    dplyr::mutate(
      x = pmin(x, xend),
      xend = pmax(x, xend)
    ) %>%
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
    tibble::tibble(
      x = i1,
      xend = i2,
      span = i2 - i1,
      p_value = p_val,
      label = sig_stars_from_p(p_val)
    )
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
    ## Even vertical spacing: y_anchor + span*span_step can bunch or cross when y_anchor differs by pair.
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
  comp_df %>%
    dplyr::mutate(
      y_low = pmax(0, y - 1.8),
      x_mid = (x + xend) / 2
    )
}

add_pairwise_signif_brackets <- function(
  plot_obj,
  ann_df,
  text_size = fsz(3),
  text_y_gap = 0.45,
  line_width = 0.35,
  line_color = "black",
  text_color = "black"
) {
  if (is.null(ann_df) || nrow(ann_df) == 0L) return(plot_obj)
  label_pt <- pmax(0.55, text_size - 2)
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
        size = label_pt,
        fontface = "bold",
        colour = text_color,
        vjust = 0
      )
  }
  plot_obj
}

format_logp_sci <- function(log_p, sig_digits = 6L) {
  if (is.na(log_p)) return("NA")
  if (!is.finite(log_p) && log_p < 0) return("0")
  log10_p <- log_p / log(10)
  exp10 <- floor(log10_p)
  mant <- 10^(log10_p - exp10)
  if (mant >= 10) {
    mant <- mant / 10
    exp10 <- exp10 + 1
  }
  paste0(formatC(mant, format = "f", digits = max(0L, as.integer(sig_digits) - 1L)), "e", ifelse(exp10 >= 0, "+", ""), exp10)
}

fmt_fdr <- function(p) {
  if (is.na(p)) return("FDR = NA")
  if (p < 0.001) return("FDR < 0.001")
  if (p < 0.01) return("FDR < 0.01")
  if (p < 0.05) return("FDR < 0.05")
  "FDR >= 0.05"
}

# Two-way radar (ggradar): mint / blue / lilac; lines + points only (no polygon fill).
# radial_zoom=TRUE (F2): metrics in [0.5,1] rescaled linearly to full ggradar [0,1] (0.5 = inner/smallest ring, 1 = outer);
# ring labels 0.50 / 0.75 / 1.00; plot.extent tightened vs default; no 0.6/0.8 rings.
# Full-scale plot: extra 0.6 and 0.8 reference circles; fonts +30% on scale annotations vs old 3.84 base.
two_way_spider_ggradar_plot <- function(ggradar_df, spider_order, RESPONSE_TYPE, radial_zoom = FALSE) {
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    return(NULL)
  }
  ggradar_df <- as.data.frame(ggradar_df)
  radar_legend_levels <- c("1 vs 7-10", "1-2 vs 6-10", "1-3 vs 5-10")
  u_grp <- unique(as.character(ggradar_df$group))
  grp_levels <- c(intersect(radar_legend_levels, u_grp), setdiff(u_grp, radar_legend_levels))
  ggradar_df$group <- factor(ggradar_df$group, levels = grp_levels)

  col_mint <- "#8FBC72"
  col_blue <- "#52A8D8"
  col_lilac <- "#a5a5ef"
  label_to_col <- c(
    "1 vs 7-10" = col_mint,
    "1-2 vs 6-10" = col_blue,
    "1-3 vs 5-10" = col_lilac
  )
  grp_cols <- unname(label_to_col[grp_levels])
  grp_cols[is.na(grp_cols)] <- "#888888"
  radar_grid_grey <- "#404040"

  fmt_ring_lbl <- function(x) {
    format(round(x, 3), trim = TRUE, nsmall = 2, digits = 3)
  }

  if (isTRUE(radial_zoom)) {
    # Map true metric band [0.5, 1] linearly to display [0, 1]: inner ring = 0.5 (small), outer = 1.0.
    plot_df_for_ggradar <- ggradar_df
    for (j in seq_len(ncol(plot_df_for_ggradar))[-1L]) {
      x <- as.numeric(plot_df_for_ggradar[[j]])
      plot_df_for_ggradar[[j]] <- pmax(0, pmin(1, (x - 0.5) / 0.5))
    }
    grid_min <- 0
    grid_max <- 1
    grid_mid <- 0.5
    values_radar_lbl <- c(fmt_ring_lbl(0.5), fmt_ring_lbl(0.75), fmt_ring_lbl(1))
  } else {
    plot_df_for_ggradar <- ggradar_df
    grid_min <- 0
    grid_max <- 1
    grid_mid <- 0.5
    values_radar_lbl <- c("0", "0.5", "1")
  }

  centre_y <- grid_min - (1 / 9) * (grid_max - grid_min)
  offset <- abs(centre_y)
  fn_circle <- getFromNamespace("funcCircleCoords", "ggradar")
  grid_lbl_x <- -0.1 * (grid_max - centre_y)

  # ggradar draws ring labels at grid.label.size * 0.8; same base for 0.6 / 0.8 annotations (+30% vs prior 3.84).
  radar_grid_label_pt <- 3.84 * 1.3
  radar_scale_ann_pt <- radar_grid_label_pt * 0.8

  p <- ggradar::ggradar(
    plot.data = plot_df_for_ggradar,
    axis.labels = spider_order,
    grid.min = grid_min,
    grid.mid = grid_mid,
    grid.max = grid_max,
    values.radar = values_radar_lbl,
    gridline.min.linetype = "solid",
    gridline.mid.linetype = "solid",
    gridline.max.linetype = "solid",
    grid.line.width = 0.4,
    gridline.min.colour = radar_grid_grey,
    gridline.mid.colour = radar_grid_grey,
    gridline.max.colour = radar_grid_grey,
    grid.label.size = radar_grid_label_pt,
    axis.label.size = if (isTRUE(radial_zoom)) 5.8 else 6,
    axis.label.offset = if (isTRUE(radial_zoom)) 1.28 else 1.15,
    group.colours = grp_cols,
    fill = FALSE,
    line.alpha = 0.72,
    draw.points = TRUE,
    point.alpha = 0.82,
    group.point.size = 2.5,
    group.line.width = 0.95,
    background.circle.colour = "white",
    background.circle.transparency = 0,
    axis.line.colour = radar_grid_grey,
    plot.legend = TRUE,
    legend.title = "Ranks compared",
    legend.position = "right",
    legend.text.size = 16.8,
    plot.title = paste0("Radar plot (", RESPONSE_TYPE, ")"),
    base.size = 14.4,
    plot.extent.x.sf = if (isTRUE(radial_zoom)) 1.48 else 1,
    plot.extent.y.sf = if (isTRUE(radial_zoom)) 1.42 else 1.2,
    x.centre.range = (if (isTRUE(radial_zoom)) 0.035 else 0.02) * (grid_max - centre_y)
  )

  p_out <- p
  if (!isTRUE(radial_zoom)) {
    r06 <- 0.6 + offset
    r08 <- 0.8 + offset
    path06 <- fn_circle(c(0, 0), r06, npoints = 360)
    path08 <- fn_circle(c(0, 0), r08, npoints = 360)
    lab_ref <- data.frame(
      x = grid_lbl_x,
      y = c(r06, r08),
      text = c("0.6", "0.8")
    )
    p_out <- p_out +
      ggplot2::geom_path(
        data = path06,
        ggplot2::aes(x = .data$x, y = .data$y),
        inherit.aes = FALSE,
        colour = radar_grid_grey,
        linewidth = 0.38,
        linetype = "solid"
      ) +
      ggplot2::geom_path(
        data = path08,
        ggplot2::aes(x = .data$x, y = .data$y),
        inherit.aes = FALSE,
        colour = radar_grid_grey,
        linewidth = 0.38,
        linetype = "solid"
      ) +
      ggplot2::geom_text(
        data = lab_ref,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$text),
        inherit.aes = FALSE,
        size = radar_scale_ann_pt,
        colour = radar_grid_grey,
        hjust = 1
      )
  }

  if (isTRUE(radial_zoom)) {
    p_out <- p_out + ggplot2::labs(subtitle = NULL)
  }
  p_out +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(hjust = 0.5, colour = "black"),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(colour = "black")
    )
}

# Radar PNG as soon as two_way_metrics exists (before two-way Excel / line plots).
# Also writes F2_two_way_radar_<RESPONSE_TYPE>.png in the same folder as F1_two_way_radar_*.png (zoomed radial window on [0,1]).
save_two_way_radar_plot <- function(two_way_metrics, two_way_dir, RESPONSE_TYPE, dpi_radar) {
  spider_metrics <- c(
    "PR AUC" = "Precision_Recall_AUC",
    "PPV" = "PPV",
    "ROC AUC" = "ROC_AUC",
    "Accuracy" = "Accuracy",
    "F1" = "F1",
    "specificity" = "Specificity",
    "sensitivity" = "Sensitivity",
    "NPV" = "NPV"
  )
  spider_order <- names(spider_metrics)
  spider_cols <- unname(spider_metrics)
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    log_note("Two-way spider plot skipped: package 'ggradar' not available (install via remotes::install_github('ricardo-bion/ggradar')).")
    return(invisible(NULL))
  }
  if (!all(spider_cols %in% colnames(two_way_metrics)) || nrow(two_way_metrics) == 0) {
    return(invisible(NULL))
  }
  ggradar_df <- two_way_metrics %>%
    dplyr::select(two_way_label, dplyr::all_of(spider_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(spider_cols), ~ pmax(0, pmin(1, as.numeric(.x))))) %>%
    dplyr::rename(group = two_way_label)
  names(ggradar_df) <- c("group", make.names(spider_order))
  p_spider <- two_way_spider_ggradar_plot(ggradar_df, spider_order, RESPONSE_TYPE, radial_zoom = FALSE)
  if (!is.null(p_spider)) {
    out_png <- file.path(two_way_dir, paste0("F1_two_way_radar_", RESPONSE_TYPE, ".png"))
    ggplot2::ggsave(out_png, p_spider, width = 10, height = 10, dpi = dpi_radar, bg = "white")
    log_note(paste0("Saved radar plot (early): ", out_png))
  }
  p_zoom <- two_way_spider_ggradar_plot(ggradar_df, spider_order, RESPONSE_TYPE, radial_zoom = TRUE)
  if (!is.null(p_zoom)) {
    f2_png <- file.path(two_way_dir, paste0("F2_two_way_radar_", RESPONSE_TYPE, ".png"))
    ggplot2::ggsave(f2_png, p_zoom, width = 10, height = 10, dpi = dpi_radar, bg = "white")
    log_note(paste0("Saved zoomed radial radar (F2): ", f2_png))
  }
  invisible(NULL)
}

n_workers <- parallel::detectCores(logical = TRUE)
if (is.na(n_workers) || n_workers < 1) n_workers <- 1L

run_step("Step 1: Load packages and parallel config", function() {
  required_pkgs <- c(
    "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "tibble",
    "lme4", "lmerTest", "scam", "mgcv", "writexl", "parallel",
    "patchwork", "gridExtra", "RColorBrewer", "survival", "survminer"
  )
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  # ggradar: not on CRAN for some R versions; install from GitHub when missing
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github("ricardo-bion/ggradar", upgrade = "never", quiet = TRUE)
  }
  log_note(paste0("R version: ", R.version.string))
  tryCatch(log_note(paste0("sessionInfo: ", paste(capture.output(sessionInfo()), collapse = "\n"))), error = function(e) log_note("sessionInfo capture skipped"))

  for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")) {
    tryCatch(do.call(Sys.setenv, setNames(list(as.character(n_workers)), v)), error = function(e) NULL)
  }
  log_note(paste0("Parallel config: using ", n_workers, " workers (all logical cores). BLAS/OpenMP threads set to ", n_workers, "."))
})

if (!isTRUE(PLOT_ONLY_MODE)) {
df <- run_step("Step 2: Read and standardize data", function() {
  df_raw <- readxl::read_xlsx(input_path)
  df_local <- df_raw %>%
    rename_with(~ stringr::str_replace_all(., "\\s+", "_")) %>%
    rename_with(~ toupper(.))

  required_candidates <- list(
    REPORT_ID = c("REPORT_ID", "PATIENT", "PATIENT_ID", "SAMPLE_ID", "MODEL", "MODEL_ID"),
    COMPOUND = c("COMPOUND", "DRUG", "COMPOUND_NAME"),
    SCORE = c("LEVEL", "SCORE", "PREDICTIVE_SCORE", "PRED_SCORE"),
    TIME_TO_DOUBLE = c(RESPONSE_COL_NAME, "TIMETODOUBLE", "TIME_TO_DOUBLE", "TIME_TO_DOUBLING", "TTD")
  )

  resolve_col <- function(candidates, df_cols) {
    hit <- candidates[candidates %in% df_cols]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  col_map <- sapply(required_candidates, resolve_col, df_cols = colnames(df_local))
  if (any(is.na(col_map))) {
    stop(paste("Missing required columns:", paste(names(col_map)[is.na(col_map)], collapse = ", ")))
  }

  has_response <- "RESPONSE_BESTAVGRESPONSE" %in% colnames(df_local)
  ## After rename (spaces→underscore, toupper), e.g. "Response_Best_Response" → RESPONSE_BEST_RESPONSE
  has_response_best <- "RESPONSE_BEST_RESPONSE" %in% colnames(df_local)
  df_local <- df_local %>%
    transmute(
      REPORT_ID = as.character(.data[[col_map["REPORT_ID"]]]),
      COMPOUND = as.character(.data[[col_map["COMPOUND"]]]),
      SCORE = as.numeric(.data[[col_map["SCORE"]]]),
      TIME_TO_DOUBLE = as.numeric(.data[[col_map["TIME_TO_DOUBLE"]]]),
      RESPONSE_BESTAVGRESPONSE = if (has_response) as.character(.data[["RESPONSE_BESTAVGRESPONSE"]]) else NA_character_,
      RESPONSE_BEST_RESPONSE = if (has_response_best) as.character(.data[["RESPONSE_BEST_RESPONSE"]]) else NA_character_
    ) %>%
    filter(!is.na(REPORT_ID), !is.na(COMPOUND), !is.na(SCORE), !is.na(TIME_TO_DOUBLE))

  if (has_response_best) {
    log_note("Optional column RESPONSE_BEST_RESPONSE (e.g. Excel: Response_Best_Response) retained in df / df_all_ranks for mRECIST Best_Response panels.")
  }
  log_note(paste0("Using all non-missing TimeToDouble rows: ", nrow(df_local), " rows, ", dplyr::n_distinct(df_local$REPORT_ID), " patients."))
  df_local
})
log_note(paste0("Post Step 2: ", nrow(df), " rows, ", dplyr::n_distinct(df$REPORT_ID), " patients."))


df_all_ranks <- run_step("Step 3: Per-patient normal ranks based on LEVEL (all data)", function() {
  ## Per patient: rank compounds by LEVEL descending, where highest LEVEL gets rank 1.
  ## Ties: equal LEVEL (e.g. same compound, different dose) share rank via dense_rank.
  ## Ranks >= RANK_POOL_FROM are pooled into one category (labeled RANK_POOL_LABEL, e.g. "7-10").
  ## Keep patients with >=2 compounds so rank-based comparisons remain meaningful.
  all_ranks <- df %>%
    group_by(REPORT_ID) %>%
    filter(n() >= 2L) %>%
    mutate(
      rank_bin = dplyr::dense_rank(-SCORE),
      rank_pooled = pmin(rank_bin, RANK_POOL_FROM),
      rank_scaled = 1 - (rank_pooled - 1) / (RANK_POOL_FROM - 1)
    ) %>%
    ungroup()

  # Defensive guard: ensure rank columns exist for downstream summaries/plots.
  if (!"rank_bin" %in% colnames(all_ranks)) {
    stop("Step 3 failed to create rank_bin (normal rank).")
  }
  if (!"rank_pooled" %in% colnames(all_ranks)) {
    stop("Step 3 failed to create rank_pooled.")
  }
  if (!"rank_scaled" %in% colnames(all_ranks)) {
    all_ranks <- all_ranks %>%
      mutate(rank_scaled = 1 - (rank_pooled - 1) / (RANK_POOL_FROM - 1))
  }

  coverage <- all_ranks %>%
    group_by(REPORT_ID) %>%
    summarise(n_ranks = n(), .groups = "drop")
  write.csv(coverage, file.path(analysis_dir, "rank_coverage_by_patient.csv"), row.names = FALSE)

  write.csv(all_ranks, file.path(analysis_dir, "data_all_ranks_normal_ranks.csv"), row.names = FALSE)

  if (nrow(all_ranks) == 0) {
    stop("No data after computing normal ranks. Cannot proceed.")
  }
  n_excl_single <- dplyr::n_distinct(df$REPORT_ID) - dplyr::n_distinct(all_ranks$REPORT_ID)
  if (n_excl_single > 0) {
    log_note(paste0("Excluded ", n_excl_single, " patients with only 1 compound (normal rank comparison undefined)."))
  }
  log_note(paste0("Normal ranks: ", nrow(all_ranks), " rows, ", dplyr::n_distinct(all_ranks$REPORT_ID), " patients (all other data retained)."))
  all_ranks
})


## rank_bin = raw per-patient rank; rank_pooled = ranks 1..6 plus pooled "7-10"; rank_scaled in [0,1] for smooth models.
rank_levels <- 1L:RANK_POOL_FROM
rank_labels <- c(as.character(1L:(RANK_POOL_FROM - 1L)), RANK_POOL_LABEL)
max_rank <- RANK_POOL_FROM
n_unique_ranks <- length(rank_levels)
log_note(paste0("Normal-rank scope: pooled ranks 1 to ", max_rank, " (1=top LEVEL; ", max_rank, " = ", RANK_POOL_LABEL, ")."))

scam_out <- run_step("Step 4: SCAM monotonicity test on normal ranks", function() {
  # rank_scaled in [0,1]: higher = better LEVEL (rank 1), lower = worse rank.
  scam_k <- 7L
  model_df <- df_all_ranks %>%
    mutate(REPORT_ID = factor(REPORT_ID)) %>%
    as.data.frame()

  # Guard against missing/invalid rank fields in incoming data.
  if (!"rank_pooled" %in% colnames(model_df) || all(is.na(model_df$rank_pooled))) {
    stop(paste0("Step 4 input is missing rank_pooled. Columns: ", paste(colnames(model_df), collapse = ", ")))
  }
  if (!"rank_scaled" %in% colnames(model_df) || all(is.na(model_df$rank_scaled))) {
    model_df <- model_df %>%
      dplyr::mutate(rank_scaled = 1 - (rank_pooled - 1) / (RANK_POOL_FROM - 1))
    log_note("Step 4: rank_scaled was missing/invalid and was reconstructed from rank_pooled.")
  }

  model <- scam::scam(
    TIME_TO_DOUBLE ~ s(rank_scaled, bs = "mpi", k = scam_k) + s(REPORT_ID, bs = "re"),
    data = model_df
  )
  sm <- summary(model)
  capture.output(sm, file = file.path(analysis_dir, "model_scam_all_ranks_summary.txt"))

  st_raw <- sm$s.table
  if (is.null(st_raw) || (is.matrix(st_raw) && (nrow(st_raw) == 0L || ncol(st_raw) == 0L))) {
    edf <- NA_real_
    p_val <- NA_real_
    s_tab <- tibble::tibble(edf = edf, p.value = p_val, term = "s(rank_scaled)")
  } else {
    s_tab <- as.data.frame(st_raw)
    s_tab$term <- rownames(s_tab)
    rownames(s_tab) <- NULL
    nms <- names(s_tab)
    if (!is.null(nms) && length(nms) == ncol(s_tab)) {
      names(s_tab) <- make.names(nms)
    }
    edf <- if ("edf" %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab$edf[1] else NA_real_
    pcol <- grep("^p", names(s_tab), ignore.case = TRUE, value = TRUE)[1]
    p_val <- if (!is.na(pcol) && pcol %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab[[pcol]][1] else NA_real_
  }

  pred_seq <- seq(0, 1, length.out = 50)
  pred_grid <- tibble::tibble(
    rank_scaled = pred_seq,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred <- predict(model, newdata = pred_grid, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)")
  pred_grid <- pred_grid %>%
    mutate(
      pred_ttd = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  rank_scaled_mid <- if (max_rank > 1L) 1 - (rank_levels - 1) / (max_rank - 1) else rep(1, length(rank_levels))
  pred_grid_ranks <- tibble::tibble(
    rank_scaled = rank_scaled_mid,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred_dec <- predict(model, newdata = pred_grid_ranks, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)")
  pred_grid_ranks <- pred_grid_ranks %>%
    mutate(
      rank_int = rank_levels,
      pred_ttd = as.numeric(pred_dec$fit),
      se = as.numeric(pred_dec$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  obs_rank <- model_df %>%
    group_by(.data$rank_pooled) %>%
    summarise(
      mean_ttd = mean(TIME_TO_DOUBLE, na.rm = TRUE),
      se_ttd = sd(TIME_TO_DOUBLE, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    rename(rank_int = rank_pooled)

  fig_w <- FIG_SCAM_MONO_WIDTH
  fig_h <- FIG_SCAM_MONO_HEIGHT
  plot_mono <- ggplot2::ggplot(obs_rank, ggplot2::aes(factor(rank_int, levels = rank_levels, labels = rank_labels), mean_ttd)) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd, ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd),
      data = obs_rank, width = 0.15
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd, label = paste0("(n = ", n, ")")),
      data = obs_rank, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(2.8)
    ) +
    ggplot2::geom_line(data = pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd), color = "#1f78b4", linewidth = 1, group = 1) +
    ggplot2::geom_ribbon(data = pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high),
                         inherit.aes = FALSE, alpha = 0.2, fill = "#1f78b4", group = 1) +
    ggplot2::theme_minimal(base_size = fsz(12)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 30, r = 5, b = 5, l = 5)) +
    ggplot2::labs(
      title = as_ggplot_label(paste0("SCAM monotonicity (normal ranks, all data): EDF = ", round(edf, 3), ", ", fmt_p(p_val))),
      x = "Rank within patient (1 = highest LEVEL)",
      y = paste("Avg", RESPONSE_COL_LABEL)
    )

  list(
    model = model,
    s_table = s_tab,
    edf = edf,
    p_value = p_val,
    pred_grid = pred_grid,
    pred_grid_ranks = pred_grid_ranks,
    obs_rank = obs_rank
  )
})

wiggle_out <- run_step("Step 4a: Unconstrained GAM visualization (wiggly reference)", function() {
  model_df <- df_all_ranks %>%
    mutate(REPORT_ID = factor(REPORT_ID)) %>%
    as.data.frame()

  gam_k <- min(20L, max(4L, max_rank))  # k limited by number of distinct ranks (e.g. 7 when pooled)
  gam_model <- mgcv::gam(
    TIME_TO_DOUBLE ~ s(rank_scaled, bs = "tp", k = gam_k) + s(REPORT_ID, bs = "re"),
    data = model_df,
    method = "REML"
  )
  gam_sm <- summary(gam_model)
  capture.output(gam_sm, file = file.path(analysis_dir, "model_gam_unconstrained_all_ranks_summary.txt"))

  gam_st <- gam_sm$s.table
  gam_edf <- NA_real_
  gam_p <- NA_real_
  if (!is.null(gam_st) && nrow(gam_st) >= 1L && ncol(gam_st) >= 1L) {
    gam_edf <- if ("edf" %in% colnames(gam_st)) gam_st[1, "edf"] else NA_real_
    pcol <- grep("^p", colnames(gam_st), ignore.case = TRUE, value = TRUE)[1]
    gam_p <- if (!is.na(pcol) && pcol %in% colnames(gam_st)) gam_st[1, pcol] else NA_real_
  }

  gam_result_label <- if (is.na(gam_p) || is.na(gam_edf)) {
    "Inconclusive: GAM summary statistics unavailable."
  } else if (gam_p < 0.05 && gam_edf > 1.1) {
    "Significant nonlinearity: unconstrained GAM indicates a wiggly rank-response pattern."
  } else if (gam_p < 0.05 && gam_edf <= 1.1) {
    "Significant but near-linear trend: unconstrained GAM supports rank association with limited wiggle."
  } else {
    "No strong evidence of nonlinear rank effect in the unconstrained GAM."
  }

  rank_scaled_mid <- if (max_rank > 1L) 1 - (rank_levels - 1) / (max_rank - 1) else rep(1, length(rank_levels))
  pred_grid_gam <- tibble::tibble(
    rank_scaled = rank_scaled_mid,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred_gam <- predict(gam_model, newdata = pred_grid_gam, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)")
  pred_grid_gam <- pred_grid_gam %>%
    mutate(
      rank_int = rank_levels,
      pred_ttd = as.numeric(pred_gam$fit),
      se = as.numeric(pred_gam$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )

  fig_w <- FIG_F22_JPG_WIDTH
  fig_h <- FIG_F22_JPG_HEIGHT
  tryCatch({
    obs_rank_4a <- scam_out$obs_rank
    y_min_4a <- min(obs_rank_4a$mean_ttd - 1.96 * obs_rank_4a$se_ttd, scam_out$pred_grid_ranks$ci_low, pred_grid_gam$ci_low, na.rm = TRUE)
    y_max_4a <- max(obs_rank_4a$mean_ttd + 1.96 * obs_rank_4a$se_ttd, scam_out$pred_grid_ranks$ci_high, pred_grid_gam$ci_high, na.rm = TRUE)
    y_pad_4a <- 0.16 * (y_max_4a - y_min_4a + 1e-8)
    y_lab_off_4a <- 0.04 * (y_max_4a - y_min_4a + 1e-8)
    p_wiggle <- ggplot2::ggplot(obs_rank_4a, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
      ggplot2::geom_point(size = 1.8, color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd), width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd + y_lab_off_4a, label = paste0("(n = ", n, ")")),
        data = obs_rank_4a, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(3.7), color = "black"
      ) +
      ggplot2::geom_ribbon(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
        inherit.aes = FALSE, alpha = 0.12) +
      ggplot2::geom_line(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"), linewidth = 1.05) +
      ggplot2::geom_ribbon(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
        inherit.aes = FALSE, alpha = 0.08) +
      ggplot2::geom_line(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"),
        linewidth = 1.05, linetype = "22") +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = "#2166ac", "GAM" = "#b2182b")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = "#2166ac", "GAM" = "#b2182b"), guide = "none") +
      f22_scam_gam_caption_gam_layer(gam_edf, gam_p, corner = "top_right", size = fsz(3.8)) +
      f22_scam_gam_caption_scam_layer(scam_out$edf, scam_out$p_value, corner = "top_right", size = fsz(3.8)) +
      ggplot2::theme_minimal(base_size = fsz(15)) +
      ggplot2::theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.text = ggplot2::element_text(color = "black", size = fsz(12)),
        plot.title = ggplot2::element_text(color = "black", face = "bold", size = fsz(16)),
        plot.subtitle = ggplot2::element_text(color = "black", size = fsz(13)),
        axis.title = ggplot2::element_text(color = "black", size = fsz(14)),
        axis.text = ggplot2::element_text(color = "black", size = fsz(12)),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.ticks.length = grid::unit(0.2, "cm"),
        plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 10)),
        plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 12)),
        plot.margin = ggplot2::margin(t = 30, r = 8, b = 8, l = 8)
      ) +
      ggplot2::labs(
        title = NULL,
        subtitle = NULL,
        x = "DDA rank groups",
        y = "Avg PFS (days)"
      ) +
      ggplot2::coord_cartesian(ylim = c(y_min_4a, y_max_4a + y_pad_4a), clip = "off")

    ## F22b JPEG: SCAM curve linewidth = half of F22; y-axis label PFS; oxford colors; ylim 0–60
    scam_lw_f22b_jpg <- 1.05 / 2
    p_wiggle_ox <- ggplot2::ggplot(obs_rank_4a, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
      ggplot2::geom_point(size = 1.8, color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd), width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd + y_lab_off_4a, label = paste0("(n = ", n, ")")),
        data = obs_rank_4a, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(3.7), color = "black"
      ) +
      ggplot2::geom_ribbon(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
        inherit.aes = FALSE, alpha = 0.12) +
      ggplot2::geom_line(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"), linewidth = scam_lw_f22b_jpg) +
      ggplot2::geom_ribbon(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
        inherit.aes = FALSE, alpha = 0.08) +
      ggplot2::geom_line(data = pred_grid_gam, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"),
        linewidth = 1.05, linetype = "dashed") +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4"), breaks = c("SCAM", "GAM")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = "#002147", "GAM" = "#C9B6E4"), guide = "none") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed"), linewidth = c(1.2, 1.4)), reverse = FALSE)) +
      f22_scam_gam_caption_gam_layer(gam_edf, gam_p, corner = "f22b_ranks", size = fsz(3.8)) +
      f22_scam_gam_caption_scam_layer(scam_out$edf, scam_out$p_value, corner = "f22b_ranks", size = fsz(3.8)) +
      ggplot2::theme_minimal(base_size = fsz(15)) +
      ggplot2::theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.key.width = grid::unit(1.3, "cm"),
        legend.text = ggplot2::element_text(color = "black", size = fsz(12)),
        plot.title = ggplot2::element_text(color = "black", face = "bold", size = fsz(16)),
        plot.subtitle = ggplot2::element_text(color = "black", size = fsz(13)),
        axis.title = ggplot2::element_text(color = "black", size = fsz(14)),
        axis.text = ggplot2::element_text(color = "black", size = fsz(12)),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.ticks.length = grid::unit(0.2, "cm"),
        plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 10)),
        plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 12)),
        plot.margin = ggplot2::margin(t = 30, r = 8, b = 8, l = 8)
      ) +
      ggplot2::labs(
        title = NULL,
        subtitle = NULL,
        x = "DDA rank groups",
        y = "Avg PFS (days)"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 66), clip = "off")
  }, error = function(e) {
    log_note(paste0("Step 4a: Skipped SCAM vs GAM plot (plot/ggsave error): ", conditionMessage(e)))
  })
  write.csv(pred_grid_gam, file.path(analysis_dir, "gam_unconstrained_pred_by_rank.csv"), row.names = FALSE)
  gam_p_txt <- if (is.na(gam_p)) {
    "NA"
  } else if (gam_p < 1e-100) {
    "p < 1e-100"
  } else if (gam_p > 0) {
    paste0("p = ", format(gam_p, scientific = TRUE, digits = 16))
  } else {
    # p-value can underflow to exact 0; recover a scientific string from log-scale tail probability.
    log_p_alt <- NA_real_
    stat_col <- c("Chi.sq", "Chisq", "F")
    stat_name <- stat_col[stat_col %in% colnames(gam_st)][1]
    ref_df <- if ("Ref.df" %in% colnames(gam_st)) as.numeric(gam_st[1, "Ref.df"]) else NA_real_
    stat_val <- if (!is.na(stat_name)) as.numeric(gam_st[1, stat_name]) else NA_real_
    if (!is.na(stat_val) && !is.na(ref_df) && ref_df > 0) {
      if (identical(stat_name, "F")) {
        df2 <- if (!is.null(gam_sm$residual.df)) as.numeric(gam_sm$residual.df) else NA_real_
        if (!is.na(df2) && df2 > 0) {
          log_p_alt <- stats::pf(stat_val, df1 = ref_df, df2 = df2, lower.tail = FALSE, log.p = TRUE)
        }
      } else {
        log_p_alt <- stats::pchisq(stat_val, df = ref_df, lower.tail = FALSE, log.p = TRUE)
      }
    }
    p_alt_txt <- format_logp_sci(log_p_alt, sig_digits = 6L)
    p_alt_num <- suppressWarnings(as.numeric(p_alt_txt))
    if (!is.na(p_alt_num) && p_alt_num >= 1e-100) paste0("p = ", p_alt_txt) else "p < 1e-100"
  }

  writeLines(
    c(
      paste0("Timestamp: ", timestamp_now()),
      paste0("GAM EDF: ", round(gam_edf, 4)),
      paste0("GAM p-value: ", gam_p_txt)
    ),
    file.path(analysis_dir, "gam_unconstrained_interpretation.txt")
  )

  list(
    model = gam_model,
    edf = gam_edf,
    p_value = gam_p,
    pred_grid = pred_grid_gam
  )
})

log_scam_gam_pvalues(scam_out, wiggle_out, analysis_dir, "normal_rank_data_modeling")

## Step 4b inactivated for a while (do not delete - will be needed shortly).
# scam_sensitivity <- run_step("Step 4b: SCAM k-sensitivity (multiple k in parallel)", function() {
#   model_df <- df_all_ranks %>%
#     mutate(REPORT_ID = factor(REPORT_ID)) %>%
#     as.data.frame()
#
#   fit_scam_k <- function(k_val) {
#     m <- scam::scam(
#       TIME_TO_DOUBLE ~ s(rank_scaled, bs = "mpi", k = k_val) + s(REPORT_ID, bs = "re"),
#       data = model_df
#     )
#     st <- summary(m)$s.table
#     edf_k <- NA_real_
#     pval <- NA_real_
#     if (!is.null(st) && nrow(st) >= 1L && ncol(st) >= 1L) {
#       edf_k <- if ("edf" %in% colnames(st)) st[1, "edf"] else NA_real_
#       pcol <- grep("^p", colnames(st), ignore.case = TRUE, value = TRUE)[1]
#       pval <- if (!is.na(pcol) && pcol %in% colnames(st)) st[1, pcol] else NA_real_
#     }
#     tibble::tibble(k = k_val, edf = edf_k, p_value = pval)
#   }
#
#   ## Use only k=7 for speed; lower k = more conservative (rigid) smooth = more trustworthy
#   k_vals <- c(min(7L, n_unique_ranks))
#   n_k <- min(n_workers, length(k_vals))
#   if (n_k > 1) {
#     cl <- parallel::makeCluster(n_k)
#     on.exit(parallel::stopCluster(cl), add = TRUE)
#     parallel::clusterExport(cl, "model_df", envir = environment())
#     parallel::clusterEvalQ(cl, library(scam))
#     sensitivity_list <- parallel::parLapply(cl, k_vals, fit_scam_k)
#   } else {
#     sensitivity_list <- lapply(k_vals, fit_scam_k)
#   }
#   sens_df <- dplyr::bind_rows(sensitivity_list)
#   write.csv(sens_df, file.path(analysis_dir, "scam_k_sensitivity.csv"), row.names = FALSE)
#   log_note(paste0("SCAM k-sensitivity: ", paste(sprintf("k=%d EDF=%.3f p=%s", sens_df$k, sens_df$edf, format.pval(sens_df$p_value, digits = 2)), collapse = "; ")))
#   sens_df
# })
scam_sensitivity <- tibble::tibble(k = NA_integer_, edf = NA_real_, p_value = NA_real_)

## Mixed model (publication): TIME_TO_DOUBLE ~ rank_f + (1|REPORT_ID); rank 1 = reference.
## P-values: lmerTest summary (Satterthwaite). Forest CIs: confint(profile), Wald if profiling fails.
effects_out <- run_step("Mixed-model effect sizes", function() {
  df_eff <- df_all_ranks %>% mutate(rank_f = factor(rank_pooled, levels = rank_levels))
  model <- lmerTest::lmer(TIME_TO_DOUBLE ~ rank_f + (1 | REPORT_ID), data = df_eff)
  sm <- summary(model)
  capture.output(sm, file = file.path(analysis_dir, "model_lmer_all_ranks_summary.txt"))

  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>%
    filter(grepl("^rank_f", term)) %>%
    mutate(
      rank_int = suppressWarnings(as.integer(gsub("^rank_f", "", term))),
      estimate = Estimate,
      std_error = `Std. Error`,
      p_value = `Pr(>|t|)`
    ) %>%
    dplyr::select(rank_int, estimate, std_error, p_value)

  ci <- tryCatch(
    suppressMessages(confint(model, method = "profile")),
    error = function(e) {
      log_note(paste0("Profile CIs failed (", conditionMessage(e), "); falling back to Wald."))
      suppressMessages(confint(model, method = "Wald"))
    }
  )
  ci_df <- as.data.frame(ci)
  ci_df$term <- rownames(ci_df)
  rownames(ci_df) <- NULL
  cn <- colnames(ci_df)
  cilo <- cn[grep("2\\.5|lower", cn, ignore.case = TRUE)[1]]
  cihi <- cn[grep("97\\.5|upper", cn, ignore.case = TRUE)[1]]
  if (is.na(cilo)) cilo <- cn[1]
  if (is.na(cihi)) cihi <- cn[2]
  ci_df <- ci_df %>%
    filter(grepl("^rank_f", term)) %>%
    mutate(rank_int = as.integer(gsub("^rank_f", "", term)), ci_low = .data[[cilo]], ci_high = .data[[cihi]]) %>%
    dplyr::select(rank_int, ci_low, ci_high)

  effects <- coef_df %>%
    left_join(ci_df, by = "rank_int") %>%
    bind_rows(tibble::tibble(rank_int = 1L, estimate = 0, std_error = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_)) %>%
    arrange(rank_int) %>%
    mutate(
      p_label = ifelse(is.na(p_value), "ref", fmt_p(p_value)),
      p_label_plot = ifelse(is.na(p_value), "ref", fmt_p_plot(p_value))
    )
  idx_eff <- !is.na(effects$p_value)
  effects$p_adj_fdr <- NA_real_
  effects$p_adj_bonf <- NA_real_
  if (any(idx_eff)) {
    effects$p_adj_fdr[idx_eff] <- p.adjust(effects$p_value[idx_eff], method = "BH")
    effects$p_adj_bonf[idx_eff] <- p.adjust(effects$p_value[idx_eff], method = "bonferroni")
  }

  write.csv(effects, file.path(analysis_dir, "effect_sizes_all_ranks_vs_top_rank.csv"), row.names = FALSE)

  x_max <- max(effects$ci_high, na.rm = TRUE)
  x_pad <- 0.06 * (max(effects$ci_high, na.rm = TRUE) - min(effects$ci_low, na.rm = TRUE) + 1e-8)
  eff_height <- max(6, 4 + max_rank * 0.2)
  label_size_eff <- max(2, min(3.5, 3.5 - max_rank * 0.02))
  p_fx <- ggplot2::ggplot(effects, ggplot2::aes(x = estimate, y = factor(rank_int, levels = rev(rank_levels), labels = rev(rank_labels)))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0.2, orientation = "y", color = "#1f78b4") +
    ggplot2::geom_point(ggplot2::aes(color = rank_int == 1L), size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(x = x_max + x_pad, label = p_label_plot), parse = TRUE, hjust = 0, size = label_size_eff) +
    ggplot2::scale_color_manual(values = c("FALSE" = "#1f78b4", "TRUE" = "black"), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Effect sizes by rank vs rank 1 top (mixed model)",
      x = "Estimated TimeToDouble difference vs rank 1 top",
      y = "Rank (1 = highest LEVEL)"
    ) +
    ggplot2::coord_cartesian(xlim = c(min(effects$ci_low, na.rm = TRUE), x_max + 3 * x_pad))

  effects
})

if (isTRUE(SAVE_PLOT_DATA_CACHE)) {
  run_step("Step 8a: Save cached data bundle for plotting-only reruns", function() {
    plot_cache <- list(
      df = as.data.frame(df),
      df_all_ranks = as.data.frame(df_all_ranks),
      scam_out = scam_out,
      wiggle_out = wiggle_out,
      effects_out = as.data.frame(effects_out),
      rank_levels = rank_levels,
      rank_labels = rank_labels,
      max_rank = max_rank
    )
    saveRDS(plot_cache, plot_data_cache_path)
    log_note(paste0("Saved plot data cache: ", plot_data_cache_path))
  })
}
} else {
  run_step("Step 2-8a bypass: Load cached data bundle for plotting", function() {
    if (!file.exists(plot_data_cache_path)) {
      stop(paste0("PLOT_ONLY_MODE=TRUE but cache not found: ", plot_data_cache_path,
                  ". Run once with PLOT_ONLY_MODE=FALSE and SAVE_PLOT_DATA_CACHE=TRUE."))
    }
    plot_cache <- readRDS(plot_data_cache_path)
    validate_plot_cache(plot_cache)
    assign("df", as.data.frame(plot_cache$df), envir = .GlobalEnv)
    assign("df_all_ranks", as.data.frame(plot_cache$df_all_ranks), envir = .GlobalEnv)
    assign("scam_out", plot_cache$scam_out, envir = .GlobalEnv)
    assign("wiggle_out", plot_cache$wiggle_out, envir = .GlobalEnv)
    assign("effects_out", as.data.frame(plot_cache$effects_out), envir = .GlobalEnv)
    assign("rank_levels", as.integer(plot_cache$rank_levels), envir = .GlobalEnv)
    if ("rank_labels" %in% names(plot_cache)) {
      assign("rank_labels", as.character(plot_cache$rank_labels), envir = .GlobalEnv)
    } else {
      assign("rank_labels", as.character(plot_cache$rank_levels), envir = .GlobalEnv)
    }
    assign("max_rank", as.integer(plot_cache$max_rank), envir = .GlobalEnv)
    assign("n_unique_ranks", length(as.integer(plot_cache$rank_levels)), envir = .GlobalEnv)
    log_note(paste0("Loaded plot data cache: ", plot_data_cache_path))
  })
}

run_step("Step 9: Publication-quality figures (peer review & PI perspective)", function() {
  ca_registry_reset()

  save_fig <- function(plot_obj, path, width, height, dpi = 300) {
    tryCatch(
      {
        # Save PNG (or whatever extension is in 'path')
        ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
        log_note(paste0("Saved: ", basename(path)))

        # Also save an SVG version in the same folder (figures)
        base_no_ext <- sub("\\.[^.]*$", "", basename(path))
        svg_path <- file.path(dirname(path), paste0(base_no_ext, ".svg"))
        ggplot2::ggsave(svg_path, plot_obj, width = width, height = height, dpi = 320, device = "svg")
        log_note(paste0("Saved: ", basename(svg_path)))
      },
      error = function(e) {
        log_note(paste0("SKIP figure ", basename(path), " (error): ", conditionMessage(e)))
      }
    )
  }

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

  ## mRECIST / Durable Benefit: horizontal y-axis aid lines only (no vertical grid).
  theme_pub_nogrid <- theme_pub +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "grey85")
    )

  # Keep F23 panel text at current approved sizes.
  theme_f23_lock <- ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = 13, color = "black"),
    axis.title = ggplot2::element_text(size = 14, color = "black"),
    axis.text = ggplot2::element_text(size = 12, color = "black"),
    legend.text = ggplot2::element_text(size = 12, color = "black"),
    legend.title = ggplot2::element_text(size = 12, color = "black")
  )


  ## SCAM residual diagnostics
  res_scam <- residuals(scam_out$model)
  pred_scam <- fitted(scam_out$model)
  res_df <- tibble::tibble(fitted = pred_scam, residual = res_scam)
  qq_df <- tibble::tibble(residual = res_scam)
  p6a <- ggplot2::ggplot(qq_df, ggplot2::aes(sample = residual)) +
    ggplot2::stat_qq(alpha = 0.5) +
    ggplot2::stat_qq_line(color = "#b2182b", linewidth = 0.8) +
    ggplot2::labs(title = "A. Normal Q-Q plot", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_pub
  p6b <- ggplot2::ggplot(res_df, ggplot2::aes(fitted, residual)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_smooth(se = TRUE, color = "#b2182b", alpha = 0.2) +
    ggplot2::labs(title = "B. Residuals vs fitted", x = paste("Fitted", RESPONSE_COL_LABEL), y = "Residual") +
    theme_pub
  p6 <- p6a + p6b + patchwork::plot_layout(ncol = 2)
  ## F22: SCAM vs GAM with CI (standalone, blue-red and oxfordblue-purple)
  obs_rank_f22 <- scam_out$obs_rank
  y_min_f22 <- min(obs_rank_f22$mean_ttd - 1.96 * obs_rank_f22$se_ttd, scam_out$pred_grid_ranks$ci_low, wiggle_out$pred_grid$ci_low, na.rm = TRUE)
  y_max_f22 <- max(obs_rank_f22$mean_ttd + 1.96 * obs_rank_f22$se_ttd, scam_out$pred_grid_ranks$ci_high, wiggle_out$pred_grid$ci_high, na.rm = TRUE)
  y_pad_f22 <- 0.16 * (y_max_f22 - y_min_f22 + 1e-8)
  y_lab_off_f22 <- 0.04 * (y_max_f22 - y_min_f22 + 1e-8)
  make_f22_plot <- function(scam_color, gam_color, text_mult = 1, line_width_mult = 1, scam_line_width_mult = line_width_mult, body_text_extra = 1, caption_corner = "top_right") {
    ggplot2::ggplot(obs_rank_f22, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
      ggplot2::geom_point(size = 1.8, color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd), width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd + 1.96 * se_ttd + y_lab_off_f22, label = paste0("(n = ", n, ")")),
        data = obs_rank_f22, angle = 90, hjust = -0.1, vjust = 0.5, size = fsz(3.7 * text_mult * body_text_extra), color = "black"
      ) +
      ggplot2::geom_ribbon(
        data = scam_out$pred_grid_ranks,
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
        inherit.aes = FALSE, alpha = 0.12
      ) +
      ggplot2::geom_line(data = scam_out$pred_grid_ranks, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"), linewidth = 1.05 * scam_line_width_mult) +
      ggplot2::geom_ribbon(
        data = wiggle_out$pred_grid,
        ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
        inherit.aes = FALSE, alpha = 0.08
      ) +
      ggplot2::geom_line(data = wiggle_out$pred_grid, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"), linewidth = 1.05 * line_width_mult, linetype = "22") +
      ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = scam_color, "GAM" = gam_color), breaks = c("SCAM", "GAM")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = scam_color, "GAM" = gam_color), guide = "none") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
      f22_scam_gam_caption_gam_layer(
        wiggle_out$edf, wiggle_out$p_value,
        corner = caption_corner, size = fsz(3.8 * text_mult * body_text_extra),
        line_gap = 3.5 * text_mult * body_text_extra
      ) +
      f22_scam_gam_caption_scam_layer(
        scam_out$edf, scam_out$p_value,
        corner = caption_corner, size = fsz(3.8 * text_mult * body_text_extra),
        line_gap = 3.5 * text_mult * body_text_extra
      ) +
      ggplot2::labs(
        title = NULL,
        subtitle = NULL,
        x = "DDA rank groups",
        y = "Avg PFS (days)"
      ) +
      theme_pub +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = fsz(16 * text_mult), color = "black"),
        plot.subtitle = ggplot2::element_text(size = fsz(13 * text_mult), color = "black"),
        axis.title = ggplot2::element_text(size = fsz(14 * text_mult * body_text_extra), color = "black"),
        axis.text = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.text = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.title = ggplot2::element_text(size = fsz(12 * text_mult * body_text_extra), color = "black"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      ggplot2::coord_cartesian(ylim = c(y_min_f22, y_max_f22 + y_pad_f22), clip = "off")
  }
  f22b_text_mult <- 1.5 * 1.2
  f22b_n_top <- max(
    obs_rank_f22$mean_ttd + 1.96 * obs_rank_f22$se_ttd + y_lab_off_f22 * f22b_text_mult,
    na.rm = TRUE
  )
  f22b_y_hi <- max(66, f22b_n_top + F22B_NLABEL_Y_PAD)
  p22_blue_red <- make_f22_plot("#2166ac", "#b2182b")

  p22_ox_purple <- make_f22_plot("#002147", "#006400", text_mult = 1.5, line_width_mult = 2, scam_line_width_mult = 1, body_text_extra = 1.2, caption_corner = "f22b_ranks") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 3)),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.width = grid::unit(1.3, "cm"),
      plot.margin = FIG_F22B_MARGIN
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed"), linewidth = c(1.2, 1.4)))) +
    ggplot2::coord_cartesian(ylim = c(0, f22b_y_hi), clip = "off", expand = FALSE)

  ## Effect forest — FDR-only labels (reference = rank 1 top)
  effects_plot <- effects_out %>%
    mutate(
      label_f8 = ifelse(is.na(p_adj_fdr), "ref", paste0("FDR = ", formatC(p_adj_fdr, format = "e", digits = 2)))
    )
  x_eff_min <- min(effects_plot$ci_low, effects_plot$estimate, na.rm = TRUE)
  x_eff_max <- max(effects_plot$ci_high, effects_plot$estimate, na.rm = TRUE)
  x_span <- max(x_eff_max - x_eff_min, 0.1)
  x_pad <- 0.25 * x_span
  x_label_pos <- x_eff_max + x_pad
  x_plot_min <- x_eff_min - 0.05 * x_span
  x_plot_max <- x_eff_max + x_pad * 3
  label_size_f8 <- fsz(4 * 1.5 * 1.2)
  p8 <- ggplot2::ggplot(effects_plot, ggplot2::aes(x = estimate, y = factor(rank_int, levels = rev(rank_levels), labels = rev(rank_labels)))) +
    ggplot2::geom_vline(xintercept = 0, color = "gray70", linewidth = 0.9) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0.252, linewidth = 1.44, color = "black", orientation = "y") +
    ggplot2::geom_point(size = 7.21, color = "black") +
    ggplot2::geom_text(ggplot2::aes(x = x_label_pos, label = label_f8), hjust = 0, size = label_size_f8, color = "black") +
    ggplot2::labs(
      title = NULL,
      x = "Effect size (avg PFS) (mixed model)",
      y = "Rank bins"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_plot_min, x_plot_max), clip = "off") +
    theme_pub
  fig8_h <- 6
  p8 <- p8 + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = fsz(16 * 1.5), color = "black"),
    axis.title = ggplot2::element_text(size = fsz(14 * 1.5 * 1.2), color = "black"),
    axis.text = ggplot2::element_text(size = fsz(12 * 1.5 * 1.2), color = "black"),
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.8),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = grid::unit(0.24, "cm"),
    legend.position = "none"
  )

  ## F88: Disease control rate (DCR) by rank – responder = CR/PR/SD, non-responder = PD (Response_BestAvgResponse)
  if ("RESPONSE_BESTAVGRESPONSE" %in% colnames(df_all_ranks)) {
    dcr_by_rank <- df_all_ranks %>%
      dplyr::mutate(resp = trimws(toupper(as.character(RESPONSE_BESTAVGRESPONSE)))) %>%
      dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
      dplyr::mutate(is_dcr = resp %in% c("CR", "PR", "SD")) %>%
      dplyr::group_by(rank_pooled) %>%
      dplyr::summarise(n_dcr = sum(is_dcr, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        dcr_pct = 100 * n_dcr / n,
        se_pct = dplyr::if_else(n > 0, 100 * sqrt((dcr_pct / 100) * (1 - dcr_pct / 100) / n), 0)
      )
    n_ranks_f88 <- length(rank_levels)
    f88_colors <- setNames(
      rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(max(n_ranks_f88, 1))),
      as.character(rank_levels)
    )
    ## F99C: ORR by rank (ORR = CR + PR over evaluable CR/PR/SD/PD), stacked CR/PR composition
    orr_by_rank <- df_all_ranks %>%
      dplyr::mutate(resp = trimws(toupper(as.character(RESPONSE_BESTAVGRESPONSE)))) %>%
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
    if (nrow(orr_by_rank) > 0) {
      n_tot_orr <- sum(orr_by_rank$n)
      r_tot_orr <- sum(orr_by_rank$n_orr)
      p_overall_orr <- if (n_tot_orr > 0) r_tot_orr / n_tot_orr else 0
      ca_ann_f99c <- ca_trend_annotate(
        scores = orr_by_rank$rank_pooled, n = orr_by_rank$n, success = orr_by_rank$n_orr,
        metric = NULL, panel_id = "F99C_ORR_by_rank", analysis = "normal_rank", stratification = "rank_pooled"
      )
      trend_lab_orr <- ca_ann_f99c$label
      # Logistic regression (binomial) fit for ORR trend line in F99C.
      # Response uses counts: success=n_orr, failure=(n - n_orr).
      lm_fit_orr <- stats::lm(orr_pct ~ rank_pooled, data = orr_by_rank)  # fallback
      logit_fit_orr <- tryCatch(
        stats::glm(cbind(n_orr, n - n_orr) ~ rank_pooled, family = stats::binomial(), data = orr_by_rank),
        error = function(e) NULL
      )
      trend_df_orr <- data.frame(
        rank_pooled = rank_levels,
        trend_fit = if (!is.null(logit_fit_orr)) {
          100 * stats::predict(logit_fit_orr, newdata = data.frame(rank_pooled = rank_levels), type = "response")
        } else {
          stats::predict(lm_fit_orr, newdata = data.frame(rank_pooled = rank_levels))
        }
      )

      orr_components <- df_all_ranks %>%
        dplyr::mutate(resp = trimws(toupper(as.character(RESPONSE_BESTAVGRESPONSE)))) %>%
        dplyr::filter(resp %in% c("CR", "PR")) %>%
        dplyr::group_by(rank_pooled, resp) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::right_join(
          expand.grid(rank_pooled = rank_levels, resp = c("CR", "PR"), stringsAsFactors = FALSE),
          by = c("rank_pooled", "resp")
        ) %>%
        dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
        dplyr::group_by(rank_pooled) %>%
        dplyr::mutate(n_orr_rank = sum(n)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(orr_by_rank %>% dplyr::select(rank_pooled, orr_pct), by = "rank_pooled") %>%
        dplyr::mutate(
          height = dplyr::if_else(n_orr_rank > 0, orr_pct * n / n_orr_rank, 0)
        )
      orr_colors <- c("CR" = "#81c784", "PR" = "#90caf9")
      orr_components <- orr_components %>%
        dplyr::mutate(resp = factor(resp, levels = c("CR", "PR")))
      rank_totals_f99c <- orr_by_rank %>%
        dplyr::select(rank_pooled, n_orr, n, orr_pct, se_pct) %>%
        dplyr::mutate(
          ci_low = pmax(0, orr_pct - 1.96 * se_pct),
          ci_high = pmin(100, orr_pct + 1.96 * se_pct),
          label_resp_total = paste0("(", n_orr, "/", n, ")")
        )
      y_n_label_f99c <- min(100, max(rank_totals_f99c$ci_high, na.rm = TRUE) + 4)
        ## F99D: mRECIST by rank (stacked CR/PR/SD/PD) + ORR/DCR logistic regression curves
      mrecist_by_rank <- df_all_ranks %>%
        dplyr::mutate(resp = trimws(toupper(as.character(RESPONSE_BESTAVGRESPONSE)))) %>%
        dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
        dplyr::group_by(rank_pooled, resp) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(rank_pooled) %>%
        dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
        dplyr::ungroup() %>%
        # Stacking bottom→top: PD, SD, PR, CR (PR below CR in the bar).
        dplyr::mutate(resp = factor(resp, levels = c("PD", "SD", "PR", "CR")))

      # DCR (CR+PR+SD) trend from binomial logistic regression
      dcr_ordered_f99d <- dcr_by_rank %>% dplyr::arrange(rank_pooled)
      lm_fit_dcr_f99d <- stats::lm(dcr_pct ~ rank_pooled, data = dcr_ordered_f99d)  # fallback
      # Use a linear trend for the visual curve (CA test is still computed separately).
      trend_df_dcr_f99d <- data.frame(
        rank_pooled = rank_levels,
        trend_fit_dcr = stats::predict(lm_fit_dcr_f99d, newdata = data.frame(rank_pooled = rank_levels))
      )

      # ORR (CR+PR) trend from binomial logistic regression
      orr_ordered_f99d <- orr_by_rank %>% dplyr::arrange(rank_pooled)
      lm_fit_orr_f99d <- stats::lm(orr_pct ~ rank_pooled, data = orr_ordered_f99d)  # fallback
      # Use a linear trend for the visual curve (CA test is still computed separately).
      trend_df_orr_f99d <- data.frame(
        rank_pooled = rank_levels,
        trend_fit_orr = stats::predict(lm_fit_orr_f99d, newdata = data.frame(rank_pooled = rank_levels))
      )

      # Cochran–Armitage trend test annotations for ORR and DCR (stats::prop.trend.test)
      p_overall_dcr_f99d <- sum(dcr_ordered_f99d$n_dcr) / sum(dcr_ordered_f99d$n)
      ca_ann_f99d_dcr <- ca_trend_annotate(
        scores = dcr_ordered_f99d$rank_pooled, n = dcr_ordered_f99d$n, success = dcr_ordered_f99d$n_dcr,
        metric = "DCR", panel_id = "F99D_mRECIST_by_rank_ORR_DCR", analysis = "normal_rank", stratification = "rank_pooled"
      )
      trend_lab_dcr_f99d <- ca_ann_f99d_dcr$label

      p_overall_orr_f99d <- sum(orr_ordered_f99d$n_orr) / sum(orr_ordered_f99d$n)
      ca_ann_f99d_orr <- ca_trend_annotate(
        scores = orr_ordered_f99d$rank_pooled, n = orr_ordered_f99d$n, success = orr_ordered_f99d$n_orr,
        metric = "ORR", panel_id = "F99D_mRECIST_by_rank_ORR_DCR", analysis = "normal_rank", stratification = "rank_pooled"
      )
      trend_lab_orr_f99d <- ca_ann_f99d_orr$label

      mrecist_fill_title_avg <- "mRECIST\n(BestAvgResponse)"
      mrecist_title_avg <- "Best Average Response by\nDrug Rank Within Tumor"
      mrecist_colors <- c(
        "CR" = "#1565c0",
        "PR" = "#bbdefb",
        "SD" = "#d4d0f5",
        "PD" = "#d9d9d9"
      )

      overall_dcr_pct_f99d <- 100 * p_overall_dcr_f99d
      overall_orr_pct_f99d <- 100 * p_overall_orr_f99d

      p99d <- ggplot2::ggplot(
        mrecist_by_rank,
        ggplot2::aes(
          x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
          y = pct,
          fill = resp
        )
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_stack(),
          width = 0.6,
          colour = "black",
          linewidth = 0.25
        ) +
        ggplot2::geom_line(
          data = trend_df_dcr_f99d,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_dcr, group = 1, color = "CA (DCR%)"),
          linewidth = 1,
          linetype = "solid",
          inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data = trend_df_orr_f99d,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_orr, group = 1, color = "CA (ORR%)"),
          linewidth = 1,
          linetype = "solid",
          inherit.aes = FALSE
        ) +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_orr_f99d, parse = TRUE, hjust = 1.05, vjust = 2.10, size = fsz(3), colour = "black") +
        ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_dcr_f99d, parse = TRUE, hjust = 1.05, vjust = 1.00, size = fsz(3), colour = "black") +
        ggplot2::scale_fill_manual(values = mrecist_colors, name = mrecist_fill_title_avg) +
        ggplot2::scale_color_manual(
          values = c("CA (DCR%)" = "#008000", "CA (ORR%)" = "#FFA500"),
          breaks = c("CA (DCR%)", "CA (ORR%)"),
          name = NULL
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(order = 1),
          color = ggplot2::guide_legend(order = 2)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
          x = "Drug Rank Within Tumor",
          y = "Best Response (%)",
          title = mrecist_title_avg
        ) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.box = "vertical",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )
      ann_f99d_orr <- build_adjacent_prop_annotations(
        summary_df = orr_ordered_f99d,
        rank_col = "rank_pooled",
        success_col = "n_orr",
        total_col = "n",
        trend_df = trend_df_orr_f99d,
        trend_rank_col = "rank_pooled",
        trend_value_col = "trend_fit_orr",
        layer_offset = 2.2,
        span_step = 2.8
      )
      ann_f99d_dcr <- build_adjacent_prop_annotations(
        summary_df = dcr_ordered_f99d,
        rank_col = "rank_pooled",
        success_col = "n_dcr",
        total_col = "n",
        trend_df = trend_df_dcr_f99d,
        trend_rank_col = "rank_pooled",
        trend_value_col = "trend_fit_dcr",
        layer_offset = 6.8,
        span_step = 3.1
      )
      p99d <- add_pairwise_signif_brackets(p99d, ann_f99d_orr)
      p99d <- add_pairwise_signif_brackets(p99d, ann_f99d_dcr)

      ## F456 (BestAvgResponse): F99D / RESPONSE_BESTAVGRESPONSE with F456 title layout (complement to Best_Response export).
      p456_avg <- p99d +
        ggplot2::labs(
          title = "Best Response by DDA Score Tiers",
          subtitle = "(BestAvgResponse data)"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 6), hjust = 0.5)
        )

      ## F99L: same stacked bars as F99D, no ORR/DCR trend lines, no Cochran–Armitage text
      mrecist_colors_f99l <- mrecist_colors
      p99l <- ggplot2::ggplot(
        mrecist_by_rank,
        ggplot2::aes(
          x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
          y = pct,
          fill = resp
        )
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_stack(),
          width = 0.6,
          colour = "black",
          linewidth = 0.25
        ) +
        ggplot2::scale_fill_manual(values = mrecist_colors_f99l, name = mrecist_fill_title_avg) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
          x = "Drug Rank Within Tumor",
          y = "Best Response (%)",
          title = mrecist_title_avg
        ) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      ## F99X: same as F99L (bars only); CR = blue, PR = light blue
      mrecist_colors_f99x <- c(
        "CR" = "#1565c0",
        "PR" = "#bbdefb",
        "SD" = "#d4d0f5",
        "PD" = "#d9d9d9"
      )
      p99x <- ggplot2::ggplot(
        mrecist_by_rank,
        ggplot2::aes(
          x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
          y = pct,
          fill = resp
        )
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_stack(),
          width = 0.6,
          colour = "black",
          linewidth = 0.25
        ) +
        ggplot2::scale_fill_manual(values = mrecist_colors_f99x, name = mrecist_fill_title_avg) +
        ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
          x = "Drug Rank Within Tumor",
          y = "Best Response (%)",
          title = mrecist_title_avg
        ) +
        theme_pub_nogrid +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
          legend.position = "right",
          legend.key.height = grid::unit(0.55, "cm"),
          legend.spacing.y = grid::unit(0.15, "cm"),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
        )

      ## F456 / F461–F463: mRECIST by rank using RESPONSE_BEST_RESPONSE (mirrors F99D, F99L, F99X)
      title_mrecist_br <- "Best Response by\nDrug Rank Within Tumor"
      mrecist_fill_title_br <- "mRECIST\n(BestResponse)"
      if ("RESPONSE_BEST_RESPONSE" %in% colnames(df_all_ranks)) {
        mrecist_by_rank_f456 <- df_all_ranks %>%
          dplyr::mutate(resp = trimws(toupper(as.character(.data[["RESPONSE_BEST_RESPONSE"]])))) %>%
          dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
          dplyr::group_by(rank_pooled, resp) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(rank_pooled) %>%
          dplyr::mutate(n_total = sum(n), pct = dplyr::if_else(n_total > 0, 100 * n / n_total, 0)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(resp = factor(resp, levels = c("PD", "SD", "PR", "CR")))

        dcr_by_rank_f456 <- df_all_ranks %>%
          dplyr::mutate(resp = trimws(toupper(as.character(.data[["RESPONSE_BEST_RESPONSE"]])))) %>%
          dplyr::filter(resp %in% c("CR", "PR", "SD", "PD")) %>%
          dplyr::mutate(is_dcr = resp %in% c("CR", "PR", "SD")) %>%
          dplyr::group_by(rank_pooled) %>%
          dplyr::summarise(n_dcr = sum(is_dcr, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(
            dcr_pct = 100 * n_dcr / n,
            se_pct = dplyr::if_else(n > 0, 100 * sqrt((dcr_pct / 100) * (1 - dcr_pct / 100) / n), 0)
          )

        orr_by_rank_f456 <- df_all_ranks %>%
          dplyr::mutate(resp = trimws(toupper(as.character(.data[["RESPONSE_BEST_RESPONSE"]])))) %>%
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

        if (nrow(mrecist_by_rank_f456) > 0L && nrow(dcr_by_rank_f456) > 0L && nrow(orr_by_rank_f456) > 0L) {
          dcr_ordered_f456 <- dcr_by_rank_f456 %>% dplyr::arrange(rank_pooled)
          lm_fit_dcr_f456 <- stats::lm(dcr_pct ~ rank_pooled, data = dcr_ordered_f456)
          trend_df_dcr_f456 <- data.frame(
            rank_pooled = rank_levels,
            trend_fit_dcr = stats::predict(lm_fit_dcr_f456, newdata = data.frame(rank_pooled = rank_levels))
          )
          orr_ordered_f456 <- orr_by_rank_f456 %>% dplyr::arrange(rank_pooled)
          lm_fit_orr_f456 <- stats::lm(orr_pct ~ rank_pooled, data = orr_ordered_f456)
          trend_df_orr_f456 <- data.frame(
            rank_pooled = rank_levels,
            trend_fit_orr = stats::predict(lm_fit_orr_f456, newdata = data.frame(rank_pooled = rank_levels))
          )
          ca_ann_f456_dcr <- ca_trend_annotate(
            scores = dcr_ordered_f456$rank_pooled, n = dcr_ordered_f456$n, success = dcr_ordered_f456$n_dcr,
            metric = "DCR", panel_id = "mRECIST_Best_Response", analysis = "normal_rank", stratification = "rank_pooled"
          )
          trend_lab_dcr_f456 <- ca_ann_f456_dcr$label
          ca_ann_f456_orr <- ca_trend_annotate(
            scores = orr_ordered_f456$rank_pooled, n = orr_ordered_f456$n, success = orr_ordered_f456$n_orr,
            metric = "ORR", panel_id = "mRECIST_Best_Response", analysis = "normal_rank", stratification = "rank_pooled"
          )
          trend_lab_orr_f456 <- ca_ann_f456_orr$label
          ## F461: mirror F99D; F456 keeps same ggplot object
          p461 <- ggplot2::ggplot(
            mrecist_by_rank_f456,
            ggplot2::aes(
              x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
              y = pct,
              fill = resp
            )
          ) +
            ggplot2::geom_col(
              position = ggplot2::position_stack(),
              width = 0.6,
              colour = "black",
              linewidth = 0.25
            ) +
            ggplot2::geom_line(
              data = trend_df_dcr_f456,
              ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_dcr, group = 1, color = "CA (DCR%)"),
              linewidth = 1,
              linetype = "solid",
              inherit.aes = FALSE
            ) +
            ggplot2::geom_line(
              data = trend_df_orr_f456,
              ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_orr, group = 1, color = "CA (ORR%)"),
              linewidth = 1,
              linetype = "solid",
              inherit.aes = FALSE
            ) +
            ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_orr_f456, parse = TRUE, hjust = 1.05, vjust = 2.10, size = fsz(3), colour = "black") +
            ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_dcr_f456, parse = TRUE, hjust = 1.05, vjust = 1.00, size = fsz(3), colour = "black") +
            ggplot2::scale_fill_manual(values = mrecist_colors, name = mrecist_fill_title_br) +
            ggplot2::scale_color_manual(
              values = c("CA (DCR%)" = "#008000", "CA (ORR%)" = "#FFA500"),
              breaks = c("CA (DCR%)", "CA (ORR%)"),
              name = NULL
            ) +
            ggplot2::guides(
              fill = ggplot2::guide_legend(order = 1),
              color = ggplot2::guide_legend(order = 2)
            ) +
            ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::labs(
              x = "Drug Rank Within Tumor",
              y = "Best Response (%)",
              title = title_mrecist_br
            ) +
            theme_pub_nogrid +
            ggplot2::theme(
              axis.ticks.x = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
              legend.position = "right",
              legend.box = "vertical",
              legend.key.height = grid::unit(0.55, "cm"),
              legend.spacing.y = grid::unit(0.15, "cm"),
              legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
            )
          ann_f456_orr <- build_adjacent_prop_annotations(
            summary_df = orr_ordered_f456,
            rank_col = "rank_pooled",
            success_col = "n_orr",
            total_col = "n",
            trend_df = trend_df_orr_f456,
            trend_rank_col = "rank_pooled",
            trend_value_col = "trend_fit_orr",
            layer_offset = 2.2,
            span_step = 2.8
          )
          ann_f456_dcr <- build_adjacent_prop_annotations(
            summary_df = dcr_ordered_f456,
            rank_col = "rank_pooled",
            success_col = "n_dcr",
            total_col = "n",
            trend_df = trend_df_dcr_f456,
            trend_rank_col = "rank_pooled",
            trend_value_col = "trend_fit_dcr",
            layer_offset = 6.8,
            span_step = 3.1
          )
          p461 <- add_pairwise_signif_brackets(p461, ann_f456_orr)
          p461 <- add_pairwise_signif_brackets(p461, ann_f456_dcr)
          p456 <- p461 +
            ggplot2::labs(
              title = "Best Response by DDA Score Tiers",
              subtitle = "(BestResponse data)"
            ) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(size = 14, margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
              plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 6), hjust = 0.5)
            )
          ## F462: mirror F99L (bars only, light green PR)
          p462 <- ggplot2::ggplot(
            mrecist_by_rank_f456,
            ggplot2::aes(
              x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
              y = pct,
              fill = resp
            )
          ) +
            ggplot2::geom_col(
              position = ggplot2::position_stack(),
              width = 0.6,
              colour = "black",
              linewidth = 0.25
            ) +
            ggplot2::scale_fill_manual(values = mrecist_colors_f99l, name = mrecist_fill_title_br) +
            ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::labs(x = "Drug Rank Within Tumor", y = "Best Response (%)", title = title_mrecist_br) +
            theme_pub_nogrid +
            ggplot2::theme(
              axis.ticks.x = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
              legend.position = "right",
              legend.key.height = grid::unit(0.55, "cm"),
              legend.spacing.y = grid::unit(0.15, "cm"),
              legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
            )
          ## F463: mirror F99X (bars only; CR blue, PR light blue)
          p463 <- ggplot2::ggplot(
            mrecist_by_rank_f456,
            ggplot2::aes(
              x = factor(rank_pooled, levels = rank_levels, labels = rank_labels),
              y = pct,
              fill = resp
            )
          ) +
            ggplot2::geom_col(
              position = ggplot2::position_stack(),
              width = 0.6,
              colour = "black",
              linewidth = 0.25
            ) +
            ggplot2::scale_fill_manual(values = mrecist_colors_f99x, name = mrecist_fill_title_br) +
            ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::labs(x = "Drug Rank Within Tumor", y = "Best Response (%)", title = title_mrecist_br) +
            theme_pub_nogrid +
            ggplot2::theme(
              axis.ticks.x = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
              axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
              legend.position = "right",
              legend.key.height = grid::unit(0.55, "cm"),
              legend.spacing.y = grid::unit(0.15, "cm"),
              legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
            )
        } else {
          p456 <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = "mRECIST by ranks (Response_Best_Response) not available (insufficient data)")
          p461 <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = "F461: mirror F99D RESPONSE_BEST_RESPONSE not available (insufficient data)")
          p462 <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = "F462: mirror F99L RESPONSE_BEST_RESPONSE not available (insufficient data)")
          p463 <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = "F463: mirror F99X RESPONSE_BEST_RESPONSE not available (insufficient data)")
        }
      } else {
        p456 <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "RESPONSE_BEST_RESPONSE not in df (add Response_Best_Response to Excel; re-run Step 2, refresh plot cache)")
        p461 <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "F461: RESPONSE_BEST_RESPONSE column not found")
        p462 <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "F462: RESPONSE_BEST_RESPONSE column not found")
        p463 <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "F463: RESPONSE_BEST_RESPONSE column not found")
      }

      ## Durable Benefit (DB) stacked bars with CA trend
      ## DB definition: TIME_TO_DOUBLE > 42
      db_colors <- c("DB" = "#81c784", "NonDB" = "#d9d9d9")

      db_df_rank <- df_all_ranks %>%
        dplyr::filter(is.finite(TIME_TO_DOUBLE), !is.na(rank_pooled)) %>%
        dplyr::mutate(is_db = TIME_TO_DOUBLE > 42)

      if (nrow(db_df_rank) > 0L) {
        db_by_rank <- db_df_rank %>%
          dplyr::group_by(rank_pooled) %>%
          dplyr::summarise(
            n_db = sum(is_db, na.rm = TRUE),
            n = dplyr::n(),
            .groups = "drop"
          ) %>%
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

        lm_db_rank <- stats::lm(db_pct ~ rank_pooled, data = db_by_rank)
        trend_df_db_rank <- data.frame(
          rank_pooled = rank_levels,
          trend_fit_db = stats::predict(lm_db_rank, newdata = data.frame(rank_pooled = rank_levels))
        )

        ca_ann_f467 <- ca_trend_annotate(
          scores = db_by_rank$rank_pooled, n = db_by_rank$n, success = db_by_rank$n_db,
          metric = "DB", panel_id = "Durable_Benefit_by_rank", analysis = "normal_rank", stratification = "rank_pooled"
        )
        trend_lab_db_rank <- ca_ann_f467$label

        p467 <- ggplot2::ggplot(
          db_comp_rank,
          ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = pct, fill = db_status)
        ) +
          ggplot2::geom_col(position = ggplot2::position_stack(), width = 0.6, colour = "black", linewidth = 0.25) +
          ggplot2::geom_line(
            data = trend_df_db_rank,
            ggplot2::aes(x = factor(rank_pooled, levels = rank_levels, labels = rank_labels), y = trend_fit_db, group = 1, color = "CA (DB%)"),
            linewidth = 1,
            linetype = "solid",
            inherit.aes = FALSE
          ) +
          ggplot2::annotate("text", x = Inf, y = Inf, label = trend_lab_db_rank, parse = TRUE, hjust = 1.05, vjust = 1.5, size = fsz(3), colour = "black") +
          ggplot2::scale_fill_manual(values = db_colors, name = "Durable Benefit\n(>42 days)") +
          ggplot2::scale_color_manual(values = c("CA (DB%)" = "black"), breaks = c("CA (DB%)"), name = NULL) +
          ggplot2::guides(
            fill = ggplot2::guide_legend(order = 1),
            color = ggplot2::guide_legend(order = 2)
          ) +
          ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, 0.1))) +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::labs(
            x = "Drug Rank Within Tumor",
            y = "Durable Benefit (%)",
            title = "Durable Benefit by Rank"
          ) +
          theme_pub_nogrid +
          ggplot2::theme(
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 0), face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
            legend.position = "right",
            legend.box = "vertical",
            legend.key.height = grid::unit(0.55, "cm"),
            legend.spacing.y = grid::unit(0.15, "cm"),
            legend.title = ggplot2::element_text(size = ggplot2::rel(0.7))
          )
        ann_f467_db <- build_adjacent_prop_annotations(
          summary_df = db_by_rank,
          rank_col = "rank_pooled",
          success_col = "n_db",
          total_col = "n",
          trend_df = trend_df_db_rank,
          trend_rank_col = "rank_pooled",
          trend_value_col = "trend_fit_db",
          layer_offset = 2.2,
          span_step = 2.8
        )
        p467 <- add_pairwise_signif_brackets(p467, ann_f467_db)
      } else {
        p467 <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = "Durable Benefit by ranks not available")
      }

    } else {
      p99c <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F99C: No evaluable CR/PR data by rank")
      p99d <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F99D: mRECIST by ranks not available (insufficient ORR data)")
      p99l <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F99L: mRECIST by ranks not available (insufficient ORR data)")
      p99x <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F99X: mRECIST by ranks not available (insufficient ORR data)")
      p456 <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "mRECIST by ranks not available (insufficient ORR data)")
      p456_avg <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "mRECIST by ranks (Response_BestAvgResponse) not available (insufficient ORR data)")
      p461 <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F461: mRECIST by rank not available (insufficient ORR data)")
      p462 <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F462: mRECIST by rank not available (insufficient ORR data)")
      p463 <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "F463: mRECIST by rank not available (insufficient ORR data)")
      p467 <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Durable Benefit by ranks not available (insufficient data)")
    }
  } else {
    p88 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F88: Response_BestAvgResponse not available")
    p99 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99: Response_BestAvgResponse not available")
    p99b <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99B: Response_BestAvgResponse not available")
    p34 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F34: Response_BestAvgResponse not available")
    p99c <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99C: Response_BestAvgResponse not available")
    p99d <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99D: Response_BestAvgResponse not available")
    p99l <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99L: Response_BestAvgResponse not available")
    p99x <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F99X: Response_BestAvgResponse not available")
    p456 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Response_BestAvgResponse not available")
    p456_avg <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Response_BestAvgResponse not available")
    p461 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F461: Response_BestAvgResponse not available")
    p462 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F462: Response_BestAvgResponse not available")
    p463 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "F463: Response_BestAvgResponse not available")
    p467 <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Response_BestAvgResponse not available")
  }

  ca_cmp_export <- ca_registry_export(analysis_dir)
  log_note(paste0(
    "Cochran-Armitage registry exported (n=", nrow(ca_cmp_export),
    "; Z mismatches=", sum(!ca_cmp_export$z_match, na.rm = TRUE),
    "; p mismatches=", sum(!ca_cmp_export$p_match, na.rm = TRUE), ")."
  ))

  save_fig(p6, file.path(fig_dir, "SCAM_residual_diagnostics.png"), FIG_F6_WIDTH, FIG_F6_HEIGHT)
  save_fig(p22_ox_purple, file.path(fig_dir, "SCAM_vs_GAM_all_ranks.png"), FIG_F22B_WIDTH, FIG_F22B_HEIGHT)
  save_fig(p8, file.path(fig_dir, "Effect_forest_enhanced.png"), FIG_F8_WIDTH, FIG_F8_HEIGHT)
  save_fig(p456, file.path(fig_dir, "mRECIST_by_rank_ORR_DCR_Response_Best_Response.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)
  save_fig(p456_avg, file.path(fig_dir, "mRECIST_by_rank_ORR_DCR_Response_BestAvgResponse.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)
  save_fig(p467, file.path(fig_dir, "Durable_Benefit_by_rank_CA_trend.png"), FIG_F99_WIDTH, FIG_F99_HEIGHT)

  log_note("Step 9 complete: publication figures saved.")
  invisible(NULL)
})

log_step("All done: normal-rank analysis (all data) completed.")



## SECTION 3: Leave-one-TARGET-out SCAM vs GAM (LOU)

## Leave-one-TARGET-out analysis: SCAM vs unconstrained GAM on normal ranks.
##
## For each iteration, one TARGET is withheld (or all data for the baseline "ALL" run).
## Per-patient compound ranks are derived from DDA LEVEL; ranks 7+ are pooled as "7-10".
## Monotonic SCAM and unconstrained GAM are fit on rank_scaled with a random intercept
## per patient. Each panel overlays both smooths on mean PFS by rank; EDF and smooth
## p-values are taken from each model's summary s.table (first smooth term).
##
## Outputs (peer review):
##   - Per iteration: SCAM vs GAM PNG, SCAM/GAM model summary text files
##   - Run level: SCAM_GAM_Final_table.xlsx

input_path <- file.path(DATA_DIR, DEFAULT_INPUT_XLSX)
output_root <- file.path(DATA_MODELING_OUTPUT_ROOT, "leave_one_out_target", "normal_ranks")
date <- format(Sys.Date(), "%Y_%m_%d")

set.seed(42)

LOU_ENABLED <- TRUE
LOU_INCLUDE_BASELINE <- TRUE
LOU_TARGETS_TEST <- NULL  # e.g. c("EGFR") for a smoke test

RESPONSE_COL_NAME <- "TIMETODOUBLE"

RANK_POOL_FROM <- 7L
RANK_POOL_LABEL <- "7-10"

FIG_WIDTH <- 13
FIG_HEIGHT <- 10
FIG_F22B_WIDTH <- FIG_WIDTH
FIG_F22B_HEIGHT <- FIG_HEIGHT

FONT_SCALE_ALL <- 1.40
fsz <- function(x) x * FONT_SCALE_ALL

F22_TITLE_PT <- 16L
F22_SUBTITLE_PT <- 12L
F22_AXIS_TITLE_PT <- 14L
F22_AXIS_TEXT_PT <- 11L
F22_LEGEND_PT <- 11L
F22_N_LABEL_PT <- 3.2
F22_CORNER_STAT_PT <- 3.4
F22_SCAM_COLOR <- "#002147"
F22_GAM_ORANGE_COLOR <- "#006400"
F22B_STANDALONE_FONT_MULT <- 1.30
F22B_AXIS_N_FONT_MULT <- 1.20

f22b_sz <- function(x) fsz(x * F22B_STANDALONE_FONT_MULT)
f22b_axis_n_sz <- function(x) f22b_sz(x * F22B_AXIS_N_FONT_MULT)

F22B_LOU_LEGEND_X <- 0.83
F22B_LOU_SCAM_GAM_TEXT_PLUS_PT <- 2
F22B_LOU_STATS_Y_GAM <- 17
F22B_LOU_STATS_Y_SCAM <- 11
F22B_LOU_STATS_X <- 1.5
F22B_LOU_STATS_HJUST <- 0

lou_date_dir <- file.path(output_root, date)
dir.create(lou_date_dir, showWarnings = FALSE, recursive = TRUE)
lou_run_summary_rows <- list()
left_out_target <- "ALL"
n_rows_excluded <- 0L

.step_counter <- 0L

TARGET_MAP_LOOKUP <- c(
  "MAPK" = "MAP2K1",
  "BRAF" = "BRAF",
  "CDK4,CDK6,CDK4/6" = "CDK4/6",
  "HSP90" = "HSP90",
  "PIK3CA" = "PIK3CA",
  "PIK3CA,PIK3CB,PIK3CG,PIK3CD,panPI3K" = "panPI3K",
  "EGFR" = "EGFR",
  "MET" = "MET",
  "FGFR" = "FGFR",
  "FGFR2/4" = "FGFR2/4",
  "ERBB3" = "ERBB3",
  "MDM2" = "MDM2",
  "JAK" = "JAK",
  "ERBB2" = "ERBB2",
  "NTRK1,NTRK2,NTRK3" = "NTRK",
  "SMO" = "SMO",
  "ESR1" = "ESR1",
  "IGF1R" = "IGF1R",
  "ALK" = "ALK"
)

safe_target_dir <- function(site) {
  if (is.na(site) || !nzchar(trimws(as.character(site)))) return("ALL")
  text <- trimws(as.character(site))
  if (identical(text, "ALL")) return("ALL")
  text <- gsub("[<>:\"/\\|?*]", "_", text)
  text <- gsub("[,\\s]+", "_", text)
  text <- gsub("_+", "_", text)
  text <- gsub("^_|_$", "", text)
  if (!nzchar(text)) return("unknown_target")
  text
}

apply_target_remap <- function(df) {
  if (!"TARGET" %in% names(df)) stop("apply_target_remap: TARGET column missing.")
  df$TARGET_RAW <- as.character(df$TARGET)
  mapped <- unname(TARGET_MAP_LOOKUP[df$TARGET_RAW])
  unmapped <- is.na(mapped)
  if (any(unmapped)) {
    unknown <- unique(df$TARGET_RAW[unmapped])
    message("TARGET remap: unmapped values kept as-is: ", paste(unknown, collapse = ", "))
    mapped[unmapped] <- df$TARGET_RAW[unmapped]
  }
  df$TARGET <- mapped
  df$map_target <- mapped
  df
}

lookup_map_target <- function(target_key) {
  key <- trimws(as.character(target_key))
  if (!nzchar(key) || tolower(key) == "all") return(key)
  mapped <- unname(TARGET_MAP_LOOKUP[key])
  if (!is.na(mapped)) return(mapped)
  key
}

f22b_leave_out_filename_tag <- function(target_key) {
  if (tolower(trimws(as.character(target_key))) == "all") return("leave_out_baseline_ALL")
  paste0("leave_out_", safe_target_dir(target_key))
}

save_f22b_plot <- function(plot_obj, fig_w, fig_h, dpi = 300L) {
  tag <- f22b_leave_out_filename_tag(left_out_target)
  base <- paste0("SCAM_vs_GAM_oxfordblue_orange_", tag)
  collect_dir <- file.path(lou_date_dir, "SCAM_vs_GAM_oxfordblue_orange_all_LOU")
  dir.create(collect_dir, showWarnings = FALSE, recursive = TRUE)
  out_png <- file.path(collect_dir, paste0(base, ".png"))
  ggplot2::ggsave(out_png, plot_obj, width = fig_w, height = fig_h, dpi = dpi)
  message("Saved: ", out_png)
  invisible(out_png)
}

record_lou_run_summary <- function(status, err_msg = NA_character_) {
  row <- data.frame(
    left_out_TARGET = left_out_target,
    n_rows_excluded = as.integer(n_rows_excluded),
    n_rows_remaining = if (exists("df") && is.data.frame(df)) nrow(df) else NA_integer_,
    n_patients = if (exists("df") && is.data.frame(df)) dplyr::n_distinct(df$REPORT_ID) else NA_integer_,
    scam_p = if (exists("scam_out") && !is.null(scam_out$p_value)) scam_out$p_value else NA_real_,
    scam_edf = if (exists("scam_out") && !is.null(scam_out$edf)) scam_out$edf else NA_real_,
    gam_p = if (exists("wiggle_out") && !is.null(wiggle_out$p_value)) wiggle_out$p_value else NA_real_,
    gam_edf = if (exists("wiggle_out") && !is.null(wiggle_out$edf)) wiggle_out$edf else NA_real_,
    status = status,
    error_message = err_msg,
    stringsAsFactors = FALSE
  )
  lou_run_summary_rows[[length(lou_run_summary_rows) + 1L]] <<- row
  invisible(row)
}

init_lou_output_dirs <- function(target_label) {
  sub_name <- if (identical(target_label, "ALL")) "ALL" else paste0("leave_out_", safe_target_dir(target_label))
  analysis_dir <<- file.path(lou_date_dir, sub_name)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  invisible(analysis_dir)
}

standardize_pdx_df <- function(df_local) {
  df_local <- df_local %>%
    dplyr::rename_with(~ stringr::str_replace_all(., "\\s+", "_")) %>%
    dplyr::rename_with(~ toupper(.))
  required_candidates <- list(
    REPORT_ID = c("REPORT_ID", "PATIENT", "PATIENT_ID", "SAMPLE_ID", "MODEL", "MODEL_ID"),
    COMPOUND = c("COMPOUND", "DRUG", "COMPOUND_NAME"),
    SCORE = c("LEVEL", "SCORE", "PREDICTIVE_SCORE", "PRED_SCORE"),
    TIME_TO_DOUBLE = c(RESPONSE_COL_NAME, "TIMETODOUBLE", "TIME_TO_DOUBLE", "TIME_TO_DOUBLING", "TTD")
  )
  resolve_col <- function(candidates, df_cols) {
    hit <- candidates[candidates %in% df_cols]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  col_map <- sapply(required_candidates, resolve_col, df_cols = colnames(df_local))
  if (any(is.na(col_map))) {
    stop(paste("Missing required columns:", paste(names(col_map)[is.na(col_map)], collapse = ", ")))
  }
  df_local %>%
    dplyr::transmute(
      REPORT_ID = as.character(.data[[col_map["REPORT_ID"]]]),
      COMPOUND = as.character(.data[[col_map["COMPOUND"]]]),
      SCORE = as.numeric(.data[[col_map["SCORE"]]]),
      TIME_TO_DOUBLE = as.numeric(.data[[col_map["TIME_TO_DOUBLE"]]])
    ) %>%
    dplyr::filter(!is.na(REPORT_ID), !is.na(COMPOUND), !is.na(SCORE), !is.na(TIME_TO_DOUBLE))
}

load_raw_pdx_with_target <- function() {
  df_raw <- readxl::read_xlsx(input_path)
  df_local <- df_raw %>%
    dplyr::rename_with(~ stringr::str_replace_all(., "\\s+", "_")) %>%
    dplyr::rename_with(~ toupper(.))
  target_candidates <- c("TARGET", "TREATMENT_TARGET", "TREATMENTTARGET")
  col_target <- target_candidates[target_candidates %in% colnames(df_local)][1]
  if (is.na(col_target)) stop("Missing TARGET column in input data.")
  df_local$TARGET <- as.character(df_local[[col_target]])
  df_local
}

write_final_outputs <- function() {
  write_scam_gam_final_table()
  invisible(NULL)
}

run_step <- function(step_name, fn) {
  .step_counter <<- .step_counter + 1L
  message("Step ", .step_counter, ": ", step_name)
  out <- fn()
  message("Step ", .step_counter, " complete: ", step_name)
  out
}

theme_f22b_standalone <- function() {
  ggplot2::theme_minimal(base_size = f22b_sz(F22_AXIS_TITLE_PT)) +
    ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.text = ggplot2::element_text(color = "black", size = f22b_sz(F22_LEGEND_PT)),
      plot.title = ggplot2::element_text(color = "black", face = "bold", size = f22b_sz(F22_TITLE_PT), hjust = 0.5, margin = ggplot2::margin(b = 10)),
      plot.subtitle = ggplot2::element_text(color = "black", size = f22b_sz(F22_SUBTITLE_PT), hjust = 0.5, margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(color = "black", size = f22b_axis_n_sz(F22_AXIS_TITLE_PT)),
      axis.text = ggplot2::element_text(color = "black", size = f22b_axis_n_sz(F22_AXIS_TEXT_PT)),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.ticks.length = grid::unit(0.2, "cm"),
      plot.margin = ggplot2::margin(t = 30, r = 8, b = 8, l = 8),
      legend.key.width = grid::unit(1.3, "cm")
    )
}

lou_target_table_label <- function(target_key) {
  key <- trimws(as.character(target_key))
  if (tolower(key) == "all") return("None")
  if (!nzchar(key)) return("Unknown")
  lookup_map_target(key)
}

lou_f22b_plot_title <- function(target_key) {
  key <- trimws(as.character(target_key))
  if (tolower(key) == "all") return("All")
  if (!nzchar(key)) return("Unknown TARGET")
  lookup_map_target(key)
}

fmt_smooth_p_plot_label <- function(p, cap_below = 1e-30, sci_digits = 2L) {
  if (is.na(p) || !is.finite(as.numeric(p))) return("p = NA")
  p <- as.numeric(p)
  if (p < cap_below) {
    return(paste0("p < ", trimws(format(cap_below, scientific = TRUE))))
  }
  paste0("p = ", formatC(p, format = "e", digits = sci_digits))
}

lou_f22b_stats_label_line <- function(model, edf, p) {
  paste0(model, " EDF = ", round(edf, 3), ", ", fmt_smooth_p_plot_label(p))
}

build_f22b_scam_gam_plot <- function(
    obs_rank, pred_scam, pred_gam, scam_edf, scam_p, gam_edf, gam_p,
    target_key, x_lab, gam_color = F22_GAM_ORANGE_COLOR) {
  n_label_sz <- f22b_axis_n_sz(F22_N_LABEL_PT)
  scam_lw <- 1.05 / 2
  y_lo <- 0
  y_hi <- 60
  y_span <- y_hi - y_lo
  y_lab_off <- 0.04 * y_span
  stats_sz <- f22b_sz(F22_CORNER_STAT_PT + F22B_LOU_SCAM_GAM_TEXT_PLUS_PT)
  ggplot2::ggplot(obs_rank, ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = mean_ttd)) +
    ggplot2::geom_point(size = 1.8, color = "black") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_ttd - 1.96 * se_ttd, ymax = mean_ttd + 1.96 * se_ttd),
      width = 0.2, alpha = 0.9, color = "black", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = factor(rank_int, levels = rank_levels, labels = rank_labels),
        y = mean_ttd + 1.96 * se_ttd + y_lab_off,
        label = paste0("(n = ", n, ")")
      ),
      data = obs_rank, angle = 90, hjust = -0.1, vjust = 0.5, size = n_label_sz, color = "black"
    ) +
    ggplot2::geom_ribbon(
      data = pred_scam,
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "SCAM"),
      inherit.aes = FALSE, alpha = 0.12
    ) +
    ggplot2::geom_line(
      data = pred_scam,
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "SCAM"),
      linewidth = scam_lw
    ) +
    ggplot2::geom_ribbon(
      data = pred_gam,
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), ymin = ci_low, ymax = ci_high, group = 1, fill = "GAM"),
      inherit.aes = FALSE, alpha = 0.08
    ) +
    ggplot2::geom_line(
      data = pred_gam,
      ggplot2::aes(x = factor(rank_int, levels = rank_levels, labels = rank_labels), y = pred_ttd, group = 1, color = "GAM"),
      linewidth = 1.05, linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(name = NULL, values = c("SCAM" = F22_SCAM_COLOR, "GAM" = gam_color), breaks = c("SCAM", "GAM")) +
    ggplot2::scale_fill_manual(name = NULL, values = c("SCAM" = F22_SCAM_COLOR, "GAM" = gam_color), guide = "none") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "dashed"), linewidth = c(1.2, 1.4)), reverse = FALSE)) +
    ggplot2::annotate(
      "text", x = F22B_LOU_STATS_X, y = F22B_LOU_STATS_Y_GAM,
      label = lou_f22b_stats_label_line("GAM", gam_edf, gam_p),
      hjust = F22B_LOU_STATS_HJUST, vjust = 0.5, size = stats_sz, color = "black"
    ) +
    ggplot2::annotate(
      "text", x = F22B_LOU_STATS_X, y = F22B_LOU_STATS_Y_SCAM,
      label = lou_f22b_stats_label_line("SCAM", scam_edf, scam_p),
      hjust = F22B_LOU_STATS_HJUST, vjust = 0.5, size = stats_sz, color = "black"
    ) +
    ggplot2::labs(
      title = lou_f22b_plot_title(target_key),
      x = x_lab,
      y = "Avg PFS (days)"
    ) +
    theme_f22b_standalone() +
    ggplot2::theme(
      legend.position = c(F22B_LOU_LEGEND_X, 0.98),
      legend.justification = c(1, 1)
    ) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi), clip = "off")
}

fmt_edf_table <- function(edf) {
  if (is.na(edf) || !is.finite(as.numeric(edf))) return("NA")
  formatC(round(as.numeric(edf), 3), format = "f", digits = 3)
}

fmt_pval_sci2 <- function(p) {
  if (is.na(p)) return("NA")
  formatC(p, format = "e", digits = 2)
}

fmt_pval_sci2_table <- function(p) {
  if (is.na(p) || !is.finite(as.numeric(p))) return("NA")
  gsub("e", "E", fmt_pval_sci2(as.numeric(p)), fixed = TRUE)
}

fmt_pval_gam_table <- function(p) {
  if (is.na(p) || !is.finite(as.numeric(p))) return("NA")
  p <- as.numeric(p)
  if (p < 1e-30) return("< 1e-30")
  gsub("e", "E", fmt_pval_sci2(p), fixed = TRUE)
}

build_scam_gam_final_table <- function(run_sum) {
  if (is.null(run_sum) || nrow(run_sum) == 0L) {
    return(data.frame(
      `TARGET excluded` = character(0),
      `Treatments excluded (n)` = integer(0),
      `Remaining (n)` = integer(0),
      `SCAM EDF` = character(0),
      `SCAM p value` = character(0),
      `GAM EDF` = character(0),
      `GAM p value` = character(0),
      Pipeline = character(0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  rs <- run_sum
  keys <- unique(as.character(rs$left_out_TARGET))
  target_order <- c("ALL", sort(setdiff(keys, "ALL"), method = "shell"))
  remaining_n <- as.integer(rs$n_rows_remaining)
  order_key <- match(tolower(rs$left_out_TARGET), tolower(target_order))
  in_order <- !is.na(order_key)
  extra <- rs[!in_order, , drop = FALSE]
  extra_rem <- remaining_n[!in_order]
  rs <- rs[in_order, , drop = FALSE]
  remaining_n <- remaining_n[in_order]
  ord <- order(order_key[in_order])
  rs <- rs[ord, , drop = FALSE]
  remaining_n <- remaining_n[ord]
  if (nrow(extra) > 0) {
    extra$._order_key <- nrow(rs) + seq_len(nrow(extra))
    rs <- dplyr::bind_rows(rs, extra[order(extra$._order_key), , drop = FALSE])
    remaining_n <- c(remaining_n, extra_rem[order(extra$._order_key)])
  }
  gam_p_col <- if ("gam_p" %in% names(rs)) rs$gam_p else NA_real_
  gam_edf_col <- if ("gam_edf" %in% names(rs)) rs$gam_edf else NA_real_
  data.frame(
    `TARGET excluded` = vapply(rs$left_out_TARGET, lou_target_table_label, character(1)),
    `Treatments excluded (n)` = as.integer(rs$n_rows_excluded),
    `Remaining (n)` = remaining_n,
    `SCAM EDF` = vapply(rs$scam_edf, fmt_edf_table, character(1)),
    `SCAM p value` = vapply(rs$scam_p, fmt_pval_sci2_table, character(1)),
    `GAM EDF` = vapply(gam_edf_col, fmt_edf_table, character(1)),
    `GAM p value` = vapply(gam_p_col, fmt_pval_gam_table, character(1)),
    Pipeline = "Ranking",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_scam_gam_final_table <- function() {
  if (length(lou_run_summary_rows) == 0L) {
    message("SCAM_GAM_Final_table skipped: no completed LOU iterations.")
    return(invisible(NULL))
  }
  run_sum <- dplyr::bind_rows(lou_run_summary_rows)
  display <- build_scam_gam_final_table(run_sum)
  out_path <- file.path(lou_date_dir, "SCAM_GAM_Final_table.xlsx")
  writexl::write_xlsx(list(SCAM_GAM_Final_table = display), path = out_path)
  message("Saved: ", out_path)
  invisible(out_path)
}

extract_smooth_stats <- function(sm) {
  st_raw <- sm$s.table
  if (is.null(st_raw) || (is.matrix(st_raw) && (nrow(st_raw) == 0L || ncol(st_raw) == 0L))) {
    return(list(edf = NA_real_, p_value = NA_real_, s_table = NULL))
  }
  s_tab <- as.data.frame(st_raw)
  s_tab$term <- rownames(st_raw)
  rownames(s_tab) <- NULL
  nms <- names(s_tab)
  if (!is.null(nms) && length(nms) == ncol(s_tab)) names(s_tab) <- make.names(nms)
  edf <- if ("edf" %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab$edf[1] else NA_real_
  pcol <- grep("^p", names(s_tab), ignore.case = TRUE, value = TRUE)[1]
  p_val <- if (!is.na(pcol) && pcol %in% names(s_tab) && nrow(s_tab) >= 1L) s_tab[[pcol]][1] else NA_real_
  list(edf = edf, p_value = p_val, s_table = s_tab)
}

predict_smooth_by_rank <- function(model, model_df, rank_levels, max_rank) {
  rank_scaled_mid <- if (max_rank > 1L) 1 - (rank_levels - 1) / (max_rank - 1) else rep(1, length(rank_levels))
  pred_grid <- tibble::tibble(
    rank_scaled = rank_scaled_mid,
    REPORT_ID = model_df$REPORT_ID[1]
  )
  pred <- predict(model, newdata = pred_grid, type = "response", se.fit = TRUE, exclude = "s(REPORT_ID)")
  pred_grid %>%
    mutate(
      rank_int = rank_levels,
      pred_ttd = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      ci_low = pred_ttd - 1.96 * se,
      ci_high = pred_ttd + 1.96 * se
    )
}

observed_mean_by_rank <- function(model_df) {
  model_df %>%
    group_by(.data$rank_pooled) %>%
    summarise(
      mean_ttd = mean(TIME_TO_DOUBLE, na.rm = TRUE),
      se_ttd = sd(TIME_TO_DOUBLE, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    rename(rank_int = rank_pooled)
}

n_workers <- parallel::detectCores(logical = TRUE)
if (is.na(n_workers) || n_workers < 1) n_workers <- 1L

run_step("Load packages", function() {
  required_pkgs <- c(
    "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "tibble",
    "scam", "mgcv", "writexl", "parallel"
  )
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")) {
    tryCatch(do.call(Sys.setenv, setNames(list(as.character(n_workers)), v)), error = function(e) NULL)
  }
})

df_raw_full <- run_step("Load raw PDX data with TARGET", function() {
  out <- load_raw_pdx_with_target()
  apply_target_remap(out)
})

if (!isTRUE(LOU_ENABLED)) stop("LOU_ENABLED must be TRUE for leave-one-TARGET-out SCAM vs GAM analysis.")

targets_all <- sort(unique(na.omit(df_raw_full$map_target)))
if (!is.null(LOU_TARGETS_TEST)) targets_all <- intersect(targets_all, LOU_TARGETS_TEST)
lou_iterations <- c(if (isTRUE(LOU_INCLUDE_BASELINE)) "ALL", targets_all)
message("LOU: ", length(lou_iterations), " iterations (baseline + ", length(targets_all), " TARGETs).")

for (left_out_target in lou_iterations) {
  message("LOU iteration (left out: ", left_out_target, ")")
  .step_counter <- 0L
  init_lou_output_dirs(left_out_target)
  lou_iter_err <- NULL
  lou_iter_status <- "OK"

  tryCatch({
    df <- run_step("Filter and standardize data", function() {
      df_filtered <- if (identical(left_out_target, "ALL")) {
        df_raw_full
      } else {
        df_raw_full %>% dplyr::filter(TARGET != left_out_target)
      }
      n_rows_excluded <<- if (identical(left_out_target, "ALL")) 0L else (nrow(df_raw_full) - nrow(df_filtered))
      standardize_pdx_df(df_filtered)
    })

    df_all_ranks <- run_step("Compute per-patient normal ranks", function() {
      all_ranks <- df %>%
        group_by(REPORT_ID) %>%
        filter(n() >= 2L) %>%
        mutate(
          rank_bin = dplyr::dense_rank(-SCORE),
          rank_pooled = pmin(rank_bin, RANK_POOL_FROM),
          rank_scaled = 1 - (rank_pooled - 1) / (RANK_POOL_FROM - 1)
        ) %>%
        ungroup()
      if (nrow(all_ranks) == 0) stop("No data after computing normal ranks.")
      all_ranks
    })

    rank_levels <- 1L:RANK_POOL_FROM
    rank_labels <- c(as.character(1L:(RANK_POOL_FROM - 1L)), RANK_POOL_LABEL)
    max_rank <- RANK_POOL_FROM

    scam_out <- run_step("Fit monotonic SCAM", function() {
      model_df <- df_all_ranks %>%
        mutate(REPORT_ID = factor(REPORT_ID)) %>%
        as.data.frame()
      model <- scam::scam(
        TIME_TO_DOUBLE ~ s(rank_scaled, bs = "mpi", k = 7L) + s(REPORT_ID, bs = "re"),
        data = model_df
      )
      sm <- summary(model)
      capture.output(sm, file = file.path(analysis_dir, "model_scam_all_ranks_summary.txt"))
      stats <- extract_smooth_stats(sm)
      list(
        model = model,
        edf = stats$edf,
        p_value = stats$p_value,
        pred_grid_ranks = predict_smooth_by_rank(model, model_df, rank_levels, max_rank),
        obs_rank = observed_mean_by_rank(model_df)
      )
    })

    wiggle_out <- run_step("Fit unconstrained GAM and build SCAM vs GAM panel", function() {
      model_df <- df_all_ranks %>%
        mutate(REPORT_ID = factor(REPORT_ID)) %>%
        as.data.frame()
      gam_k <- min(20L, max(4L, max_rank))
      gam_model <- mgcv::gam(
        TIME_TO_DOUBLE ~ s(rank_scaled, bs = "tp", k = gam_k) + s(REPORT_ID, bs = "re"),
        data = model_df,
        method = "REML"
      )
      gam_sm <- summary(gam_model)
      capture.output(gam_sm, file = file.path(analysis_dir, "model_gam_unconstrained_all_ranks_summary.txt"))
      gam_stats <- extract_smooth_stats(gam_sm)
      gam_edf <- gam_stats$edf
      gam_p <- gam_stats$p_value
      pred_grid_gam <- predict_smooth_by_rank(gam_model, model_df, rank_levels, max_rank)
      obs_rank <- scam_out$obs_rank

      save_f22b_plot(
        build_f22b_scam_gam_plot(
          obs_rank, scam_out$pred_grid_ranks, pred_grid_gam,
          scam_out$edf, scam_out$p_value, gam_edf, gam_p,
          left_out_target, "DDA rank groups", F22_GAM_ORANGE_COLOR
        ),
        FIG_F22B_WIDTH, FIG_F22B_HEIGHT
      )

      list(model = gam_model, edf = gam_edf, p_value = gam_p, pred_grid = pred_grid_gam)
    })

    invisible(NULL)
  }, error = function(e) {
    lou_iter_err <<- conditionMessage(e)
    lou_iter_status <<- "FAILED"
    message("LOU iteration FAILED for ", left_out_target, ": ", lou_iter_err)
  })

  if (identical(lou_iter_status, "OK")) {
    record_lou_run_summary("OK")
  } else {
    record_lou_run_summary(lou_iter_status, lou_iter_err)
  }
}

write_final_outputs()
message("Leave-one-TARGET-out SCAM vs GAM analysis complete. Outputs in: ", lou_date_dir)
