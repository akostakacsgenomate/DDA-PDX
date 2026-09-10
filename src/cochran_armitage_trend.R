# Cochran-Armitage trend test helpers (DCR / ORR / Durable Benefit panels).
# Uses stats::prop.trend.test (Cochran-Armitage test for trend in proportions)
# and a manual normal approximation (Z and two-sided P).

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

#' Manual Cochran-Armitage Z (normal approximation).
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
