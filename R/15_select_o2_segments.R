## ============================================================
## FILE : R/select_o2_segments.R
## PURPOSE : Select segments representing O2 consumption using
##           kernel mode estimation on slope magnitudes
## AUTHOR : D. Nerini — June 2026
## ============================================================

## ---- select_o2_segments ----
#' From a slope classification table, select the subset of
#' decreasing segments (class == "dec") that best represent
#' active oxygen consumption.
#'
#' Steps:
#'   1. Restrict to decreasing segments longer than min_duration.
#'   2. Find the mode of the negative-slope distribution (kernel
#'      density estimate with bandwidth h, weighted by duration).
#'   3. Select segments within 1.5 * h of the mode.
#'   4. Iteratively expand the selection to include any segment
#'      that is not significantly different from an already
#'      selected segment (similarity threshold alpha_sim).
#'   5. Compute mean O2 consumption as the duration-weighted
#'      average slope across selected segments.
#'
#' @param slope_table   data.frame with columns segment, slope, duration, class
#' @param p_sim         Pairwise similarity p-value matrix from slope_inference
#' @param alpha_sim     Similarity threshold for segment expansion
#' @param min_duration  Minimum duration (time units) for a segment to be eligible
#' @return List with:
#'   - selected    : subset of slope_table rows that are retained
#'   - slope_mode  : modal slope (negative value)
#'   - conso_O2    : weighted mean O2 consumption rate (positive value, -slope)
select_o2_segments <- function(slope_table, p_sim, alpha_sim, min_duration) {
  # Candidate: decreasing segments that are long enough
  cand <- slope_table[slope_table$class == "dec" & slope_table$duration > min_duration, ]
  if (nrow(cand) == 0)
    return(list(selected = data.frame(), slope_mode = NA, conso_O2 = NA))

  # Kernel density mode of -slope (absolute consumption rate)
  rho <- -cand$slope
  w   <- cand$duration
  h   <- 1.06 * weighted.mean(abs(rho - weighted.mean(rho, w)), w) * length(rho)^(-0.2)
  if (!is.finite(h) || h <= 1e-8) h <- sd(rho)
  if (!is.finite(h) || h <= 1e-8) h <- 0.1 * max(abs(rho))
  if (!is.finite(h) || h <= 1e-8) h <- 1e-6

  grid <- seq(min(rho), max(rho), length.out = 500)
  dens <- sapply(grid, function(r) sum(w * dnorm((r - rho) / h) / h))
  rho_mode    <- grid[which.max(dens)]
  slope_mode  <- -rho_mode

  # Select segments close to the mode
  cand$dist_scaled <- abs(rho - rho_mode) / h
  sel_ids <- cand$segment[cand$dist_scaled <= 1.5]

  # Iteratively expand to include similar segments
  changed <- TRUE
  while (changed) {
    old <- sel_ids
    add <- cand$segment[sapply(cand$segment,
                               function(j) any(p_sim[j, sel_ids] > alpha_sim))]
    sel_ids <- sort(unique(c(sel_ids, add)))
    changed <- length(sel_ids) > length(old)
  }

  selected <- cand[cand$segment %in% sel_ids, ]
  conso_O2 <- if (nrow(selected) > 0) -weighted.mean(selected$slope, selected$duration) else NA

  list(selected  = selected,
       slope_mode = slope_mode,
       conso_O2   = conso_O2)
}