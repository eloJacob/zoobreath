## ============================================================
## FILE : R/merge_segments.R
## PURPOSE : Iterative segment merging driven by slope similarity
## AUTHOR : D. Nerini — June 2026
## ============================================================

## ---- merge_segments ----
#' Repeatedly merge adjacent segments whose slopes are not
#' significantly different (similarity threshold alpha_sim).
#' Each iteration: refit spline on current knots, run slope
#' inference, check pairwise similarity, and merge the pair with
#' the highest p-value among all similar adjacent pairs.
#' Short segments (below min_n) are repaired before each attempt.
#' Iterations stop when no similar adjacent pair remains.
#'
#' @param x           Numeric predictor vector
#' @param y           Numeric response vector
#' @param br0         Initial break points (vector, must span range of x)
#' @param df          B-spline degrees of freedom
#' @param lambda      Roughness penalty weight
#' @param alpha_slope Significance level for slope classification
#' @param alpha_sim   Similarity threshold; adjacent segments with
#'                    p-value > alpha_sim are candidates for merging
#' @param min_n       Minimum observations per segment
#' @return Final vector of break points after merging stops
merge_segments <- function(x, y, br0, df, lambda, alpha_slope, alpha_sim, min_n) {
  br_current <- br0; changed <- TRUE
  while (changed) {
    knots_c <- br_current[-c(1, length(br_current))]
    fit_c   <- tryCatch(fit_quant_pspline(x, y, df, lambda, knots_c),
                        error = function(e) NULL)
    if (is.null(fit_c)) break

    test_c <- tryCatch(slope_inference(fit_c, alpha_slope, alpha_sim),
                       error = function(e) NULL)
    if (is.null(test_c)) break

    nseg_c <- nrow(test_c$slope_table)
    if (nseg_c <= 1) break

    # Identify adjacent pairs that are not significantly different
    merge_adj <- sapply(seq_len(nseg_c - 1),
                        function(j) test_c$similar[j, j + 1])
    if (!any(merge_adj)) { changed <- FALSE; break }

    # Pick the pair with the highest p-value to merge next
    p_adj <- sapply(seq_len(nseg_c - 1),
                    function(j) test_c$p_sim[j, j + 1])
    j_remove <- which.max(ifelse(merge_adj, p_adj, -Inf))

    # Apply min_n repair to decide which break to remove
    candidate <- repair_short_segments(x, br_current[-(j_remove + 1)], min_n)
    if (length(candidate) < length(br_current)) {
      br_current <- candidate
    } else {
      changed <- FALSE
    }
  }
  repair_short_segments(x, br_current, min_n)
}