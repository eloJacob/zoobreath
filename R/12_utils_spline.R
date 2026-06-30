## ============================================================
## FILE : R/utils_spline.R
## PURPOSE : Utility functions for segment handling, spline
##           fitting, and prediction with penalised LAD B-splines
## AUTHOR : D. Nerini — June 2026
## ============================================================

library(splines2)
library(CVXR)

## ------------------------------------------------------------
## SEGMENT UTILITIES
## ------------------------------------------------------------

## ---- make_segment ----
#' Assign observations to segments based on break points
#' @param x Numeric vector of positions
#' @param breaks Numeric vector of break points (must span range of x)
#' @return Integer vector of segment labels (1-indexed)
make_segment <- function(x, breaks) {
  cut(x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
}

## ---- segment_counts ----
#' Count observations per segment
#' @param x Numeric vector of positions
#' @param breaks Numeric vector of break points
#' @return Integer vector of counts per segment
segment_counts <- function(x, breaks) {
  seg <- make_segment(x, breaks)
  as.integer(table(factor(seg, levels = seq_len(length(breaks) - 1))))
}

## ---- repair_short_segments ----
#' Iteratively merge the shortest adjacent pair of segments
#' until every segment has at least min_n observations.
#' Merges to the side (left or right) that results in the
#' larger combined neighbour segment.
#' @param x Numeric vector of positions
#' @param breaks Numeric vector of break points
#' @param min_n Minimum number of points required per segment
#' @return Updated vector of break points
repair_short_segments <- function(x, breaks, min_n = 10) {
  breaks <- sort(unique(breaks))
  repeat {
    counts <- segment_counts(x, breaks)
    if (all(counts >= min_n) || length(counts) <= 1) break
    j <- which.min(counts); nseg <- length(counts)
    if (j == 1) {
      remove_id <- 2
    } else if (j == nseg) {
      remove_id <- length(breaks) - 1
    } else {
      left_size  <- counts[j - 1] + counts[j]
      right_size <- counts[j]     + counts[j + 1]
      remove_id  <- if (left_size >= right_size) j else j + 1
    }
    breaks <- breaks[-remove_id]
  }
  breaks
}

## ------------------------------------------------------------
## SPLINE FITTING & PREDICTION
## ------------------------------------------------------------

## ---- fit_quant_pspline ----
#' Fit a penalised LAD B-spline to (x, y) data.
#' Uses CVXR to minimise  sum |y_i - f(x_i)| + lambda * ||D2 f||^2
#' where D2 is the second-difference penalty matrix.
#'
#' @param x        Numeric vector of predictor values
#' @param y        Numeric vector of response values (same length as x)
#' @param df       Degrees of freedom for the initial B-spline basis
#' @param lambda   Roughness penalty weight
#' @param knots    Optional numeric vector of internal knot locations
#' @return List with basis matrix B, coefficient vector beta, knot
#'         locations, scaled coordinates, fitted values, penalty matrix,
#'         lambda, and solver status
fit_quant_pspline <- function(x, y, df, lambda, knots = NULL) {
  ok <- is.finite(x) & is.finite(y); x0 <- x[ok]; y0 <- y[ok]
  ord <- order(x0); x0 <- x0[ord]; y0 <- y0[ord]
  xr <- range(x0); scale_x <- diff(xr)
  if (scale_x <= 0) stop("x is constant or invalid.")
  xs <- (x0 - xr[1]) / scale_x
  if (!is.null(knots)) {
    knots_s <- (knots - xr[1]) / scale_x
    knots_s <- knots_s[is.finite(knots_s) & knots_s > 0 & knots_s < 1]
    knots_s <- sort(unique(knots_s))
    if (length(knots_s) == 0) knots_s <- NULL
  } else {
    knots_s <- NULL
  }
  B <- if (is.null(knots_s)) {
    bSpline(xs, df = df, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  } else {
    bSpline(xs, knots = knots_s, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  }
  B <- as.matrix(B); P <- ncol(B)
  if (any(!is.finite(B))) stop("B-spline basis contains non-finite values.")
  D2 <- diff(diag(P), differences = 2)
  beta <- Variable(P)
  res  <- y0 - B %*% beta
  sol  <- solve(Problem(Minimize(sum(abs(res)) + lambda * sum_squares(D2 %*% beta))),
                solver = "CLARABEL", verbose = FALSE)
  if (!sol$status %in% c("optimal", "optimal_inaccurate"))
    stop(paste("CVXR solver failed:", sol$status))
  beta_hat <- as.vector(sol$getValue(beta))
  if (any(!is.finite(beta_hat))) stop("beta_hat contains non-finite values.")
  knots_out <- attr(B, "knots")
  if (length(knots_out) > 0) knots_out <- xr[1] + knots_out * scale_x
  list(B = B, beta = beta_hat, knots = knots_out, knots_s = attr(B, "knots"),
       x_range = xr, x_scaled = xs, x_used = x0, y_used = y0,
       fitted = drop(B %*% beta_hat), D2 = D2, lambda = lambda, status = sol$status)
}

## ---- predict_fit ----
#' Predict fitted values at new x locations using a fitted pspline object
#' @param fit   Fitted object returned by fit_quant_pspline
#' @param xnew  Numeric vector of new predictor values
#' @return Numeric vector of predicted values
predict_fit <- function(fit, xnew) {
  xs <- pmin(pmax((xnew - fit$x_range[1]) / diff(fit$x_range), 0), 1)
  Bnew <- bSpline(xs, knots = fit$knots_s, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  drop(as.matrix(Bnew) %*% fit$beta)
}

## ---- basis_at ----
#' Evaluate the B-spline basis matrix at new x locations
#' @param fit   Fitted object returned by fit_quant_pspline
#' @param xnew  Numeric vector of new predictor values
#' @return Matrix with one row per xnew and one column per basis function
basis_at <- function(fit, xnew) {
  xs <- pmin(pmax((xnew - fit$x_range[1]) / diff(fit$x_range), 0), 1)
  as.matrix(bSpline(xs, knots = fit$knots_s, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1)))
}