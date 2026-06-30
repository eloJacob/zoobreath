## ============================================================
## FILE : R/slope_inference.R
## PURPOSE : Asymptotic slope inference for penalised LAD spline
##           using the sandwich variance formula
## AUTHOR : D. Nerini — June 2026
## ============================================================

## ---- slope_inference ----
#' Compute asymptotic slope estimates and inference for each
#' knot-delimited segment of a fitted penalised LAD B-spline.
#'
#' Uses the sandwich variance formula:
#'   V(beta) = A^{-1} * B * A^{-1}
#' where A = 2*f(0)*B'B + 2*lambda*D2'D2
#' and B = (1/4)*B'B (LAD dispersion matrix).
#' Segment slopes and their variance are derived via the basis
#' derivative matrix L.
#'
#' Pairwise similarity tests between adjacent slopes are also
#' computed using the implied covariance matrix.
#'
#' @param fit        Fitted object returned by fit_quant_pspline
#' @param alpha_slope Significance level for classifying individual
#'                    segment slopes (NS / inc / dec)
#' @param alpha_sim  Significance threshold for declaring two slopes
#'                    as similar (used to decide whether to merge)
#' @return List with:
#'   - breaks       : break points used
#'   - slope_table  : data.frame with slope, se, z, p_value, class
#'   - L            : basis-derivative matrix (nseg x P)
#'   - Vbeta        : variance matrix of coefficients
#'   - Vslopes      : variance matrix of segment slopes
#'   - p_sim        : pairwise p-values for slope similarity
#'   - z_sim        : pairwise z-statistics for slope similarity
#'   - similar      : logical matrix; similar[i,j] = TRUE if
#'                    slopes i and j are not significantly different
slope_inference <- function(fit, alpha_slope = 0.05, alpha_sim = 0.1) {
  br <- sort(unique(c(min(fit$x_used), fit$knots, max(fit$x_used))))
  nseg <- length(br) - 1
  B <- fit$B; beta <- fit$beta; D2 <- fit$D2; lambda <- fit$lambda; P <- length(beta)

  # Basis derivative matrix L: each row j gives dB/dx evaluated
  # at the centre of segment j (linear B-spline, so constant on
  # each interval; we approximate using the two endpoint values)
  L <- matrix(NA, nseg, P)
  for (j in seq_len(nseg)) {
    Bz <- basis_at(fit, c(br[j], br[j + 1]))
    L[j, ] <- (Bz[2, ] - Bz[1, ]) / (br[j + 1] - br[j])
  }
  slopes <- drop(L %*% beta)

  # Residuals for density estimation at zero (LAD density f(0))
  resid <- fit$y_used - fit$fitted
  f0 <- approx(density(resid)$x, density(resid)$y, xout = 0, rule = 2)$y
  if (!is.finite(f0) || f0 <= 0) f0 <- 1 / (2 * mad(resid))

  # Sandwich variance: A = 2*f(0)*B'B + 2*lambda*D2'D2
  A <- 2 * f0 * crossprod(B) + 2 * lambda * crossprod(D2)
  Vbeta <- solve(A, 0.25 * crossprod(B)) %*% solve(A)
  Vslopes <- L %*% Vbeta %*% t(L)
  se <- sqrt(diag(Vslopes)); zval <- slopes / se; p_two <- 2 * pnorm(-abs(zval))

  slope_table <- data.frame(
    segment   = seq_len(nseg),
    start     = br[-length(br)],
    end       = br[-1],
    duration  = diff(br),
    n_points  = segment_counts(fit$x_used, br),
    slope     = slopes,
    se        = se,
    z         = zval,
    p_value   = p_two,
    class     = ifelse(p_two > alpha_slope, "NS",
                       ifelse(slopes > 0, "inc", "dec"))
  )

  # Pairwise similarity between segment slopes
  p_sim <- z_sim <- matrix(NA, nseg, nseg)
  for (i in seq_len(nseg)) for (j in seq_len(nseg)) {
    se_ij <- sqrt(Vslopes[i, i] + Vslopes[j, j] - 2 * Vslopes[i, j])
    z_sim[i, j] <- (slopes[i] - slopes[j]) / se_ij
    p_sim[i, j] <- 2 * pnorm(-abs(z_sim[i, j]))
  }
  diag(p_sim) <- 1

  list(breaks    = br,
       slope_table = slope_table,
       L          = L,
       Vbeta      = Vbeta,
       Vslopes    = Vslopes,
       p_sim      = p_sim,
       z_sim      = z_sim,
       similar    = p_sim > alpha_sim)
}