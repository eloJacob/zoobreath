## ============================================================
## FILE : R/slope_inference.R
## PURPOSE : Slope inference for penalised LAD spline
##           with bootstrap
## AUTHOR : D. Nerini — June 2026
## ============================================================

compute_contrast_matrix <- function(fit){
  br <- sort(unique(c(min(fit$x_used), fit$knots, max(fit$x_used))))
  nseg <- length(br) - 1
  beta <- fit$beta
  P <- length(beta)
  
  # Basis derivative matrix L: each row j gives dB/dx evaluated
  # at the centre of segment j (linear B-spline, so constant on
  # each interval; we approximate using the two endpoint values)
  L <- matrix(NA, nseg, P)
  for (j in seq_len(nseg)) {
    Bz <- basis_at(fit, c(br[j], br[j + 1]))
    L[j, ] <- (Bz[2, ] - Bz[1, ]) / (br[j + 1] - br[j])
  }
  return(list(br = br, L = L))
}

create_bootstrap_replicate <- function(fit){
  new_data <- data.frame(x = fit$x_used)
  residuals <- fit$y_used - fit$fitted
  residuals <- residuals - median(residuals)
  new_data$y <- fit$fitted + sample(residuals, length(fit$fitted), replace = TRUE)
  return(new_data)
}

compute_V_bootstrap <- function(fit, knots_final, n_boot = 500){
  coefs_boot <- matrix(nrow = n_boot, ncol = length(fit$beta))
  
  for(i in 1:n_boot){
    data_boot <- create_bootstrap_replicate(fit)
    boot_fit <- fit_quant_pspline(data_boot$x, data_boot$y, ncol(fit$B), fit$lambda, knots = knots_final)
    coefs_boot[i, ] <- boot_fit$beta
  }
  
  return(coefs_boot)
}

## ---- slope_inference ----
#' Compute asymptotic slope estimates and inference for each
#' knot-delimited segment of a fitted penalised LAD B-spline.
#'
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
slope_inference <- function(fit, knots_final, alpha_slope = 0.05, alpha_sim = 0.1) {
  contrast_mat <- compute_contrast_matrix(fit)
  br <- contrast_mat$br
  nseg <- length(br) - 1
  L <- contrast_mat$L
  
  slopes <- drop(L %*% fit$beta)
  coefs_boot <- compute_V_bootstrap(fit, knots_final)
  
  Vbeta <- cov(coefs_boot)
  
  Vslopes<- L %*% Vbeta %*% t(L)
  
  se <- sqrt(diag(Vslopes))
  zval <- slopes / se
  p_two <- 2 * pnorm(-abs(zval))

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

  return(list(breaks    = br,
       slope_table = slope_table,
       L          = L,
       Vbeta      = Vbeta,
       Vslopes    = Vslopes,
       p_sim      = p_sim,
       z_sim      = z_sim,
       similar    = p_sim > alpha_sim))
}