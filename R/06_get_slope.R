#' Estimate linear drift slope with significance check
#'
#' @param df       Dataframe
#' @param channel  Name of the channel (e.g., "umol.l.Ch3")
#' @param alpha    Significance threshold (default 0.05)
#' @return         Slope (umol/L per hour), with a warning if not significant

get_slope <- function(df, channel, alpha = 0.05) {
  
  fit     <- lm(df[[channel]] ~ df$Time.h)
  slope   <- fit$coefficients[2]
  p_value <- summary(fit)$coefficients[2, 4]
  r2      <- summary(fit)$r.squared
  
  if (p_value > alpha) {
    message(sprintf(
      "[get_slope] WARNING: slope for '%s' is NOT significant (p = %.3f, R² = %.3f, slope = %.4f umol/L/h).\n  -> Consider not applying this correction.",
      channel, p_value, r2, slope
    ))
  } else {
    message(sprintf(
      "[get_slope] Slope for '%s': %.4f umol/L/h (p = %.3f, R² = %.3f).",
      channel, slope, p_value, r2
    ))
  }
  
  return(slope)
}