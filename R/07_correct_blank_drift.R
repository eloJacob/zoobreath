#' Correct one O2 channel for instrumental blank drift
#' Subtracts the linear trend estimated from blank channel
#'
#' @param df          Dataframe
#' @param channel     Name of the O2 column to correct (e.g., "umol.l.Ch1")
#' @param blank_slope Drift slope (umol/L per hour) from blank linear regression
#' @return            Corrected O2 vector
correct_blank_drift <- function(df, channel, blank_slope) {
  df[[channel]] - blank_slope * df$Time.h
}