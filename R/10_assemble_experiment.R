#' Assemble final long-format dataframe for one experiment
#' Wrapper around df_to_long() for two devices
#'
#' @param d1  Processed & corrected dataframe from device 1
#' @param d2  Processed & corrected dataframe from device 2
#' @return    Long-format dataframe with columns: Time, O2, Condition, ID
assemble_experiment <- function(d1, d2) {
  df      <- rbind(df_to_long(d1), df_to_long(d2))
  df$Time <- as.numeric(df$Time)
  df$O2   <- as.numeric(df$O2)
  return(df)
}
