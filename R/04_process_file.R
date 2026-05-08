#' Load, clean and calibrate one pyroscience CSV file
#'
#' @param filepath    Path to CSV file
#' @param cal         Calibration coefficient list
#' @param cols        Named list of dphi column names:
#'                    list(Ch1 = "col_name", Ch2 = ..., Ch3 = ...)
#'                    If NULL, uses default column positions [7], [8], [9]
#' @param conditions  Named list describing the biological condition of each channel:
#'                    list(Ch1 = "ZooHP", Ch2 = "ZooHP", Ch3 = "BlankATM")
#'                    Accepted values: "ZooHP", "ZooATM", "BlankATM", "BlankHP"
#' @param ids         Named list of individual IDs for each channel:
#'                    list(Ch1 = "ZooHP.1", Ch2 = "ZooHP.2", Ch3 = "Blank")
#' @return            Cleaned dataframe with:
#'                      - umol.l.Ch1/Ch2/Ch3 columns (calibrated O2)
#'                      - condition.Ch1/Ch2/Ch3 columns (biological condition label)
#'                      - id.Ch1/Ch2/Ch3 columns (individual ID)
process_file <- function(filepath, cal,
                         cols       = NULL,
                         conditions = list(Ch1 = "ZooHP",   Ch2 = "ZooHP",
                                           Ch3 = "BlankATM"),
                         ids        = list(Ch1 = "Zoo.HP1", Ch2 = "Zoo.HP2",
                                           Ch3 = "Blank")) {
  
  raw <- read.csv(filepath, sep = ";")
  df  <- Clean.pyroscience.data(raw, unit = "%airsat", salinity = 38.5, dphi = TRUE)
  
  # --- Calibrate each channel ---
  if (is.null(cols)) {
    # Device 1 type: dphi in fixed column positions [7], [8], [9]
    df$umol.l.Ch1 <- calibrate_channel(df, df[, 7], cal)
    df$umol.l.Ch2 <- calibrate_channel(df, df[, 8], cal)
    df$umol.l.Ch3 <- calibrate_channel(df, df[, 9], cal)
  } else {
    # Device 2 type: dphi in named columns
    df$umol.l.Ch1 <- calibrate_channel(df, df[[cols$Ch1]], cal)
    df$umol.l.Ch2 <- calibrate_channel(df, df[[cols$Ch2]], cal)
    df$umol.l.Ch3 <- calibrate_channel(df, df[[cols$Ch3]], cal)
  }
  
  # --- Store condition labels and IDs as attributes ---
  df$condition.Ch1 <- conditions$Ch1
  df$condition.Ch2 <- conditions$Ch2
  df$condition.Ch3 <- conditions$Ch3
  
  df$id.Ch1 <- ids$Ch1
  df$id.Ch2 <- ids$Ch2
  df$id.Ch3 <- ids$Ch3
  
  return(df)
}
