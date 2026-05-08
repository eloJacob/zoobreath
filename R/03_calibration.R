
#' Modified Stern-Volmer equation (Bittig et al. 2018)
#' @param c0      Calibration coefficient 
#' @param t       Temperature (°C)
#' @param PR      Phase shifts measurements (°)
Oxy_corr <- function(t, PR, c0, c1, c2, c3, c4, c5, c6) {
  KSV <- c0 + c1 * t + c2 * t^2
  P0  <- c3 + c4 * t
  PC  <- c5 + c6 * PR
  corr <- ((P0 / PC) - 1) / KSV
  return(corr)
}


#' Step 1: Apply sensor calibration (Oxy_corr)
#' Output is in atm
#'
#' @param temp    Temperature vector (°C)
#' @param dphi    Raw dphi values from sensor
#' @param cal     Calibration coefficient list
#' @return        O2 values in mbar (as atmospheric fraction)
apply_oxy_corr <- function(temp, dphi, cal) {
  Oxy_corr(temp, dphi,
           cal$c0, cal$c1, cal$c2,
           cal$c3, cal$c4, cal$c5, cal$c6)
}

#' Step 2: Convert atmospheric O2 fraction to umol/L
#' Oxy_corr output (o2_atm) is multiplied by 1013.25 to convert to hPa,
#' then passed to o2_unit_conv
#'
#' @param o2_atm    O2 values as atmospheric fraction (output of Oxy_corr)
#' @param temp      Temperature vector (°C)
#' @param ref_pres  Reference pressure in mbar (single value, taken from row 2)
#' @return          O2 values in umol/L
convert_to_umol <- function(o2_atm, temp, ref_pres) {
  presens::o2_unit_conv(
    o2_atm * 1013.25,        # atmospheric fraction -> hPa
    from     = "hPa",
    to       = "umol_per_l",
    temp     = temp,
    air_pres = ref_pres / 1000  # mbar -> bar
  )$umol_per_l
}

#' Full calibration pipeline for one channel
#'
#' @param df        Cleaned pyroscience dataframe from Clean.pyroscience.data
#' @param dphi_col  Raw dphi values (vector)
#' @param cal       Calibration coefficient list
#' @return          O2 values in umol/L
calibrate_channel <- function(df, dphi_col, cal) {
  o2_atm <- apply_oxy_corr(df$`Temperature[C]`, dphi_col, cal)
  convert_to_umol(
    o2_atm   = o2_atm,
    temp     = df$`Temperature[C]`,
    ref_pres = df[2, "Pressure[mbar]"]
  )
}
