#' Remove the first hour of data (stabilization period)
#' Assumes 1 measurement per 30 seconds -> 120 rows = 1 hour
remove_first_hour <- function(df) df[-c(1:120), ]