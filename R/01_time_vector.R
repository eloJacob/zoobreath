#' Convert date and time columns to decimal hours from experiment start
#'
#' @param a  Date column (character or factor)
#' @param b  Time column (character or factor)
#' @param d  Start day   (e.g. 23)
#' @param m  Start month (e.g. 04)
#' @param y  Start year  (e.g. 2024)
#' @param date_format  Date format: "d/m/y" (default) or "d-m-y"
#' @return   Numeric vector of elapsed time in hours since first observation

time_to_hours <- function(a, b, d, m, y, date_format = "d/m/y") {
  
  dd <- chron::chron(
    dates  = as.character(a),
    times  = as.character(b),
    format = c(dates = date_format, times = "h:m:s"),
    origin = c(month = m, day = d, year = y)
  )
  
  p1 <- as.numeric(dd)
  h1 <- (p1 - p1[1]) * 24
  
  return(h1)
}