#' Convert one processed dataframe to long format
#' Uses condition and id columns stored by process_file()
#'
#' @param df  Processed dataframe (output of process_file)
#' @return    Long-format dataframe with columns: Time, O2, Condition, ID
df_to_long <- function(df) {
  rbind(
    data.frame(Time      = df$Time.h,
               O2        = df$umol.l.Ch1,
               Condition = df$condition.Ch1,
               ID        = df$id.Ch1),
    data.frame(Time      = df$Time.h,
               O2        = df$umol.l.Ch2,
               Condition = df$condition.Ch2,
               ID        = df$id.Ch2),
    data.frame(Time      = df$Time.h,
               O2        = df$umol.l.Ch3,
               Condition = df$condition.Ch3,
               ID        = df$id.Ch3)
  )
}
