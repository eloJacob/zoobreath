#' Plot raw O2 data for one experiment (2 devices)
#' Colors and legend are driven by the condition labels stored in the dataframes
#'
#' @param d1     Processed dataframe from device 1 (HP chambers)
#' @param d2     Processed dataframe from device 2 (ATM chambers)
#' @param title  Plot title
#' @param ylim   Y-axis limits (default c(150, 250))
#' @param type   "raw" or "corr" — determines output subfolder
#' @param filename  Name of the output file (e.g. "10MPA_exp1.png")
plot_experiment <- function(d1, d2, title, ylim = c(150, 250),
                            type = "raw", filename = NULL) {
  
  # ---- Output path ----
  if (!is.null(filename)) {
    folder <- if (type == "corr") {
      "outputs/figures/experiment_vizualization_corr"
    } else {
      "outputs/figures/experiment_vizualization_raw"
    }
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    filepath <- file.path(folder, filename)
    png(filepath, width = 1200, height = 800, res = 150)
  }
  
  # ---- Build long-format data ----
  long <- rbind(df_to_long(d1), df_to_long(d2))
  long$Color <- CONDITION_COLORS[long$Condition]
  
  # ---- Layout: 3/4 O2, 1/4 Temperature ----
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
  par(mar = c(2, 4, 3, 2))
  
  # ---- O2 panel ----
  plot(long$Time, long$O2,
       col  = long$Color,
       pch  = 16, cex = 0.5,
       xlab = "", ylab = "O2 (µmol/L)",
       ylim = ylim, main = title)
  
  present_conditions <- unique(long$Condition)
  legend("bottomleft",
         legend = CONDITION_LABELS[present_conditions],
         col    = CONDITION_COLORS[present_conditions],
         pch    = 16, cex = 0.7)
  
  # ---- Temperature panel ----
  par(mar = c(4, 4, 1, 2))
  
  temp_all <- rbind(
    data.frame(Time = d1$Time.h, Temp = d1$`Temperature[C]`),
    data.frame(Time = d2$Time.h, Temp = d2$`Temperature[C]`)
  )
  temp_all <- temp_all[order(temp_all$Time), ]
  
  plot(temp_all$Time, temp_all$Temp,
       col  = "firebrick",
       pch  = 16, cex = 0.3,
       xlab = "Time (h)", ylab = "T° (°C)",
       main = "")
  
  # ---- Close device if saving ----
  if (!is.null(filename)) dev.off()
  
  # Reset layout
  par(mfrow = c(1, 1))
}
