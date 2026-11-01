## ============================================================
## FILE : R/make_publication_plots.R
## PURPOSE : Publication-quality multi-panel summary figure
## AUTHOR : D. Nerini — June 2026
## ============================================================

## ---- make_publication_plots ----
#' Generate a single multi-panel PDF figure showing O2 time
#' series with retained slope segments highlighted for all
#' individuals (optionally filtered by group).
#'
#' Layout: 4 rows x 5 columns of panels (up to 20 IDs per page).
#' Within each panel:
#'   - Grey scatter points (centred O2).
#'   - Coloured spline segments: red-orange for ATM IDs,
#'     gold for HP IDs; thinner grey for non-retained segments.
#'   - Dashed linear model overlay.
#'   - Vertical dotted blue lines at final knot locations.
#'   - Annotation of mean O2 consumption slope in the corner.
#'
#' @param results_all   Named list; output$results from analyze_oxy
#' @param summary_table data.frame; output$summary_table from analyze_oxy
#' @param group_filter  Optional character; passed to grepl() to
#'                      select a subset of IDs (e.g. "Zoo" or "HP")
#' @param output_dir    Character; path to output directory
#' @return Invisible; saves PDF to output_dir/fig_publication.pdf
make_publication_plots <- function(results_all,
                                   summary_table,
                                   group_filter = NULL,
                                   output_dir   = "oxy_results") {

  ids_plot <- names(results_all)
  if (!is.null(group_filter))
    ids_plot <- ids_plot[grepl(group_filter, ids_plot)]

  n_col <- 5
  n_row <- 4

  col_other    <- "#AAAAAA"
  col_points   <- adjustcolor("black", alpha.f = 0.4)

  get_col_retained <- function(id) {
    if (grepl("HP", id, ignore.case = TRUE)) "#DABD61FF" else "#882255"
  }

  pw     <- 10.27
  ph     <- n_row * 1.8 + 1.2
  pdf_out <- file.path(output_dir, "fig_publication.pdf")

  pdf(pdf_out, width = pw, height = ph)

  par(mfrow  = c(n_row, n_col),
      mar    = c(2.5, 2.5, 2.8, 0.8),
      oma    = c(6, 3, 3, 1))


  for (id in ids_plot) {
    r       <- results_all[[id]]
    x       <- r$x
    y_raw   <- r$y
    y       <- y_raw - mean(y_raw, na.rm = TRUE)

    sf      <- r$slope_final
    sel     <- r$selected_segments
    y_fit_mean <- mean(predict_fit(r$fit_final, r$x), na.rm = TRUE)

    col_retained <- get_col_retained(id)

    # Linear model overlay for comparison
    lm_mod <- tryCatch(lm(y ~ x, na.action = na.omit), error = function(e) NULL)
    lm_ok  <- !is.null(lm_mod)
    if (lm_ok) {
      lm_a <- coef(lm_mod)[1]
      lm_b <- coef(lm_mod)[2]
    }

    plot(x, y,
         pch    = 16, cex = 0.35,
         col    = col_points,
         main   = id, cex.main = 0.9,
         xlab   = "", ylab = "",
         ylim   = range(y, na.rm = TRUE))
    mtext(side = 3, line = 0.15, cex = 0.50,
          text = "- Spline ret.   - - LM")

    for (j in seq_len(nrow(sf))) {
      zz    <- seq(sf$start[j], sf$end[j], length.out = 50)
      yyzz  <- predict_fit(r$fit_final, zz) - y_fit_mean
      lines(zz, yyzz,
            col = ifelse(j %in% sel$segment, col_retained, col_other),
            lwd = ifelse(j %in% sel$segment, 2.5, 1))
    }

    if (lm_ok) abline(a = lm_a, b = lm_b, col = "black", lwd = 1.5, lty = 2)
    abline(v = r$knots_final, lty = 2, col = "steelblue", lwd = 0.7)

    # Slope annotation
    st <- summary_table[summary_table$id == id, ]
    if (nrow(st) > 0 && is.finite(st$conso_O2)) {
      usr   <- par("usr")
      x_pos <- usr[1] + 0.03 * diff(usr[1:2])
      y_pos <- usr[3] + 0.06 * diff(usr[3:4])
      dy    <- diff(usr[3:4]) * 0.10   # espacement vertical entre lignes
      
      # Ligne 1 : Slope
      label1 <- "Slope (umol/L/h) = "
      val1   <- as.character(round(st$conso_O2, 4))
      text(x_pos, y_pos + 2 * dy,
           labels = label1, adj = c(0, 0), cex = 0.50, col = "grey50", font = 3)
      text(x_pos + strwidth(label1, cex = 0.50, font = 3),
           y_pos + 2 * dy,
           labels = val1, adj = c(0, 0), cex = 0.50, col = col_retained, font = 2)
      
      # Ligne 2 : DW
      if (!is.null(st$DW_mg) && is.finite(st$DW_mg)) {
        label2 <- "DW (mg) = "
        val2   <- as.character(round(st$DW_mg, 4))
        text(x_pos, y_pos + dy,
             labels = label2, adj = c(0, 0), cex = 0.50, col = "grey50", font = 3)
        text(x_pos + strwidth(label2, cex = 0.50, font = 3),
             y_pos + dy,
             labels = val2, adj = c(0, 0), cex = 0.50, col = col_retained, font = 2)
      }
      
      # Ligne 3 : Mass-normalised respiration
      if (!is.null(st$Resp_umol_mgDW_h) && is.finite(st$Resp_umol_mgDW_h)) {
        label3 <- "O2 cons. (umol/mgDW/h) = "
        val3   <- as.character(round(st$Resp_umol_mgDW_h, 5))
        text(x_pos, y_pos,
             labels = label3, adj = c(0, 0), cex = 0.50, col = "grey50", font = 3)
        text(x_pos + strwidth(label3, cex = 0.50, font = 3),
             y_pos,
             labels = val3, adj = c(0, 0), cex = 0.50, col = col_retained, font = 2)
      }
    }
    
  }

  mtext("Time", side = 1, outer = TRUE, line = 0.1, cex = 1.0)
  mtext("Centred O2 (umol/L)", side = 2, outer = TRUE, line = 1.5, cex = 1.0)
  mtext("O2 time series - retained slopes highlighted",
        side = 3, outer = TRUE, line = 1.0, cex = 1.1, font = 2)
  
  # Shared legend — dans la marge oma basse
  plot.new()
  #legend("bottom",
  #       legend = c("Knots","LM","Spline retained (ATM)", "Spline retained (HP)",
  #                  "Other Spline segment"),
   #      col    = c("steelblue","black", "#882255", "#DABD61FF", col_other),
  #       lwd    = c(2.5, 2.5, 1, 1.5, 0.7),
   #      lty    = c(1, 1, 1, 2, 2),
     #    horiz  = TRUE, bty = "n", cex = 1,
   #      inset = c(0, -0.2),xpd=T) 
  
  dev.off()
  cat("Publication figure saved:", pdf_out, "\n")

  invisible()
}