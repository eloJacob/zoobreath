## ============================================================
## FILE : R/analyze_oxy.R
## PURPOSE : Per-ID analysis loop — spline fitting, segment
##           merging, slope inference, O2 segment selection,
##           per-ID PDF diagnostics, global PDF and summary CSV
## AUTHOR : D. Nerini — June 2026
## ============================================================

## ---- analyze_oxy ----
#' Main analysis function for a full dataset of respirometry
#' time series. Processes one individual (ID) at a time.
#'
#' Pipeline per ID:
#'   1. Fit initial penalised LAD B-spline (full basis, no internal knots).
#'   2. Extract initial knot candidates from short segments (repair_short_segments).
#'   3. Refit on initial knots.
#'   4. Iteratively merge segments (merge_segments).
#'   5. Refit on final knots.
#'   6. Run slope inference and hierarchical clustering of slopes.
#'   7. Select O2-consuming segments (select_o2_segments).
#'   8. Save per-ID diagnostic PDF.
#'   9. After all IDs: save global summary PDF and summary CSV.
#'
#' @param data          data.frame with columns: time, O2, ... (any extra columns ok)
#'                      The function expects columns 1 and 2 to be time and O2,
#'                      and column 4 to be the ID grouping variable.
#'                      Adjust column indices to match your data structure.
#' @param df            Degrees of freedom for initial B-spline basis
#' @param lambda        Roughness penalty weight
#' @param alpha_slope   Significance level for individual slope classification
#' @param alpha_sim     Similarity threshold for segment merging and expansion
#' @param min_n         Minimum observations required per segment
#' @param min_duration  Minimum time duration for a segment to be eligible for
#'                      O2 consumption estimation
#' @param output_dir    Character; path to output directory
#' @return Invisible list with:
#'   - results      : named list, one element per ID; each contains all
#'                    intermediate objects and final estimates
#'   - summary_table: data.frame with id, n_segments, n_seg_selected,
#'                    slope_mode, conso_O2
analyze_oxy <- function(data,
                        df           = 40,
                        lambda       = 3,
                        alpha_slope  = 0.05,
                        alpha_sim    = 0.4,
                        min_n        = 20,
                        min_duration = 0,
                        output_dir   = "oxy_results") {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  oxy <- by(data[, 1:2], data[, 4], list)
  ids <- names(oxy)
  
  results_all   <- list()
  summary_table <- data.frame()
  
  for (id in ids) {
    
    cat("\n========== Processing:", id, "==========\n")
    
    dat <- oxy[[id]]
    x   <- dat[[1]]
    y   <- dat[[2]]
    
    if (sum(is.finite(x) & is.finite(y)) < min_n * 2) {
      cat("  Not enough data — skipped.\n"); next
    }
    
    ## ---- initial fit ------------------------------------------------
    fit0 <- tryCatch(
      fit_quant_pspline(x, y, df, lambda),
      error = function(e) { cat("  fit0 failed:", e$message, "\n"); NULL })
    if (is.null(fit0)) next
    
    br0    <- repair_short_segments(
      fit0$x_used,
      sort(unique(c(min(fit0$x_used), fit0$knots, max(fit0$x_used)))),
      min_n)
    knots0 <- br0[-c(1, length(br0))]
    
    fit0 <- tryCatch(
      fit_quant_pspline(x, y, df, lambda, knots0),
      error = function(e) { cat("  fit0 refit failed:", e$message, "\n"); NULL })
    if (is.null(fit0)) next
    
    ## ---- iterative merging ------------------------------------------
    br_final    <- merge_segments(x, y, br0,
                                  df, lambda, alpha_slope, alpha_sim, min_n)
    knots_final <- br_final[-c(1, length(br_final))]
    
    ## ---- final fit --------------------------------------------------
    fit_final <- tryCatch(
      fit_quant_pspline(x, y, df, lambda, knots_final),
      error = function(e) { cat("  fit_final failed:", e$message, "\n"); NULL })
    if (is.null(fit_final)) next
    
    test_final <- tryCatch(
      slope_inference(fit_final, alpha_slope, alpha_sim),
      error = function(e) { cat("  inference failed:", e$message, "\n"); NULL })
    if (is.null(test_final)) next
    
    slope_final <- test_final$slope_table
    
    ## ---- clustering -------------------------------------------------
    D        <- abs(test_final$z_sim); diag(D) <- 0
    hc       <- hclust(as.dist(D), method = "average")
    grp      <- cutree(hc, h = 2)
    slope_final$group <- grp
    
    ## ---- O2 selection -----------------------------------------------
    o2 <- select_o2_segments(slope_final, test_final$p_sim,
                             alpha_sim, min_duration)
    
    cat("  Modal slope    :", o2$slope_mode, "\n")
    cat("  O2 consumption :", o2$conso_O2,   "\n")
    
    ## ---- colours for classification plot ----------------------------
    cols        <- ifelse(slope_final$class == "NS",  "black",
                          ifelse(slope_final$class == "inc", "red", "blue"))
    resid_final <- fit_final$y_used - fit_final$fitted
    outliers    <- abs(resid_final) > 5 * mad(resid_final)
    nseg_final  <- nrow(slope_final)
    xx          <- seq(min(x), max(x), length.out = 1000)
    yy0         <- predict_fit(fit0,      xx)
    yyf         <- predict_fit(fit_final, xx)
    
    ## ---- per-ID PDF -------------------------------------------------
    pdf_file <- file.path(output_dir, paste0(id, "_analysis.pdf"))
    pdf(pdf_file, width = 10, height = 7)
    
    ## page 1 : dendrogram
    if (nseg_final > 2) {
      par(mar = c(5, 4, 4, 2))
      plot(hc,
           main   = paste0("Dendrogram of slope similarity \u2014 ", id),
           xlab   = "Segment", ylab = "Distance (|z|)",
           sub    = "", hang = -1,
           labels = paste0("Seg ", seq_len(nseg_final)))
      abline(h = 2, lty = 2, col = "red", lwd = 2)
      legend("topright", legend = "Merge threshold (h = 2)",
             lty = 2, col = "red", lwd = 2, bty = "n")
    }
    
    ## page 2 : slope classification
    par(mar = c(5, 4, 4, 2))
    plot(x, y, pch = 16, cex = 0.5,
         main = paste0("Slope classification \u2014 ", id),
         xlab = "Time", ylab = "O2")
    lines(xx, yy0, lwd = 2, lty = 2, col = "grey40")
    lines(xx, yyf, lwd = 3, col = "black")
    
    for (j in seq_len(nseg_final)) {
      zz   <- seq(slope_final$start[j], slope_final$end[j], length.out = 50)
      yyzz <- predict_fit(fit_final, zz)
      lines(zz, yyzz, col = cols[j], lwd = 4)
      xmid <- mean(zz)
      text(xmid, predict_fit(fit_final, xmid),
           labels = paste0("s=",  round(slope_final$slope[j],   4),
                           "\np=", round(slope_final$p_value[j], 3)),
           cex = 0.6, col = cols[j], pos = 3)
    }
    abline(v = br0[-c(1, length(br0))], lty = 3, col = "grey70")
    abline(v = knots_final, lty = 2, lwd = 2, col = "blue")
    legend("topright",
           legend = c("Initial fit", "Final fit", "NS", "Increasing",
                      "Decreasing", "Initial knots", "Final knots"),
           lwd = c(2, 3, 4, 4, 4, 1, 2), lty = c(2, 1, 1, 1, 1, 3, 2),
           col = c("grey40", "black", "black", "red", "blue", "grey70", "blue"),
           bty = "n", cex = 0.8)
    
    dev.off()
    cat("  PDF saved:", pdf_file, "\n")
    
    ## ---- store results ----------------------------------------------
    results_all[[id]] <- list(
      id                = id,
      x                 = x,
      y                 = y,
      fit0              = fit0,
      fit_final         = fit_final,
      test_final        = test_final,
      slope_final       = slope_final,
      selected_segments = o2$selected,
      conso_O2          = o2$conso_O2,
      slope_mode        = o2$slope_mode,
      hc                = hc,
      knots_final       = knots_final,
      br0               = br0,
      outliers          = outliers
    )
    
    summary_table <- rbind(summary_table, data.frame(
      id             = id,
      n_segments     = nseg_final,
      n_seg_selected = nrow(o2$selected),
      slope_mode     = o2$slope_mode,
      conso_O2       = o2$conso_O2,
      stringsAsFactors = FALSE
    ))
  }
  
  ## ---- global PDF ---------------------------------------------------
  pdf_global <- file.path(output_dir, "all_ids_retained_slopes.pdf")
  pdf(pdf_global, width = 10, height = 7)
  
  for (id in names(results_all)) {
    r   <- results_all[[id]]
    x   <- r$x;  y <- r$y
    xx  <- seq(min(x), max(x), length.out = 1000)
    sf  <- r$slope_final
    sel <- r$selected_segments
    kf  <- r$knots_final
    
    conso_label <- if (is.finite(r$conso_O2))
      paste0("Mean O2 consumption = ", round(r$conso_O2, 5), " O2/time unit")
    else "No retained segment"
    
    par(mar = c(5, 4, 5, 2))
    plot(x, y, pch = 16, cex = 0.5,
         main = paste0(id, "\n", conso_label),
         xlab = "Time", ylab = "O2")
    
    for (j in seq_len(nrow(sf))) {
      zz   <- seq(sf$start[j], sf$end[j], length.out = 50)
      yyzz <- predict_fit(r$fit_final, zz)
      lines(zz, yyzz,
            col = ifelse(j %in% sel$segment, "red", "grey60"),
            lwd = ifelse(j %in% sel$segment, 5, 3))
    }
    
    abline(v = kf, lty = 2, col = "blue", lwd = 2)
    if (length(kf) > 0) {
      usr <- par("usr")
      points(kf, rep(usr[3] + 0.02 * diff(usr[3:4]), length(kf)),
             pch = 16, cex = 2, col = "blue")
    }
    
    points(r$fit_final$x_used[r$outliers],
           r$fit_final$y_used[r$outliers],
           pch = 21, bg = "darkgrey", col = "black", cex = 1.8, lwd = 1.5)
    
    if (nrow(sel) > 0 && is.finite(r$conso_O2)) {
      x0    <- min(sel$start)
      y_ref <- predict_fit(r$fit_final, x0)
      abline(a = y_ref - (-r$conso_O2) * x0,
             b = -r$conso_O2, col = "darkred", lwd = 2, lty = 2)
    }
    
    legend("topright",
           legend = c("Retained segments", "Other segments",
                      "Final knots", "Outliers", "Mean slope"),
           col    = c("red", "grey60", "blue", "black", "darkred"),
           lwd    = c(5, 3, 2, NA, 2),
           lty    = c(1, 1, 2, NA, 2),
           pch    = c(NA, NA, 16, 21, NA),
           pt.bg  = c(NA, NA, "blue", "darkgrey", NA),
           bty    = "n", cex = 0.8)
  }
  dev.off()
  cat("\nGlobal PDF saved:", pdf_global, "\n")
  
  ## ---- summary CSV --------------------------------------------------
  csv_file <- file.path(output_dir, "summary_conso_O2.csv")
  write.csv(summary_table, csv_file, row.names = FALSE)
  cat("Summary CSV saved:", csv_file, "\n")
  
  invisible(list(results = results_all, summary_table = summary_table))
}