## ============================================================
## 03 — O2 CONSUMPTION ANALYSIS FROM TIME SERIES
# ================================================
# Description : Penalised LAD B-spline + asymptotic slope inference
# Author      :  David Nerini
# Date        : 2026
# ================================================

#This script estimates oxygen consumption rates from respirometry time series.
# Rather than fitting a single linear regression over the entire recording, 
#it allows the O2 signal to be decomposed into segments with statistically distinct
#slopes, which is more appropriate for noisy biological signals that may include
# artefacts, behavioural transitions, or instrumental drift.
#Method overview
#1. Penalised LAD spline fitting
#Each time series is fitted using a B-spline with a roughness penalty,
#optimised by least absolute deviations (LAD/median regression) rather 
#than least squares. This makes the fit robust to outliers and sudden jumps in the signal.
#2. Slope inference per segment
#For each spline segment, a slope (O2 change per time unit) is 
#estimated with its standard error using an asymptotic sandwich variance
#estimator appropriate for LAD regression. Each slope is tested against 
#zero and classified as decreasing (dec), increasing (inc), or non-significant (NS).
#3. Iterative segment merging
#Adjacent segments with statistically similar slopes are iteratively merged, 
#reducing over-segmentation. The merging criterion is a pairwise z-test on slope 
#differences, controlled by ALPHA_SIM.
#4. O2 consumption selection
#Among all decreasing segments, the script identifies the cluster of segments 
#closest to the modal slope using a weighted kernel density estimate. The final 
#O2 consumption estimate is the duration-weighted mean slope of the retained 
#segments.

# ------------------------------------------------------------
# USAGE
# ------------------------------------------------------------
# This script expects a data.frame called 'zoo_all' in the
# environment, with at least four columns:
#   [1] time  (numeric, hours or seconds)
#   [2] O2    (numeric, dissolved oxygen concentration)
#   [3] any additional covariates (ignored)
#   [4] ID    (character, grouping variable — e.g. "Zoo.1", "Zoo.2 HP")
# ------------------------------------------------------------

# Load all helper functions
source("R/12_utils_spline.R")
source("R/13_slope_inference.R")
source("R/14_merge_segments.R")
source("R/15_select_o2_segments.R")
source("R/16_analyze_oxy.R")
source("R/17_make_publication_plots.R")

# ------------------------------------------------------------
# 1. PARAMETERS
# ------------------------------------------------------------

OUTPUT_DIR   <- "oxy_results"  # output directory for PDFs and CSV

DF           <- 40             # degrees of freedom for initial B-spline basis
LAMBDA       <- 3              # roughness penalty weight
ALPHA_SLOPE  <- 0.05           # significance level for slope classification
ALPHA_SIM    <- 0.4            # similarity threshold for segment merging
MIN_N        <- 20             # minimum points required per segment
MIN_DURATION <- 0              # minimum time duration for a segment to be eligible

CHAMBER_VOL_mL <- 5            # respiration chamber volume (mL)

# Allometric relationship for dry weight estimation:
# DW (mg) = ALLO_A × Length (mm) ^ ALLO_B
# Source: doi:10.6620/ZS.2017.56-13  (Neocalanus flemingeri / plumchrus)
ALLO_A <- 0.01841
ALLO_B <- 2.457

# ------------------------------------------------------------
# 2. SPLINE-BASED O2 CONSUMPTION ANALYSIS
# ------------------------------------------------------------
out <- analyze_oxy(
  data         = zoo_all,
  df           = DF,
  lambda       = LAMBDA,
  alpha_slope  = ALPHA_SLOPE,
  alpha_sim    = ALPHA_SIM,
  min_n        = MIN_N,
  min_duration = MIN_DURATION,
  output_dir   = OUTPUT_DIR
)


# ------------------------------------------------------------
# 3. QUICK INSPECTION
# ------------------------------------------------------------
cat("\n====== SUMMARY TABLE ======\n")
print(out$summary_table)

cat("\n====== EXAMPLE SLOPE TABLE : Zoo.3 ======\n")
print(out$results[["Zoo.3"]]$slope_final)

# ------------------------------------------------------------
# 4. BIOMASS ESTIMATION FROM PROSOME LENGTH
# ------------------------------------------------------------
# Dry weight (DW) estimated from prosome length using the
# allometric relationship:
#   DW (mg) = ALLO_A × Length (mm) ^ ALLO_B
# See PARAMETERS section for coefficient values and source.

length_data <- read.csv(
  "data/raw/length_respi_Neocalanus_gracilis_final.csv",
  sep = ";"
)

length_data$DW_mg <- ALLO_A * (length_data[["Length..mm."]])^ALLO_B

# ------------------------------------------------------------
# 5. MASS-NORMALISED RESPIRATION RATE
# ------------------------------------------------------------
# Units: µmol O2 / mg DW / h
#
# Formula:
#   R = Slope (µmol/L/h) × chamber volume (L) / DW (mg)
#
# Slope is negative by convention (O2 is consumed).
# The sign is preserved here; reverse for display if needed.

CHAMBER_VOL_L <- CHAMBER_VOL_mL / 1000   # 5 mL → 0.005 L

# Utiliser un nom distinct pour le tableau résumé
zoo_summary <- dplyr::left_join(out$summary_table, length_data, by = c("id" = "ID"))

# Calcul
zoo_summary$Resp_umol_mgDW_h <- (zoo_summary$conso_O2 * CHAMBER_VOL_L) / zoo_summary$DW_mg

# Aperçu
zoo_summary[, c("id", "conso_O2", "DW_mg", "Resp_umol_mgDW_h")]

# ------------------------------------------------------------
# 6. PUBLICATION FIGURE 
# ------------------------------------------------------------
make_publication_plots(
  results_all   = out$results,
  summary_table = zoo_summary, 
  group_filter  = "Zoo",
  output_dir    = OUTPUT_DIR
)

