# ================================================
# 01 - Import, Cleaning and Calibration
# ================================================
# Description : Load raw Pyroscience files, apply O2 calibration,
#               remove stabilisation periods, correct for blank drift,
#               and export cleaned per-experiment datasets.
# Author      : Élodie M.A. Jacob
# Date        : 2026
# ================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(presens)

source("R/01_time_vector.R")
source("R/02_clean_pyroscience_data.R")
source("R/03_calibration.R")
source("R/04_process_file.R")
source("R/05_remove_first_hours.R")
source("R/06_get_slope.R")
source("R/07_correct_blank_drift.R")
source("R/08_df_to_long.R")
source("R/09_plot_experiment.R")
source("R/10_assemble_experiment.R")

# ============================================================
# 1. GLOBAL PARAMETERS
# ============================================================

# Calibration coefficients (35-point calibration; Chirugien et al., 2025)
cal <- list(
  c0 = 6.55,    c1 = 0.0863,   c2 = 0.000626,
  c3 = 0.00923, c4 = 9.28e-06, c5 = -0.00106, c6 = 0.000196
)

# Color and label mapping for biological conditions
CONDITION_COLORS <- c(
  "ZooHP"    = "dodgerblue4",
  "ZooATM"   = "gold",
  "BlankATM" = "grey60",
  "BlankHP"  = "black"
)

CONDITION_LABELS <- c(
  "ZooHP"    = "Zoo under pressure",
  "ZooATM"   = "Zoo at ATM pressure",
  "BlankATM" = "Blank ATM",
  "BlankHP"  = "Blank HP"
)

# Output directory
dir.create("outputs/data_corrected", recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. EXPERIMENT 1 – 10 MPa
# ============================================================
# Device 1 (HP) : Ch1 = ZooHP.1 | Ch2 = ZooHP.2 | Ch3 = BlankATM
# Device 2 (ATM): Ch1 = BlankATM | Ch2 = ZooATM  | Ch3 = BlankHP

d_10MPA.1 <- process_file(
  "data/raw/copepoda_respi_hp_10mpa_bis.csv", cal,
  conditions = list(Ch1 = "ZooHP",   Ch2 = "ZooHP",   Ch3 = "BlankATM"),
  ids        = list(Ch1 = "ZooHP.1", Ch2 = "ZooHP.2", Ch3 = "Blank")
)

d_10MPA.2 <- process_file(
  "data/raw/2024-09-19_103252_copepoda_respi_hp_10mpa_bis copie.csv", cal,
  cols       = list(Ch1 = "dphi......A.Ch.2.Main.",
                    Ch2 = "dphi......A.Ch.3.Main.",
                    Ch3 = "dphi......A.Ch.4.Main."),
  conditions = list(Ch1 = "BlankATM", Ch2 = "ZooATM", Ch3 = "BlankHP"),
  ids        = list(Ch1 = "Blank.2",  Ch2 = "Zoo.1",  Ch3 = "BlankHP")
)

# Raw visualization
plot_experiment(d_10MPA.1, d_10MPA.2,
                title    = "Raw data – 10 MPa (Experiment 1)",
                ylim     = c(180, 240),
                type     = "raw",
                filename = "10MPA_exp1.png")

# Remove stabilisation period
d_10MPA.1 <- remove_first_hour(d_10MPA.1)
d_10MPA.2 <- remove_first_hour(d_10MPA.2)

# Blank drift correction
slope_blank_atm_10      <- get_slope(d_10MPA.2, "umol.l.Ch1")  # ATM blank
slope_blank_hp_10       <- get_slope(d_10MPA.2, "umol.l.Ch3")  # HP blank
d_10MPA.1$umol.l.Ch1   <- correct_blank_drift(d_10MPA.1, "umol.l.Ch1", slope_blank_hp_10)
d_10MPA.1$umol.l.Ch2   <- correct_blank_drift(d_10MPA.1, "umol.l.Ch2", slope_blank_hp_10)

# Corrected visualization
plot_experiment(d_10MPA.1, d_10MPA.2,
                title    = "Corrected data – 10 MPa (Experiment 1)",
                ylim     = c(180, 240),
                type     = "corr",
                filename = "10MPA_exp1_corr.png")

# Assemble long-format dataset
d_10MPA <- assemble_experiment(d_10MPA.1, d_10MPA.2)


# ============================================================
# 3. EXPERIMENT 2 – 10 MPa bis
# ============================================================
# Device 1 (HP) : Ch1 = ZooHP.3 | Ch2 = Excluded | Ch3 = ZooATM
# Device 2 (ATM): Ch1 = BlankATM | Ch2 = ZooATM  | Ch3 = BlankHP

d_10MPA_bis.1 <- process_file(
  "data/raw/copepoda_respi_hp_10mpa_tris.csv", cal,
  conditions = list(Ch1 = "ZooHP",    Ch2 = "Blank",    Ch3 = "ZooATM"),
  ids        = list(Ch1 = "ZooHP.3",  Ch2 = "Excluded", Ch3 = "Zoo.2")
)

d_10MPA_bis.2 <- process_file(
  "data/raw/2024-09-23_102817_copepoda_10mpa_tris.csv", cal,
  cols       = list(Ch1 = "dphi......A.Ch.2.Main.",
                    Ch2 = "dphi......A.Ch.3.Main.",
                    Ch3 = "dphi......A.Ch.4.Main."),
  conditions = list(Ch1 = "BlankATM", Ch2 = "ZooATM", Ch3 = "BlankHP"),
  ids        = list(Ch1 = "BlankATM", Ch2 = "Zoo.3",  Ch3 = "BlankHP")
)

# Raw visualization
plot_experiment(d_10MPA_bis.1, d_10MPA_bis.2,
                title    = "Raw data – 10 MPa (Experiment 2)",
                ylim     = c(180, 240),
                type     = "raw",
                filename = "10MPA_exp2.png")

# Remove stabilisation period
d_10MPA_bis.1 <- remove_first_hour(d_10MPA_bis.1)
d_10MPA_bis.2 <- remove_first_hour(d_10MPA_bis.2)

# Blank drift correction
slope_blank_atm_10bis     <- get_slope(d_10MPA_bis.2, "umol.l.Ch1")
d_10MPA_bis.1$umol.l.Ch3  <- correct_blank_drift(d_10MPA_bis.1, "umol.l.Ch3", slope_blank_atm_10bis)
d_10MPA_bis.2$umol.l.Ch2  <- correct_blank_drift(d_10MPA_bis.2, "umol.l.Ch2", slope_blank_atm_10bis)

slope_blank_hp_10bis      <- get_slope(d_10MPA_bis.2, "umol.l.Ch3")
d_10MPA_bis.1$umol.l.Ch1  <- correct_blank_drift(d_10MPA_bis.1, "umol.l.Ch1", slope_blank_hp_10bis)
# Mask artefact before 3.6 h on Ch1
d_10MPA_bis.1$umol.l.Ch1[d_10MPA_bis.1$Time.h > 3.6] <- NA

# Corrected visualization
plot_experiment(d_10MPA_bis.1, d_10MPA_bis.2,
                title    = "Corrected data – 10 MPa (Experiment 2)",
                ylim     = c(180, 250),
                type     = "corr",
                filename = "10MPA_exp2_corr.png")

# Assemble and remove excluded channel
d_10MPA_bis <- assemble_experiment(d_10MPA_bis.1, d_10MPA_bis.2)
d_10MPA_bis <- d_10MPA_bis[d_10MPA_bis$Condition != "Excluded", ]


# ============================================================
# 4. EXPERIMENT 3 – 10 MPa tris
# ============================================================
# Device 1 (HP) : Ch1 = ZooHP.4 | Ch2 = ZooHP.5 | Ch3 = BlankATM
# Device 2 (ATM): Ch1 = ZooATM  | Ch2 = ZooATM  | Ch3 = BlankHP
# Note: HP blank slope not significant; both HP channels still corrected

d_10MPA_tris.1 <- process_file(
  "data/raw/copepoda_respi_hp_10mpa_4.csv", cal,
  conditions = list(Ch1 = "ZooHP",   Ch2 = "ZooHP",   Ch3 = "BlankATM"),
  ids        = list(Ch1 = "ZooHP.4", Ch2 = "ZooHP.5", Ch3 = "Blank")
)

d_10MPA_tris.2 <- process_file(
  "data/raw/2024-09-24_090213_copepoda_10mpa_4.csv", cal,
  cols       = list(Ch1 = "dphi......A.Ch.2.Main.",
                    Ch2 = "dphi......A.Ch.3.Main.",
                    Ch3 = "dphi......A.Ch.4.Main."),
  conditions = list(Ch1 = "ZooATM", Ch2 = "ZooATM", Ch3 = "BlankHP"),
  ids        = list(Ch1 = "Zoo.4",  Ch2 = "Zoo.5",  Ch3 = "BlankHP")
)

# Raw visualization
plot_experiment(d_10MPA_tris.1, d_10MPA_tris.2,
                title    = "Raw data – 10 MPa (Experiment 3)",
                ylim     = c(180, 250),
                type     = "raw",
                filename = "10MPA_exp3.png")

# Remove stabilisation period
d_10MPA_tris.1 <- remove_first_hour(d_10MPA_tris.1)
d_10MPA_tris.2 <- remove_first_hour(d_10MPA_tris.2)

# Blank drift correction
slope_blank_atm_10tris     <- get_slope(d_10MPA_tris.1, "umol.l.Ch3")
d_10MPA_tris.2$umol.l.Ch1  <- correct_blank_drift(d_10MPA_tris.2, "umol.l.Ch1", slope_blank_atm_10tris)
d_10MPA_tris.2$umol.l.Ch2  <- correct_blank_drift(d_10MPA_tris.2, "umol.l.Ch2", slope_blank_atm_10tris)

slope_blank_hp_10tris      <- get_slope(d_10MPA_tris.2, "umol.l.Ch3")
d_10MPA_tris.1$umol.l.Ch2  <- correct_blank_drift(d_10MPA_tris.1, "umol.l.Ch2", slope_blank_hp_10tris)

# Corrected visualization
plot_experiment(d_10MPA_tris.1, d_10MPA_tris.2,
                title    = "Corrected data – 10 MPa (Experiment 3)",
                ylim     = c(180, 250),
                type     = "corr",
                filename = "10MPA_exp3_corr.png")

# Assemble
d_10MPA_tris <- assemble_experiment(d_10MPA_tris.1, d_10MPA_tris.2)


# ============================================================
# 5. EXPERIMENT 4 – 5 MPa
# ============================================================
# Device 1 (HP) : Ch1 = ZooHP.6 | Ch2 = ZooHP.7 | Ch3 = BlankATM
# Device 2 (ATM): Ch1 = ZooATM  | Ch2 = ZooATM  | Ch3 = BlankHP
# Note: truncated at 4 h (artefacts observed beyond this point)

d_5MPA.1 <- process_file(
  "data/raw/copepoda_respi_hp_5mpa_bis.csv", cal,
  conditions = list(Ch1 = "ZooHP",   Ch2 = "ZooHP",   Ch3 = "BlankATM"),
  ids        = list(Ch1 = "ZooHP.6", Ch2 = "ZooHP.7", Ch3 = "Blank")
)

d_5MPA.2 <- process_file(
  "data/raw/2024-09-24_131308_copepoda_5mpa_bis copie.csv", cal,
  cols       = list(Ch1 = "dphi......A.Ch.2.Main.",
                    Ch2 = "dphi......A.Ch.3.Main.",
                    Ch3 = "dphi......A.Ch.4.Main."),
  conditions = list(Ch1 = "ZooATM", Ch2 = "ZooATM", Ch3 = "BlankHP"),
  ids        = list(Ch1 = "Zoo.6",  Ch2 = "Zoo.7",  Ch3 = "BlankHP")
)

# Raw visualization
plot_experiment(d_5MPA.1, d_5MPA.2,
                title    = "Raw data – 5 MPa (Experiment 1)",
                ylim     = c(180, 250),
                type     = "raw",
                filename = "5MPA_exp1.png")

# Remove stabilisation period
d_5MPA.1 <- remove_first_hour(d_5MPA.1)
d_5MPA.2 <- remove_first_hour(d_5MPA.2)

# Truncate at 4 h
d_5MPA.1 <- d_5MPA.1[d_5MPA.1$Time.h <= 4, ]
d_5MPA.2 <- d_5MPA.2[d_5MPA.2$Time.h <= 4, ]

# Blank drift correction
slope_blank_atm_5       <- get_slope(d_5MPA.1, "umol.l.Ch3")
d_5MPA.2$umol.l.Ch1     <- correct_blank_drift(d_5MPA.2, "umol.l.Ch1", slope_blank_atm_5)
d_5MPA.2$umol.l.Ch2     <- correct_blank_drift(d_5MPA.2, "umol.l.Ch2", slope_blank_atm_5)

slope_blank_hp_5        <- get_slope(d_5MPA.2, "umol.l.Ch3")
d_5MPA.1$umol.l.Ch1     <- correct_blank_drift(d_5MPA.1, "umol.l.Ch1", slope_blank_hp_5)
d_5MPA.1$umol.l.Ch2     <- correct_blank_drift(d_5MPA.1, "umol.l.Ch2", slope_blank_hp_5)

# Corrected visualization
plot_experiment(d_5MPA.1, d_5MPA.2,
                title    = "Corrected data – 5 MPa (Experiment 1)",
                ylim     = c(180, 250),
                type     = "corr",
                filename = "5MPA_exp1_corr.png")

# Assemble
d_5MPA <- assemble_experiment(d_5MPA.1, d_5MPA.2)


# ============================================================
# 6. EXPERIMENT 5 – 5 MPa bis
# ============================================================
# Device 1 (HP) : Ch1 = ZooHP.8 (DEAD) | Ch2 = ZooHP.9 | Ch3 = BlankATM
# Device 2 (ATM): Ch1 = ZooATM         | Ch2 = ZooATM  | Ch3 = BlankHP
# Note: truncated at 4 h; ZooHP.8 excluded (animal found dead)
#       Only HP blank available; Ch1 ATM (Zoo.8) slope not significant

d_5MPA_bis.1 <- process_file(
  "data/raw/copepoda_respi_hp_5mpa_tris.csv", cal,
  conditions = list(Ch1 = "ZooHP",   Ch2 = "ZooHP",   Ch3 = "BlankATM"),
  ids        = list(Ch1 = "ZooHP.8", Ch2 = "ZooHP.9", Ch3 = "Blank")
)

d_5MPA_bis.2 <- process_file(
  "data/raw/2024-09-25_093922_copepoda_5mpa_tris copie.csv", cal,
  cols       = list(Ch1 = "dphi......A.Ch.2.Main.",
                    Ch2 = "dphi......A.Ch.3.Main.",
                    Ch3 = "dphi......A.Ch.4.Main."),
  conditions = list(Ch1 = "ZooATM", Ch2 = "ZooATM", Ch3 = "BlankHP"),
  ids        = list(Ch1 = "Zoo.8",  Ch2 = "Zoo.9",  Ch3 = "BlankHP")
)

# Raw visualization
plot_experiment(d_5MPA_bis.1, d_5MPA_bis.2,
                title    = "Raw data – 5 MPa (Experiment 2)",
                ylim     = c(120, 250),
                type     = "raw",
                filename = "5MPA_exp2.png")

# Remove stabilisation period
d_5MPA_bis.1 <- remove_first_hour(d_5MPA_bis.1)
d_5MPA_bis.2 <- remove_first_hour(d_5MPA_bis.2)

# Truncate at 4 h
d_5MPA_bis.1 <- d_5MPA_bis.1[d_5MPA_bis.1$Time.h <= 4, ]
d_5MPA_bis.2 <- d_5MPA_bis.2[d_5MPA_bis.2$Time.h <= 4, ]

# Blank drift correction (HP blank only; ATM blank slope not significant)
slope_blank_hp_5bis      <- get_slope(d_5MPA_bis.2, "umol.l.Ch3")
d_5MPA_bis.2$umol.l.Ch2  <- correct_blank_drift(d_5MPA_bis.2, "umol.l.Ch2", slope_blank_hp_5bis)
d_5MPA_bis.1$umol.l.Ch1  <- correct_blank_drift(d_5MPA_bis.1, "umol.l.Ch1", slope_blank_hp_5bis)
d_5MPA_bis.1$umol.l.Ch2  <- correct_blank_drift(d_5MPA_bis.1, "umol.l.Ch2", slope_blank_hp_5bis)

# Corrected visualization
plot_experiment(d_5MPA_bis.1, d_5MPA_bis.2,
                title    = "Corrected data – 5 MPa (Experiment 2)",
                ylim     = c(180, 250),
                type     = "corr",
                filename = "5MPA_exp2_corr.png")

# Flag dead animal then assemble
d_5MPA_bis.1$condition.Ch1 <- "Excluded"
d_5MPA_bis.1$id.Ch1        <- "ZooHP.8 (dead)"

d_5MPA_bis <- assemble_experiment(d_5MPA_bis.1, d_5MPA_bis.2)
d_5MPA_bis <- d_5MPA_bis[d_5MPA_bis$Condition != "Excluded", ]


# ============================================================
# 7. SAVE OUTPUTS
# ============================================================

write.csv(d_10MPA,      "outputs/data_corrected/d_10MPA_exp1.csv", row.names = FALSE)
write.csv(d_10MPA_bis,  "outputs/data_corrected/d_10MPA_exp2.csv", row.names = FALSE)
write.csv(d_10MPA_tris, "outputs/data_corrected/d_10MPA_exp3.csv", row.names = FALSE)
write.csv(d_5MPA,       "outputs/data_corrected/d_5MPA_exp1.csv",  row.names = FALSE)
write.csv(d_5MPA_bis,   "outputs/data_corrected/d_5MPA_exp2.csv",  row.names = FALSE)
