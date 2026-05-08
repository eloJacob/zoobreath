# ================================================
# 02 - Merge & explore oxygen time series
# ================================================
# Description : Load corrected experiments, merge into a single dataset,
#               and produce exploratory figures of O2 time series.
# Author      : Élodie M.A. Jacob
# Date        : 2026
# ================================================

library(dplyr)
library(ggplot2)

# ---- 1. Load corrected experiments ----

d_10MPA      <- read.csv("outputs/data_corrected/d_10MPA_exp1.csv")
d_10MPA_bis  <- read.csv("outputs/data_corrected/d_10MPA_exp2.csv")
d_10MPA_tris <- read.csv("outputs/data_corrected/d_10MPA_exp3.csv")
d_5MPA       <- read.csv("outputs/data_corrected/d_5MPA_exp1.csv")
d_5MPA_bis   <- read.csv("outputs/data_corrected/d_5MPA_exp2.csv")

# ---- 2. Combine all experiments ----

all_experiments <- bind_rows(
  mutate(d_10MPA,      Experiment = "10MPA_exp1", Pressure_MPa = 10),
  mutate(d_10MPA_bis,  Experiment = "10MPA_exp2", Pressure_MPa = 10),
  mutate(d_10MPA_tris, Experiment = "10MPA_exp3", Pressure_MPa = 10),
  mutate(d_5MPA,       Experiment = "5MPA_exp1",  Pressure_MPa = 5),
  mutate(d_5MPA_bis,   Experiment = "5MPA_exp2",  Pressure_MPa = 5)
)

# ATM controls are at 0.1 MPa
all_experiments$Pressure_MPa[all_experiments$Condition == "ZooATM"] <- 0.1

# Save merged dataset
write.csv(all_experiments, "outputs/data_corrected/all_experiments.csv", row.names = FALSE)

# ---- 3. Subset zooplankton chambers only ----

zoo_all <- all_experiments[all_experiments$Condition %in% c("ZooHP", "ZooATM"), ]

# ---- 4. Figures ----

# -- 4a. O2 time series by pressure (3 panels) --
p_by_pressure <- ggplot(zoo_all, aes(x = Time, y = O2, color = ID, group = ID)) +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ Pressure_MPa, labeller = label_both, nrow = 3) +
  labs(title = "O2 consumption – 0.1 / 5 / 10 MPa",
       x = "Time (h)", y = "O2 (µmol/L)", color = "Individual") +
  theme_bw(base_size = 11)

ggsave("outputs/figures/exploration/all_series_by_pressure.png",
       p_by_pressure, width = 8, height = 10, dpi = 300)

# -- 4b. O2 time series by pressure × experiment --
p_by_exp <- ggplot(zoo_all, aes(x = Time, y = O2, color = ID, group = ID)) +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_grid(Pressure_MPa ~ Experiment, labeller = label_both) +
  labs(title = "O2 consumption – by pressure & experiment",
       x = "Time (h)", y = "O2 (µmol/L)", color = "Individual") +
  theme_bw(base_size = 11) +
  theme(strip.text.y = element_text(face = "bold", size = 12),
        strip.text.x = element_text(size = 10))

ggsave("outputs/figures/exploration/all_series_by_exp.png",
       p_by_exp, width = 14, height = 10, dpi = 300)

