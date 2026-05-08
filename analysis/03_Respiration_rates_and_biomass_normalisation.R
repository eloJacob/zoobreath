# ================================================
# 03 - Respiration rates & biomass normalisation
# ================================================
# Description : Compute per-individual O2 consumption slopes (µmol O2/L/ind/h),
#               estimate dry weight from body length, and calculate
#               respiration rates normalised to dry weight (µmol O2/mg DW/h).
# Author      : Élodie M.A. Jacob
# Date        : 2026
# ================================================

library(dplyr)
library(ggplot2)
library(rstatix)


# ============================================================
# 0. LOAD DATA
# ============================================================

# Merged dataset produced by script 02
all_experiments <- read.csv("outputs/data_corrected/all_experiments.csv")

# Keep zooplankton chambers only (exclude blanks)
zoo_all <- all_experiments[all_experiments$Condition %in% c("ZooHP", "ZooATM"), ]


# ============================================================
# 1. INDIVIDUAL RESPIRATION SLOPES (µmol O2 / L / h)
# ============================================================
# One linear regression per individual (O2 ~ Time).
# The slope (coefficient of Time) gives the O2 change rate in µmol/L/h.
# Slopes are expected to be negative (O2 decreases over time).

zoo_all$Slope <- NA

for (i in unique(zoo_all$ID)) {
  idx               <- zoo_all$ID == i
  fit               <- lm(O2 ~ Time, data = zoo_all[idx, ])
  zoo_all$Slope[idx] <- coef(fit)["Time"]
}


# ============================================================
# 2. BIOMASS ESTIMATION FROM PROSOME LENGTH
# ============================================================
# Dry weight (DW) estimated from prosome length using the allometric
# relationship for Neocalanus flemingeri / plumchrus (Pacific):
#   DW (mg) = 0.01841 × Length (mm) ^ 2.457
# Source: doi:10.6620/ZS.2017.56-13

length_data       <- read.csv("data/raw/length_respi_Neocalanus_gracilis_final.csv", sep = ";")
length_data$DW_mg <- 0.01841 * (length_data[, "Length..mm."])^2.457

zoo_all <- left_join(zoo_all, length_data, by = "ID")


# ============================================================
# 3. MASS-NORMALISED RESPIRATION RATE (µmol O2 / mg DW / h)
# ============================================================
# Slope (µmol/L/h) × chamber volume (L) / dry weight (mg).
# Result is negative by convention (O2 consumption); sign is
# reversed only for display purposes in section 5.

CHAMBER_VOL_mL <- 5   # volume of the respiration chamber in mL

zoo_all$Resp_umol_mgDW_h <- (zoo_all$Slope * (CHAMBER_VOL_mL / 1000)) / zoo_all$DW_mg


# ============================================================
# 4. SUMMARY TABLE (one row per individual)
# ============================================================

resp_summary <- zoo_all %>%
  group_by(ID, Condition, Pressure_MPa, Experiment, DW_mg) %>%
  summarise(
    Slope_umol_L_h   = unique(Slope),
    Resp_umol_mgDW_h = unique(Resp_umol_mgDW_h),
    .groups = "drop"
  )

# Ordered pressure factor for consistent plotting
resp_summary$Pressure_char <- factor(
  as.character(resp_summary$Pressure_MPa),
  levels = c("0.1", "5", "10")
)

print(resp_summary)


# ============================================================
# 5. STATISTICAL TEST
# ============================================================
# Kruskal-Wallis test (non-parametric, small and unequal group sizes).
# O2 consumption expressed as positive values for the test.

kw       <- kruskal.test(-Resp_umol_mgDW_h ~ Pressure_char, data = resp_summary)
kw_label <- paste0("Kruskal-Wallis: p = ", round(kw$p.value, 3))

cat(kw_label, "\n")


# ============================================================
# 6. FIGURE – Individual respiration rates by pressure
# ============================================================
# Display convention: O2 consumption shown as positive values (-Slope).
# Each point = one individual; large point = group mean; bars = ± 1 SD.
# Kruskal-Wallis p-value annotated in the top-right corner.


pressure_colors <- c("0.1" = "#DABD61FF",
                     "5"   = "#117733",
                     "10"  = "#882255")

df_summary <- resp_summary %>%
  group_by(Pressure_char) %>%
  summarise(
    mean_resp   = mean(  -Resp_umol_mgDW_h * 100, na.rm = TRUE),
    sd_resp     = sd(    -Resp_umol_mgDW_h * 100, na.rm = TRUE),
    median_resp = median(-Resp_umol_mgDW_h * 100, na.rm = TRUE),
    .groups = "drop"
  )

p_resp <- ggplot(
  resp_summary,
  aes(
    x     = Pressure_char,
    y     = -Resp_umol_mgDW_h * 100,
    fill  = Pressure_char,
    color = Pressure_char
  )
) +
  
  # Individual data points
  geom_point(
    alpha    = 0.7,
    size     = 3,
    shape    = 16,
    color="black",
    position = position_jitter(width = 0.08, seed = 42)
  ) +
  
  # Mean ± sd error bars
  geom_errorbar(
    data = df_summary,
    aes(
      x    = Pressure_char,
      y    = mean_resp,
      ymin = mean_resp - sd_resp,
      ymax = mean_resp + sd_resp
    ),
    width     = 0.2,
    linewidth = 1.2
  ) +
  
  # Group mean points
  geom_point(
    data  = df_summary,
    aes(x = Pressure_char, y = mean_resp),
    shape = 16,
    size  = 6
  ) +
  
  # Median point
  geom_point(
    data  = df_summary,
    aes(x = Pressure_char, y = median_resp),
    shape = 17,
    size  = 6
  ) +
  
  # Kruskal-Wallis p-value
  annotate(
    "text",
    x = Inf, y = Inf,
    label = kw_label,
    hjust = 1.05, vjust = 1.5,
    size  = 3.5, color = "grey40"
  ) +
  
  scale_fill_manual(values  = pressure_colors) +
  scale_color_manual(values = pressure_colors) +
  
  scale_x_discrete(
    labels = c("0.1" = "0.1 MPa", "5" = "5 MPa", "10" = "10 MPa")
  ) +
  
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  labs(
    x = "",
    y = expression(O[2]~consumption~"("*µmol~O[2]~mg~DW^{-1}~h^{-1}~"\u00d7"~10^{-2}*")")
  ) +
  
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

print(p_resp)


# ============================================================
# 7. SAVE OUTPUTS
# ============================================================

write.csv(resp_summary,
          "outputs/results/respiration_rates.csv",
          row.names = FALSE)

ggsave("outputs/figures/resp_rates_by_pressure.png",
       p_resp, width = 5, height = 5, dpi = 300)
