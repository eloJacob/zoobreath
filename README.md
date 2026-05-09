# zoobreath
Analysis pipeline to quantify oxygen consumption rates of Neocalanus gracilis under different hydrostatic pressure conditions (0.1, 5, 10 MPa).
This repository contains scripts to process raw Pyroscience oxygen data, apply calibration and corrections, compute respiration rates, and generate publication-ready figures.

## Repository structure
├── data/
│   ├── raw/                  # Raw Pyroscience files + length data
│   └── corrected/            # Cleaned datasets (generated)
│
├── outputs/
│   ├── data_corrected/       # Per-experiment cleaned data
│   ├── figures/              # Figures (exploration + final)
│   └── results/              # Final respiration rates table
│
├── R/
│   ├── 01_time_vector.R
│   ├── 02_clean_pyroscience_data.R
│   ├── 03_calibration.R
│   ├── 04_process_file.R
│   ├── 05_remove_first_hours.R
│   ├── 06_get_slope.R
│   ├── 07_correct_blank_drift.R
│   ├── 08_df_to_long.R
│   ├── 09_plot_experiment.R
│   └── 10_assemble_experiment.R
│
├── analysis/
│   ├── 01_Clean_importe_visualize.R
│   ├── 02_Merge_and_explore_oxygen_time_series.R
│   ├── 3_respiration_rates_and_biomass_normalisation.R
│
└── README.md

## Workflow overview
### 1. Import, cleaning & calibration
Script: 01_Clean_importe_visualize.R

Load raw Pyroscience files with phase shifts mesurement across time and oxygen concentration deduced per Chanel activated
Clean the pyroscience data and asigned chanel number to its id. 
Apply oxygen calibration coefficients from a 35-point calibration decribed in Chirugien et al. (2025)
Vizualize oxygen concentration across times (in hours) per chanel and temperature
Remove the first hour
Correct zooplankton chamber for blank drift
Export cleaned datasets per experiment

Output: 
outputs/data_corrected/*.csv
Raw & corrected diagnostic plots

### 2. Merge & exploratory analysis
Script: 02_Merge_and_explore_oxygen_time_series.R

Merge all experiments into a single dataset
Assign pressure conditions (0.1, 5, 10 MPa)
Filter zooplankton chambers (i.e. without controls chambers)
Visualize O₂ time series:
By pressure
By pressure × experiment

Output: 
all_experiments.csv
Exploratory figures

### 3. Respiration rates calculation
Script: 03_respiration_rates_and_biomass_normalisation.R

Calculate individual slopes using a linear model: O2 ~ Time
Extract slope (µmol O₂ L⁻¹ h⁻¹ ind⁻¹)
Estimate Dry weight from total length copepods using: 
DW = 0.01841 × Length²·⁴⁵⁷(Source: Ikeda et al., 2017)
Mass-normalised respiration rate(µmol O₂ mg DW⁻¹ h⁻¹)
Kruskal–Wallis test across pressure treatments

Output:
respiration_rates.csv
Final figure: respiration rate vs pressure


## Experimental design

Species: Neocalanus gracilis
Conditions:
0.1 MPa (ATM control)
5 MPa (2 experiments)
or 10 MPa (3 experiments)

Multiple independent experiments  (1 blank at 0.1 MPa + 1 under pressure (5 or 10 MPa) + 2 copepods at 0.1 Mpa + 2 copepods at 10 or 5 MPa)
Individuals tracked separately (ID-based analysis)

## Dependencies
R packages used:
dplyr
ggplot2
patchwork
stringr
presens
rstatix


## Notes
Blank correction differs depending on experiment (ATM vs HP blank availability)
Some channels are excluded due to artefacts or mortality

## Author
Élodie M.A. Jacob
2026
