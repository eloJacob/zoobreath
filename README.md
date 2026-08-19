# zoobreath
Analysis pipeline to quantify oxygen consumption rates of Neocalanus gracilis under different hydrostatic pressure conditions (0.1, 5, 10 MPa).
This repository contains scripts to process raw Pyroscience oxygen data, apply calibration and corrections, compute respiration rates.

## Repository structure
├── data/
│   ├── raw/                # Raw Pyroscience files + length data
│   └── corrected/          # Cleaned datasets (generated)
│
├── outputs/
│   ├── data_corrected/     # Per-experiment cleaned data
│   ├── figures/            # Figures (exploration + final)
│   └── results/            # Final respiration rates table
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
│   ├── 10_assemble_experiment.R
│   ├── 12_utils_spline.R
│   ├── 13_slope_inference.R
│   ├── 14_merge_segments.R
│   ├── 15_select_o2_segments.R
│   ├── 16_analyze_oxy.R
│   └── 17_make_publication_plots.R
│
├── analysis/
│   ├── 01_Clean_importe_visualize.R
│   ├── 02_Merge_and_explore_oxygen_time_series.R
│   ├── 03_respiration_rates_and_biomass_normalisation.R
│   └── README.md

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

This script estimates oxygen consumption rates from respirometry time series 
using a **robust penalised spline segmentation approach**, rather than a single 
linear regression over the full recording. This method is more appropriate for 
noisy biological signals that may include artefacts, behavioural transitions, 
or instrumental drift.


## Method overview

1. **Penalised LAD spline fitting**
   Each time series (`zoo_all`) is fitted with a linear (degree 1) B-spline 
   basis (`splines2::bSpline()`), with `df = 10 × duration (h)` degrees of 
   freedom, and a second-order difference penalty (`λ = 0.008`) on the 
   coefficients. Coefficients are estimated by **L1 (LAD) loss minimisation** 
   using `CVXR::solve()` (CLARABEL solver), making the fit robust to outliers 
   and abrupt jumps.

2. **Segmentation & slope inference**
   Initial segment boundaries correspond to the internal knots of the spline 
   basis. Segments with fewer than `MIN_N = 20` points are merged with their 
   most similar neighbour (`repair_short_segment()`). For each segment, a 
   slope (µmol O₂ L⁻¹ h⁻¹) and its variance are estimated via a linear 
   contrast of the spline coefficients, with uncertainty obtained through 
   **500 residual bootstrap replicates** (median-centred residuals, preserving 
   the LAD zero-median assumption).

3. **Iterative segment merging**
   Adjacent segments with statistically indistinguishable slopes are merged 
   iteratively using a pairwise Z-test on slope differences 
   (`merge_segments()`), controlled by `ALPHA_SIM = 0.4` — a relaxed threshold 
   chosen to avoid under-segmentation (see Figure S4).

4. **Classification & selection of representative segments**
   Each final segment's slope is tested against zero (`ALPHA_SLOPE = 0.05`) 
   and classified as:
   - `"dec"` — significantly decreasing
   - `"inc"` — significantly increasing
   - `"NS"`  — non-significant

   Among all decreasing segments, a **duration-weighted kernel density 
   estimate** (Sheather–Jones bandwidth) identifies the modal slope. Segments 
   statistically similar to this mode (same Z-test, p > 0.30) are 
   progressively added to form the final selected cluster.

5. **Final oxygen consumption rate**
   The volumetric oxygen consumption rate (V_O₂, µmol O₂ L⁻¹ h⁻¹) is computed 
   as the **duration-weighted mean slope** of all selected segments.
   
7. **Biomass normalization**
   Individual **dry weight (DW, mg)** is estimated from prosome length (L, mm) 
  using the allometric relationship for *Neocalanus* sp. (Yang et al., 2017; 
  doi:10.6620/ZS.2017.56-13):
  DW = 0.01841 × L^2.457
 The volumetric rate is converted into a **mass-specific respiration rate** 
(µmol O₂ mg DW⁻¹ h⁻¹):
  R = (V_O₂ × v_chamber) / DW
where `v_chamber = 0.005 L` is the respiration chamber volume.

#### Key parameters

| Parameter        | Value    | Description                                   |
|------------------|----------|------------------------------------------------|
| `DF`             | 40       | Degrees of freedom for B-spline basis          |
| `LAMBDA`         | 0.008    | Roughness penalty weight                       |
| `ALPHA_SLOPE`    | 0.05     | Significance threshold for slope classification|
| `ALPHA_SIM`      | 0.4      | Similarity threshold for segment merging       |
| `MIN_N`          | 20       | Minimum points required per segment            |
| `MIN_DURATION`   | 0        | Minimum segment duration                       |
| `CHAMBER_VOL_mL` | 5        | Respiration chamber volume (mL)                |
| `ALLO_A`         | 0.01841  | Allometric coefficient (dry weight)            |
| `ALLO_B`         | 2.457    | Allometric exponent (dry weight)               |


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

## Authors
Élodie M.A. Jacob
Mathilde Couteyen-Carpaye
2026
