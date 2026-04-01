# Independent Validation of Test-Adjusted COVID-19 Incidence Estimates 
# Using Wastewater Surveillance Data in Ontario, Canada

## Overview

This repository contains analysis code and data for:

Fisman DN, Tuite AR. Independent validation of test-adjusted COVID-19 
incidence estimates using wastewater surveillance data in Ontario, Canada. 
[Journal, Year, DOI to be added on publication]

## Description

We compared wastewater-based SARS-CoV-2 surveillance signals with crude 
reported and test-adjusted COVID-19 case counts across 111 weeks in Ontario, 
Canada (July 2020 – August 2022). Test adjustment methodology is described 
in Bosco et al. (BMC Infectious Diseases, 2025; https://doi.org/10.1186/s12879-025-10968-6).

## Repository Contents

- `wws_analysis_clean.R` — Main analysis script (Spearman correlations, 
  linear regression, negative binomial DLNM, stratified analyses)
- `figure1.R` — Time series and adjustment ratio figure
- `figure2.R` — Scatter plot regression figures  
- `figure3.R` — Lagged Spearman correlations and DLNM lag-response figure
- `merged_wastewater_cases.csv` — Analysis dataset

## Data

Wastewater data were obtained from provincial data files made available 
by the Ontario COVID-19 Science Table. Case and testing data were obtained 
from Ontario's Case and Contact Management System (CCM).

## Requirements

R version 4.5.0 or later. Required packages:

install.packages(c("tidyverse", "dlnm", "MASS", "tseries", 
                   "lmtest", "vars", "patchwork", "ggplot2",
                   "plotly", "htmlwidgets", "viridis"))

## Usage

Run wws_analysis_clean.R first to fit all models, then run 
figure1.R, figure2.R, and figure3.R to generate figures.

## License

MIT License

## Contact

David Fisman: david.fisman@utoronto.ca
