# Code for Online Updating Proportional Hazards Test

This repository contains `R` code for *Modified Kaplan-Meier estimator and Nelson-Aalen
estimator with geographical weighting for survival data*.

1. `Hest.R` implements geographically weighted N-A estimator.

2. `plotband.R` implements plot of estimator of `Hest.R`.

3. `HT.R` implements the local test.

4. `covtemp.R` implements the global test.

5. `mod_Hcov.R` implements covariance estimator. 

6. `simulated_survdata.csv` contains a total of 1,000,000 simulated observations
with three covariates, `Age`, `Sex` and `Black`. The average censoring rate
is 50.2%.
