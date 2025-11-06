# Statistical Inferences and Predictions for Areal Data and Spatial Data Fusion with Hausdorff–Gaussian Processes (JABE-D-25-00205)

This repository contains the code and data to reproduce the results in the paper
to be published in the Journal of Agricultural, Biological, and Environmental
Statistics (JABES)

---

## Repository Structure

* `/code`: Contains all `R` scripts to run the analyses, including helper
  functions in `/code/utils`.

* `/data`: Contains raw and processed data for the PM2.5 application.

* `/src`: Contains `C++` source code (`haus-dist.cpp`) used for efficient
  computation of the Hausdorff distance on the sphere (used only in the
  supplementary material).

* `/stan_models`: Contains all `Stan` models for Bayesian inference, including
  modular utility functions.

* `run-cv.sh`: A shell script for executing the cross-validation analysis.

> **Note:** Scripts prefixed with `*-00-*` are primarily for data pre-processing
> and supplemental checks. They are provided for completeness, but the main
> analyses are contained in the `*-01-*` (and higher) scripts.

### General Scripts

* `all-psd-check.R`: `R` script to generate Figure 1 from the paper (Positive
  Semi-Definite checks).

### Areal Data Analysis

* `areal-01-analyses.R`: Main `R` script for running the areal data analyses
  (e.g., comparing BYM2 and DAGAR).

* `areal-01-nu.R`: R script for the sensitivity analysis on the $\nu$ parameter.

* `areal-00-new-par.R`: (Processing) Script for parameter setup.

* `areal-00-supplemental-psd.R`: (Processing) Script for supplemental PSD checks
  related to the areal models. In this script, we check the PSD property for
  different areal datasets.

### Spatial Data Fusion Analysis

* `fusion-01-prep-cv.R`: Prepares the data folds for cross-validation.

* `fusion-02-cv-hgp.R`: Runs the 10-fold cross-validation using the
  Hausdorff–Gaussian Process (HGP) model.

* `fusion-03-cv-inla.R`: Runs the 10-fold cross-validation using the INLA model
  for data fusion.

* `fusion-04-cv-summary.R`: Summarizes and compares the cross-validation results
  (e.g., generates tables/figures).

* `fusion-05-analyses.R`: Runs the final, main analyses for the full data fusion
  application.

* `fusion-00-preproc-sat.R`: (Processing) Pre-processes the raw satellite data.

* `fusion-00-preproc-stations.R`: (Processing) Pre-processes the raw monitoring
  station data.

* `fusion-00-preproc-join.R`: (Processing) Joins satellite and station data.

* `fusion-00-supplemental-psd.R`: (Processing) Script for supplemental PSD
  checks related to the fusion models.

### Utility Scripts (`/code/utils`)

* `aux_dagar.R`: Helper functions for the DAGAR model.

* `aux_icar.R`: Helper functions for the ICAR model.

* `cv-task-funs.R`: Helper functions for managing cross-validation tasks.

* `summary_inla.R`: Helper functions for summarizing R-INLA model outputs.

---

## Stan Models

The `/stan_models` directory contains the `Stan` code for all Bayesian models.

* `/stan_models/areal`:
    * `bym2.stan`: Implements the BYM2 model.
    * `dagar.stan`: Implements the DAGAR model (see acknowledgments).

* `hgp-gaussian.stan`: The main Hausdorff-Gaussian Process (HGP) model with a
  Gaussian likelihood.

* `hgp-lognormal.stan`: HGP model with a Lognormal likelihood.

* `hgp-poisson.stan`: HGP model with a Poisson likelihood.

* `gaus-pred.stan`, `ln-pred.stan`: Stan models for generating predictions.

* `gaussian-gq.stan`, `ln-gq.stan`: Stan models for generating quantities from
  Gaussian/Lognormal models. That is, these are generated to generate samples
  from the "fitted values" and compute the log-likelihood for each sample from
  the MCMC algorithm. These are useful for model comparison purposes.

* `/stan_models/utils`: Contains modular `Stan` functions for code re-use (e.g.,
  `gp_lpdf.stan`, `pexp_corr.stan`).

---

## Data

The data for the PM2.5 application is located in `/data/pm25`.

* `stations2010_2012.csv`: Ground-level PM2.5 monitoring station data.
* `satellite.tif`: Raster file~(.tif) of satellite-derived PM2.5 data.
* `counties.rds`: `.rds` file containing the geometry of Ventura and Los Angeles
  counties.
* `dataset.rds`: The main processed dataset used in the fusion analysis.
* `/cluster`: Contains processed `.rds` files used for the cross-validation,
  such as the adjacency matrix (`W.rds`), distance matrices (`hmat.rds`), and
  saved cross-validation results (`cv-hgp.rds`, `cv-inla.rds`).

---

## 🚀 How to Reproduce

Ensure you have a working installation of `R`, with (at least) the packages
`rstan`, `sf`, `spdep`, `dplyr`, and `ggplot2` installed. A `C++` compiler is
also necessary.

Notice that, we expect the code to be replicable but not necessarily
reproducible since `Stan` results depend on the OS and architecture where the
code is executed.

---

## Acknowledgments

* The `dagar.stan` model in `/stan_models/areal` is an adaptation of `Stan` code
  kindly provided by Dr. Abhi Datta.
  
* **arXiv:** https://arxiv.org/html/2208.07900v3#S5  
