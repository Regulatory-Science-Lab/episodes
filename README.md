 <!-- badges: start -->
  [![R-CMD-check](https://github.com/Regulatory-Science-Lab/episodes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Regulatory-Science-Lab/episodes/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# episodes

`episodes` is an R package for estimating transition probabilities from prescription data. The package has specific functions catered to the Flatiron oncogenomic database but can be used on other cancer prescription datasets as well to estimate transition health states from the data. 


This R package provides utilities to generate:
- **Death tables**
- **Summary tables of state counts**
- **Kaplan-Meier survival curves**
- **Weibull shape and scale parameter estimates**

It is designed to integrate with the [`{targets}`](https://docs.ropensci.org/targets/) pipeline system for reproducible and modular workflows.

## 📦 Installation

You can install the development version of the package directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Regulatory-Science-Lab/episodes")
```

## 🚀 Usage


The episodes pipeline can be run in a single line with the `run_episodes_pipeline` function. It takes in a `tibble` with the following columns:

```r
run_episodes_pipeline(
  tibble::tibble(
    tumour = c("nsclc", "breast"),
    treatment = c("cisplatin|carboplatin", "capecitabine"),
    exact = c(FALSE, TRUE)
  ), dir = "YOUR DIRECTORY TO SAVE RESULTS"
)
```
where `exact` represents whether we want to search for a combination of drugs (`exact = FALSE`) or only the specific treatment given by itself (`exact = TRUE`).

This will:

- Prepare drug episode data

- Construct state transitions

- Estimate Weibull parameters per state

- Output results to a timestamped Excel file

For the Excel file to be saved, you MUST have a template Excel file that is compatible to be populated with the health state estimates saved in the `dir` (directory) that you specify (see `h:\PREDiCText\nirupama\weibull_estimates\01_public_parameters_updatedNT.xlsx`)
