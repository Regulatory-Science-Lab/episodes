# episodes

`episodes` is an R package for estimating transition probabilities from prescription data. The package has specific functions catered to the Flatiron oncogenomic database but can be used on other cancer prescription datasets as well to estimate transition health states from the data. 

The package can be installed from github with the following command. 

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
devtools::install_github("Regulatory-Sciene-Lab/episodes")
```

This package is intended to be used with a _targets.R workflow. Once your _targets.R file is set up, run the pipeline with:

```r
targets::tar_make()
```

This will:

Prepare drug episode data

Construct state transitions

Estimate Weibull parameters per state

Output results to a timestamped Excel file

To modify the tumour types and treatments being modeled, edit the tumour_list object in your `_targets.R` file:

```r
tibble::tibble(
  tumour = c("nsclc", "breast"),
  treatment = c("cisplatin|carboplatin", "capecitabine"),
  exact = c(FALSE, TRUE)
)
```