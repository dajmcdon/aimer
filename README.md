# aimer: Amplified Initially Marginal Eigenvector Regression

The aimer package implements the aimer algorithm as described in [Ding and McDonald 2017]( https://doi.org/10.1093/bioinformatics/btx265).
The main goal is to use marginal regression to select a subset of coefficients, use
matrix approximation to estimate the principal components, and then extend those estimates
to the original predictor space before thresholding.

The package uses fast C++ routines to quickly perform matrix computations and cross validation for
selection of tuning parameters.

The package provides functions for estimating the model, choosing tuning parameters, and
generating simulated data as well as methods for prediction, plotting, and extraction.

## Installation

The easiest way to install is to write
```
devtools::install_github('dajmcdon/aimer')
```
at the `R` command prompt. If you want to also access the vignette (takes some time to build ~1 minute), you can use
```
devtools::install_github('dajmcdon/aimer', build_vignettes = TRUE)
```

Alternatively, you can clone the repository.

