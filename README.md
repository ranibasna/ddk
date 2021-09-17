
# DDK <a href='https://github.com/ranibasna/ddk/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Last-changedate](https://img.shields.io/badge/last%20change-2021--09--17-brightgreen.svg)](/commits/master)
[![R-CMD-check](https://github.com/ranibasna/ddk/workflows/R-CMD-check/badge.svg)](https://github.com/ranibasna/ddk/actions)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ranibasna/ddk?branch=master&svg=true)](https://ci.appveyor.com/project/ranibasna/ddk)
[![Codecov test
coverage](https://codecov.io/gh/ranibasna/ddk/branch/master/graph/badge.svg)](https://codecov.io/gh/ranibasna/ddk?branch=master)
<!-- badges: end -->

<!-- [![Codecov test coverage](https://codecov.io/gh/ranibasna/NumericalTransformation/branch/master/graph/badge.svg)](https://codecov.io/gh/ranibasna/NumericalTransformation?branch=master) -->

# Data driven orthogonal basis selection for functional data analysis

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->

This is an R implementation for the data driven orthogonal basis
selection for functional data analysis paper. For more details see [the
arxiv link](https://arxiv.org/pdf/2103.07453.pdf) for the paper, Machine
Learning Assisted Orthonormal Basis Selection for Functional Data
Analysis.

In this work, we propose machine learning style techniques for the
placement of the knots. The chosen knots are used to build orthogonal
splines basis functions *f*<sub>*k*</sub>(*t*), *k* = 1, ..., *K* that
are used in basis function expansion to convert the data from discrete
recorded data into a functional one. The method address the problem of
the choice of the initial functional basis selection for functional data
analysis.

# DKK algorithm

For a more detailed description see the paper.

**Input** : knot search interval, *θ* - validation hyperparameter,
𝒳<sub>*t**r**a**i**n*</sub> - training data, 𝒳<sub>*v**a**l**i**d*</sub>
- validation data - Find the first global knot *ξ* using equation on
𝒳<sub>*t**r**a**i**n*</sub>.

-   Add split to *I* at the location of *ξ* and add *ξ* to
    𝒦<sup>(0)</sup>, the initial set of knots.

-   Set *s* = 1.

**While** stopping condition (see the stopping condition in the paper)
on 𝒳<sub>*v**a**l**i**d*</sub> is satisfied

1.  Find the new optimal placement using equation ;
2.  Add split to *I* at the location of the new selected knot
    *ξ*<sub>*s*</sub>.
3.  *s* = *s* + 1

## Installation

You can install the development version from
[GitHub](https://github.com/ranibasna/ddk) with:

``` r
# install.packages("devtools")
devtools::install_github("ranibasna/ddk")
```

## Reproducing the results in the paper

For reproducing the results presented in the paper for the Wine data,
see [this
link](https://ranibasna.github.io/ddk/articles/DKK_and_Functional_analysis_on_wine_data.html)
