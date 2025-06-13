
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stocBM

A set of functions to get the partitions of a stochastic blockmodel for
mobility networks; that is, for a cross-classification matrix in which
rowsums differ between rows AND rowsums and columnsums differ. In other
words, a stochastic blockmodel for directed, weighted networks, which is
conditional on the in- and out-degree distribution.

## Installation

You can install the current version of stocBM from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("perblock/stocBM")
```

## Small example

First, we create a matrix of counts as an example and add some
structure:

``` r
test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_mat[1:15, 1:15] <- test_mat[1:15, 1:15] + 2  # Add some structure
test_mat[25:40, 16:24] <- test_mat[25:40, 16:24] + 5  # Add some structure

# Look at the first rows and columns
test_mat[1:10, 1:10]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]   14   13   13   10    9   12   10   15   15    14
#>  [2,]    7    5   12    9    9   12   13   15   14     7
#>  [3,]   12   11   10   11    8   14   15   11    9    10
#>  [4,]   16    9   20   16   17   17   15   12   14    12
#>  [5,]   17   14   14   14   10   13   12   13   16    15
#>  [6,]    7   10   11   11   17   16    7   16    9    12
#>  [7,]   13   12   14   17   13   13   12   15    8     8
#>  [8,]   19   15   16   13   13   15    9    9   11     6
#>  [9,]   12    8   11   12   13   12   14   10   11    14
#> [10,]   20   18   12   17   18    7    9   14   17    18
```
