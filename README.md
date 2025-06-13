
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
library(stocBM)

test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_mat[1:15, 1:15] <- test_mat[1:15, 1:15] + 2  # Add some structure
test_mat[25:40, 16:24] <- test_mat[25:40, 16:24] + 5  # Add some structure

# Look at the first rows and columns
test_mat[1:10, 1:10]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]    9   10   17   11   13    9    6    9    7    11
#>  [2,]   10   13    9   11    9    8   15   10   14    10
#>  [3,]    9   16   10   13   12    9   14   12    8    17
#>  [4,]    8    9   12    8    7   13   14    8   13     9
#>  [5,]   10   19   14    9   11   11   11   17    8     7
#>  [6,]   16   12   11    7   15   14   14    7   10    11
#>  [7,]   12    7   12   10   16    8   10   17   11     8
#>  [8,]   10    7   10   16   12   13   14   12   17    10
#>  [9,]   18   11   16   12   11    8   13   11   10    14
#> [10,]   11   17   12    9   13   11    9   13    8     8
```

Now we can run a simple first test to get a (probably) preliminary
partitioning. In case the partitions are very clear, this might already
be pretty good.

``` r
init_res <- repeat_until_no_improvement(test_mat, n_blocks = 3, max_iter = 1000)
#> Iteration: 1 Fit (AIC): 330975.3 
#> Iteration: 2 Fit (AIC): 330972.2 
#> Iteration: 3 Fit (AIC): 330967.8 
#> Iteration: 4 Fit (AIC): 330964.5 
#> Iteration: 5 Fit (AIC): 330962 
#> Iteration: 6 Fit (AIC): 330959.6 
#> Iteration: 7 Fit (AIC): 330956.2 
#> Iteration: 8 Fit (AIC): 330951.5 
#> Iteration: 9 Fit (AIC): 330947.5 
#> Iteration: 10 Fit (AIC): 330943.7 
#> Iteration: 11 Fit (AIC): 330939.3 
#> Iteration: 12 Fit (AIC): 330935.7 
#> Iteration: 13 Fit (AIC): 330931.6 
#> Iteration: 14 Fit (AIC): 330925.7 
#> Iteration: 15 Fit (AIC): 330922 
#> Iteration: 16 Fit (AIC): 330916.9 
#> Iteration: 17 Fit (AIC): 330912.1 
#> Iteration: 18 Fit (AIC): 330908.9 
#> Iteration: 19 Fit (AIC): 330905.8 
#> Iteration: 20 Fit (AIC): 330904.4 
#> Iteration: 21 Fit (AIC): 330903.6 
#> Iteration: 22 Fit (AIC): 330901.1 
#> Iteration: 23 Fit (AIC): 330899 
#> Iteration: 24 Fit (AIC): 330898.3 
#> Iteration: 25 Fit (AIC): 330896.6 
#> Iteration: 26 Fit (AIC): 330894.2 
#> Iteration: 27 Fit (AIC): 330889.6 
#> Iteration: 28 Fit (AIC): 330882.4 
#> Iteration: 29 Fit (AIC): 330875.9 
#> Iteration: 30 Fit (AIC): 330867.7 
#> Iteration: 31 Fit (AIC): 330861.7 
#> Iteration: 32 Fit (AIC): 330856.8 
#> Iteration: 33 Fit (AIC): 330850.9 
#> Iteration: 34 Fit (AIC): 330844.7 
#> Iteration: 35 Fit (AIC): 330840.4 
#> Iteration: 36 Fit (AIC): 330831.7 
#> Iteration: 37 Fit (AIC): 330818.6 
#> Iteration: 38 Fit (AIC): 330810.1 
#> Iteration: 39 Fit (AIC): 330793.1 
#> Iteration: 40 Fit (AIC): 330773.8 
#> Iteration: 41 Fit (AIC): 330765.8
init_res$partition
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1
```

Ok, that looks not ideal, there are some misclassified nodes. Next, we
do a more extensive sampling of different good partitions and a plot of
average same block membership:

``` r
repeat_result <- repeat_sample_likely_partition(test_mat, n_blocks = 3, n_runs = 10, n_iter = 1000)
#> Run: 1 
#> Iteration: 100 Fit (AIC): 330935.5 
#> Iteration: 200 Fit (AIC): 330886.8 
#> Iteration: 300 Fit (AIC): 330865.3 
#> Iteration: 400 Fit (AIC): 330803.5 
#> Iteration: 500 Fit (AIC): 330803.5 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Iteration: 1100 Fit (AIC): 330765.8 
#> Iteration: 1200 Fit (AIC): 330765.8 
#> Iteration: 1300 Fit (AIC): 330765.8 
#> Iteration: 1400 Fit (AIC): 330765.8 
#> Iteration: 1500 Fit (AIC): 330765.8 
#> Iteration: 1600 Fit (AIC): 330765.8 
#> Iteration: 1700 Fit (AIC): 330765.8 
#> Iteration: 1800 Fit (AIC): 330765.8 
#> Iteration: 1900 Fit (AIC): 330765.8 
#> Iteration: 2000 Fit (AIC): 330765.8 
#> Iteration: 2100 Fit (AIC): 330765.8 
#> Iteration: 2200 Fit (AIC): 330765.8 
#> Iteration: 2300 Fit (AIC): 330765.8 
#> Iteration: 2400 Fit (AIC): 330765.8 
#> Iteration: 2500 Fit (AIC): 330765.8 
#> Iteration: 2600 Fit (AIC): 330765.8 
#> Iteration: 2700 Fit (AIC): 330765.8 
#> Iteration: 2800 Fit (AIC): 330765.8 
#> Iteration: 2900 Fit (AIC): 330765.8 
#> Iteration: 3000 Fit (AIC): 330765.8 
#> Iteration: 3100 Fit (AIC): 330765.8 
#> Iteration: 3200 Fit (AIC): 330765.8 
#> Iteration: 3300 Fit (AIC): 330765.8 
#> Iteration: 3400 Fit (AIC): 330765.8 
#> Iteration: 3500 Fit (AIC): 330765.8 
#> Iteration: 3600 Fit (AIC): 330765.8 
#> Iteration: 3700 Fit (AIC): 330765.8 
#> Iteration: 3800 Fit (AIC): 330765.8 
#> Iteration: 3900 Fit (AIC): 330765.8 
#> Iteration: 4000 Fit (AIC): 330765.8 
#> Iteration: 4100 Fit (AIC): 330765.8 
#> Iteration: 4200 Fit (AIC): 330765.8 
#> Iteration: 4300 Fit (AIC): 330765.8 
#> Iteration: 4400 Fit (AIC): 330765.8 
#> Iteration: 4500 Fit (AIC): 330773.8 
#> Iteration: 4600 Fit (AIC): 330765.8 
#> Iteration: 4700 Fit (AIC): 330765.8 
#> Iteration: 4800 Fit (AIC): 330765.8 
#> Iteration: 4900 Fit (AIC): 330765.8 
#> Iteration: 5000 Fit (AIC): 330765.8 
#> Run: 2 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 3 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 4 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 5 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 6 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 7 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 8 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 9 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330765.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8 
#> Run: 10 
#> Iteration: 100 Fit (AIC): 330765.8 
#> Iteration: 200 Fit (AIC): 330765.8 
#> Iteration: 300 Fit (AIC): 330765.8 
#> Iteration: 400 Fit (AIC): 330765.8 
#> Iteration: 500 Fit (AIC): 330765.8 
#> Iteration: 600 Fit (AIC): 330765.8 
#> Iteration: 700 Fit (AIC): 330773.8 
#> Iteration: 800 Fit (AIC): 330765.8 
#> Iteration: 900 Fit (AIC): 330765.8 
#> Iteration: 1000 Fit (AIC): 330765.8
repeat_result$partitions
#> [[1]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[2]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[3]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[4]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[5]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[6]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[7]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[8]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[9]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> [[10]]
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2

average_partition_matrix(repeat_result)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

    #>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    #>  [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [4,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [5,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [6,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [7,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [8,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #>  [9,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [10,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [11,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [12,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [13,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [14,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [15,]    1    1    1    1    1    1    1    1    1     1     1     1     1
    #> [16,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [18,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [19,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [21,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [22,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [23,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [24,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [25,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [26,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [27,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [28,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [29,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [30,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [31,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [32,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [33,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [34,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [35,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [36,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [37,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [38,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [39,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #> [40,]    0    0    0    0    0    0    0    0    0     0     0     0     0
    #>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
    #>  [1,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [2,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [3,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [4,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [5,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [6,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [7,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [8,]     1     1     0     0     0     0     0     0     0     0     0     0
    #>  [9,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [10,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [11,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [12,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [13,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [14,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [15,]     1     1     0     0     0     0     0     0     0     0     0     0
    #> [16,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [17,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [18,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [19,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [20,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [21,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [22,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [23,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [24,]     0     0     1     1     1     1     1     1     1     1     1     0
    #> [25,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [26,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [27,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [28,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [29,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [30,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [31,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [32,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [33,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [34,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [35,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [36,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [37,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [38,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [39,]     0     0     0     0     0     0     0     0     0     0     0     1
    #> [40,]     0     0     0     0     0     0     0     0     0     0     0     1
    #>       [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
    #>  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
    #>  [9,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [11,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [12,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [13,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [14,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [15,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [16,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [17,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [18,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [19,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [20,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [21,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [22,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [23,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [24,]     0     0     0     0     0     0     0     0     0     0     0     0
    #> [25,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [26,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [27,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [28,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [29,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [30,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [31,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [32,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [33,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [34,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [35,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [36,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [37,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [38,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [39,]     1     1     1     1     1     1     1     1     1     1     1     1
    #> [40,]     1     1     1     1     1     1     1     1     1     1     1     1
    #>       [,38] [,39] [,40]
    #>  [1,]     0     0     0
    #>  [2,]     0     0     0
    #>  [3,]     0     0     0
    #>  [4,]     0     0     0
    #>  [5,]     0     0     0
    #>  [6,]     0     0     0
    #>  [7,]     0     0     0
    #>  [8,]     0     0     0
    #>  [9,]     0     0     0
    #> [10,]     0     0     0
    #> [11,]     0     0     0
    #> [12,]     0     0     0
    #> [13,]     0     0     0
    #> [14,]     0     0     0
    #> [15,]     0     0     0
    #> [16,]     0     0     0
    #> [17,]     0     0     0
    #> [18,]     0     0     0
    #> [19,]     0     0     0
    #> [20,]     0     0     0
    #> [21,]     0     0     0
    #> [22,]     0     0     0
    #> [23,]     0     0     0
    #> [24,]     0     0     0
    #> [25,]     1     1     1
    #> [26,]     1     1     1
    #> [27,]     1     1     1
    #> [28,]     1     1     1
    #> [29,]     1     1     1
    #> [30,]     1     1     1
    #> [31,]     1     1     1
    #> [32,]     1     1     1
    #> [33,]     1     1     1
    #> [34,]     1     1     1
    #> [35,]     1     1     1
    #> [36,]     1     1     1
    #> [37,]     1     1     1
    #> [38,]     1     1     1
    #> [39,]     1     1     1
    #> [40,]     1     1     1

Now that looks better. We can also use this output to optimise the
different sampled partitions to find the best (from our current
partitioning).

``` r
get_good_partitions(repeat_result)
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8 
#> Iteration: 1 Fit (AIC): 330765.8
#> $n_unique_partitions
#> [1] 1
#> 
#> $all_fits
#>      fit      fit      fit      fit      fit      fit      fit      fit 
#> 330765.8 330765.8 330765.8 330765.8 330765.8 330765.8 330765.8 330765.8 
#>      fit      fit 
#> 330765.8 330765.8 
#> 
#> $best_partition
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
#> 
#> $best_fit
#> $best_fit$coefficients
#> NULL
#> 
#> $best_fit$fit
#> [1] 330765.8
```

This tells us how many different solutions are found (fewer = better)
and what the best partition is.
