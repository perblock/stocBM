% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimating_blocks_functions.R
\name{multipartite_repeat_until_no_improvement}
\alias{multipartite_repeat_until_no_improvement}
\title{multipartite_repeat_until_no_improvement}
\usage{
multipartite_repeat_until_no_improvement(
  table_aa,
  table_ab,
  table_bb,
  init_partition = NULL,
  n_blocks = c(2, 2),
  max_iter = 100
)
}
\arguments{
\item{table_aa}{A square matrix representing the mobility table
in which rows and columns are partitioned into the same blocks (a).}

\item{table_ab}{A square matrix representing the mobility table
in which rows are partitioned into blocks (a) and columns are partitioned
into different blocks (b).}

\item{table_bb}{A square matrix representing the mobility table
in which rows and columns are partitioned
into the same blocks (b).}

\item{init_partition}{An optional initial partition list.}

\item{n_blocks}{A vector of two numbers representing the number of
blocks for (a) and (b).}

\item{max_iter}{The maximum number of iterations to perform.}
}
\value{
A list containing the final partitions (a) and (b)
and the summed fit of the log-linear model for the three matrices.
}
\description{
Repeatedly update the multipartite partition of three square matrices
until no further improvement in fit is possible.
This function uses the `update_multipartite_partition_node_swap` and
`update_multipartite_partition_node_move` functions to
improve the partition iteratively.
}
\examples{
# Create a square matrix with 40 rows and 40 columns
mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
# Create a square matrix with 40 rows and 40 columns
mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
# Create a square matrix with 40 rows and 40 columns
mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
# Update multipartite partition and repeat until no improvement in fit is possible
multipartite_repeat_until_no_improvement(mat_aa, mat_ab, mat_bb,
n_blocks = c(2, 3), max_iter = 2)
}
