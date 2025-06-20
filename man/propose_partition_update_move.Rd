% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition_change_functions.R
\name{propose_partition_update_move}
\alias{propose_partition_update_move}
\title{propose_partition_update_move}
\usage{
propose_partition_update_move(mobility_table, partition)
}
\arguments{
\item{mobility_table}{A square matrix representing the mobility table.}

\item{partition}{A vector representing the current partition of the
square matrix into blocks.}
}
\value{
A list containing the updated partition and the fit of the log-linear model.
#' @examples
# Create a square matrix with 40 rows and 40 columns
test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
# Partition the matrix into 2 blocks
init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
}
\description{
Propose an update of the calculated partition
by moving the block membership of one node.
Accept the update with a probability that is proportional to the
exp() of the difference in fit of the log-linear model of the
current partition and the proposed partition.
}
