% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition_change_functions.R
\name{update_partition_node_move}
\alias{update_partition_node_move}
\title{update_partition_node_move}
\usage{
update_partition_node_move(mobility_table, partition)
}
\arguments{
\item{mobility_table}{A square matrix representing the mobility table.}

\item{partition}{A vector representing the current partition of the
square matrix into blocks.}
}
\value{
A list containing the updated partition and the fit of the log-linear model.
}
\description{
Update the partition of a square matrix by moving the block membership of one node.
This function iteratively checks all nodes and moves a node to a different block
that results in the best fit of the log-linear model that predicts cell counts
}
\examples{
# Create a square matrix with 40 rows and 40 columns
test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
# Partition the matrix into 2 blocks
init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
# Update the partition by moving a node
update_partition_node_move(test_mat, init_partition)
}
