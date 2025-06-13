###### Functions to estimate Blocks ######

#' average_partition_matrix
#'
#' Calculate the average partition matrix from multiple runs of the
#' `repeat_sample_likely_partition` function.
#' This function takes the results from `repeat_sample_likely_partition`
#' and computes a matrix indicating the average likelihood of each pair of nodes
#' being in the same partition.
#' @param repeat_result A list containing the results from `repeat_sample_likely_partition`.
#' @return A matrix representing the average partition likelihood between nodes.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Repeat sampling a likely partition of the matrix with 2 blocks for 10 runs
#' repeat_result <- repeat_sample_likely_partition(test_mat, n_blocks = 2, n_runs = 10, n_iter = 1000)
#' avg_matrix <- average_partition_matrix(repeat_result)
#' @export
average_partition_matrix <- function(repeat_result) {
  partitions <- repeat_result$partitions
  n_runs <- length(partitions)
  n_nodes <- length(partitions[[1]])

  # Initialize an empty matrix to accumulate the results
  avg_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)

  for (partition in partitions) {
    # Create a matrix indicating whether each pair of nodes is in the same partition
    same_partition_matrix <- outer(partition, partition, FUN = "==")
    avg_matrix <- avg_matrix + same_partition_matrix
  }

  # Average the matrix by dividing by the number of runs
  avg_matrix <- avg_matrix / n_runs

  # Plot the average partition matrix as a heatmap
  heatmap(avg_matrix, Rowv = NA, Colv = NA, scale = "none",
          main = "Average Partition Matrix Heatmap")

  return(avg_matrix)
}


#' repeat_until_no_improvement
#'
#' Repeatedly update the partition of a square matrix
#' until no further improvement in fit is possible.
#' This function uses the `update_partition_node_swap` and
#' `update_partition_node_move` functions to
#' improve the partition iteratively.
#' @param mobility_table A square matrix representing the mobility table.
#' @param init_partition An optional initial partition vector.
#' @param n_blocks The number of blocks to partition the matrix into.
#' @param max_iter The maximum number of iterations to perform.
#' @return A list containing the final partition and the fit of the model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Update partition and repeat until no improvement in fit is possible
#' repeat_until_no_improvement(test_mat, n_blocks = 2, max_iter = 100)
#' @export
repeat_until_no_improvement <- function(mobility_table, init_partition = NULL,
                                        n_blocks = 2, max_iter = 100) {
  if(is.null(init_partition)){
    # If no initial partition is provided, create a random one
    partition <- partition_square_matrix(mobility_table, n_blocks = n_blocks)
  } else {
    # Use the provided initial partition
    partition <- init_partition
  }
  fit <- calculate_log_linear_model(mobility_table, partition)

  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1

    # Print the current iteration number and fit
    cat("Iteration:", iter, "Fit (AIC):", fit$fit, "\n")

    # Update partition by moving nodes
    result_move <- update_partition_node_move(mobility_table, partition)
    new_partition <- result_move$partition
    new_fit <- result_move$fit

    # If the fit improved, update the partition and fit
    if (new_fit$fit < fit$fit) {
      partition <- new_partition
      fit <- new_fit
      next
    }

    # Update partition by swapping nodes
    result_swap <- update_partition_node_swap(mobility_table, partition)
    new_partition <- result_swap$partition
    new_fit <- result_swap$fit

    # If the fit improved, update the partition and fit
    if (new_fit$fit < fit$fit) {
      partition <- new_partition
      fit <- new_fit
      next
    }

    # If no improvement was found, break the loop
    break
  }

  return(list(partition = partition, fit = fit))
}

