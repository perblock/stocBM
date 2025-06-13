




# take all the partitions from the repeat_sample_likely_partition function
# for each partition create a matrix indicating whether
# each pair of nodes is in the same partition
# then sum over all these matrices and divide by the number of runs
# of the repeat_sample_likely_partition function
# return the final matrix and plot it as a heatmap
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




# First use the update_partition_node_swap function
# and then the update_partition_node_move function
# until no further improvement in fit is possible
# and return the partition and the fit
# print the number of iterations and allow to set a maximum number of iterations
# use the partition_square_matrix function so that
# the only arguments are mobility_table and max_iterations
# only writing the function not the whole code
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

