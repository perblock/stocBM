###### Functions for changing partitions ######

#' update_partition_node_swap
#'
#' Update the partition of a square matrix by swapping the block membership of two nodes.
#' This function iteratively checks all pairs of nodes and swaps the block membership
#' of the pair that results in the best fit of the log-linear model that predicts cell counts
#' of the square matrix based on block membership.
#'
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the square matrix into blocks.
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
#' # Update the partition by swapping nodes
#' update_partition_node_swap(test_mat, init_partition)
#' @export
update_partition_node_swap <- function(mobility_table, partition){
  n <- nrow(mobility_table)
  n_blocks <- length(unique(partition))
  best_fit <- calculate_log_linear_model(mobility_table, partition)
  best_partition <- partition

  for (i in seq_len(n - 1)) {
    for (j in seq((i + 1), n)) {
      # Swap the partition of nodes i and j
      new_partition <- partition
      new_partition[c(i, j)] <- new_partition[c(j, i)]

      # Check if the new partition is valid
      if (length(unique(new_partition)) < n_blocks ||
          any(table(new_partition) < 4)) {
        next
      }

      # Calculate the fit of the new partition
      new_fit <- calculate_log_linear_model(mobility_table, new_partition)

      # if the fit of the new partition cannot be calculated, skip
      if (is.null(new_fit) || is.na(new_fit$fit)) {
        next
      }

      # If the new fit is better, update the best fit and partition
      if (new_fit$fit < best_fit$fit) {
        best_fit <- new_fit
        best_partition <- new_partition
      }
    }
  }

  return(list(partition = best_partition, fit = best_fit))
}



#' update_partition_node_move
#'
#' Update the partition of a square matrix by moving the block membership of one node.
#' This function iteratively checks all nodes and moves a node to a different block
#' that results in the best fit of the log-linear model that predicts cell counts
# of the square matrix based on block membership.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the
#' square matrix into blocks.
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
#' # Update the partition by moving a node
#' update_partition_node_move(test_mat, init_partition)
#' @export
update_partition_node_move <- function(mobility_table, partition) {
  n <- nrow(mobility_table)
  n_blocks <- length(unique(partition))
  best_fit <- calculate_log_linear_model(mobility_table, partition)
  best_partition <- partition

  for (i in seq_len(n)) {
    current_block <- partition[i]
    other_blocks <- setdiff(unique(partition), current_block)

    for (new_block in other_blocks) {
      # Move node i to the new block
      new_partition <- partition
      new_partition[i] <- new_block

      # Check if the new partition is valid
      if (length(unique(new_partition)) < n_blocks ||
          any(table(new_partition) < 4)) {
        next
      }

      # Calculate the fit of the new partition
      new_fit <- calculate_log_linear_model(mobility_table, new_partition)

      # if the fit of the new partition cannot be calculated, skip
      if (is.null(new_fit) || is.na(new_fit$fit)) {
        next
      }

      # skip step if update cannot be calculated

      # If the new fit is better, update the best fit and partition
      if (new_fit$fit < best_fit$fit) {
        best_fit <- new_fit
        best_partition <- new_partition
      }
    }
  }

  return(list(partition = best_partition, fit = best_fit))
}

#' propose_partition_update_swap
#'
#' Propose an update of the calculated partition
#' by swapping the block membership of two nodes.
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in fit of the log-linear model of the
#' current partition and the proposed partition.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the
#' square matrix into blocks.
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
#' # Propose an update of the partition by swapping nodes
#' propose_partition_update_swap(test_mat, init_partition)
#' @export
propose_partition_update_swap <- function(mobility_table, partition) {
  n <- nrow(mobility_table)
  current_fit <- calculate_log_linear_model(mobility_table, partition)

  # Randomly select two nodes to swap
  i <- sample(seq_len(n - 1), 1)
  j <- sample((i + 1):n, 1)

  # Swap the partition of nodes i and j
  new_partition <- partition
  new_partition[c(i, j)] <- new_partition[c(j, i)]

  # Calculate the fit of the new partition
  new_fit <- calculate_log_linear_model(mobility_table, new_partition)

  # if the fit of the new partition cannot be calculated, reject the update
  if (is.null(new_fit) || is.na(new_fit$fit)) {
    return(list(partition = partition, fit = current_fit))
  }

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}


#' propose_partition_update_move
#'
#' Propose an update of the calculated partition
#' by moving the block membership of one node.
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in fit of the log-linear model of the
#' current partition and the proposed partition.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the
#' square matrix into blocks.
#' @return A list containing the updated partition and the fit of the log-linear model.
#' #' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = 2)
propose_partition_update_move <- function(mobility_table, partition) {
  n <- nrow(mobility_table)
  current_fit <- calculate_log_linear_model(mobility_table, partition)

  # Randomly select a node to move
  i <- sample(seq_len(n), 1)
  current_block <- partition[i]
  other_blocks <- setdiff(unique(partition), current_block)

  # Randomly select a new block for the node
  new_block <- sample(other_blocks, 1)

  # Move node i to the new block
  new_partition <- partition
  new_partition[i] <- new_block

  # Calculate the fit of the new partition
  new_fit <- calculate_log_linear_model(mobility_table, new_partition)

  # if the fit of the new partition cannot be calculated, reject the update
  if (is.null(new_fit) || is.na(new_fit$fit)) {
    return(list(partition = partition, fit = current_fit))
  }

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}

#' sample_likely_partition
#'
#' Sample a likely partition of the square matrix
#' Use either propose_partition_update_move or propose_partition_update_swap,
#' randomly for a specified number of iterations
#' and return the partition and the fit
#' @param mobility_table A square matrix representing the mobility table.
#' @param initial_partition A vector representing the initial partition of the square matrix into blocks.
#' @param n_blocks The number of blocks to partition the matrix into.
#' @param n_iter The number of iterations to run the sampling.
#' @return A list containing the sampled partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Sample a likely partition of the matrix with 2 blocks
#' sample_likely_partition(test_mat, n_blocks = 2, n_iter = 1000)
#' @export
sample_likely_partition <- function(mobility_table, initial_partition = NULL,
                                    n_blocks = 2, n_iter = 1000) {
  if(is.null(initial_partition)) {
    # If no initial partition is provided, create a random one
    partition <- partition_square_matrix(mobility_table, n_blocks = n_blocks)
  } else {
    # Use the provided initial partition
    partition <- initial_partition
  }
  fit <- calculate_log_linear_model(mobility_table, partition)

  iter <- 0
  while (iter < n_iter) {
    iter <- iter + 1

    # Print the current iteration number and fit
    if (iter %% 100 == 0) {  # Print every 100 iterations
      cat("Iteration:", iter, "Fit (AIC):", fit$fit, "\n")
    }

    # Randomly choose to either move or swap nodes
    if (runif(1) < 0.5) {
      result_move <- propose_partition_update_move(mobility_table, partition)
    } else {
      result_move <- propose_partition_update_swap(mobility_table, partition)
    }

    partition <- result_move$partition
    fit <- result_move$fit

  }

  return(list(partition = partition, fit = fit))
}


#' repeat_sample_likely_partition
#'
#' Repeat the function sample_likely_partition for a specified number of runs
#' the first time with a random partition
#' and then with the partition from the previous run
#' return the partition and the fit of each run as two lists
#' @param mobility_table A square matrix representing the mobility table.
#' @param n_blocks The number of blocks to partition the matrix into.
#' @param n_runs The number of runs to perform.
#' @param n_iter The number of iterations to run the sampling in each run.
#' @return A list containing two lists (partitions and fits), where each list contains the
#' partition and fit from each run; as well as the mobility table.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Repeat sampling a likely partition of the matrix with 2 blocks for 10 runs
#' repeat_sample_likely_partition(test_mat, n_blocks = 2, n_runs = 10, n_iter = 1000)
#' @export
repeat_sample_likely_partition <- function(mobility_table, n_blocks = 2, n_runs = 10, n_iter = 100) {
  partitions <- vector("list", n_runs)
  fits <- vector("list", n_runs)

  for (run in seq_len(n_runs)) {
    cat("Run:", run, "\n")
    if (run == 1) {
      # First run with a random partition
      result <- sample_likely_partition(mobility_table, n_blocks = n_blocks, n_iter = 5 * n_iter)
    } else {
      # Subsequent runs with the partition from the previous run
      result <- sample_likely_partition(mobility_table,
                                        initial_partition = partitions[[run - 1]],
                                        n_blocks = n_blocks,
                                        n_iter = n_iter)
    }

    partitions[[run]] <- result$partition
    fits[[run]] <- result$fit
  }

  return(list(partitions = partitions, fits = fits,
              mobility_table = mobility_table))
}

