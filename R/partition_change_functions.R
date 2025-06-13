# update the calculated partition by swapping the block membership of two nodes
# so that the fit of the log linear model
# in which block membership predicts cell counts of the
# square matrix is maximized
# only writing the function not the whole code
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

      # If the new fit is better, update the best fit and partition
      if (new_fit$fit < best_fit$fit) {
        best_fit <- new_fit
        best_partition <- new_partition
      }
    }
  }

  return(list(partition = best_partition, fit = best_fit))
}

# update the calculated partition by changing the block membership of one nodes
# so that the fit of the log linear model
# in which block membership predicts cell counts of the
# square matrix is maximized
# do not allow moving a node so that one partition has fewer than 4 nodes
# only writing the function not the whole code
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

      # If the new fit is better, update the best fit and partition
      if (new_fit$fit < best_fit$fit) {
        best_fit <- new_fit
        best_partition <- new_partition
      }
    }
  }

  return(list(partition = best_partition, fit = best_fit))
}


# Propose an update of the calculated partition
# by swapping the block membership of two nodes
# Accept the update with a probability that is proportional to the
# exp() of the difference in fit of the log linear model of the
# current partition and the proposed partition.
# only writing the function not the whole code
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

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}


# Propose an update of the calculated partition
# by moving the block membership of one nodes
# Accept the update with a probability that is proportional to the
# exp() of the difference in fit of the log linear model of the
# current partition and the proposed partition.
# only writing the function not the whole code
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

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}


# Use either propose_partition_update_move or propose_partition_update_swap, randomly
# for a specified number of iterations
# and return the partition and the fit
# print the number of iterations
# only writing the function not the whole code
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
    cat("Iteration:", iter, "Fit (AIC):", fit$fit, "\n")

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


# Repeat the function sample_likely_partition for a specified number of runs
# the first time with a random partition
# and then with the partition from the previous run
# return the partition and the fit of each run as two lists
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

  return(list(partitions = partitions, fits = fits))
}
