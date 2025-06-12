### some basic functions to calculate a stochastic blockmodel by
# (1) randomly partition a square matrix (mobility table)
# (2) calculate the parameters of a log-linear model of
#     this partition
# (3) update the partition by moving one node between blocks
#     or swapping the block membership of two nodes
#     so that the fit of the log linear model is maximized
# (4) repeat (3) until no further improvement is possible
#     and return the partition and the fit
# (5) build a sampler that samples likely partitions

# Partition a square matrix into blocks
# if the number of nodes does not divide into blocks,
# allow blocks of approximately same sizes
# only writing the function not the whole code
partition_square_matrix <- function(mobility_table, n_blocks) {
  n <- nrow(mobility_table)
  if (n %% n_blocks != 0) {
    block_size <- floor(n / n_blocks)
    remainder <- n %% n_blocks
    sizes <- rep(block_size, n_blocks)
    sizes[1:remainder] <- sizes[1:remainder] + 1
  } else {
    sizes <- rep(n / n_blocks, n_blocks)
  }

  partition <- rep(1:n_blocks, times = sizes)
  partition <- sample(partition)  # Shuffle the partition

  return(partition)
}



# Calculate the parameters of a log-linear model
# in which block membership predicts cell counts of the
# square matrix;
# additionally, add a parameter for each row and column of the matrix.
# then return the parameters for the log-linear model
# AND the fit of the model in terms of aic
# only writing the function not the whole code
calculate_log_linear_model <- function(mobility_table, partition, include_coeffs = FALSE) {
  n <- nrow(mobility_table)
  n_blocks <- length(unique(partition))

  if(include_coeffs) {
    # Create a data frame for the model
    df <- data.frame(
      count = as.vector(mobility_table),
      block_row = rep(partition, times = n),
      block_col = rep(partition, each = n),
      row_no = rep(1:n, times = n),
      col_no = rep(1:n, each = n)
    )

    # Fit the log-linear model
    model <- glm(count ~ 1 + factor(block_row) * factor(block_col) +
                   factor(row_no) + factor(col_no),
                 data = df, family = poisson())

    coeffs <- coef(model)
    # Calculate AIC
    fit_aic <- AIC(model)
  } else {
    LL <- 0
    for(r in 1:n_blocks){
      for(s in 1:n_blocks){
        m_rs <- sum(mobility_table[partition == r, partition == s])
        k_r <- sum(mobility_table[partition == r,])
        k_s <- sum(mobility_table[,partition == s])

        LL <- LL + m_rs * log( m_rs / (k_r * k_s) )
      }
    }
    fit_aic <- -2*LL
    coeffs <- NULL
  }

  # Return the model coefficients and fit
  return(list(coefficients = coeffs, fit = fit_aic))
}



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










