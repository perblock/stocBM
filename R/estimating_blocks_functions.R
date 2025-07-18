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


#' multipartite_repeat_until_no_improvement
#'
#' Repeatedly update the multipartite partition of three square matrices
#' until no further improvement in fit is possible.
#' This function uses the `update_multipartite_partition_node_swap` and
#' `update_multipartite_partition_node_move` functions to
#' improve the partition iteratively.
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned
#' into the same blocks (b).
#' @param init_partition An optional initial partition list.
#' @param n_blocks A vector of two numbers representing the number of
#' blocks for (a) and (b).
#' @param max_iter The maximum number of iterations to perform.
#' @return A list containing the final partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Update multipartite partition and repeat until no improvement in fit is possible
#' multipartite_repeat_until_no_improvement(mat_aa, mat_ab, mat_bb,
#' n_blocks = c(2, 3), max_iter = 2)
#' @export
multipartite_repeat_until_no_improvement <- function(table_aa, table_ab, table_bb,
                                                     init_partition = NULL,
                                                     n_blocks = c(2, 2), max_iter = 100) {
  if(is.null(init_partition)){
    # If no initial partition is provided, create a random one
    partition <- partition_square_matrix(table_ab, n_blocks = n_blocks)
  } else {
    # Use the provided initial partition
    partition <- init_partition
  }

  fit_aa <- calculate_log_linear_model(table_aa, partition$rows)
  fit_ab <- calculate_log_linear_model(table_ab, partition)
  fit_bb <- calculate_log_linear_model(table_bb, partition$cols)

  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1

    # Print the current iteration number and fit
    cat("Iteration:", iter, "Fit (AIC):", fit_aa$fit + fit_ab$fit + fit_bb$fit, "\n")

    # Update multipartite partition by moving nodes
    result_move <- update_multipartite_partition_node_move(table_aa, table_ab, table_bb,
                                                           partition, swap_a = TRUE)
    new_partition <- result_move$partition
    new_fit_aa <- calculate_log_linear_model(table_aa, new_partition$rows)
    new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
    new_fit_bb <- calculate_log_linear_model(table_bb, new_partition$cols)

    # If the fit improved, update the partition and fit
    if (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit <
        fit_aa$fit + fit_ab$fit + fit_bb$fit) {
      partition <- new_partition
      fit_aa <- new_fit_aa
      fit_ab <- new_fit_ab
      fit_bb <- new_fit_bb
      next
    }

    result_move <- update_multipartite_partition_node_move(table_aa, table_ab, table_bb,
                                                           partition, swap_a = FALSE)
    new_partition <- result_move$partition
    new_fit_aa <- calculate_log_linear_model(table_aa, new_partition$rows)
    new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
    new_fit_bb <- calculate_log_linear_model(table_bb, new_partition$cols)

    # If the fit improved, update the partition and fit
    if (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit <
        fit_aa$fit + fit_ab$fit + fit_bb$fit) {
      partition <- new_partition
      fit_aa <- new_fit_aa
      fit_ab <- new_fit_ab
      fit_bb <- new_fit_bb
      next
    }

    # Update multipartite partition by swapping nodes
    result_swap <- update_multipartite_partition_node_swap(table_aa, table_ab, table_bb,
                                                           partition, swap_a = TRUE)
    new_partition <- result_swap$partition
    new_fit_aa <- calculate_log_linear_model(table_aa, new_partition$rows)
    new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
    new_fit_bb <- calculate_log_linear_model(table_bb, new_partition$cols)

    # If the fit improved, update the partition and fit
    if (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit <
        fit_aa$fit + fit_ab$fit + fit_bb$fit) {
      partition <- new_partition
      fit_aa <- new_fit_aa
      fit_ab <- new_fit_ab
      fit_bb <- new_fit_bb
      next
    }

    result_swap <- update_multipartite_partition_node_swap(table_aa, table_ab, table_bb,
                                                           partition, swap_a = FALSE)
    new_partition <- result_swap$partition
    new_fit_aa <- calculate_log_linear_model(table_aa, new_partition$rows)
    new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
    new_fit_bb <- calculate_log_linear_model(table_bb, new_partition$cols)

    # If the fit improved, update the partition and fit
    if (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit <
        fit_aa$fit + fit_ab$fit + fit_bb$fit) {
      partition <- new_partition
      fit_aa <- new_fit_aa
      fit_ab <- new_fit_ab
      fit_bb <- new_fit_bb
      next
    }
    # If no improvement was found, break the loop
    break
  }
  return(list(partition = partition,
              fit = fit_aa$fit + fit_ab$fit + fit_bb$fit))
}



#' get_good_partitions
#'
#' Get good partitions from the results of `repeat_sample_likely_partition`.
#' This function calculates the log-linear model for each partition,
#' repeats the partition update until no improvement is found,
#' and returns the best partition and its fit.
#' @param repeat_result A list containing the results from `repeat_sample_likely_partition`.
#' @param max_iter The maximum number of iterations to perform for each partition.
#' @return A list containing the number of unique partitions, all fits, the best partition,
#' and the best fit.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' test_mat[1:15, 1:15] <- test_mat[1:15, 1:15] + 2  # Add some structure
#' test_mat[25:40, 16:24] <- test_mat[25:40, 16:24] + 5  # Add some structure
#' # Repeat sampling a likely partition of the matrix with 3 blocks for 10 runs
#' repeat_result <- repeat_sample_likely_partition(test_mat, n_blocks = 3, n_runs = 10, n_iter = 1000)
#' # Get good partitions from the repeat result
#' get_good_partitions(repeat_result, max_iter = 100)
#' @export
get_good_partitions <- function(repeat_result,
                                max_iter = 100) {
  n_blocks <- length(unique(repeat_result$partitions[[1]]))
  partitions <- repeat_result$partitions
  unique_partitions <- list()
  fits <- list()

  for (i in seq_along(partitions)) {
    partition <- partitions[[i]]
    fit <- calculate_log_linear_model(repeat_result$mobility_table,
                                      partition)

    result <- repeat_until_no_improvement(repeat_result$mobility_table,
                                          partition, n_blocks, max_iter)

    unique_partitions[[i]] <- result$partition
    fits[[i]] <- result$fit
  }

  unique_partitions <- unique(unique_partitions)
  n_unique_partitions <- length(unique_partitions)

  best_partition_index <- which.min(sapply(fits, function(x) x$fit))
  best_partition <- unique_partitions[[best_partition_index]]
  best_fit <- fits[[best_partition_index]]

  return(list(n_unique_partitions = n_unique_partitions,
              all_fits = unlist(fits),
              best_partition = best_partition,
              best_fit = best_fit))
}


#' get_good_multipartite_partitions
#'
#' Get good multipartite partitions from the results of `repeat_sample_likely_multipartite_partition`.
#' This function calculates the log-linear model for each partition,
#' repeats the multipartite partition update until no improvement is found,
#' and returns the best partition and its fit.
#' @param repeat_result A list containing the results from `repeat_sample_likely_multipartite_partition`.
#' @param max_iter The maximum number of iterations to perform for each partition.
#' @return A list containing the number of unique partitions, all fits, the best partition,
#' and the best fit.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Repeat sampling a likely multipartite partition of the matrices with 2 blocks
#' # for (a) and 3 blocks for (b) for 10 runs
#' repeat_result <- repeat_sample_likely_multipartite_partition(mat_aa, mat_ab, mat_bb,
#' n_blocks = c(2, 3), n_runs = 10, n_iter = 1000)
#' # Get good multipartite partitions from the repeat result
#' get_good_multipartite_partitions(repeat_result, max_iter = 100)
#' @export
get_good_multipartite_partitions <- function(repeat_result,
                                             max_iter = 100) {
  n_blocks <- length(unique(repeat_result$partitions[[1]]$rows))
  partitions <- repeat_result$partitions
  unique_partitions <- list()
  fits <- list()

  for (i in seq_along(partitions)) {
    partition <- partitions[[i]]
    fit_aa <- calculate_log_linear_model(repeat_result$table_aa, partition$rows)
    fit_ab <- calculate_log_linear_model(repeat_result$table_ab, partition)
    fit_bb <- calculate_log_linear_model(repeat_result$table_bb, partition$cols)

    result <- multipartite_repeat_until_no_improvement(repeat_result$table_aa,
                                                       repeat_result$table_ab,
                                                       repeat_result$table_bb,
                                                       partition, n_blocks, max_iter)

    unique_partitions[[i]] <- result$partition
    fits[[i]] <- result$fit
  }

  unique_partitions <- unique(unique_partitions)
  n_unique_partitions <- length(unique_partitions)

  best_partition_index <- which.min(sapply(fits, function(x) x))
  best_partition <- unique_partitions[[best_partition_index]]
  best_fit <- fits[[best_partition_index]]

  return(list(n_unique_partitions = n_unique_partitions,
              all_fits = unlist(fits),
              best_partition = best_partition,
              best_fit = best_fit))
}


