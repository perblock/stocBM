###### Functions for changing partitions ######

#' update_partition_node_swap
#'
#' Update the partition of a square matrix by swapping the block membership of two nodes.
#' This function iteratively checks all pairs of nodes and swaps the block membership
#' of the pair that results in the best fit of the log-linear model that predicts cell counts
#' of the square matrix based on block membership.
#' In case the rows and columns of the square matrix are partitioned into different blocks,
#' this function will either swap the block membership of rows or columns,
#' depending on a parameter or randomly.
#'
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the square matrix into blocks.
#' if the rows and columns are partitioned into different blocks,
#' the partition should be a list with two vectors: one for the rows and one for the columns.
#' @param swap_rows A boolean indicating whether to swap the block membership of rows (TRUE)
#' or columns (FALSE).
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = c(2,3))
#' # Update the partition by swapping nodes
#' update_partition_node_swap(test_mat, init_partition, swap_rows = TRUE)
#' @export
update_partition_node_swap <- function(mobility_table, partition, swap_rows = NULL) {
  n <- nrow(mobility_table)
  swap_rows <- ifelse(is.null(swap_rows), sample(c(TRUE, FALSE), 1), swap_rows)

  if (is.list(partition)) {
    partition_rows <- partition$rows
    n_blocks_rows <- length(unique(partition$rows))
    partition_cols <- partition$cols
    n_blocks_cols <- length(unique(partition$cols))
  } else {
    n_blocks <- length(unique(partition))
  }

  best_fit <- calculate_log_linear_model(mobility_table, partition)
  best_partition <- partition

  for (i in seq_len(n - 1)) {
    for (j in seq((i + 1), n)) {
      # Swap the block membership of nodes i and j
      if (is.list(partition)){
        if (swap_rows) {
          new_partition <- partition
          new_partition$rows[c(i, j)] <- new_partition$rows[c(j, i)]
        } else {
          new_partition <- partition
          new_partition$cols[c(i, j)] <- new_partition$cols[c(j, i)]
        }
      } else {
        new_partition <- partition
        new_partition[c(i, j)] <- new_partition[c(j, i)]
      }

      # # Check if the new partition is valid
      # if (length(unique(new_partition$rows)) < n_blocks ||
      #     length(unique(new_partition$cols)) < n_blocks ||
      #     any(table(new_partition$rows) < 4) ||
      #     any(table(new_partition$cols) < 4)) {
      #   next
      # }

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


#' update_multipartite_partition_node_swap
#'
#' Update one of two partitions that determine cell counts
#' of three square matrices by swapping the block membership of two nodes.
#' The partitions guide three matrices:
#' 1. A sqare matrix in which rows and columns are partitioned
#' into the same blocks (a),
#' 2. A square matrix in which rows are partitioned into blocks (a)
#' and columns are partitioned into different blocks (b),
#' 3. A square matrix in which rows and columns are partitioned
#' into the same blocks (b).
#' This function iteratively checks all pairs of nodes and swaps the block membership
#' of the pair that results in the best fit of the
#' log-linear model that predicts cell counts.
#' In each updating step, this function will either swap the
#' block membership of blocks (a) or blocks (b),
#' depending on a parameter or randomly.
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (b).
#' @param partition A list containing two vectors:
#' one for the partitioning of (a)
#' and one for the partitioning of (b).
#' @param swap_a A boolean indicating whether to swap the block membership of blocks (a) (TRUE)
#' or blocks (b) (FALSE). Default is NULL, which means it will be chosen randomly.
#' @return A list containing the updated partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrices into 2 blocks for (a) and 3 blocks for (b)
#' init_partition <- partition_square_matrix(mat_ab, n_blocks = c(2, 3))
#' # Update the partition by swapping nodes
#' update_multipartite_partition_node_swap(mat_aa, mat_ab, mat_bb,
#'                                         init_partition, swap_a = TRUE)
#' @export
update_multipartite_partition_node_swap <- function(table_aa, table_ab, table_bb,
                                                    partition, swap_a = NULL) {
  n <- nrow(table_aa)
  swap_a <- ifelse(is.null(swap_a), sample(c(TRUE, FALSE), 1), swap_a)

  if (is.list(partition)) {
    partition_a <- partition$rows
    n_blocks_a <- length(unique(partition$rows))
    partition_b <- partition$cols
    n_blocks_b <- length(unique(partition$cols))
  } else {
    stop("partition should be a list with two vectors:
         one for the rows/a and one for the columns/b.")
  }

  best_fit_aa <- calculate_log_linear_model(table_aa, partition_a)
  best_fit_ab <- calculate_log_linear_model(table_ab, partition)
  best_fit_bb <- calculate_log_linear_model(table_bb, partition_b)

  best_partition <- partition

  for (i in seq_len(n - 1)) {
    for (j in seq((i + 1), n)) {
      # Swap the block membership of nodes i and j
      if (swap_a) {
        new_partition_a <- partition_a
        new_partition_a[c(i, j)] <- new_partition_a[c(j, i)]
        new_partition_b <- partition_b
        new_partition <- list(rows = new_partition_a, cols = new_partition_b)
      } else {
        new_partition_b <- partition_b
        new_partition_b[c(i, j)] <- new_partition_b[c(j, i)]
        new_partition_a <- partition_a
        new_partition <- list(rows = new_partition_a, cols = new_partition_b)
      }

      # Calculate the fit of the new partitions
      new_fit_aa <- calculate_log_linear_model(table_aa, new_partition_a)
      new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
      new_fit_bb <- calculate_log_linear_model(table_bb, new_partition_b)

      # if the fit of the new partitions cannot be calculated, skip
      if (is.null(new_fit_aa) || is.na(new_fit_aa$fit) ||
          is.null(new_fit_ab) || is.na(new_fit_ab$fit) ||
          is.null(new_fit_bb) || is.na(new_fit_bb$fit)) {
        next
      }

      # If the new fit is better, update the best fit and partition
      if ((new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit) <
          (best_fit_aa$fit + best_fit_ab$fit + best_fit_bb$fit)) {
        best_fit_aa <- new_fit_aa
        best_fit_ab <- new_fit_ab
        best_fit_bb <- new_fit_bb
        best_partition <- list(rows = new_partition_a, cols = new_partition_b)
      }
    }
  }
  return(list(partition = best_partition,
              fit = best_fit_aa$fit + best_fit_ab$fit + best_fit_bb$fit))
}




#' update_partition_node_move
#'
#' Update the partition of a square matrix by moving the block membership of one node.
#' This function iteratively checks all nodes and moves a node to a different block
#' that results in the best fit of the log-linear model that predicts cell counts
#' of the square matrix based on block membership.
#' In case the rows and columns of the square matrix are partitioned into different blocks,
#' this function will either move the block membership of one node in rows or columns,
#' depending on a parameter or randomly.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the current partition of the square matrix into blocks.
#' if the rows and columns are partitioned into different blocks,
#' the partition should be a list with two vectors: one for the rows and one for the columns.
#' @param move_rows A boolean indicating whether to swap the block membership of rows (TRUE)
#' or columns (FALSE). Default is NULL, which means it will be chosen randomly.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' init_partition <- partition_square_matrix(test_mat, n_blocks = c(2,3))
#' # Update the partition by moving a node
#' update_partition_node_move(test_mat, init_partition)
#' @export
update_partition_node_move <- function(mobility_table, partition, move_rows = NULL) {
  n <- nrow(mobility_table)
  move_rows <- ifelse(is.null(move_rows), sample(c(TRUE, FALSE), 1), move_rows)

  if (is.list(partition)) {
    partition_rows <- partition$rows
    n_blocks_rows <- length(unique(partition$rows))
    partition_cols <- partition$cols
    n_blocks_cols <- length(unique(partition$cols))
  } else {
    n_blocks <- length(unique(partition))
  }

  best_fit <- calculate_log_linear_model(mobility_table, partition)
  best_partition <- partition

  for (i in seq_len(n)) {
    if (is.list(partition)) {
      current_block <- ifelse(move_rows, partition_rows[i], partition_cols[i])
      if(move_rows) {
        other_blocks <- setdiff(unique(partition_rows), current_block)
      } else {
        other_blocks <- setdiff(unique(partition_cols), current_block)
      }
      for (new_block in other_blocks) {
        # Move node i to the new block
        if (move_rows) {
          new_partition <- partition
          new_partition$rows[i] <- new_block
        } else {
          new_partition <- partition
          new_partition$cols[i] <- new_block
        }
      }
    } else {
      current_block <- partition[i]
      other_blocks <- setdiff(unique(partition), current_block)
      for (new_block in other_blocks) {
        # Move node i to the new block
        new_partition <- partition
        new_partition[i] <- new_block
      }
    }

    # # Check if the new partition is valid
    # if (length(unique(new_partition$rows)) < n_blocks_rows ||
    #     length(unique(new_partition$cols)) < n_blocks_cols ||
    #     any(table(new_partition$rows) < 4) ||
    #     any(table(new_partition$cols) < 4)) {
    #   next
    # }

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

  return(list(partition = best_partition, fit = best_fit))
}


#' update_multipartite_partition_node_move
#'
#' Update one of two partitions that determine cell counts
#' of three square matrices by moving the block membership of one node.
#' The partitions guide three matrices:
#' 1. A sqare matrix in which rows and columns are partitioned
#' into the same blocks (a),
#' 2. A square matrix in which rows are partitioned into blocks (a)
#' and columns are partitioned into different blocks (b),
#' 3. A square matrix in which rows and columns are partitioned
#' into the same blocks (b).
#' This function iteratively checks all nodes and moves a node
#' to a different block that results in the best fit of the
#' log-linear model that predicts cell counts.
#' In each updating step, this function will either move the
#' block membership of blocks (a) or blocks (b),
#' depending on a parameter or randomly.
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned
#' into the same blocks (b).
#' @param partition A list containing two vectors:
#' one for the partitioning of (a)
#' and one for the partitioning of (b).
#' @param swap_a A boolean indicating whether to move the block membership of blocks (a) (TRUE)
#' or blocks (b) (FALSE). Default is NULL, which means it will be chosen randomly.
#' @return A list containing the updated partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrices into 2 blocks for (a) and 3 blocks for (b)
#' init_partition <- partition_square_matrix(mat_ab, n_blocks = c(2, 3))
#' # Update the partition by moving a node
#' update_multipartite_partition_node_move(mat_aa, mat_ab, mat_bb,
#'                                         init_partition, swap_a = TRUE)
#' @export
update_multipartite_partition_node_move <- function(table_aa, table_ab, table_bb,
                                                    partition, swap_a = NULL) {
  n <- nrow(table_aa)
  swap_a <- ifelse(is.null(swap_a), sample(c(TRUE, FALSE), 1), swap_a)

  if (is.list(partition)) {
    partition_a <- partition$rows
    n_blocks_a <- length(unique(partition$rows))
    partition_b <- partition$cols
    n_blocks_b <- length(unique(partition$cols))
  } else {
    stop("partition should be a list with two vectors:
         one for the rows/a and one for the columns/b.")
  }

  best_fit_aa <- calculate_log_linear_model(table_aa, partition_a)
  best_fit_ab <- calculate_log_linear_model(table_ab, partition)
  best_fit_bb <- calculate_log_linear_model(table_bb, partition_b)

  best_partition <- partition

  for (i in seq_len(n)) {
    if (swap_a) {
      current_block <- partition_a[i]
      other_blocks <- setdiff(unique(partition_a), current_block)
      for (new_block in other_blocks) {
        # Move node i to the new block
        new_partition_a <- partition_a
        new_partition_a[i] <- new_block
        new_partition_b <- partition_b
        new_partition <- list(rows = new_partition_a, cols = new_partition_b)

        # Calculate the fit of the new partitions
        new_fit_aa <- calculate_log_linear_model(table_aa, new_partition_a)
        new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
        new_fit_bb <- calculate_log_linear_model(table_bb, new_partition_b)

        # if the fit of the new partitions cannot be calculated, skip
        if (is.null(new_fit_aa) || is.na(new_fit_aa$fit) ||
            is.null(new_fit_ab) || is.na(new_fit_ab$fit) ||
            is.null(new_fit_bb) || is.na(new_fit_bb$fit)) {
          next
        }

        # If the new fit is better, update the best fit and partition
        if ((new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit) <
            (best_fit_aa$fit + best_fit_ab$fit + best_fit_bb$fit)) {
          best_fit_aa <- new_fit_aa
          best_fit_ab <- new_fit_ab
          best_fit_bb <- new_fit_bb
          best_partition <- list(rows = new_partition_a, cols = new_partition_b)
        }
      }
    } else {
      current_block <- partition_b[i]
      other_blocks <- setdiff(unique(partition_b), current_block)
      for (new_block in other_blocks) {
        # Move node i to the new block
        new_partition_b <- partition_b
        new_partition_b[i] <- new_block
        new_partition_a <- partition_a
        new_partition <- list(rows = new_partition_a, cols = new_partition_b)

        # Calculate the fit of the new partitions
        new_fit_aa <- calculate_log_linear_model(table_aa, new_partition_a)
        new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
        new_fit_bb <- calculate_log_linear_model(table_bb, new_partition_b)

        # if the fit of the new partitions cannot be calculated, skip
        if (is.null(new_fit_aa) || is.na(new_fit_aa$fit) ||
            is.null(new_fit_ab) || is.na(new_fit_ab$fit) ||
            is.null(new_fit_bb) || is.na(new_fit_bb$fit)) {
          next
        }

        # If the new fit is better, update the best fit and partition
        if ((new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit) <
            (best_fit_aa$fit + best_fit_ab$fit + best_fit_bb$fit)) {
          best_fit_aa <- new_fit_aa
          best_fit_ab <- new_fit_ab
          best_fit_bb <- new_fit_bb
          best_partition <- list(rows = new_partition_a, cols = new_partition_b)
        }
      }
    }
  }
  return(list(partition = best_partition,
              fit = best_fit_aa$fit + best_fit_ab$fit + best_fit_bb$fit))
}



#' propose_partition_update_swap
#'
#' Propose an update of the calculated partition
#' by swapping the block membership of two nodes.
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in fit of the log-linear model of the
#' current partition and the proposed partition.
#' In case the rows and columns of the square matrix are partitioned into different blocks,
#' this function will either swap the block membership of two nodes within rows or columns,
#' depending on a parameter or randomly.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A list containing two vectors: one for the row partition
#' and one for the column partition. If only one vector is provided,
#' the rows and columns are partitioned into the same blocks.
#' @param swap_rows A boolean indicating whether to swap the block membership of rows (TRUE)
#' or columns (FALSE).
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks for rows and 3 blocks for columns
#' init_partition <- partition_square_matrix(test_mat, n_blocks = c(2, 3))
#' # Propose an update of the partition by swapping nodes
#' propose_partition_update_swap(test_mat, init_partition, swap_rows = TRUE)
#' @export
propose_partition_update_swap <- function(mobility_table, partition, swap_rows = NULL) {
  n <- nrow(mobility_table)
  swap_rows <- ifelse(is.null(swap_rows), sample(c(TRUE, FALSE), 1), swap_rows)

  if (is.list(partition)) {
    partition_rows <- partition$rows
    n_blocks_rows <- length(unique(partition$rows))
    partition_cols <- partition$cols
    n_blocks_cols <- length(unique(partition$cols))
  } else {
    n_blocks <- length(unique(partition))
  }

  current_fit <- calculate_log_linear_model(mobility_table, partition)

  # Randomly select two nodes to swap
  i <- sample(seq_len(n - 1), 1)
  j <- sample((i + 1):n, 1)

  # Swap the block membership of nodes i and j
  if (is.list(partition)){
    if (swap_rows) {
      new_partition <- partition
      new_partition$rows[c(i, j)] <- new_partition$rows[c(j, i)]
    } else {
      new_partition <- partition
      new_partition$cols[c(i, j)] <- new_partition$cols[c(j, i)]
    }
  } else {
    new_partition <- partition
    new_partition[c(i, j)] <- new_partition[c(j, i)]
  }

  # Calculate the fit of the new partition
  new_fit <- calculate_log_linear_model(mobility_table, new_partition)

  # if the fit of the new partition cannot be calculated, reject the update
  if (is.null(new_fit) || is.na(new_fit$fit)) {
    return(list(partition = partition, fit = current_fit))
  }

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < (exp(fit_diff))) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}


#' propose_multipartite_partition_update_swap
#'
#' Propose an update of one out of two given partitions
#' by swapping the block membership of two nodes.
#' The partitions guide three matrices:
#' 1. A sqare matrix in which rows and columns are partitioned
#' into the same blocks (a),
#' 2. A square matrix in which rows are partitioned into blocks (a)
#' and columns are partitioned into different blocks (b),
#' 3. A square matrix in which rows and columns are partitioned
#' into the same blocks (b).
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in summed fit of the log-linear model
#' of the three matrices comparing the current partitioning
#' and the proposed partition.
#' In each updating step, this function will either swap the
#' block membership of blocks (a) or blocks (b),
#' depending on a parameter or randomly.
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (b).
#' @param partition A list containing two vectors:
#' one for the partitioning of (a)
#' and one for the partitioning of (b).
#' @param swap_a A boolean indicating whether to swap the block membership of blocks (a) (TRUE)
#' or blocks (b) (FALSE).
#' @return A list containing the updated partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrices into 2 blocks for (a) and 3 blocks for (b)
#' init_partition <- partition_square_matrix(mat_ab, n_blocks = c(2, 3))
#' # Propose an update of the partition by swapping nodes
#' propose_multipartite_partition_update_swap(mat_aa, mat_ab, mat_bb,
#'  init_partition, swap_a = TRUE)
#' @export
propose_multipartite_partition_update_swap <- function(table_aa, table_ab, table_bb,
                                                    partition, swap_a = NULL) {
  n <- nrow(table_aa)
  swap_a <- ifelse(is.null(swap_a), sample(c(TRUE, FALSE), 1), swap_a)

  if (is.list(partition)) {
    partition_a <- partition$rows
    n_blocks_a <- length(unique(partition$rows))
    partition_b <- partition$cols
    n_blocks_b <- length(unique(partition$cols))
  } else {
    stop("partition should be a list with two vectors:
         one for the rows/a and one for the columns/b.")
  }

  current_fit_aa <- calculate_log_linear_model(table_aa, partition_a)
  current_fit_ab <- calculate_log_linear_model(table_ab, partition)
  current_fit_bb <- calculate_log_linear_model(table_bb, partition_b)

  # Randomly select two nodes to swap
  i <- sample(seq_len(n - 1), 1)
  j <- sample((i + 1):n, 1)

  # Swap the block membership of nodes i and j
  new_partition <- partition
  if (swap_a) {
    new_partition_a <- partition_a
    new_partition_a[c(i, j)] <- new_partition_a[c(j, i)]
    new_partition_b <- partition_b
    new_partition$rows <- new_partition_a
  } else {
    new_partition_b <- partition_b
    new_partition_b[c(i, j)] <- new_partition_b[c(j, i)]
    new_partition_a <- partition_a
    new_partition$cols <- new_partition_b
  }

  # Calculate the fit of the new partitions
  new_fit_aa <- calculate_log_linear_model(table_aa, new_partition_a)
  new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
  new_fit_bb <- calculate_log_linear_model(table_bb, new_partition_b)

  # if the fit of the new partitions cannot be calculated, reject the update
  if (is.null(new_fit_aa) || is.na(new_fit_aa$fit) ||
      is.null(new_fit_ab) || is.na(new_fit_ab$fit) ||
      is.null(new_fit_bb) || is.na(new_fit_bb$fit)) {
    return(list(partition = partition,
                fit = current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit))
  }

  # Calculate the difference in fit
  fit_diff <- (current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit) / 2 -
    (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit) / 2
  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition,
                fit = new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit))
  } else {
    return(list(partition = partition,
                fit = current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit))
  }
}



#' propose_partition_update_move
#'
#' Propose an update of the calculated partition
#' by moving the block membership of one node.
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in fit of the log-linear model of the
#' current partition and the proposed partition.
#' In case the rows and columns of the square matrix are partitioned into different blocks,
#' this function will either move the block membership of one node wihtin rows or columns,
#' depending on a parameter or randomly.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A list containing two vectors: one for the row partition
#' and one for the column partition. If only one vector is provided,
#' the rows and columns are partitioned into the same blocks.
#' @param swap_rows A boolean indicating whether to move the block membership of rows (TRUE)
#' or columns (FALSE).
#' @return A list containing the updated partition and the fit of the log-linear model.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks for rows and 3 blocks for columns
#' init_partition <- partition_square_matrix(test_mat, n_blocks = c(2, 3))
#' # Propose an update of the partition by moving a node
#' propose_partition_update_move(test_mat, init_partition, swap_rows = TRUE)
#' @export
propose_partition_update_move <- function(mobility_table, partition, swap_rows = NULL) {
  n <- nrow(mobility_table)
  swap_rows <- ifelse(is.null(swap_rows), sample(c(TRUE, FALSE), 1), swap_rows)

  if (is.list(partition)) {
    partition_rows <- partition$rows
    n_blocks_rows <- length(unique(partition$rows))
    partition_cols <- partition$cols
    n_blocks_cols <- length(unique(partition$cols))
  } else {
    n_blocks <- length(unique(partition))
  }

  current_fit <- calculate_log_linear_model(mobility_table, partition)

  # Randomly select a node to move
  i <- sample(seq_len(n), 1)

  # Get the current block of the node
  if (is.list(partition)) {
    current_block <- if (swap_rows) partition_rows[i] else partition_cols[i]
    other_blocks <- setdiff(if (swap_rows) unique(partition_rows) else unique(partition_cols),
                            current_block)
  } else {
    current_block <- partition[i]
    other_blocks <- setdiff(unique(partition), current_block)
  }

  # Randomly select a new block for the node
  new_block <- sample(other_blocks, 1)

  # Move node i to the new block
  new_partition <- partition
  if (is.list(partition)) {
    if (swap_rows) {
      new_partition$rows[i] <- new_block
    } else {
      new_partition$cols[i] <- new_block
    }
  } else {
    new_partition[i] <- new_block
  }

  # Calculate the fit of the new partition
  new_fit <- calculate_log_linear_model(mobility_table, new_partition)

  # if the fit of the new partition cannot be calculated, reject the update
  if (is.null(new_fit) || is.na(new_fit$fit)) {
    return(list(partition = partition, fit = current_fit))
  }

  # Calculate the difference in fit
  fit_diff <- current_fit$fit/2 - new_fit$fit/2

  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < (exp(fit_diff))) {
    return(list(partition = new_partition, fit = new_fit))
  } else {
    return(list(partition = partition, fit = current_fit))
  }
}

#' propose_multipartite_partition_update_move
#'
#' Propose an update of one out of two given partitions
#' by moving the block membership of one node.
#' The partitions guide three matrices:
#' 1. A sqare matrix in which rows and columns are partitioned
#' into the same blocks (a),
#' 2. A square matrix in which rows are partitioned into blocks (a)
#' and columns are partitioned into different blocks (b),
#' 3. A square matrix in which rows and columns are partitioned
#' into the same blocks (b).
#' Accept the update with a probability that is proportional to the
#' exp() of the difference in summed fit of the log-linear model
#' of the three matrices comparing the current partitioning
#' and the proposed partition.
#' In each updating step, this function will either move the
#' block membership of blocks (a) or blocks (b),
#' depending on a parameter or randomly.
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned
#' into the same blocks (b).
#' @param partition A list containing two vectors:
#' one for the partitioning of (a)
#' and one for the partitioning of (b).
#' @param swap_a A boolean indicating whether to move the block membership of blocks (a) (TRUE)
#' or blocks (b) (FALSE).
#' @return A list containing the updated partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrices into 2 blocks for (a) and 3 blocks for (b)
#' init_partition <- partition_square_matrix(mat_ab, n_blocks = c(2, 3))
#' # Propose an update of the partition by moving a node
#' propose_multipartite_partition_update_move(mat_aa, mat_ab, mat_bb,
#' init_partition, swap_a = TRUE)
#' @export
propose_multipartite_partition_update_move <- function(table_aa, table_ab, table_bb,
                                                    partition, swap_a = NULL) {
  n <- nrow(table_aa)
  swap_a <- ifelse(is.null(swap_a), sample(c(TRUE, FALSE), 1), swap_a)

  if (is.list(partition)) {
    partition_a <- partition$rows
    n_blocks_a <- length(unique(partition$rows))
    partition_b <- partition$cols
    n_blocks_b <- length(unique(partition$cols))
  } else {
    stop("partition should be a list with two vectors:
         one for the rows/a and one for the columns/b.")
  }

  current_fit_aa <- calculate_log_linear_model(table_aa, partition_a)
  current_fit_ab <- calculate_log_linear_model(table_ab, partition)
  current_fit_bb <- calculate_log_linear_model(table_bb, partition_b)

  # Randomly select a node to move
  i <- sample(seq_len(n), 1)

  # Get the current block of the node
  if (swap_a) {
    current_block <- partition_a[i]
    other_blocks <- setdiff(unique(partition_a), current_block)
  } else {
    current_block <- partition_b[i]
    other_blocks <- setdiff(unique(partition_b), current_block)
  }

  # Randomly select a new block for the node
  new_block <- sample(other_blocks, 1)

  # Move node i to the new block
  new_partition <- partition
  if (swap_a) {
    new_partition_a <- partition_a
    new_partition_a[i] <- new_block
    new_partition$rows <- new_partition_a
    new_partition_b <- partition_b
  } else {
    new_partition_b <- partition_b
    new_partition_b[i] <- new_block
    new_partition$cols <- new_partition_b
    new_partition_a <- partition_a
  }

  # Calculate the fit of the new partitions
  new_fit_aa <- calculate_log_linear_model(table_aa, new_partition_a)
  new_fit_ab <- calculate_log_linear_model(table_ab, new_partition)
  new_fit_bb <- calculate_log_linear_model(table_bb, new_partition_b)

  # if the fit of the new partitions cannot be calculated, reject the update
  if (is.null(new_fit_aa) || is.na(new_fit_aa$fit) ||
      is.null(new_fit_ab) || is.na(new_fit_ab$fit) ||
      is.null(new_fit_bb) || is.na(new_fit_bb$fit)) {
    return(list(partition = partition,
                fit = current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit))
  }
  # Calculate the difference in fit
  fit_diff <- (current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit) / 2 -
    (new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit) / 2
  # Accept the update with a probability proportional to exp(fit_diff)
  if (runif(1) < exp(fit_diff)) {
    return(list(partition = new_partition,
                fit = new_fit_aa$fit + new_fit_ab$fit + new_fit_bb$fit))
  } else {
    return(list(partition = partition,
                fit = current_fit_aa$fit + current_fit_ab$fit + current_fit_bb$fit))
  }
}

#' sample_likely_partition
#'
#' Sample a likely partition of the square matrix
#' Use either propose_partition_update_move or propose_partition_update_swap,
#' randomly for a specified number of iterations
#' and return the partition and the fit
#' @param mobility_table A square matrix representing the mobility table.
#' @param initial_partition A vector representing the initial partition of the
#' square matrix into blocks. If NULL, a random partition will be created.
#' @param n_blocks The number of blocks to partition the matrix into.
#' If two numbers are provided, the first number is the number of blocks for rows
#' and the second number is the number of blocks for columns.
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


#' sample_likely_multipartite_partition
#'
#' Sample a likely multipartite partition of three square matrices
#' Use either propose_multipartite_partition_update_move or
#' propose_multipartite_partition_update_swap,
#' randomly for a specified number of iterations
#' and return the partition and the fit
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned
#' into the same blocks (b).
#' @param initial_partition A list containing two vectors:
#' one for the partitioning of (a)
#' and one for the partitioning of (b).
#' If NULL, a random partition will be created.
#' @param n_blocks A vector of two numbers representing the number of
#' blocks for (a) and (b).
#' @param n_iter The number of iterations to run the sampling.
#' @return A list containing the sampled partitions (a) and (b)
#' and the summed fit of the log-linear model for the three matrices.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Sample a likely multipartite partition of the matrices with 2 blocks for (a) and 3 blocks for (b)
#' sample_likely_multipartite_partition(mat_aa, mat_ab, mat_bb,
#' n_blocks = c(2, 3), n_iter = 1000)
#' @export
sample_likely_multipartite_partition <- function(table_aa, table_ab, table_bb,
                                                  initial_partition = NULL,
                                                  n_blocks = c(2, 2), n_iter = 1000) {
  if(is.null(initial_partition)) {
    # If no initial partition is provided, create a random one
    partition <- partition_square_matrix(table_ab, n_blocks = n_blocks)
  } else {
    # Use the provided initial partition
    partition <- initial_partition
  }

  fit_aa <- calculate_log_linear_model(table_aa, partition$rows)
  fit_ab <- calculate_log_linear_model(table_ab, partition)
  fit_bb <- calculate_log_linear_model(table_bb, partition$cols)

  iter <- 0
  while (iter < n_iter) {
    iter <- iter + 1

    # Print the current iteration number and fit
    if (iter %% 100 == 0) {  # Print every 100 iterations
      cat("Iteration:", iter, "Fit (AIC):", fit_aa$fit + fit_ab$fit + fit_bb$fit, "\n")
    }

    # Randomly choose to either move or swap nodes
    if (runif(1) < 0.5) {
      result_move <- propose_multipartite_partition_update_move(table_aa, table_ab, table_bb,
                                                                partition)
    } else {
      result_move <- propose_multipartite_partition_update_swap(table_aa, table_ab, table_bb,
                                                                partition)
    }

    partition <- result_move$partition
    fit_aa <- calculate_log_linear_model(table_aa, partition$rows)
    fit_ab <- calculate_log_linear_model(table_ab, partition)
    fit_bb <- calculate_log_linear_model(table_bb, partition$cols)

  }

  return(list(partition = partition, fit = fit_aa$fit + fit_ab$fit + fit_bb$fit))
}


#' repeat_sample_likely_partition
#'
#' Repeat the function sample_likely_partition for a specified number of runs
#' the first time with a random partition
#' and then with the partition from the previous run
#' return the partition and the fit of each run as two lists
#' @param mobility_table A square matrix representing the mobility table.
#' @param n_blocks The number of blocks to partition the matrix into.
#' If two numbers are provided, the first number is the number of blocks for rows
#' and the second number is the number of blocks for columns.
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


#' repeat_sample_likely_multipartite_partition
#'
#' Repeat the function sample_likely_multipartite_partition for a specified number of runs
#' the first time with a random partition
#' and then with the partition from the previous run
#' return the partition and the fit of each run as two lists
#' @param table_aa A square matrix representing the mobility table
#' in which rows and columns are partitioned into the same blocks (a).
#' @param table_ab A square matrix representing the mobility table
#' in which rows are partitioned into blocks (a) and columns are partitioned
#' into different blocks (b).
#' @param table_bb A square matrix representing the mobility table
#' in which rows and columns are partitioned
#' into the same blocks (b).
#' @param n_blocks A vector of two numbers representing the number of
#' blocks for (a) and (b).
#' @param n_runs The number of runs to perform.
#' @param n_iter The number of iterations to run the sampling in each run.
#' @return A list containing two lists (partitions and fits), where each list contains the
#' partition and fit from each run; as well as the three mobility tables.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' mat_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Create a square matrix with 40 rows and 40 columns
#' mat_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Repeat sampling a likely multipartite partition of the matrices with 2 blocks for
#' # (a) and 3 blocks for (b) for 10 runs
#' repeat_sample_likely_multipartite_partition(mat_aa, mat_ab, mat_bb,
#' n_blocks = c(2, 3), n_runs = 10, n_iter = 1000)
#' @export
repeat_sample_likely_multipartite_partition <- function(table_aa, table_ab, table_bb,
                                                        n_blocks = c(2, 2), n_runs = 10, n_iter = 100) {
  partitions <- vector("list", n_runs)
  fits <- vector("list", n_runs)

  for (run in seq_len(n_runs)) {
    cat("Run:", run, "\n")
    if (run == 1) {
      # First run with a random partition
      result <- sample_likely_multipartite_partition(table_aa, table_ab, table_bb,
                                                      n_blocks = n_blocks, n_iter = 5 * n_iter)
    } else {
      # Subsequent runs with the partition from the previous run
      result <- sample_likely_multipartite_partition(table_aa, table_ab, table_bb,
                                                      initial_partition = partitions[[run - 1]],
                                                      n_blocks = n_blocks,
                                                      n_iter = n_iter)
    }

    partitions[[run]] <- result$partition
    fits[[run]] <- result$fit
  }

  return(list(partitions = partitions, fits = fits,
              table_aa = table_aa, table_ab = table_ab, table_bb = table_bb))
}


