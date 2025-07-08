###### Basic functions ######

#' partition_square_matrix
#'
#' randomly partition a square matrix into a specified number of
#' (approximately) equal sized blocks;
#' rows and columns can be partitioned into the same
#' number of blocks, or different numbers of blocks.
#'
#' @param mobility_table A square matrix representing the mobility table.
#' @param n_blocks The number of blocks to partition the matrix into. If
#' this is a vector of two numbers, the first number is the number of
#' blocks for the rows and the second number is the number of blocks
#' for the columns.
#' @return A vector representing the partition(s) into blocks.
#'
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' partition_square_matrix(test_mat, n_blocks = 2)
#' @export
partition_square_matrix <- function(mobility_table, n_blocks){

  if(length(n_blocks) == 1){
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

  } else if(length(n_blocks) == 2){
    n_blocks_rows <- n_blocks[1]
    n_blocks_cols <- n_blocks[2]

    n_rows <- n_cols <- ncol(mobility_table)

    if (n_rows %% n_blocks_rows != 0) {
      block_size_rows <- floor(n_rows / n_blocks_rows)
      remainder_rows <- n_rows %% n_blocks_rows
      sizes_rows <- rep(block_size_rows, n_blocks_rows)
      sizes_rows[1:remainder_rows] <- sizes_rows[1:remainder_rows] + 1
    } else {
      sizes_rows <- rep(n_rows / n_blocks_rows, n_blocks_rows)
    }
    if (n_cols %% n_blocks_cols != 0) {
      block_size_cols <- floor(n_cols / n_blocks_cols)
      remainder_cols <- n_cols %% n_blocks_cols
      sizes_cols <- rep(block_size_cols, n_blocks_cols)
      sizes_cols[1:remainder_cols] <- sizes_cols[1:remainder_cols] + 1
    } else {
      sizes_cols <- rep(n_cols / n_blocks_cols, n_blocks_cols)
    }
    partition_rows <- rep(1:n_blocks_rows, times = sizes_rows)
    partition_cols <- rep(1:n_blocks_cols, times = sizes_cols)
    partition_rows <- sample(partition_rows)  # Shuffle the partition
    partition_cols <- sample(partition_cols)  # Shuffle the partition
    partition <- list(rows = partition_rows, cols = partition_cols)
    return(partition)
  } else {
    stop("n_blocks should be a single number or a vector of two numbers.")
  }
}


#' calculate_log_linear_model
#'
#' Calculate the parameters of a log-linear model
#' in which block membership predicts cell counts of the
#' square matrix;
#' the rows and columns can be members of different blocks.
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A list containing two vectors: one for the row partition
#' and one for the column partition. If only one vector is provided,
#' the rows and columns are partitioned into the same blocks.
#' @param include_coeffs A boolean indicating whether to return the
#' coefficients of the model.
#' @return A list containing the coefficients of the model (if `include_coeffs` is TRUE)
#' and the fit of the model (AIC) - if `include_coeffs` is FALSE,
#' the returned fit is includes an unknown constant.
#' # @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks for rows and 3 blocks for columns
#' partition <- partition_square_matrix(test_mat, n_blocks = c(2, 3))
#' # Calculate the log-linear model
#' calculate_log_linear_model(test_mat, partition, include_coeffs = TRUE)
#' @export
calculate_log_linear_model <- function(mobility_table, partition,
                                       include_coeffs = FALSE) {
  n <- nrow(mobility_table)

  if(is.list(partition)) {
    partition_rows <- partition$rows
    partition_cols <- partition$cols
  } else {
    partition_rows <- partition
    partition_cols <- partition
  }

  n_blocks_rows <- length(unique(partition_rows))
  n_blocks_cols <- length(unique(partition_cols))

  if(include_coeffs) {
    # Create a data frame for the model
    df <- data.frame(
      count = as.vector(mobility_table),
      block_row = rep(partition_rows, times = n),
      block_col = rep(partition_cols, each = n),
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
    for(r in 1:n_blocks_rows){
      for(s in 1:n_blocks_cols){
        m_rs <- sum(mobility_table[partition_rows == r, partition_cols == s])
        k_r <- sum(mobility_table[partition_rows == r,])
        k_s <- sum(mobility_table[,partition_cols == s])

        LL <- LL + m_rs * log( m_rs / (k_r * k_s) )
      }
    }
    fit_aic <- -2*LL
    coeffs <- NULL
  }

  # Return the model coefficients and fit
  return(list(coefficients = coeffs, fit = fit_aic))
}
