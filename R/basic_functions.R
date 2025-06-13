###### Basic functions ######

#' partition_square_matrix
#'
#' Randomly partition a square matrix into a specified number of
#' (approximately) equal sized blocks.
#'
#' @param mobility_table A square matrix representing the mobility table.
#' @param n_blocks The number of blocks to partition the matrix into.
#' @return A vector representing the partition of the square matrix into blocks.
#'
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' partition_square_matrix(test_mat, n_blocks = 2)
#' @export
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


#' calculate_log_linear_model
#'
#' Calculate the parameters of a log-linear model
#' in which block membership predicts cell counts of the
#' square matrix;
#' optionally return the coefficients of the model
#' @param mobility_table A square matrix representing the mobility table.
#' @param partition A vector representing the partition of the square matrix
#' into blocks.
#' @param include_coeffs A boolean indicating whether to return the
#' coefficients of the model.
#' @return A list containing the coefficients of the model (if `include_coeffs` is TRUE)
#' and the fit of the model (AIC) - if `include_coeffs` is FALSE,
#' the returned fit is includes an unknown constant.
#' @examples
#' # Create a square matrix with 40 rows and 40 columns
#' test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
#' # Partition the matrix into 2 blocks
#' partition <- partition_square_matrix(test_mat, n_blocks = 2)
#' # Calculate the log-linear model
#' calculate_log_linear_model(test_mat, partition, include_coeffs = TRUE)
#' @export
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
