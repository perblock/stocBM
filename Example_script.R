##### WARNING: This is an example script for testing purposes only. ######

# create a matrix with 40 rows and 40 columns;
# numbers in the matrix are integers that follow a poisson distribution

test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)

test_mat[1:15, 1:15] <- test_mat[1:15, 1:15] + 2  # Add some structure
test_mat[25:40, 16:24] <- test_mat[25:40, 16:24] + 5  # Add some structure

# use the function 'partition_square_matrix'
# from the R file 'basic_functions.R' on test_mat
# with 2 blocks

init_partition <- partition_square_matrix(test_mat, 2)
init_partition2 <- partition_square_matrix(test_mat, 2)

calculate_log_linear_model(test_mat, init_partition, include_coeffs = TRUE)$fit
calculate_log_linear_model(test_mat, init_partition)$fit

calculate_log_linear_model(test_mat, init_partition2, include_coeffs = TRUE)$fit
calculate_log_linear_model(test_mat, init_partition2)$fit

calculate_log_linear_model(test_mat, init_partition)$fit - calculate_log_linear_model(test_mat, init_partition2)$fit
calculate_log_linear_model(test_mat, init_partition, include_coeffs = TRUE)$fit - calculate_log_linear_model(test_mat, init_partition2, include_coeffs = TRUE)$fit


# use the function 'calculate_fit' from the R file 'basic_functions.R'
calculate_log_linear_model(test_mat, init_partition)
calculate_log_linear_model(test_mat, c(rep(1, 20), rep(2, 20)))

update_partition_node_swap(test_mat, init_partition)

repeat_until_no_improvement(test_mat, n_blocks = 3, max_iter = 1000)

all(init_partition == propose_partition_update_move(test_mat, init_partition)$partition)

all(init_partition == propose_partition_update_swap(test_mat, init_partition)$partition)

sample_likely_partition(test_mat, n_blocks = 3)

repeat_sample_likely_partition(test_mat, n_blocks = 3, n_runs = 10)


