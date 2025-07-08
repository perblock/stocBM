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

calculate_log_linear_model(test_mat, init_partition)$fit -
  calculate_log_linear_model(test_mat, init_partition2)$fit
calculate_log_linear_model(test_mat, init_partition, include_coeffs = TRUE)$fit - calculate_log_linear_model(test_mat, init_partition2, include_coeffs = TRUE)$fit


init_partition <- partition_square_matrix(test_mat, c(2,3))
init_partition2 <- partition_square_matrix(test_mat, c(2,3))

calculate_log_linear_model_2(test_mat, init_partition, include_coeffs = TRUE)$fit
calculate_log_linear_model_2(test_mat, init_partition)$fit

calculate_log_linear_model_2(test_mat, init_partition2, include_coeffs = TRUE)$fit
calculate_log_linear_model_2(test_mat, init_partition2)$fit

calculate_log_linear_model_2(test_mat, init_partition)$fit -
  calculate_log_linear_model_2(test_mat, init_partition2)$fit
calculate_log_linear_model_2(test_mat, init_partition, include_coeffs = TRUE)$fit -
  calculate_log_linear_model_2(test_mat, init_partition2, include_coeffs = TRUE)$fit


init_partition <- partition_square_matrix(test_mat, 2)

update_partition_node_swap_2(test_mat, init_partition)
update_partition_node_swap(test_mat, init_partition)


init_partition <- partition_square_matrix(test_mat, 2)

update_partition_node_move_2(test_mat, init_partition)
update_partition_node_move(test_mat, init_partition)


# use the function 'calculate_fit' from the R file 'basic_functions.R'
calculate_log_linear_model(test_mat, init_partition)
calculate_log_linear_model(test_mat, c(rep(1, 20), rep(2, 20)))

update_partition_node_swap(test_mat, init_partition)

repeat_until_no_improvement(test_mat, n_blocks = 3, max_iter = 1000)

all(init_partition == propose_partition_update_move_2(test_mat, init_partition)$partition)

all(unlist(init_partition2) == unlist(propose_partition_update_move_2(test_mat, init_partition2)$partition))


all(init_partition == propose_partition_update_swap(test_mat, init_partition)$partition)

all(unlist(init_partition) == unlist(propose_partition_update_swap_2(test_mat, init_partition)$partition))


propose_partition_update_swap_2(test_mat, init_partition)

sample_likely_partition(test_mat, n_blocks = 3)

repeat_sample_likely_partition(test_mat, n_blocks = 3, n_runs = 10)



# Create a square matrix with 40 rows and 40 columns
test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_mat[1:15, 1:15] <- test_mat[1:15, 1:15] + 3  # Add some structure
test_mat[25:40, 16:24] <- test_mat[25:40, 16:24] + 4  # Add some structure
# Repeat sampling a likely partition of the matrix with 2 blocks for 10 runs
repeat_result <- repeat_sample_likely_partition(test_mat, n_blocks = 3,
                                                n_runs = 10, n_iter = 1000)
# Get good partitions from the repeat result
res <- get_good_partitions(repeat_result, max_iter = 100)
res

# now with different row/column partitions
test_mat <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_mat[1:15, 1:20] <- test_mat[1:15, 1:20] + 4  # Add some structure
test_mat[25:40, 21:40] <- test_mat[25:40, 21:40] + 4  # Add some structure
# Repeat sampling a likely partition of the matrix with 2 blocks for 10 runs
repeat_result <- repeat_sample_likely_partition(test_mat, n_blocks = c(3,2),
                                                n_runs = 30, n_iter = 1000)
# Get good partitions from the repeat result
res <- get_good_partitions(repeat_result, max_iter = 100)
res


# now multipartite partitions
test_aa <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_aa[1:15, 1:15] <- test_aa[1:15, 1:15] + 1  # Add some structure
test_aa[25:40, 16:24] <- test_aa[25:40, 16:24] + 1  # Add some structure
test_aa[16:24, 25:40] <- test_aa[16:24, 25:40] + 1  # Add some structure

test_ab <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_ab[1:15, 1:20] <- test_aa[1:15, 1:20] + 2  # Add some structure
test_ab[16:24, 1:20] <- test_aa[16:24, 1:20] + 2  # Add some structure
test_ab[25:40, 21:40] <- test_aa[25:40, 21:40] + 2  # Add some structure

test_bb <- matrix(rpois(1600, lambda = 10), nrow = 40, ncol = 40)
test_bb[1:20, 1:20] <- test_bb[1:20, 1:20] + 1  # Add some structure
test_bb[21:40, 21:40] <- test_bb[21:40, 21:40] + 1  # Add some structure

repeat_result <- repeat_sample_likely_multipartite_partition(test_aa, test_ab, test_bb,
                                                             n_blocks = c(3,2),
                                                             n_runs = 30, n_iter = 1000)
multipartite_repeat_until_no_improvement(test_aa, test_ab, test_bb,
                                         max_iter = 20, init_partition = repeat_result$partitions[[28]])

res <- get_good_multipartite_partitions(repeat_result)
res


