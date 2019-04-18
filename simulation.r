library(MASS)

setwd("/users/matt/desktop/qmss/multivariate statistical inference/group")

############ Define parameters for all simulations ############

sample_size <- 10000 # e.g. number of observations
x_dimension <- 10     # number of X covariates
y_dimension <- 3 # desired number of Y columns

############ Define functions ############

# This function takes the intended dimension of the covariate matrix
# as an argument. It then randomly selects standard deviations for each
# dimension from the specified distribution (in this case rnorm(0,3))
generate_sd <- function(dimension) {
  sd_vector <- abs(rnorm(dimension, 0, 15))
  return(sd_vector)
}

# Takes a scalar (dimension) and a matrix of 2-tuples e.g. rbind(c(1,5),c(8,9) etc.)
# The matrix of 2-tuples (correlated_pairs) is user-specified. The matrix should have
# two columns. Each column contains the indices of two columns in the simulated
# dataset that we defined to be correlated, e.g. rbind(c(1,5),c(8,9) specifies
# that in our simulated dataset, columns 1 and 5 will be correlated and columns
# 8 and 9 will be correlated. When this argument is blank, all columns are independent.
generate_correlation_matrix <- function(dimension, correlated_pairs=rbind()) {
  matrix_ <- matrix(0, nrow=dimension, ncol=dimension)
  
  # everything is correlated with itself of course
  for (i in 1:dimension) {
    matrix_[i,i] <- 1
  }
  
  if(!is.null(nrow(correlated_pairs))) {
    for (j in 1:nrow(correlated_pairs)) {
      tuple <- correlated_pairs[j,]
      correlation <- runif(1,0.7,1) * sample(c(-1,1),1) # generate a random correlation between 0.5 and 1 and then make it either positive or negative (randomly)
      matrix_[tuple[1],tuple[2]] <- correlation
      matrix_[tuple[2],tuple[1]] <- correlation
    }
  }
  return(matrix_)
}

# This function takes in a vector of standard deviations of dimension n
# and a correlation matrix of size n*n. If the dimension of correlation_matrix
# and the length of sd_vector are different, this will throw an error.
# The logic here is that element i,j of a covariance matrix is given
# by Covar(Xi,Xj).  We recall that cor(Xi, Xj) = Covar(Xi,Xj) / Sd(Xi)Sd(Xj) = 
# Covar(Xi,Xj) / Sd(Xi)Sd(Xj). So we note that that Covar(Xi,Xj) is the
# element-wise product of Cor(Xi,Xj) and Sd(Xi)Sd(Xj), two i*j matrices
# given by (obviously) the correlation matrix of X and the matrix given by
# the vector of standard deviations times its transpose.
generate_covariance_matrix <- function(sd_vector, correlation_matrix) {
  return(sd_vector %*% t(sd_vector) * correlation_matrix)
}

# Given a vector of possible values and a dimension,
# this function returns a vector of means by sampling
# from the provided value vector with replacement.
generate_means <- function(range, dimension) {
  means <- sample(range, dimension, replace=TRUE)
  return(means)
}

generate_coefficient_matrix <- function(covariate_dimension, desired_output_dimension, zeroes = FALSE) {
  if (zeroes == TRUE) {
    B <- replicate(desired_output_dimension, rep(0,covariate_dimension))
    return(B)
  } else {
    B <- replicate(desired_output_dimension, runif(covariate_dimension,0.7,1) * sample(c(-1,1),covariate_dimension, replace=TRUE)) # generate sandom correlations between 0.2 and 1 and then make it either positive or negative (randomly)
    return(B)
  }
}

# given an input covariate matrix, generate an output matrix
# with the specified output_dimension.
# by default, columns are correlated.
# use generate_output_matrix_w_uncorrelated to do the same thing but
# to add the specified number of uncorrelated columns
generate_output_matrix <- function(input_covariates, output_dimension, coefficient_matrix_B) {
  X <- data.matrix(input_covariates)
  
  input_dimension <- dim(X)[2]
  obs <- dim(X)[1]
  sd_input <- sd(X) / 2 # this is just a benchmark to make sure noise doesn't swamp signal here

  B <- coefficient_matrix_B
  e <- replicate(output_dimension, rnorm(obs,0,sd_input))

  Y <- (X %*% B) + e
  return(Y)
}

############ Case 1 - Independent X - Independent Y ############

simulated_correlations_1 <- generate_correlation_matrix(x_dimension, c()) # correlations of elements of X
simulated_sd_1 <- generate_sd(x_dimension)
simulated_Sigma_1 <- generate_covariance_matrix(simulated_sd_1, simulated_correlations_1)
simulated_means_1 <- generate_means(1:100, x_dimension)

X_1 <- mvrnorm(n = sample_size, simulated_means_1, simulated_Sigma_1)
B_1 = generate_coefficient_matrix(x_dimension, y_dimension, zeroes=TRUE)
Y_1 <- generate_output_matrix(X_1, 3, B_1)

save(X_1, B_1, Y_1, file="case1.rdata")



############ Case 2 - Independent X - Correlated Y ############

X_2 <- X_1

B_2 <- generate_coefficient_matrix(x_dimension, y_dimension)
Y_2 <- generate_output_matrix(X_2, y_dimension, B_2)

save(X_2, B_2, Y_2, file="case2.rdata")

############ Case 3 - Correlated X - Correlated Y ############

simulated_correlations_3 <- generate_correlation_matrix(x_dimension, rbind(c(1,2),c(1,3),c(2,3),c(3,4),c(4,5),c(2,4))) # correlations of elements of X
simulated_sd_3 <- generate_sd(x_dimension)
simulated_Sigma_3 <- generate_covariance_matrix(simulated_sd_3, simulated_correlations_3)
simulated_means_3 <- generate_means(1:100, x_dimension)

X_3 <- mvrnorm(n = sample_size, simulated_means_3, simulated_Sigma_3, tol=1)
B_3 <- B_2
Y_3 <- generate_output_matrix(X_3, y_dimension, B_3)

save(X_3, B_3, Y_3, simulated_correlations_3, file="case3.rdata")

