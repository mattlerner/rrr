library(MASS)

############ Define parameters for simulation ############

sample_size <- 10000 # e.g. number of observations
dimension <- 50     # number of covariates
correlated_columns <- rbind(c(2,3), c(15,22)) #specify which pairs of columns are correlated. leave as rbind() if all columns are independent

############ Define functions ############

# This function takes the intended dimension of the covariate matrix
# as an argument. It then randomly selects standard deviations for each
# dimension from the specified distribution (in this case rnorm(0,3))
generate_sd <- function(dimension) {
  sd_vector <- abs(rnorm(dimension, 0, 3))
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
      correlation <- runif(1,-1,1) # generate a random correlation between 1 and -1
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

############ Run simulation ############

simulated_correlations <- generate_correlation_matrix(dimension, correlated_columns)
simulated_sd <- generate_sd(dimension)
simulated_Sigma <- generate_covariance_matrix(simulated_sd, simulated_correlations)
simulated_means <- generate_means(1:20, dimension)

simulation <- mvrnorm(n = sample_size, simulated_means, simulated_Sigma)
simulation_frame <- as.data.frame(simulation) # this can be sent to csv

############ Output ############

# Identify columns and associated correlations and means
for (i in 1:dim(correlated_columns)[1]) {
  current_pair <- correlated_columns[i,]
  current_correlation <- simulated_correlations[current_pair[1], current_pair[2]]
  print(paste("Columns",current_pair[1],"and",current_pair[2],"are correlated at",round(current_correlation,2)))
}
