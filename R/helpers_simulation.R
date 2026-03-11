#' @noRd
# Bounds for Weibull distribution
# Updated 08.03.2026
check_bounds <- function(shape, scale)
{
  # The bounds were derived from the 222 empirical datasets
  # from Huth et al. (2025), representing most "standard"
  # conditions that occur in empirical data (round up and down)
  # (see `?weibull_descriptives` for selection criterion)
  # These ranges will be updated as more empirical data
  # are aggregated
  (shape < 0.70) | (shape > 1.50) | (scale < 0.03) | (scale > 0.19)
}

#' @noRd
# Generate from Weibull distribution
# Updated 08.03.2026
weibull_xoshiro <- function(n, shape, scale)
{
  return(scale * (-log(runif_xoshiro(n)))^(1 / shape))
}

#' @noRd
# Generate edge values
# Updated 08.03.2026
generate_edges <- function(nonzero, n, p)
{

  # Check for empirical parameter bounds and maximal partial correlation
  outside_bounds <- greater_than <- TRUE

  # Generate until no edges are greater than one
  while(greater_than | outside_bounds){

    # Generate Weibull parameters
    params <- weibull_parameters(nodes = p, sample_size = n, bootstrap = TRUE)

    # Generate edge weights
    weights <- weibull_xoshiro(nonzero, shape = params[["shape"]], scale = params[["scale"]])

    # Update greater than
    greater_than <- any(weights >= 0.85)
    # No weights greater than max weight = 0.8452469

    # Check outside of bounds
    outside_bounds <- check_bounds(params[["shape"]], params[["scale"]])

  }

  # Attach attributes to weights
  attr(weights, "params") <- params

  # Return weights
  return(weights)

}

#' @noRd
# Minimally condition the network
# Updated 11.03.2026
condition_network <- function(network, target_condition)
{

  # Obtain pre-inversion matrix
  Omega <- -network
  diag(Omega) <- 1

  # Compute one plus minimum eigenvalue
  min_eigen <- 1 + abs(min(matrix_eigenvalues(Omega)))

  # Condition network
  opt <- uniroot(
    f = function(lambda){

      # Create precision
      diag(Omega) <- min_eigen + lambda

      # Compute the precision matrix
      K <- try(solve(Omega), silent = TRUE)

      # Catch bad matrices
      if(is(K, "try-error")){
        return(target_condition)
      }else{
        return(fast_kappa(cov2cor(K)) - target_condition)
      }

    }, interval = c(0.001, 1)
  )

  # Set diagonal on precision
  diag(Omega) <- min_eigen + opt$root

  # Return results
  return(list(R = cov2cor(solve(Omega)), lambda = opt$root))

}

#' @noRd
# Simulate data
# Updated 09.03.2026
simulate_data <- function(n, R, skew, skew_range)
{

  # Check for skew range
  if(!is.null(skew_range)){
    typeof_error(skew_range, "numeric") # object type error
    length_error(skew_range, 2) # object length error
    range_error(skew_range, c(-2, 2)) # object range error
    possible_skews <- seq(-2, 2, 0.05) # possible skews
    skew_range <- round(skew_range, 2) # get to hundredths digit
    min_range <- abs(min(skew_range) - possible_skews) # difference for minimum
    min_skew <- possible_skews[which.min(min_range)] # get minimum skew
    max_range <- abs(max(skew_range) - possible_skews) # difference for maximum
    max_skew <- possible_skews[which.min(max_range)] # get maximum skew
    skew <- seq(min_skew, max_skew, 0.05) # obtain skews
  }

  # Generate data
  output <- mvrnorm_precompute(n, R)
  data <- MASS_mvrnorm_quick(p = output$p, np = output$np, coV = output$coV)
  data <- data %*% chol(R)

  # Set skew
  n_skew <- length(skew)
  if(n_skew == 1){
    skew <- rep(skew, output$p)
  }else if(n_skew != output$p){
    skew <- shuffle_replace(skew, output$p)
  }

  # Loop through columns
  node_sequence <- seq_len(output$p)
  for(i in node_sequence){
    data[,i] <- skew_continuous(skewness = skew[i], data = data[,i])
  }

  # Add column names to data
  names(skew) <- colnames(data) <- paste0(
    "V", format_integer(node_sequence, digits(output$p) - 1)
  )

  # Return data and skew
  return(list(data = data, skew = skew))

}