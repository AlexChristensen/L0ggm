#' @noRd
# Bounds for Weibull distribution
# Updated 08.03.2026
check_bounds <- function(shape, scale)
{
  # The bounds were derived from the 203 empirical datasets
  # from Huth et al. (2025), representing most "standard"
  # conditions that occur in empirical data
  # (see `?weibull_descriptives` for selection criterion)
  # These ranges will be updated as more empirical data
  # are aggregated
  (shape < 0.7177596) | (shape > 1.6276798) | (scale < 0.03275058) | (scale > 0.18792034)
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
generate_edges <- function(weibull_weights, nonzero, n, p)
{

  # Check for empirical parameter bounds and maximal partial correlation
  outside_bounds <- greater_than <- TRUE

  # Generate until no edges are greater than one
  while(greater_than | outside_bounds){

    # Generate shape value based on normal distribution of empirical
    shape <- 1.074521 + rnorm_ziggurat(1) * 0.1275993

    # Generate scale value
    scale <- predict(
      weibull_weights, interval = "none",
      newdata = data.frame(
        shape = shape, scaling = 1 / sqrt(n),
        rlp = 1 / log(p), ppo = (p * (p - 1) / 2) / n
      )
    )

    # Use residual bootstrapping
    scale <- as.numeric(scale + shuffle(weibull_weights$residuals, size = 1))

    # Generate edge weights
    weights <- weibull_xoshiro(nonzero, shape = shape, scale = scale)

    # Update greater than
    greater_than <- any(weights >= 0.85)
    # No weights greater than max weight = 0.8452469

    # Check outside of bounds
    outside_bounds <- check_bounds(shape, scale)

  }

  # Attach attributes to weights
  attr(weights, "params") <- c(shape = shape, scale = scale)

  # Return weights
  return(weights)

}

#' @noRd
# Minimally condition the network
# Updated 08.03.2026
condition_network <- function(network)
{

  # Obtain pre-inversion matrix
  Omega <- -network
  diag(Omega) <- 1

  # Compute minimum eigenvalue
  min_eigen <- min(matrix_eigenvalues(Omega))

  # Set maximum correlation
  max_R <- max(abs(network))

  # Pre-compute lower triangle
  lower_triangle <- lower.tri(network)

  # Condition network
  opt <- uniroot(
    f = function(lambda){

      # Create precision
      diag(Omega) <- 1 + abs(min_eigen) + lambda

      # Compute correlation matrix
      R <- try(cov2cor(solve(Omega)), silent = TRUE)

      # Catch bad matrices
      if(is(R, "try-error")){
        return(2)
      }else{
        return(max(abs(R[lower_triangle])) - max_R)
      }

    }, interval = c(0.001, 0.300)
  )

  # Set diagonal on precision
  diag(Omega) <- 1 + abs(min_eigen) + opt$root

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