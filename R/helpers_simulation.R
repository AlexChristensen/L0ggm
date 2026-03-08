#' @noRd
# Bounds for Weibull distribution
# Updated 08.03.2026
check_bounds <- function(shape, scale)
{
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
generate_edges <- function(nonzero, n, p)
{

  # Get Weibull model
  data("weibull_weights", package = "L0ggm", envir = environment())

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
# Get number of neighbors based on density
# Updated 08.03.2026
get_neighbors <- function(nodes, density)
{

  # Degrees of freedom
  half_df <- (nodes - 1) / 2

  # Total edges with density
  edges <- nodes * half_df * density

  # Return with minimum neighbor check
  return(
    min(
      max(round(edges / nodes), 1), # maximum neighbors
      floor(half_df) # minimum degrees of freedom
    )
  )

}
