#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Maximum Likelihood Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# MLE Gumbel Scale Parameter ----
# Updated 05.02.2026
gumbel_mle <- function(x)
{

  # Set up MLE for scale
  scale_mle <- function(scale, x, n)
  {

    # Pre-compute reused values
    x_scale <- x / scale

    # Return log-likelihood
    return(-n * log(scale) - sum(x_scale) - sum(exp(-x_scale)))

  }

  # Return parameters
  return(
    optimize(
      f = scale_mle, interval = c(1e-04, 1),
      x = x, n = length(x), maximum = TRUE
    )$maximum
  )

}

#' @noRd
# MLE Weibull Parameters ----
# Updated 10.01.2026
weibull_mle <- function(x)
{

  # Set up MLE for shape
  shape_mle <- function(k, x, n)
  {

    # Pre-compute reused values
    x_k <- x^k
    log_x <- log(x)

    # Return log-likelihood
    return(sum(x_k * log_x) / sum(x_k) - 1 / k - sum(log_x) / n)

  }

  # Obtain MLE estimate
  shape <- uniroot(
    f = shape_mle, interval = c(0.01, 20),
    x = x, n = length(x)
  )$root

  # Return parameters
  return(c(shape = shape, scale = mean(x^shape)^(1 / shape)))

}