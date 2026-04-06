#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Maximum Likelihood Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# MLE Gumbel Scale Parameter ----
# Updated 06.04.2026
gumbel_mle <- function(x)
{

  # Set up MLE for scale
  scale_mle <- function(scale, x, n)
  {

    # Pre-compute reused values
    x_scale <- x / scale

    # Return log-likelihood
    return(-n * scale + sum(x) - sum(x * exp(-x_scale)))

  }

  # Return parameters
  return(uniroot(f = scale_mle, interval = c(1e-04, max(x)), x = x, n = length(x))$root)

}

#' @noRd
# MLE Weibull Parameters ----
# Updated 06.04.2026
weibull_mle <- function(x)
{

  # Ensure no exact zeros to avoid log(0)
  x <- pmax.int(x, .Machine$double.eps)

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
    f = shape_mle, interval = c(0.01, 10),
    x = x, n = length(x), extendInt = "upX"
  )$root

  # Return parameters
  return(c(shape = shape, scale = mean(x^shape)^(1 / shape)))

}