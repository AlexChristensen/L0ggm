#%%%%%%%%%%%%%%%%%#
#### Penalties ####
#%%%%%%%%%%%%%%%%%#

#' @noRd
# Updated 25.07.2025
atan_penalty <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (gamma + 2 / pi) * atan(abs(x) / gamma))
}

#' @noRd
# Updated 10.01.2026
exp_penalty <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(lambda * (1 - exp(-(x / gamma))))

}

#' @noRd
# Updated 09.02.2026
gumbel_penalty <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute
  exp_1 <- exp(-1)

  return((lambda / (1 - exp_1)) * (exp(-exp(-abs(x) / gamma)) - exp_1))
  # theoretically, the `- exp(-1)` is necessary for the
  # penalty to converge at zero for the sparsity condition
  # in practice, this addition does not change the derivative,
  # which is used in the LLA
  # `1 - exp(-1)` is to scale lambda
}

#' @noRd
# Updated 08.03.2026
log_penalty <- function(x, lambda, gamma = 0.10, ...)
{
  return(lambda * log(1 + abs(x) / gamma) / log(1 + 1 / gamma))
}

#' @noRd
# Updated 15.03.2026
weibull_penalty <- function(x, lambda, gamma = 0.01, shape, ...)
{

  # Pre-compute
  x <- abs(x)

  # Check for shape greater than one
  if(shape > 1){

    # Mirror derivative
    shape_one  <- shape - 1
    peak       <- gamma * ((shape - 1) / shape)^(1 / shape)
    peak_gamma <- peak / gamma
    peak_value <- (shape / gamma) * peak_gamma^shape_one * exp(-peak_gamma^shape)

    # Value of Weibull penalty at peak
    at_peak <- 1 - exp(-peak_gamma^shape)

    # Shift so pieces meet at peak, anchored at zero
    shift <- peak_value * peak - at_peak

    # Return penalty
    return(lambda * swiftelse(x <= peak, peak_value * x, 1 - exp(-(x / gamma)^shape) + shift))

  }

  # Otherwise, return standard penalty
  return(lambda * (1 - exp(-(x / gamma)^shape)))

}
