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
# Updated 05.02.2026
weibull_penalty <- function(x, lambda, gamma = 0.01, shape, ...)
{

  # Pre-compute components
  x <- abs(x)

  # Return penalty
  return(lambda * (1 - exp(-(x / gamma)^shape)))

}