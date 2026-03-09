#%%%%%%%%%%%%%%%%%%%#
#### Derivatives ####
#%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Updated 27.02.2026
atan_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (gamma * (gamma + 2 / pi)) / (gamma^2 + x^2))
}

#' @noRd
# Updated 27.02.2026
exp_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  return(lambda * (1 / gamma) * exp(-(abs(x) / gamma)))
}

#' @noRd
# Updated 27.02.2026
gumbel_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # Pre-compute values
  gamma_x <- abs(x) / gamma

  # Return derivative
  return((lambda / (1 - exp(-1))) * (1 / gamma) * exp(-gamma_x - exp(-gamma_x)))

}

#' @noRd
# Updated 08.03.2026
log_derivative <- function(x, lambda, gamma = 0.10, ...)
{
  return(lambda / ((gamma + abs(x)) * log(1 + 1 / gamma)))
}

#' @noRd
# Updated 03.03.2026
weibull_derivative <- function(x, lambda, gamma = 0.01, shape, ...)
{
  # Pre-compute components
  abs_x <- pmax.int(abs(x), .Machine$double.eps)
  x_gamma <- abs_x / gamma

  # Return derivative
  return(lambda * (shape / gamma) * x_gamma^(shape - 1) * exp(-x_gamma^shape))

}