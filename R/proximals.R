#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Proximal Operators ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Updated 13.01.2026
atan_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, atan_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 05.02.2026
exp_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, exp_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 05.02.2026
gumbel_proximal <- function(x, lambda, gamma = 0.01, ...)
{
  return(l1_proximal(x, gumbel_derivative(x, lambda, gamma)))
}

#' @noRd
# Updated 25.07.2025
l1_proximal <- function(x, lambda, ...)
{
  return(sign(x) * pmax(abs(x) - lambda, 0))
}

#' @noRd
# Updated 04.03.2026
weibull_proximal <- function(x, lambda, gamma = 0.01, shape, ...)
{
  return(l1_proximal(x, weibull_derivative(x, lambda, gamma, shape)))
}