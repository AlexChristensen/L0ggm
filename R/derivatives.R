#%%%%%%%%%%%%%%%%%%%#
#### Derivatives ####
#%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Updated 15.03.2026
atan_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  # return(lambda * (gamma * (gamma + 2 / pi)) / (gamma^2 + x^2))

  # Return derivative
  return(
    .Call(
      "atan_derivative_c",
      as.numeric(x), lambda_ = lambda, gamma_ = gamma,
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Updated 15.03.2026
exp_derivative <- function(x, lambda, gamma = 0.01, ...)
{
  # return(lambda * (1 / gamma) * exp(-(abs(x) / gamma)))

  # Return derivative
  return(
    .Call(
      "exp_derivative_c",
      as.numeric(x), lambda_ = lambda, gamma_ = gamma,
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Updated 15.03.2026
gumbel_derivative <- function(x, lambda, gamma = 0.01, ...)
{

  # # Pre-compute values
  # gamma_x <- abs(x) / gamma
  #
  # # Return derivative
  # return((lambda / gamma) * (exp(-gamma_x - exp(-gamma_x)) + exp(gamma_x - exp(gamma_x))))

  # Return derivative
  return(
    .Call(
      "gumbel_derivative_c",
      as.numeric(x), lambda_ = lambda, gamma_ = gamma,
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Updated 15.03.2026
log_derivative <- function(x, lambda, gamma = 0.10, ...)
{
  # return(lambda / ((gamma + abs(x)) * log(1 + 1 / gamma)))

  # Return derivative
  return(
    .Call(
      "log_derivative_c",
      as.numeric(x), lambda_ = lambda, gamma_ = gamma,
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Updated 15.03.2026
weibull_derivative <- function(x, lambda, gamma = 0.01, shape, ...)
{

  # # Pre-compute components
  # abs_x <- pmax.int(abs(x), .Machine$double.eps)
  # x_gamma <- abs_x / gamma
  #
  # # Check for shape greater than 1
  # if(shape > 1){
  #
  #   # Pre-compute shape
  #   shape_one <- shape - 1
  #
  #   # Obtain peak
  #   peak <- gamma * ((shape - 1) / shape)^(1 / shape)
  #
  #   # Compute peak gamma
  #   peak_gamma <- peak / gamma
  #
  #   # Compute values
  #   peak_value <- (shape / gamma) * peak_gamma^shape_one * exp(-peak_gamma^shape)
  #
  #   # Return
  #   return(
  #     lambda * swiftelse(
  #       abs_x <= peak, peak_value,
  #       (shape / gamma) * x_gamma^shape_one * exp(-x_gamma^shape)
  #     )
  #   )
  #
  # }
  #
  # # Return derivative
  # return(lambda * (shape / gamma) * x_gamma^(shape - 1) * exp(-x_gamma^shape))

  # Return derivative
  return(
    .Call(
      "weibull_derivative_c",
      as.numeric(x), lambda_ = lambda,
      gamma_ = gamma, shape_ = shape,
      PACKAGE = "L0ggm"
    )
  )

}