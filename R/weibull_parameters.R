#' Predict Weibull Parameters for Edge Weight Distributions
#'
#' @description
#' Predicts the shape and scale parameters of a Weibull distribution that
#' characterizes the absolute partial correlation edge weights of a psychometric
#' network, given its number of nodes. Parameter estimates are derived from a
#' Seemingly Unrelated Regression (SUR) model fitted to empirical network data
#' from Huth et al. (2025), where absolute partial correlations were found to
#' follow a Weibull distribution more consistently than Beta, Gamma, or
#' log-normal alternatives.
#'
#' @param nodes Integer. The number of nodes (variables) in the network.
#' Must be between 8 and 40, reflecting the range of the empirical networks
#' used to fit the underlying SUR model.
#'
#' @param bootstrap Logical. If \code{TRUE}, a randomly sampled residual from
#' the SUR model fit is added to each predicted parameter, introducing
#' empirically grounded variability suitable for use in simulation or
#' bootstrapping contexts. Defaults to \code{FALSE}.
#'
#' @details
#' Both shape and scale are predicted from a single network descriptor:
#'
#' \describe{
#'   \item{\code{rlp}}{The reciprocal log of the number of nodes,
#'     \eqn{1 / \log(p)}, capturing the diminishing marginal effect of network
#'     size on edge weight distributions.}
#' }
#'
#' This predictor enters two SUR equations — one for shape and one for
#' log-transformed scale — whose coefficients are stored in the internal
#' \code{weibull_weights} dataset. SUR was used to account for the correlated
#' residuals between the shape and scale equations across networks (residual
#' correlation = 0.605). The scale equation was fitted on \code{log(scale)} to
#' satisfy residual normality; predictions are back-transformed via \code{exp()}
#' before being returned. When \code{bootstrap = TRUE}, residuals are added in
#' log space prior to back-transformation, preserving the correct error
#' structure.
#'
#' Population edge weights are treated as independent of sample size by design:
#' \code{nodes} is the sole predictor because the Weibull parameters describe
#' a fixed population partial correlation structure, not the estimation
#' characteristics of any particular sample drawn from it.
#'
#' Empirically, shape values ranged from approximately 0.72 to 1.63
#' (M = 1.08, SD = 0.13) and scale values from approximately 0.03 to 0.19
#' (M = 0.10, SD = 0.03) across the 190 networks used to fit the model. Shape
#' values near 1 indicate approximately exponential edge weight distributions;
#' values above 1 indicate a rising hazard (mode-bearing distribution).
#'
#' When \code{bootstrap = TRUE}, residuals are drawn via \code{shuffle()} --
#' a random sampling without replacement — from the empirical SUR residuals,
#' preserving the observed marginal residual distribution.
#'
#' @return A named numeric vector of length 2:
#' \describe{
#'   \item{\code{shape}}{The predicted Weibull shape parameter (\eqn{k > 0}).}
#'   \item{\code{scale}}{The predicted Weibull scale parameter (\eqn{\lambda > 0}),
#'   back-transformed from log space.}
#' }
#'
#' @examples
#' # Predict parameters for a 10-node network
#' weibull_parameters(nodes = 10)
#'
#' # With bootstrapped residuals for use in simulation
#' weibull_parameters(nodes = 10, bootstrap = TRUE)
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., &
#' Marsman, M. (2025). Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#'
#' @export
#'
# Predict Weibull parameters ----
# Updated 10.03.2026
weibull_parameters <- function(nodes, bootstrap = FALSE)
{

  # Get Weibull model
  weibull_weights <- get(data("weibull_weights", package = "L0ggm", envir = environment()))

  # Compute descriptive parameters
  parameters <- c( rlp = 1 / log(nodes))

  # Compute Weibull parameters
  shape <- unname(
    weibull_weights$shape$coefficients[1] +
    sum(weibull_weights$shape$coefficients[-1] * parameters)
  )
  scale <- unname(
    weibull_weights$scale$coefficients[1] +
    sum(weibull_weights$scale$coefficients[-1] * parameters)
  )

  # Bootstrap residuals
  if(bootstrap){

    # Draw index
    index <- shuffle(seq_along(weibull_weights$shape$residuals), 1)

    shape <- shape + weibull_weights$shape$residuals[index]
    scale <- scale + weibull_weights$scale$residuals[index]
  }

  # Return parameters (scale needs a back-transform)
  return(c(shape = shape, scale = exp(scale)))

}