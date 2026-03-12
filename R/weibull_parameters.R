#' Predict Weibull Parameters for Edge Weight Distributions
#'
#' @description
#' Predicts the shape and scale parameters of a Weibull distribution that
#' characterizes the absolute partial correlation edge weights of a psychometric
#' network, given its number of nodes and sample size. Parameter estimates are
#' derived from a Seemingly Unrelated Regression (SUR) model fitted to empirical
#' network data from Huth et al. (2025), where absolute partial correlations were
#' found to follow a Weibull distribution more consistently than Beta, Gamma, or
#' log-normal alternatives.
#'
#' @param nodes Integer. The number of nodes (variables) in the network.
#' Must be between 8 and 40, reflecting the range of the empirical networks
#' used to fit the underlying SUR model.
#'
#' @param sample_size Integer. The sample size of the dataset from which the
#' network is estimated.
#'
#' @param bootstrap Logical. If \code{TRUE}, a randomly sampled residual from
#' the SUR model fit is added to each predicted parameter, introducing
#' empirically grounded variability suitable for use in simulation or
#' bootstrapping contexts. Defaults to \code{FALSE}.
#'
#' @details
#' The shape and scale parameters are predicted from derived network descriptors
#' that differ between the two SUR equations:
#'
#' \describe{
#'   \item{\code{rlp}}{The reciprocal log of the number of nodes,
#'     \eqn{1 / \log(p)}, capturing the diminishing marginal effect of network
#'     size on edge weight distributions. Used in both the shape and scale
#'     equations.}
#'   \item{\code{beta_min}}{A noise-to-signal proxy defined as
#'     \eqn{\sqrt{\log(p) / n}}, where \eqn{p} is the number of nodes and
#'     \eqn{n} is the sample size. Larger values indicate sparser, noisier
#'     estimation conditions. Used in the shape equation only.}
#' }
#'
#' These predictors enter two SUR equations whose coefficients are stored in
#' the internal \code{weibull_weights} dataset. SUR was used to account for the
#' correlated residuals between the shape and scale equations across networks
#' (residual correlation = 0.498). Model fit was acceptable: the shape equation
#' achieved \eqn{R^2 = 0.209} (RMSE = 0.119) and the scale equation achieved
#' \eqn{R^2 = 0.686} (RMSE = 0.017). Residuals were normally distributed for
#' both equations (Shapiro-Wilk p = 0.865 and p = 0.143, respectively).
#' Heteroskedasticity was detected in the shape equation (Breusch-Pagan
#' p < 0.001) but not the scale equation (p = 0.344); robust standard errors
#' (HC3) were used for inference. Multicollinearity among shape predictors was
#' negligible (VIF < 1.001, condition number = 42.0).
#'
#' Notably, the scale equation includes only \code{rlp} and no
#' sample-size-dependent term, making predicted scale invariant to \code{n}.
#' Shape predictions vary with both \code{nodes} and \code{sample_size} via
#' \code{beta_min}.
#'
#' Empirically, shape values ranged from approximately 0.72 to 1.48
#' (M = 1.07, SD = 0.13) and scale values from approximately 0.03 to 0.16
#' (M = 0.10, SD = 0.03) across the 189 networks used to fit the model. Shape
#' values near 1 indicate approximately exponential edge weight distributions;
#' values above 1 indicate a rising hazard (mode-bearing distribution).
#'
#' When \code{bootstrap = TRUE}, residuals are drawn via \code{shuffle()} --
#' a random sampling without replacement -- from the empirical SUR residuals,
#' preserving the observed marginal residual distribution.
#'
#' @return A named numeric vector of length 2:
#' \describe{
#'   \item{\code{shape}}{The predicted Weibull shape parameter (\eqn{k > 0}).}
#'   \item{\code{scale}}{The predicted Weibull scale parameter (\eqn{\lambda > 0}).}
#' }
#'
#' @examples
#' # Predict parameters for a 10-node network with n = 500
#' weibull_parameters(nodes = 10, sample_size = 500)
#'
#' # With bootstrapped residuals for use in simulation
#' weibull_parameters(nodes = 10, sample_size = 500, bootstrap = TRUE)
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
# Updated 12.03.2026
weibull_parameters <- function(nodes, sample_size, bootstrap = FALSE)
{

  # Get Weibull model
  weibull_weights <- get(data("weibull_weights", package = "L0ggm", envir = environment()))

  # Compute descriptive parameters
  scale_parameters <- c(rlp = 1 / log(nodes))
  shape_parameters <- c(scale_parameters, beta_min = sqrt(log(nodes) / sample_size))

  # Compute Weibull parameters
  shape <- unname(
    weibull_weights$shape$coefficients[1] +
    sum(weibull_weights$shape$coefficients[-1] * shape_parameters)
  )
  scale <- unname(
    weibull_weights$scale$coefficients[1] +
    sum(weibull_weights$scale$coefficients[-1] * scale_parameters)
  )

  # Bootstrap residuals
  if(bootstrap){

    # Draw index
    index <- shuffle(seq_along(weibull_weights$shape$residuals), 1)

    shape <- shape + weibull_weights$shape$residuals[index]
    scale <- scale + weibull_weights$scale$residuals[index]
  }

  # Return parameters
  return(c(shape = shape, scale = scale))

}