#' SUR Model Coefficients and Residuals for Weibull Parameter Prediction
#'
#' A named list encoding the Seemingly Unrelated Regression (SUR) model fitted
#' to Weibull shape and scale parameters derived from 194
#' empirical networks. The object is consumed internally by \code{\link{weibull_parameters}}
#' to generate data-driven Weibull parameter estimates given a network's
#' number of nodes and sample size.
#'
#' @name weibull_weights
#'
#' @docType data
#'
#' @usage data(weibull_weights)
#'
#' @format A named list with two elements, \code{shape} and \code{scale},
#' each containing:
#' \describe{
#'   \item{\code{coefficients}}{A named numeric vector of regression
#'   coefficients from the SUR model. Both equations include an intercept
#'   and three predictors (\code{snr}, \code{rlp}, \code{beta_min}).
#'   Predictors are defined as follows:
#'   \describe{
#'     \item{\code{snr}}{Signal-to-noise ratio of the absolute partial
#'     correlations, computed as the mean divided by the standard deviation
#'     of the absolute edge weights.}
#'     \item{\code{rlp}}{Reciprocal log of node count, \eqn{1 / \log(p)},
#'     where \eqn{p} is the number of nodes.}
#'     \item{\code{beta_min}}{A noise-to-signal proxy, \eqn{\sqrt{\log(p)/n}},
#'     where \eqn{p} is the number of nodes and \eqn{n} is the sample size.}
#'   }}
#'   \item{\code{residuals}}{A numeric vector of residuals from the fitted SUR
#'   equation, used by \code{\link{weibull_parameters}} to introduce
#'   empirically grounded variability when \code{bootstrap = TRUE}.}
#' }
#'
#' @details
#' Absolute partial correlations from 222 deduplicated empirical networks
#' (Huth et al., 2025) were fitted to Beta, Gamma, log-normal, and Weibull
#' distributions via maximum likelihood. Weibull provided the best fit most
#' consistently: it outperformed each alternative by more than 2 log-likelihood
#' units far more often than the reverse (vs. Beta: 56--0; vs. Gamma: 15--2;
#' vs. log-normal: 155--13; vs. Exponential: 54--0).
#'
#' The resulting Weibull shape and scale parameters were then jointly modelled
#' as a function of network descriptors using Seemingly Unrelated Regression
#' (\code{systemfit}), which accounts for correlated residuals between the
#' shape and scale equations across networks. Prior to fitting, networks with
#' fewer than eight nodes (\eqn{p < 8}), more than 300,000 observations, or
#' fewer than one observation per edge (\eqn{\text{ope} \leq 1}) were excluded,
#' as Huth et al. (2025) demonstrated that networks in this regime show the
#' weakest statistical evidence for edge presence or absence, yielding unstable
#' parameter estimates (n = 28 excluded). This left n = 194 networks for
#' analysis. Shape and scale parameters were modelled on their original scales.
#'
#' Variance inflation factors (\eqn{\leq} 1.62) and the condition number of the
#' design matrix (kappa = 61.1) confirmed the absence of problematic
#' multicollinearity. Breusch-Pagan tests indicated statistically significant
#' heteroskedasticity in both equations; however, all predictors remained
#' significant under HC3-robust standard errors, indicating no material effect
#' on inference. Shape residuals were normally distributed (Shapiro-Wilk
#' W = 0.995, p = 0.814). Scale residuals showed a modest departure from
#' normality (Shapiro-Wilk W = 0.973, p < 0.001), consistent with slight right
#' skew in the scale outcome and test sensitivity at n = 194 rather than a
#' substantive violation, as confirmed by visual inspection of the residual
#' histogram.
#'
#' @keywords datasets
#'
#' @references
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., &
#' Marsman, M. (2025). Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#'
#' @examples
#' data("weibull_weights")
#'
#' # Inspect SUR coefficients for each equation
#' weibull_weights$shape$coefficients
#' weibull_weights$scale$coefficients
#'
#' # Predict Weibull parameters for a new network
#' weibull_parameters(nodes = 12, sample_size = 500)
#'
NULL
#----