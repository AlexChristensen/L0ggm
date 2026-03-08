#' Linear Model for Weibull Scale Parameter Prediction
#'
#' An \code{lm} object encoding the relationship between the Weibull scale
#' parameter and network-level descriptives derived from 223 empirical datasets
#' in \code{\link{weibull_descriptives}}. The model is intended for use in
#' generating data-driven Weibull scale estimates given a network's shape
#' parameter and structural properties.
#'
#' @name weibull_weights
#'
#' @docType data
#'
#' @usage data(weibull_weights)
#'
#' @format An object of class \code{lm} with the following model specification:
#' \preformatted{scale ~ shape + scaling + rlp + ppo}
#' Predictors are as defined in \code{\link{weibull_descriptives}}:
#' \describe{
#'   \item{shape}{Weibull shape parameter (\eqn{k}) estimated from the
#'   absolute partial correlations of each empirical network.}
#'   \item{scaling}{Sample-size scaling index, \eqn{1 / \sqrt{n}}.}
#'   \item{rlp}{Node-count scaling index, \eqn{1 / \ln(p)}.}
#'   \item{ppo}{Parameters-per-observation ratio,
#'   \eqn{p(p - 1) / 2\,/\,n}.}
#' }
#'
#' @details
#' The model was fit in two stages. In the first stage, an ordinary least
#' squares regression of the Weibull scale parameter on \code{shape},
#' \code{scaling}, \code{rlp}, and \code{ppo} was estimated on the full
#' \code{\link{weibull_descriptives}} sample. Influential observations were
#' then identified using Cook's distance with a conventional threshold of
#' \eqn{D_i > 4 / n}, and those observations were removed. The model was
#' subsequently re-estimated on the cleaned sample, which is the object
#' stored here.
#'
#' Variance inflation factors were examined at both stages to confirm the
#' absence of problematic multicollinearity among predictors.
#'
#' @keywords datasets
#'
#' @references
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., & Marsman, M. (2025).
#' Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#'
#' @examples
#' data("weibull_weights")
#'
#' # Inspect model summary
#' summary(weibull_weights)
#'
#' # Predict Weibull scale for a new network
#' # (e.g., p = 12 nodes, n = 200 observations, shape = 1)
#' p <- 12; n <- 200
#' new_data <- data.frame(
#'   shape   = 1,
#'   scaling = 1 / sqrt(n),
#'   rlp     = 1 / log(p),
#'   ppo     = (p * (p - 1) / 2) / n
#' )
#' predict(weibull_weights, newdata = new_data)
#'
NULL
#----