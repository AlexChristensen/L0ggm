#' Weibull Descriptive Statistics from Huth et al. (2025)
#'
#' A data frame of descriptive statistics derived from 223 empirical datasets
#' collected by Huth et al. (2025). Partial correlation matrices from each
#' dataset were fit to four candidate distributions (Beta, Gamma, Log-Normal,
#' and Weibull) via maximum likelihood estimation. The Weibull distribution
#' provided the best or near-best fit across the majority of networks (see
#' Details), and its shape and scale parameters are retained here alongside
#' network-level descriptives and derived scaling indices.
#'
#' @name weibull_descriptives
#'
#' @docType data
#'
#' @usage data(weibull_descriptives)
#'
#' @format A 223 \eqn{\times} 10 data frame with the following columns:
#' \describe{
#'   \item{number}{Integer. The original index of the dataset within the full
#'   Huth et al. (2025) data collection, after removal of exact duplicates and
#'   networks outside the node range \eqn{8 \le p \le 40} or with
#'   \eqn{n \ge 300{,}000}.}
#'   \item{shape}{Numeric. The Weibull shape parameter (\eqn{k}), estimated via
#'   maximum likelihood on the absolute values of the lower-triangle partial
#'   correlations for each network.}
#'   \item{scale}{Numeric. The Weibull scale parameter (\eqn{\lambda}),
#'   estimated via maximum likelihood on the absolute values of the
#'   lower-triangle partial correlations for each network.}
#'   \item{n}{Integer. Sample size of the empirical dataset.}
#'   \item{p}{Integer. Number of nodes (variables) in the network.}
#'   \item{mean_values}{Numeric. Mean of the absolute partial correlations in
#'   the lower triangle of the partial correlation matrix.}
#'   \item{snr_values}{Numeric. Signal-to-noise ratio of the absolute partial
#'   correlations, computed as \code{mean_values / sd(|partial correlations|)}.}
#'   \item{scaling}{Numeric. A sample-size-based scaling index equivalent to
#'   the standard error unit, computed as \eqn{1 / \sqrt{n}}.}
#'   \item{rlp}{Numeric. A node-count scaling index computed as
#'   \eqn{1 / \ln(p)}, following the approach of Christensen et al. (2025).}
#'   \item{ppo}{Numeric. Parameters-per-observation ratio, computed as the
#'   number of unique edges divided by sample size:
#'   \eqn{\frac{p(p-1)/2}{n}}. Reflects the relative complexity of the
#'   network model given the available data.}
#' }
#'
#' @details
#' Empirical partial correlation matrices were obtained from Huth et al. (2025),
#' which aggregated results across a large number of published psychological
#' network studies. Prior to analysis, exact duplicate networks (identified by
#' shared citation and zero mean absolute error between partial correlation
#' matrices) were removed. Networks were further filtered to retain only those
#' with \eqn{8 \le p \le 40} nodes and \eqn{n < 300{,}000} observations.
#'
#' Absolute partial correlations from each network's lower triangle were fit
#' to four distributions via maximum likelihood: Beta, Gamma, Log-Normal, and
#' Weibull. Fits were compared by log-likelihood. The Weibull distribution
#' provided a substantially better fit (log-likelihood difference \eqn{> 2})
#' over Beta in 60 networks, over Gamma in 16, and over Log-Normal in 184,
#' while fewer than 13 networks showed a substantially better fit for any
#' alternative distribution over Weibull. Accordingly, Weibull parameters are
#' retained as the best general-purpose characterization of empirical partial
#' correlation edge weight distributions.
#'
#' Sign structure of the partial correlations was also examined. On average,
#' fewer than 0.1% of nodes required a sign flip during estimation
#' (\eqn{\bar{x} = 0.0009}, \eqn{SD = 0.008}), while approximately 35% of
#' remaining edges were negative (\eqn{\bar{x} = 0.351}, \eqn{SD = 0.085}),
#' consistent with mixed-sign partial correlation structures commonly observed
#' in psychological network data.
#'
#' @keywords datasets
#'
#' @references
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., & Marsman, M. (2025).
#' Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#'
#' @examples
#' data("weibull_descriptives")
#'
#' # Inspect structure
#' str(weibull_descriptives)
#'
#' # Distribution of Weibull shape parameters across empirical networks
#' hist(weibull_descriptives$shape, breaks = "FD",
#'      xlab = "Shape (k)", main = "Weibull Shape Parameters")
#'
#' # Relationship between scale and parameters-per-observation
#' plot(weibull_descriptives$scale, weibull_descriptives$ppo,
#'      xlab = "Number of Nodes", ylab = "Parameters per Observation")
#'
NULL
#----