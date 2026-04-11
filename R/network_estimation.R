#' @title L0 Norm Regularized Network Estimation
#'
#' @description A general function to estimate Gaussian graphical models using
#' L0 penalty approximations. All penalties are implemented using
#' either single-pass or full Local Linear Approximation (LLA: Fan & Li, 2001; Zou & Li, 2008)
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param n Numeric (length = 1).
#' Sample size \strong{must} be provided if \code{data} provided is a correlation matrix
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal_categories}
#' (see \code{\link[L0ggm]{polychoric_matrix}} for more details)
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
#'
#' }
#'
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#'
#' @param na_data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param penalty Character (length = 1).
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"atan"} --- Arctangent (Wang & Zhu, 2016)
#' \deqn{\lambda \cdot (\gamma + \frac{2}{\pi}) \cdot \arctan\left(\frac{|x|}{\gamma}\right)}
#'
#' \item \code{"exp"} --- EXP (Wang, Fan, & Zhu, 2018)
#' \deqn{\lambda \cdot (1 - e^{-\frac{|x|}{\gamma}})}
#'
#' \item \code{"gumbel"} --- Gumbel
#' \deqn{\frac{\lambda}{1 - e^{-1}} \cdot \left(e^{-e^{-\frac{|x|}{\gamma}}} - e^{-1}\right)}
#'
#' \item \code{"log"} --- Log (Candes, Wakin, & Boyd, 2008)
#' \deqn{\frac{\lambda \cdot \log\left(1 + \frac{|x|}{\gamma}\right)}{\log\left(1 + \frac{1}{\gamma}\right)}}
#'
#' \item \code{"weibull"} --- Weibull
#' \deqn{\lambda \cdot \left(1 - e^{-\left(\frac{|x|}{\gamma}\right)^k}\right)}
#'
#' }
#'
#' @param gamma Numeric (length = 1).
#' Adjusts the shape of the penalty.
#' Defaults:
#'
#' \itemize{
#'
#' \item \code{"atan"} = 0.01
#'
#' \item \code{"exp"} = 0.01
#'
#' \item \code{"gumbel"} = 0.01
#'
#' \item \code{"log"} = 0.10
#'
#' \item \code{"weibull"} = 0.01
#'
#' }
#'
#' @param adaptive Boolean (length = 1).
#' Whether data-adaptive (gamma) parameters should be used.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to apply default gamma parameters for
#' adaptive penalties.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"exp"}
#'
#' \item \code{"gumbel"}
#'
#' \item \code{"weibull"}
#'
#' }
#'
#' All adaptive penalties use the 95\% confidence interval from zero
#'
#' @param nlambda Numeric (length = 1).
#' Number of lambda values to test.
#' Defaults to \code{50}
#'
#' @param lambda_min_ratio Numeric (length = 1).
#' Ratio of lowest lambda value compared to maximal lambda.
#' Defaults to \code{0.01}
#'
#' @param penalize_diagonal Boolean (length = 1).
#' Should the diagonal be penalized?
#' Defaults to \code{TRUE}
#'
#' @param ic Character (length = 1).
#' What information criterion should be used for model selection?
#' Available options include:
#'
#' \itemize{
#'
#' \item \code{"AIC"} --- Akaike's information criterion: \eqn{-2L + 2E}
#'
#' \item \code{"AICc"} --- AIC corrected: \eqn{AIC + \frac{2E^2 + 2E}{n - E - 1}}
#'
#' \item \code{"BIC"} --- Bayesian information criterion: \eqn{-2L + E \cdot \log{(n)}}
#'
#' \item \code{"BIC0"} --- Bayesian information criterion not (Dicker et al., 2013): \eqn{\log{\large(\frac{D}{n - E}\large)} + \large(\frac{\log{(n)}}{n}\large) \cdot E}
#'
#' \item \code{"EBIC"} --- Extended BIC: \eqn{BIC + 4E \cdot \gamma \cdot \log{(E)}}
#'
#' \item \code{"MBIC"} --- Modified Bayesian information criterion (Wang et al., 2018):  \eqn{\log{\large(\frac{D}{n - E}\large)} + \large(\frac{\log{(n)} \cdot E}{n}\large) \cdot \log{(\log{(p)}})}
#'
#' }
#'
#' Term definitions:
#'
#' \itemize{
#'
#' \item \eqn{n} --- sample size
#'
#' \item \eqn{p} --- number of variables
#'
#' \item \eqn{E} --- edges
#'
#' \item \eqn{S} --- empirical correlation matrix
#'
#' \item \eqn{K} --- estimated inverse covariance matrix (network)
#'
#' \item \eqn{L = \frac{n}{2} \cdot \log \text{det} K - \sum_{i=1}^p (SK)_{ii}}
#'
#' \item \eqn{D = n \cdot \sum_{i=1}^p (SK)_{ii} - \log \text{det} K}
#'
#' }
#'
#' Defaults to \code{"BIC"}
#'
#' @param ebic_gamma Numeric (length = 1)
#' Value to set gamma parameter in EBIC (see above).
#' Defaults to \code{0.50}
#'
#' \emph{Only used if \code{ic = "EBIC"}}
#'
#' @param fast Boolean (length = 1).
#' Whether the \code{\link[glassoFast]{glassoFast}} version should be used
#' to estimate the GLASSO.
#' Defaults to \code{TRUE}.
#'
#' The fast results \emph{may} differ by less than floating point of the original
#' GLASSO implemented by \code{\link[glasso]{glasso}} and should not impact reproducibility much (set to \code{FALSE} if concerned)
#'
#' @param LLA Boolean (length = 1).
#' Should Local Linear Approximation be used to find optimal minimum?
#' Defaults to \code{FALSE} or a single-pass approximation, which can be
#' significantly faster (Zou & Li, 2008).
#' Set to \code{TRUE} to find global minimum based on convergence (\code{LLA_threshold})
#'
#' @param LLA_threshold Numeric (length = 1).
#' When performing the Local Linear Approximation, the maximum threshold
#' until convergence is met.
#' Defaults to \code{1e-03}
#'
#' @param LLA_iter Numeric (length = 1).
#' Maximum number of iterations to perform to reach convergence.
#' Defaults to \code{10000}
#'
#' @param network_only Boolean (length = 1).
#' Whether the network only should be output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to obtain all output for the
#' network estimation method
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to \code{auto.correlate}
#'
#' @return When \code{network_only = TRUE} (default), returns a \eqn{p \times p}
#' numeric matrix of partial correlations representing the estimated Gaussian
#' graphical model. Off-diagonal entry \eqn{[i,j]} is the partial correlation
#' between variables \eqn{i} and \eqn{j} controlling for all other variables,
#' with values in \eqn{[-1, 1]}; a value of zero indicates the absence of an
#' edge. Diagonal entries are zero. Row and column names are inherited from
#' \code{data}.
#'
#' When \code{network_only = FALSE}, returns a named list with the following
#' elements:
#'
#' \describe{
#'
#' \item{\code{network}}{The \eqn{p \times p} partial correlation matrix
#' described above}
#'
#' \item{\code{K}}{The \eqn{p \times p} estimated inverse covariance
#' (precision) matrix at the optimal lambda. Diagonal entries are the
#' conditional precisions; off-diagonal entries are proportional to
#' partial covariances}
#'
#' \item{\code{R}}{The \eqn{p \times p} regularized covariance matrix
#' returned by GLASSO at the optimal lambda (the \code{w} component of the
#' GLASSO solution)}
#'
#' \item{\code{penalty}}{Character string naming the penalty function used
#' (one of \code{"atan"}, \code{"exp"}, \code{"gumbel"}, \code{"log"},
#' \code{"weibull"})}
#'
#' \item{\code{lambda}}{Numeric scalar giving the regularization parameter
#' value selected by the information criterion. Larger values correspond
#' to sparser networks}
#'
#' \item{\code{gamma}}{Numeric scalar giving the shape parameter of the
#' penalty actually used. For adaptive penalties (\code{adaptive = TRUE}),
#' this is the data-derived value; otherwise it is the default or
#' user-supplied value}
#'
#' \item{\code{correlation}}{The \eqn{p \times p} empirical correlation
#' matrix computed from \code{data}, used as input to the GLASSO}
#'
#' \item{\code{criterion}}{Character string naming the information criterion
#' used for model selection (e.g., \code{"bic"})}
#'
#' \item{\code{IC}}{Numeric scalar giving the value of the information
#' criterion at the optimal lambda}
#'
#' }
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#'
#' \strong{Log penalty} \cr
#' Candes, E. J., Wakin, M. B., & Boyd, S. P. (2008).
#' Enhancing sparsity by reweighted l1 minimization.
#' \emph{Journal of Fourier Analysis and Applications}, \emph{14}(5), 877--905.
#'
#' \strong{BIC0} \cr
#' Dicker, L., Huang, B., & Lin, X. (2013).
#' Variable selection and estimation with the seamless-L0 penalty.
#' \emph{Statistica Sinica}, \emph{23}(2), 929--962.
#'
#' \strong{Local Linear Approximation} \cr
#' Fan, J., & Li, R. (2001).
#' Variable selection via nonconcave penalized likelihood and its oracle properties.
#' \emph{Journal of the American Statistical Association}, \emph{96}(456), 1348--1360.
#'
#' \strong{EXP penalty} \cr
#' Wang, Y., Fan, Q., & Zhu, L. (2018).
#' Variable selection and estimation using a continuous approximation to the L0 penalty.
#' \emph{Annals of the Institute of Statistical Mathematics}, \emph{70}(1), 191--214.
#'
#' \strong{Atan penalty} \cr
#' Wang, Y., & Zhu, L. (2016).
#' Variable selection and parameter estimation with the Atan regularization method.
#' \emph{Journal of Probability and Statistics}, \emph{2016}, 1--12.
#'
#' \strong{Seminal simulation in network psychometrics} \cr
#' Williams, D. R. (2020).
#' Beyond lasso: A survey of nonconvex regularization in Gaussian graphical models.
#' \emph{PsyArXiv}.
#'
#' \strong{One-step Local Linear Approximation} \cr
#' Zou, H., & Li, R. (2008).
#' One-step sparse estimates in nonconcave penalized likelihood models.
#' \emph{Annals of Statistics}, \emph{36}(4), 1509--1533.
#'
#' @examples
#' # Obtain default estimator (adaptive Weibull)
#' weibull_network <- network_estimation(
#'   data = basic_smallworld, LLA = TRUE
#' )
#'
#' # Obtain Atan network
#' atan_network <- network_estimation(
#'   data = basic_smallworld, penalty = "atan", LLA = TRUE
#' )
#'
#' # Obtain static EXP network
#' exp_network <- network_estimation(
#'   data = basic_smallworld, penalty = "exp",
#'   adaptive = FALSE, LLA = TRUE
#' )
#'
#' @export
#'
# Apply non-convex regularization ----
# Updated 08.03.2026
network_estimation <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na_data = c("pairwise", "listwise"),
    penalty = c("atan", "exp", "gumbel", "log", "weibull"),
    gamma = NULL, adaptive = TRUE,
    nlambda = 50, lambda_min_ratio = 0.01, penalize_diagonal = TRUE,
    ic = c("AIC", "AICc", "BIC", "BIC0", "EBIC", "MBIC"), ebic_gamma = 0.50,
    fast = TRUE, LLA = FALSE, LLA_threshold = 1e-03, LLA_iter = 10000,
    network_only = TRUE, verbose = FALSE, ...
)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  # (keeping non-function choices for `cor_auto`)
  corr <- set_default(corr, "auto", network_estimation)
  na_data <- set_default(na_data, "pairwise", network_estimation)
  penalty <- set_default(penalty, "weibull", network_estimation)
  ic <- set_default(ic, "bic", network_estimation)

  # Argument errors (return data in case of tibble)
  data <- network_regularization_errors(
    data, n, gamma, adaptive, nlambda,
    lambda_min_ratio, penalize_diagonal, ebic_gamma,
    fast, LLA, LLA_threshold, network_only, verbose, ...
  )

  # Check whether penalty is adaptive option
  adaptive_option <- c("exp", "gumbel", "weibull")
  adaptive_flag <- adaptive & (penalty %in% adaptive_option)

  # Get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = n,
    corr = corr, na_data = na_data,
    verbose = verbose, needs_usable = FALSE, # skips usable data check
    ...
  )

  # Get outputs
  data <- output$data; n <- output$n
  S <- output$correlation_matrix

  # Get number of variables
  nodes <- dim(S)[2]

  # Obtain precision matrix
  K <- solve(S)

  # Obtain GLASSO function
  glasso_FUN <- swiftelse(fast, glassoFast::glassoFast, glasso::glasso)

  # Get function arguments
  glasso_ARGS <- obtain_arguments(glasso_FUN, FUN.args = list(...))

  # Supply correlation matrix
  glasso_ARGS[[1]] <- S

  # Get derivative function
  derivative_FUN <- switch(
    penalty,
    "atan" = atan_derivative,
    "exp" = exp_derivative,
    "gumbel" = gumbel_derivative,
    "log" = log_derivative,
    "weibull" = weibull_derivative
  )

  # Set shape (only used for Weibull)
  shape <- 1

  # Set gamma (set ahead of time for messaging on adaptive)
  if(is.null(gamma)){

    # Set defaults
    gamma <- switch(
      penalty,
      "atan" = 0.01,
      "exp" = 0.01,
      "gumbel" = 0.01,
      "log" = 0.10,
      "weibull" = 0.01
    )

  }

  # Check for adaptive and penalty is adaptive option
  if(adaptive_flag){

    # Set lower triangle
    lower_triangle <- lower.tri(S)

    # Set partial correlations
    P <- inv2pcor(K); lower_P <- abs(P[lower_triangle])

    if(penalty == "exp"){

      # Set 10th percentile
      gamma <- -log(0.90) * mean(lower_P)

    }else if(penalty == "gumbel"){

      # Set 0.4065697 percentile
      # (makes equivalent to EXP and Weibull when scale = 0.10)
      gamma <- -gumbel_mle(lower_P) * log(-log(exp(-0.90)))

    }else if(penalty == "weibull"){

      # Obtain Weibull estimates
      estimates <- weibull_mle(lower_P)

      # Set 10th percentile
      gamma <- estimates[["scale"]] * (-log(0.90))^(1 / estimates[["shape"]])

    }

  }

  # Initialize lambda matrix
  lambda_matrix <- matrix(0, nrow = nodes, ncol = nodes)

  # Simplify source for fewer computations (minimal improvement)
  S_zero_diagonal <- S - diag(nodes) # makes diagonal zero
  lambda_max <- max(abs(S_zero_diagonal)) # uses absolute rather than inverse
  lambda_min <- lambda_max * lambda_min_ratio
  lambda_sequence <- exp(seq.int(log(lambda_min), log(lambda_max), length.out = nlambda))

  # Get GLASSO output
  glasso_list <- lapply(lambda_sequence, function(lambda){

    # Set lambda matrix
    glasso_ARGS$rho <- lambda

    # Obtain estimate
    estimate <- do.call(what = glasso_FUN, args = glasso_ARGS)

    # Check for LLA
    if(LLA){

      # Obtain new K
      new_K <- estimate$wi

      # Set convergence and iterations
      convergence <- Inf; iterations <- 0

      # Loop over to convergence
      while((convergence > LLA_threshold) & (iterations < LLA_iter)){

        # Set old K
        old_K <- new_K

        # Obtain lambda matrix
        lambda_matrix[] <- derivative_FUN(x = old_K, lambda = lambda, gamma = gamma, shape = shape)

        # Check for diagonal penalization
        if(!penalize_diagonal){
          diag(lambda_matrix) <- 0
        }

        # Set lambda matrix
        glasso_ARGS$rho <- lambda_matrix

        # Obtain estimate
        estimate <- do.call(what = glasso_FUN, args = glasso_ARGS)

        # Obtain new K
        new_K <- estimate$wi

        # Increase iterations
        iterations <- iterations + 1

        # Compute convergence
        convergence <- max(abs(new_K - old_K))

      }

    }else{

      # Obtain lambda matrix
      lambda_matrix[] <- derivative_FUN(x = estimate$wi, lambda = lambda, gamma = gamma, shape = shape)

      # Check for diagonal penalization
      if(!penalize_diagonal){
        diag(lambda_matrix) <- 0
      }

      # Set lambda matrix
      glasso_ARGS$rho <- lambda_matrix

      # Obtain estimate
      estimate <- do.call(what = glasso_FUN, args = glasso_ARGS)

    }

    # Estimate
    return(estimate)

  })

  # Compute ICs
  ICs <- nvapply(glasso_list, function(element){
    information_criterion(
      S = S, K = element$wi, n = n, nodes = nodes,
      ic = ic, ebic_gamma = ebic_gamma
    )
  })

  # Optimal value
  optimal <- which.min(ICs)

  # Get R
  R <- glasso_list[[optimal]]$w

  # Get W
  W <- inv2pcor(glasso_list[[optimal]]$wi)
  dimnames(R) <- dimnames(W) <- dimnames(S)

  # Return results
  if(network_only){
    return(W)
  }else{
    return(
      list(
        network = W, K = glasso_list[[optimal]]$wi, R = R,
        penalty = penalty, lambda = lambda_sequence[[optimal]], gamma = gamma,
        correlation = S, criterion = ic, IC = ICs[[optimal]]
      )
    )
  }

}

# Bug checking ----
# data = basic_smallworld; n = NULL; corr = "auto"
# na_data = "pairwise"; penalty = "weibull"; adaptive = TRUE
# gamma = NULL; lambda = NULL; nlambda = 50
# lambda_min_ratio = 0.01; penalize_diagonal = TRUE
# optimize.lambda = FALSE; ic = "BIC"; network_only = TRUE
# ebic_gamma = 0.5; fast = TRUE; verbose = FALSE
# LLA = TRUE; LLA_threshold = 1e-04; LLA_iter = 10000

#' @noRd
# Errors ----
# Updated 19.01.2026
network_regularization_errors <- function(
    data, n, gamma, adaptive, nlambda,
    lambda_min_ratio, penalize_diagonal, ebic_gamma,
    fast, LLA, LLA_threshold, network_only, verbose, ...
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "network_estimation")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "network_estimation")
    typeof_error(n, "numeric", "network_estimation")
  }

  # 'gamma' errors
  if(!is.null(gamma)){
    length_error(gamma, 1, "network_estimation")
    typeof_error(gamma, "numeric", "network_estimation")
    range_error(gamma, c(0, Inf), "network_estimation")
  }

  # 'adaptive' errors
  length_error(adaptive, 1, "network_estimation")
  typeof_error(adaptive, "logical", "network_estimation")

  # 'nlambda' errors
  length_error(nlambda, 1, "network_estimation")
  typeof_error(nlambda, "numeric", "network_estimation")
  range_error(nlambda, c(1, Inf), "network_estimation")

  # 'lambda_min_ratio' errors
  length_error(lambda_min_ratio, 1, "network_estimation")
  typeof_error(lambda_min_ratio, "numeric", "network_estimation")
  range_error(lambda_min_ratio, c(0, 1), "network_estimation")

  # 'penalize_diagonal' errors
  length_error(penalize_diagonal, 1, "network_estimation")
  typeof_error(penalize_diagonal, "logical", "network_estimation")

  # 'ebic_gamma' errors
  length_error(ebic_gamma, 1, "network_estimation")
  typeof_error(ebic_gamma, "numeric", "network_estimation")
  range_error(ebic_gamma, c(0, Inf), "network_estimation")

  # 'fast' errors
  length_error(fast, 1, "network_estimation")
  typeof_error(fast, "logical", "network_estimation")

  # 'LLA' errors
  length_error(LLA, 1, "network_estimation")
  typeof_error(LLA, "logical", "network_estimation")

  # 'LLA_threshold' errors
  if(LLA){
    length_error(LLA_threshold, 1, "network_estimation")
    typeof_error(LLA_threshold, "numeric", "network_estimation")
    range_error(LLA_threshold, c(-Inf, 0.10), "network_estimation")
  }

  # 'network_only' errors
  length_error(network_only, 1, "network_estimation")
  typeof_error(network_only, "logical", "network_estimation")

  # 'verbose' errors
  length_error(verbose, 1, "network_estimation")
  typeof_error(verbose, "logical", "network_estimation")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return data in case of tibble
  return(data)

}

#' @noRd
# Information criterion ----
# Updated 09.03.2026
information_criterion <- function(S, K, n, nodes, ic, ebic_gamma)
{

  # Compute Gaussian likelihood (minus two for convenience)
  trSK <- sum(S * K)
  log_det_K <- log(det(K))
  L <- swiftelse(
    ic %in% c("bic0", "mbic"),
    n * trSK - log_det_K, # use deviance
    -2 * (n / 2) * (log_det_K - trSK) # use log-likelihood
  )

  # Set diagonal of K to zero
  diag(K) <- 0

  # Get parameters (edges)
  E <- edge_count(K, nodes)

  # Ensure that there is enough degrees of freedom
  df <- n - E

  # Return information criterion
  return(
    switch(
      ic,
      "aic" = L + 2 * E,
      "aicc" = swiftelse(n - E - 1 > 0, L + 2 * E + (2 * E^2 + 2 * E) / (n - E - 1), Inf),
      "bic" = L + E * log(n),
      "ebic" = L + E * log(n) + 4 * E * ebic_gamma * log(nodes),
      "bic0" = swiftelse( # see https://www.jstor.org/stable/24310368
        df > 0, log(L / df) + (log(n) * E / n), -Inf
      ),
      "mbic" = swiftelse( # see https://doi.org/10.1007%2Fs10463-016-0588-3
        df > 0, log(L / df) + (log(n) * E / n) * log(log(nodes)), -Inf
      )
    )
  )

}
