#' Simulates Stochastic Block Model Data
#'
#' @description
#' Simulates data from a Gaussian Graphical Model (GGM) with a stochastic
#' block model (SBM) structure. Nodes are partitioned into communities, with
#' edge density controlled separately within and between communities via a
#' \code{blocks x blocks} density matrix. Edge weights are derived empirically
#' from a Weibull distribution fit to real psychometric network data, and the
#' resulting network is used to generate multivariate normal data. A
#' \code{diffusion} parameter controls the proportion of the strongest edges
#' that are reserved for within-community positions versus allowed to appear
#' between communities.
#' Parameters do not have default values (except \code{negative_proportion},
#' \code{diffusion}, \code{target_condition}, and \code{max_iterations}) and
#' must each be set. See Details and Examples to get started.
#'
#' @param nodes Numeric (length = 1 or \code{blocks}).
#' Number of nodes per community block.
#' Can be a single value (applied to all blocks) or as many values as there
#' are blocks (corresponding to each block in order).
#' Minimum of three nodes per block.
#'
#' @param blocks Numeric (length = 1).
#' Number of community blocks.
#'
#' @param density_matrix Matrix (dim = \code{blocks x blocks}).
#' A symmetric numeric matrix specifying edge probabilities within and between
#' community blocks. Diagonal entries give the within-block edge probability
#' for each community; off-diagonal entries give the between-block edge
#' probability for each pair of communities. All values must be between 0
#' and 1, and the matrix must be symmetric (i.e., \code{density_matrix[i,j] ==
#' density_matrix[j,i]}). See Details for construction guidance.
#'
#' @param snr Numeric (length = 1).
#' Signal-to-noise ratio of partial correlations
#' (\eqn{\bar{|w|} / \mathrm{SD}(|w|)}).
#' Values less than 1 indicate a wider spread of partial correlation weights
#' whereas values greater than 1 indicate a narrower spread.
#' Defaults to \code{1} where the mean equals the standard deviation of the
#' absolute partial correlations.
#'
#' @param negative_proportion Numeric (length = 1).
#' Proportion of between-community edges that are negative (inhibitory).
#' Must be between 0 and 1.
#' If not provided, a value is sampled from an empirical distribution with
#' mean = 0.35 and SD = 0.09, bounded between 0.08 and 0.55.
#'
#' @param diffusion Numeric (length = 1).
#' Proportion of the strongest edges that are randomly assigned to either
#' within-block or between-block positions rather than being reserved for
#' within-block positions only. At \code{diffusion = 0.50} (default), 50\%
#' of the strongest edges are placed within communities. Higher values
#' redistribute more top-weighted edges across the network, weakening the
#' community contrast. Must be between 0 and 1.
#'
#' @param diffusion_range Numeric (length = 2).
#' If provided, the diffusion proportion is drawn uniformly from this range
#' on each call, overriding \code{diffusion}. Useful for introducing
#' variability across simulation replications. For example,
#' \code{diffusion_range = c(0.05, 0.20)} samples a diffusion value between
#' 5\% and 20\% on each draw. Both values must be between 0 and 1.
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases to generate from the population multivariate normal
#' distribution.
#'
#' @param skew Numeric (length = 1 or \code{sum(nodes)}).
#' Skew applied to each variable. Can be a single value (applied to all
#' variables) or one value per variable. Values are rounded to the nearest
#' increment of 0.05 in the range [-2, 2].
#'
#' @param skew_range Numeric (length = 2).
#' If provided, a skew value is drawn uniformly from this range for each
#' variable, overriding \code{skew}. Both values must be between -2 and 2.
#'
#' @param target_condition Numeric (length = 1).
#' Target condition number (using \code{\link{kappa}} with \code{exact = TRUE})
#' applied when ridge regularization is needed to recover a positive definite
#' precision matrix. Lower values produce better-conditioned matrices.
#' Defaults to \code{30}. Values up to \code{100} are accepted but not
#' recommended.
#'
#' @param max_correlation Numeric (length = 1).
#' Maximum allowed absolute correlation between any pair of nodes in the
#' population correlation matrix \code{R}. Draws exceeding this threshold
#' are rejected and a new attempt is made. Must be between 0 and 1.
#' Defaults to \code{0.80}.
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of attempts to find a connected network structure with a
#' valid set of edge weights. Defaults to \code{100}.
#'
#' @details
#' \strong{Constructing \code{density_matrix}}
#'
#' The \code{density_matrix} is a \code{blocks x blocks} symmetric matrix
#' where entry \eqn{[i,i]} gives the within-block edge probability for
#' community \eqn{i}, and entry \eqn{[i,j]} (for \eqn{i \neq j}) gives the
#' between-block edge probability for communities \eqn{i} and \eqn{j}. The
#' simplest construction uses a uniform off-diagonal density with
#' block-specific diagonals:
#'
#' \preformatted{
#' # Uniform within (0.90) and between (0.20) density for 3 blocks
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- 0.90
#' }
#'
#' For asymmetric community structure, each diagonal entry can differ:
#'
#' \preformatted{
#' # Varying within-block density per community
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- c(0.85, 0.90, 0.95)
#' }
#'
#' For full pairwise control over between-block densities, specify the
#' complete symmetric matrix directly:
#'
#' \preformatted{
#' dm <- matrix(c(
#'   0.90, 0.20, 0.10,
#'   0.20, 0.85, 0.30,
#'   0.10, 0.30, 0.95
#' ), nrow = 3, ncol = 3)
#' }
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data matrix (\code{sample_size x sum(nodes)}) drawn
#' from the population GGM.}
#'
#' \item{parameters}{
#' A list of input, derived, and estimated parameters:
#' \itemize{
#'   \item \code{nodes} --- Number of nodes per block (expanded to length
#'   \code{blocks})
#'   \item \code{blocks} --- Number of community blocks
#'   \item \code{sample_size} --- Number of simulated cases
#'   \item \code{skew} --- Named numeric vector of per-variable skew values
#'   actually applied
#'   \item \code{density_matrix} --- The \code{blocks x blocks} density matrix
#'   as supplied
#'   \item \code{negative_proportion} --- Proportion of negative between-block
#'   edges used
#'   \item \code{weibull} --- Weibull \code{shape} and \code{scale} parameters
#'   of the absolute edge weight distribution
#'   \item \code{diffusion} --- Numeric vector of length 2 giving the
#'   within-block reservation range. Both values equal \code{1 - diffusion}
#'   when \code{diffusion} is scalar; equals \code{1 - diffusion_range} when
#'   \code{diffusion_range} is provided
#'   \item \code{Q} --- Newman-Girvan modularity of the population network
#'   with respect to the block membership, computed via
#'   \code{igraph::modularity}
#' }
#' }
#'
#' \item{population}{
#' Population-level network parameters:
#' \itemize{
#'   \item \code{R} --- Population correlation matrix derived from the GGM
#'   \item \code{Omega} --- Population partial correlation matrix (GGM edge
#'   weights)
#'   \item \code{membership} --- Named integer vector of community block
#'   assignments for each node
#' }
#' }
#'
#' \item{convergence}{
#' Iteration and conditioning diagnostics:
#' \itemize{
#'   \item \code{iterations} --- Number of iterations needed to find a valid
#'   network
#'   \item \code{rejections} --- Character vector of rejection reasons across
#'   iterations, useful for diagnosing convergence failures
#'   \item \code{lambda} --- Ridge regularization parameter applied to
#'   condition the precision matrix; \code{NA} if no conditioning was required
#'   \item \code{condition} --- Condition number of the population correlation
#'   matrix \code{R}
#' }
#' }
#'
#' @examples
#' # Construct density matrix for 3 blocks with uniform densities
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- 0.90
#'
#' # Basic 3-block simulation with equal-sized communities
#' result <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block = 18 total
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm
#' )
#'
#' # Unequal block sizes
#' result <- simulate_sbm(
#'   nodes = c(4, 6, 8), # 18 total nodes
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm
#' )
#'
#' # Varying within-block density per community
#' dm_varying <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm_varying) <- c(0.85, 0.90, 0.95)
#'
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm_varying
#' )
#'
#' # Full pairwise between-block density control
#' dm_pairwise <- matrix(c(
#'   0.90, 0.20, 0.10,
#'   0.20, 0.85, 0.30,
#'   0.10, 0.30, 0.95
#' ), nrow = 3, ncol = 3)
#'
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm_pairwise
#' )
#'
#' # Fix the proportion of negative edges
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm,
#'   negative_proportion = 0.20
#' )
#'
#' # Introduce variability in diffusion across replications
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm,
#'   diffusion_range = c(0.30, 0.70)
#' )
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{Seminal introduction to Stochastic Block Models} \cr
#' Holland, P. W., Laskey, K. B., & Leinhardt, S. (1983).
#' Stochastic blockmodels: First steps.
#' \emph{Social Networks}, \emph{5}(2), 109--137.
#'
#' @export
#'
# Simulate SBM GGM data ----
# Updated 14.03.2026
simulate_sbm <- function(
    nodes, blocks, density_matrix,
    snr = 1, diffusion = 0.50, diffusion_range = NULL,
    negative_proportion, sample_size,
    skew = 0, skew_range = NULL,
    target_condition = 30, max_correlation = 0.80,
    max_iterations = 100
)
{

  # Check for missing negative proportion
  if(missing(negative_proportion)){

    # Compute proportion of negative edges based on empirical values
    negative_proportion <- pmin(
      pmax(
        0.34 + rnorm_ziggurat(1) * 0.086, # empirical mean +/- 1 SD
        0.083 # empirical minimum
      ), 0.55 # empirical maximum
    )

  }

  # Check for input errors
  ## Returns updated nodes to ensure same length as blocks
  nodes <- simulate_sbm_errors(
    nodes, blocks, density_matrix,
    snr, negative_proportion, sample_size, skew, diffusion,
    target_condition, max_correlation, max_iterations
  )

  # Determine total number of nodes
  total_nodes <- sum(nodes)

  # Get community sequence
  community_sequence <- seq_len(blocks)

  # Set membership
  membership <- sort(rep(community_sequence, times = nodes))

  # Set up membership matrix
  membership_matrix <-  outer(membership, membership, "==")

  # Obtain lower triangle
  lower_triangle <- lower.tri(membership_matrix)

  # Initialize iterations
  iter <- 0

  # Initialize condition
  search <- TRUE

  # Collect rejections
  rejections <- character(length = max_iterations + 1)

  # Generation loop
  while(search){

    # Increase iterations
    iter <- iter + 1

    # Stop on greater than max iterations
    if(iter > max_iterations){

      # Collect rejections
      rejection_table <- fast_table(rejections)

      # Remove NULL
      rejection_table <- rejection_table[names(rejection_table) != ""]

      # Return error
      stop(
        paste0(
          "Reached maximum iterations: Could not find weights that for the SBM structure. ",
          "Resulting rejections were due to:\n\n",
          paste0(names(rejection_table), " = ", rejection_table, collapse = "\n")
        )
      )

    }

    # Get block matrix
    block_matrix <- generate_sbm(
      membership_matrix, blocks, membership,
      membership_matrix, density_matrix, lower_triangle
    )

    # Check block matrix
    if(!igraph::is_connected(convert2igraph(block_matrix))){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where every node was connected."

      # Move to next iteration
      next

    }

    # Try to get good weights
    output <- try(
      sbm_weights(
        block_matrix, membership, membership_matrix, sample_size,
        snr, total_nodes, negative_proportion,
        diffusion, diffusion_range, target_condition,
        max_correlation, lower_triangle
      ), silent = TRUE
    )

    # Check for error
    if(is(output, "try-error")){

      # Add rejection reason
      rejections[iter] <- trimws(strsplit(output[1], split = "\n")[[1]][2])

      # Move to next iteration
      next

    }

    # Update search
    search <- FALSE

  }

  # Generate data
  data_output <- simulate_data(
    n = sample_size, p = total_nodes, R = output$R, skew = skew, skew_range = skew_range
  )

  # Return parameters
  return(
    list(
      data = data_output$data,
      parameters = list(
        nodes = nodes,
        blocks = blocks,
        sample_size = sample_size,
        skew = data_output$skew,
        density_matrix = density_matrix,
        negative_proportion = negative_proportion,
        weibull = output$params,
        diffusion = output$diffusion,
        Q = igraph::modularity(convert2igraph(abs(output$network)), output$membership)
      ),
      population = list(
        R = output$R,
        Omega = output$network,
        membership = output$membership
      ),
      convergence = list(
        iterations = iter,
        rejections = rejections,
        lambda = output$lambda,
        condition = output$condition
      )
    )
  )

}

# Bug checking ----
# nodes = 6; blocks = 3; sample_size = 1000
# density_matrix = matrix(0.20, nrow = blocks, ncol = blocks)
# diag(density_matrix) <- 0.90
# snr = 1; negative_proportion = 0.35
# skew = 0; skew_range = NULL
# diffusion = 0.05; diffusion_range = NULL
# target_condition = 10;
# max_correlation = 0.80; max_iterations = 100

#' @noRd
# Errors ----
# Updated 14.03.2026
simulate_sbm_errors <- function(
    nodes, blocks, density_matrix,
    snr, negative_proportion, sample_size, skew, diffusion,
    target_condition, max_correlation, max_iterations
)
{

  # Errors for 'blocks'
  typeof_error(blocks, "numeric", "simulate_sbm")
  length_error(blocks, 1, "simulate_sbm")
  range_error(blocks, c(1, Inf), "simulate_sbm")

  # Errors for 'nodes'
  typeof_error(nodes, "numeric", "simulate_sbm")
  length_error(nodes, c(1, blocks), "simulate_sbm")
  range_error(nodes, c(3, Inf), "simulate_sbm")

  # Ensure 'nodes' is the length of blocks
  if(length(nodes) != blocks){
    nodes <- rep(nodes, blocks)
  }

  # Object error for 'density_matrix'
  object_error(density_matrix, c("matrix", "data.frame", "tibble"), "simulate_sbm")

  # Check for matrix
  if(get_object_type(density_matrix) != "matrix"){
    density_matrix <- as.matrix(as.data.frame(density_matrix))
  }

  # Errors for 'density_matrix'
  typeof_error(density_matrix, "numeric", "simulate_sbm")
  length_error(density_matrix, blocks * blocks, "simulate_sbm")
  range_error(density_matrix, c(0, 1), "simulate_sbm")

  # Symmetric error for 'density_matrix'
  if(!is_symmetric(density_matrix)){
    stop("`simulate_sbm` can only handle symmetric 'density_matrix' input.")
  }

  # Error for 'snr'
  typeof_error(snr, "numeric", "simulate_sbm")
  length_error(snr, 1, "simulate_sbm")
  range_error(snr, c(0.01, Inf), "simulate_sbm")

  # Send warning if beyond bounds
  if((snr < 0.64) | (snr > 1.72)){
    warning(paste0(
      "The signal-to-ratio specified (", round(snr, 2), ") is beyond ",
      "the bounds of cataloged empirical data (0.64-1.72). "
    ))
  }

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric", "simulate_sbm")
  length_error(negative_proportion, 1, "simulate_sbm")
  range_error(negative_proportion, c(0, 1), "simulate_sbm")

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric", "simulate_sbm")
  length_error(sample_size, 1, "simulate_sbm")
  range_error(sample_size, c(1, Inf), "simulate_sbm")

  # Errors for 'skew'
  typeof_error(skew, "numeric", "simulate_sbm")
  length_error(skew, c(1, sum(nodes)), "simulate_sbm")
  range_error(skew, c(-2, 2), "simulate_sbm")

  # Errors for 'diffusion'
  typeof_error(diffusion, "numeric", "simulate_sbm")
  length_error(diffusion, 1, "simulate_sbm")
  range_error(diffusion, c(0, 1), "simulate_sbm")

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric", "simulate_sbm")
  length_error(target_condition, 1, "simulate_sbm")
  range_error(target_condition, c(1, 100), "simulate_sbm")

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric", "simulate_sbm")
  length_error(max_correlation, 1, "simulate_sbm")
  range_error(max_correlation, c(0, 1), "simulate_sbm")

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric", "simulate_sbm")
  length_error(max_iterations, 1, "simulate_sbm")
  range_error(max_iterations, c(1, Inf), "simulate_sbm")

  # Return updated nodes
  return(nodes)

}

#' @noRd
# Generate SBM ----
# Updated 14.03.2026
generate_sbm <- function(
    block_matrix, blocks, membership,
    membership_matrix, density_matrix, lower_triangle
)
{

  # Set upper triangles to FALSE
  block_matrix[!lower_triangle] <- FALSE

  # Set up lower block matrix
  for(i in seq_len(blocks)){

    # Get community index
    community_index <- membership == i

    # Get block edges
    community_block <- block_matrix[community_index, community_index]

    # Compute probabilities
    community_block[community_block] <- runif_xoshiro(sum(community_block)) < density_matrix[i,i]

    # Insert block back into block matrix
    block_matrix[community_index, community_index] <- community_block

    # Check for off-blocks
    if(i != blocks){

      # Set off-block
      for(j in (i + 1):blocks){

        # Get off-block index
        off_index <- membership == j

        # Get block edges
        off_block <- block_matrix[off_index, community_index]

        # Set index
        not_index <- !off_block

        # Compute probabilities
        off_block[not_index] <- runif_xoshiro(sum(not_index)) < density_matrix[i,j]

        # Insert block back into block matrix
        block_matrix[off_index, community_index] <- off_block

      }

    }

  }

  # Return block matrix
  return(block_matrix + t(block_matrix))

}

#' @noRd
# SBM weight generation ----
# Updated 12.03.2026
sbm_weights <- function(
    block_matrix, membership, membership_matrix, sample_size,
    snr, total_nodes, negative_proportion,
    diffusion, diffusion_range, target_condition,
    max_correlation, lower_triangle
)
{

  # Check for diffusion range
  if(!is.null(diffusion_range)){
    typeof_error(diffusion_range, "numeric") # object type error
    length_error(diffusion_range, 2) # object length error
    range_error(diffusion_range, c(0, 1)) # object range error
    diffusion <- diffusion_range
  }else{
    diffusion <- c(diffusion, diffusion)
  }

  # Reverse diffusion
  diffusion <- range(1 - diffusion)

  # Set lower triangles
  block_lower <- block_matrix[lower_triangle]
  membership_lower <- membership_matrix[lower_triangle]

  # Generate edges
  total_edges <- sum(block_lower)
  edges <- generate_edges(nonzero = total_edges, n = sample_size, p = total_nodes, snr = snr)

  # Sort edges
  sorted_edges <- sort(abs(edges), decreasing = TRUE)

  # Obtain edges
  within_index <- membership_lower & block_lower
  total_within <- sum(within_index)

  # Set diffusion parameter
  within_reserved <- round(total_within * runif_xoshiro(1, min = diffusion[1], max = diffusion[2]))
  within_random <- total_within - within_reserved

  # Top edges guaranteed to go into communities
  reserved_index <- seq_len(within_reserved)

  # Remaining edges within-block edges
  remaining_index <- seq_len(within_random)
  between_index <- (within_random + 1):(total_edges - within_reserved)

  # Remaining edge weights
  remaining_weights <- shuffle(sorted_edges[-reserved_index])

  # Get signs
  signs <- swiftelse(runif_xoshiro(length(between_index)) < negative_proportion, -1, 1)

  # Set edges
  block_matrix[lower_triangle][!membership_lower & block_lower] <- shuffle(remaining_weights[between_index]) * signs
  block_matrix[lower_triangle][within_index] <- shuffle(
    c(sorted_edges[reserved_index], remaining_weights[remaining_index])
  )

  # Set upper triangle to zero
  block_matrix[!lower_triangle] <- 0

  # Make symmetric
  network <- block_matrix + t(block_matrix)

  # Add names to network
  names(membership) <- row.names(network) <- colnames(network) <- paste0(
    "V", format_integer(seq_len(total_nodes), digits(total_nodes) - 1)
  )

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set condition and lambda
  lambda <- NA

  # Check for positive definite
  if(anyNA(R) || !is_positive_definite(R)){

    # Try to condition network
    output <- try(condition_network(network, target_condition), silent = TRUE)

    # Check for whether correlation matrix is conditioned
    if(is(output, "try-error")){
      stop("Ill-conditioned")
    }

    # Store outputs
    lambda <- output$lambda

    # Update network and correlations
    P <- cor2pcor(output$R)
    P[network == 0] <- 0
    network <- P
    R <- pcor2cor(network)

    # Check for positive definite again
    if(anyNA(R) || !is_positive_definite(R)){
      stop("Not positive definite")
    }

    # Update parameters
    edges <- network[lower_triangle]
    edges <- edges[edges != 0]
    params <- weibull_mle(abs(edges))

    # Attach attributes to MLE Weibull to weights
    attr(edges, "params") <- params

    # Check bounds
    if(any(check_bounds(params[["shape"]], params[["scale"]]))){
      stop("Weibull edge parameteres were not in empirical bounds")
    }

  }

  # Check for condition
  condition <- fast_kappa(R)
  if(condition > (target_condition + 1)){
    stop("Condition number exceeds target")
  }

  # Check for maximum correlation
  if(max(abs(R[lower_triangle])) > max_correlation){
    stop("Max correlation exceeds target")
  }

  # Return results
  return(
    list(
      R = R,
      network = network,
      membership = membership,
      params = attr(edges, "params"),
      diffusion = diffusion,
      lambda = lambda,
      condition = condition
    )
  )

}

