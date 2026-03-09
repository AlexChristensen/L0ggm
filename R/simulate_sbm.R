#' Simulates Stochastic Block Model Data
#'
#' Simulates data from a Gaussian Graphical Model (GGM) with a stochastic
#' block model (SBM) structure. Nodes are partitioned into communities, with
#' edge density controlled separately within and between communities. Edge
#' weights are derived empirically and the resulting network is used to
#' generate multivariate normal data. Parameters do not have default values
#' (except \code{negative_proportion} and \code{max_iterations}) and must
#' each be set. See examples to get started
#'
#' @param nodes Numeric (length = 1 or \code{blocks}).
#' Number of nodes per community block.
#' Can be a single value (applied to all blocks) or as many values as there
#' are blocks (corresponding to each block).
#' Minimum three nodes per block
#'
#' @param blocks Numeric (length = 1).
#' Number of community blocks
#'
#' @param within_density Numeric (length = 1 or \code{blocks}).
#' Probability that an edge exists between any two nodes within the same
#' community block.
#' Can be a single value (applied to all blocks) or as many values as there
#' are blocks.
#' Must be between 0 and 1
#'
#' @param between_density Numeric (length = 1 or \code{blocks}).
#' Probability that an edge exists between any two nodes in different
#' community blocks.
#' Can be a single value (applied to all block pairs) or as many values as
#' there are blocks.
#' Must be between 0 and 1
#'
#' @param negative_proportion Numeric (length = 1).
#' Proportion of edges that are negative (inhibitory).
#' Must be between 0 and 1.
#' If not provided, a value is sampled from an empirical distribution with
#' mean = 0.35 and SD = 0.09, bounded between 0.08 and 0.55
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases to generate from a random multivariate normal distribution
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of attempts to find (1) a connected network structure and
#' (2) a valid set of edge weights.
#' Defaults to \code{100}
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified SBM network model}
#'
#' \item{population}{
#' A list containing the population-level network parameters:
#'
#' \itemize{
#'
#' \item \code{R} --- Population correlation matrix derived from the GGM
#'
#' \item \code{Omega} --- Population partial correlation network (GGM)
#'
#' \item \code{membership} --- Named integer vector of community block
#' assignments for each node
#'
#' }
#'
#' }
#'
#' \item{weight_parameters}{
#' MLE parameters of the Weibull distribution fit to the absolute edge weights,
#' containing \code{shape} and \code{scale}
#' }
#'
#' \item{convergence}{
#' A list containing iteration counts and conditioning information:
#'
#' \itemize{
#'
#' \item \code{connected} --- Number of iterations needed to obtain a
#' connected network
#'
#' \item \code{weights} --- Number of iterations needed to obtain valid
#' edge weights
#'
#' \item \code{condition} --- Ridge regularization parameter lambda used
#' to condition the network if it was not initially positive definite;
#' \code{NA} if no conditioning was required
#'
#' }
#'
#' }
#'
#' @examples
#' # Generate SBM data with equal-sized blocks
#' three_block <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block
#'   blocks = 3, # 3 community blocks = 18 total nodes
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = 0.90, # 90% edge probability within blocks
#'   between_density = 0.20 # 20% edge probability between blocks
#' )
#'
#' # Generate SBM data with unequal block sizes
#' unequal_blocks <- simulate_sbm(
#'   nodes = c(4, 6, 8), # 4, 6, and 8 nodes per block = 18 total nodes
#'   blocks = 3, # 3 community blocks
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = 0.90, # 90% edge probability within blocks
#'   between_density = 0.20 # 20% edge probability between blocks
#' )
#'
#' # Generate SBM data with varying density per block
#' varying_density <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block
#'   blocks = 3, # 3 community blocks
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = c(0.85, 0.90, 0.95), # density varies by block
#'   between_density = 0.15 # 15% edge probability between blocks
#' )
#'
#' # Control the proportion of negative edges
#' with_negatives <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block
#'   blocks = 3, # 3 community blocks
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = 0.90, # 90% edge probability within blocks
#'   between_density = 0.20, # 20% edge probability between blocks
#'   negative_proportion = 0.20 # 20% of edges are negative
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
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
# Updated 08.03.2026
simulate_sbm <- function(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size, max_iterations = 100
)
{

  # Check for missing negative proportion
  if(missing(negative_proportion)){

    # Compute proportion of negative edges based on empirical values
    negative_proportion <- pmin(
      pmax(
        0.3519157 + rnorm_ziggurat(1) * 0.08531577, # empirical mean +/- 1 SD
        0.08333333 # empirical minimum
      ),
      0.54945055 # empirical maximum
    )

  }

  # Check for input errors
  simulate_sbm_errors(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size
  )

  # Determine total number of nodes
  total_nodes <- swiftelse(length(nodes) == 1, nodes * blocks, sum(nodes))

  # Get community sequence
  community_sequence <- seq_len(blocks)

  # Set membership
  membership <- sort(rep(community_sequence, times = nodes))

  # Set up membership matrix
  membership_matrix <- block_matrix <- outer(membership, membership, "==")

  # Obtain lower triangle
  lower_triangle <- lower.tri(membership_matrix)

  # Initialize generation checks
  connected <- FALSE
  bad_weights <- TRUE

  # Initialize iterations
  weights_iter <- connected_iter <- 0

  # Ensure all nodes are connected to one other node
  while(!connected){

    # Get block matrix
    block_matrix <- generate_sbm(
      block_matrix, blocks, membership, membership_matrix,
      within_density, between_density, lower_triangle
    )

    # Check block matrix
    connected <- igraph::is_connected(convert2igraph(block_matrix + t(block_matrix)))

    # Increase count
    connected_iter <- connected_iter + 1

    # Stop on greater than max iterations
    if(connected_iter > max_iterations){
      stop("Reached maximum iterations: Could not find solution where every node was connected.")
    }

  }

  # Get Weibull model
  weibull_weights <- get(data("weibull_weights", package = "L0ggm", envir = environment()))

  # Collect errors
  errors <- character(length = max_iterations + 1)

  # Generate edges for the network
  while(bad_weights){

    # Try to get good weights
    output <- try(
      sbm_weights(
        weibull_weights, block_matrix, membership, membership_matrix, sample_size,
        total_nodes, negative_proportion, lower_triangle
      ), silent = TRUE
    )

    # Set good weights
    bad_weights <- is(output, "try-error")

    # Increase count
    weights_iter <- weights_iter + 1

    # If an error, collect it
    if(bad_weights){
      errors[weights_iter] <- trimws(strsplit(output[1], split = "\n")[[1]][2])
    }

    # Stop on greater than max iterations
    if(weights_iter > max_iterations){

      # Collect errors
      error_table <- fast_table(errors)

      # Remove NULL
      error_table <- error_table[names(error_table) != ""]

      # Return error
      stop(
        paste0(
          "Reached maximum iterations: Could not find weights that for the SBM solution. ",
          "Resulting errors were due to:\n\n",
          paste0(names(error_table), " = ", error_table, collapse = "\n")
        )
      )

    }

  }

  # Generate data
  data <- MASS_mvrnorm(n = sample_size, mu = rep(0, total_nodes), Sigma = diag(total_nodes))

  # Return parameters
  return(
    list(
      data = data %*% chol(output$R),
      population = list(
        R = output$R, Omega = output$network, membership = output$membership
      ),
      weight_parameters = output$params,
      convergence = list(
        connected = connected_iter,
        weights = weights_iter,
        condition = output$lambda
      )
    )
  )

}

# Bug checking ----
# nodes = c(4, 6, 8); blocks = 3; sample_size = 1000
# within_density = 0.90; between_density = 0.20
# negative_proportion = 0.35; max_iterations = 100

#' @noRd
# Errors ----
# Updated 08.03.2026
simulate_sbm_errors <- function(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size
)
{

  # Errors for 'blocks'
  typeof_error(blocks, "numeric")
  length_error(blocks, 1)
  range_error(blocks, c(1, Inf))

  # Errors for 'nodes'
  typeof_error(nodes, "numeric")
  length_error(nodes, c(1, blocks))
  range_error(nodes, c(3, Inf))

  # Errors for 'within_density'
  typeof_error(within_density, "numeric")
  length_error(within_density, c(1, blocks))
  range_error(within_density, c(0, 1))

  # Errors for 'between_density'
  typeof_error(between_density, "numeric")
  length_error(between_density, c(1, blocks))
  range_error(between_density, c(0, 1))

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric")
  length_error(negative_proportion, 1)
  range_error(negative_proportion, c(0, 1))

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric")
  length_error(sample_size, 1)
  range_error(sample_size, c(1, Inf))

}

#' @noRd
# Generate SBM ----
# Updated 08.03.2026
generate_sbm <- function(
    block_matrix, blocks, membership, membership_matrix,
    within_density, between_density, lower_triangle
)
{

  # Set upper triangles to FALSE
  block_matrix[!lower_triangle] <- FALSE

  # Set up lower block matrix
  for(i in 1:(blocks - 1)){

    # Get community index
    community_index <- membership == i

    # Get block edges
    community_block <- block_matrix[community_index, community_index]

    # Compute probabilities
    community_block[community_block] <- runif_xoshiro(sum(community_block)) < within_density

    # Insert block back into block matrix
    block_matrix[community_index, community_index] <- community_block

    # Set off-block
    for(j in (i + 1):blocks){

      # Get off-block index
      off_index <- membership == j

      # Get block edges
      off_block <- block_matrix[off_index, community_index]

      # Set index
      not_index <- !off_block

      # Compute probabilities
      off_block[not_index] <- runif_xoshiro(sum(not_index)) < (between_density)

      # Insert block back into block matrix
      block_matrix[off_index, community_index] <- off_block

    }

  }

  # Return block matrix
  return(block_matrix + t(block_matrix))

}

#' @noRd
# SBM weight generation ----
# Updated 08.03.2026
sbm_weights <- function(
    weibull_weights, block_matrix, membership, membership_matrix,
    sample_size, total_nodes, negative_proportion, lower_triangle
)
{

  # Generate edges
  total_edges <- sum(block_matrix)
  edges <- generate_edges(weibull_weights, nonzero = total_edges, n = sample_size, p = total_nodes)

  # Sort edges
  sorted_edges <- sort(abs(edges), decreasing = TRUE)

  # Obtain edges
  within_index <- membership_matrix & block_matrix
  total_within <- sum(within_index)

  # Set mixing parameter = (U(0.15, 0.35))
  within_reserved <- round(total_within * runif_xoshiro(1, min = 0.65, max = 0.85))
  within_random <- total_within - within_reserved

  # Top edges guaranteed to go into communities
  reserved_index <- seq_len(within_reserved)

  # Remaining edges within-block edges
  remaining_index <- seq_len(within_random)
  between_index <- (within_random + 1):(total_edges - within_reserved)

  # Remaining edge weights
  remaining_weights <- shuffle(sorted_edges[-reserved_index])

  # Divide by number of between edges
  n_between <- length(between_index)
  negative_proportion <- min(negative_proportion / (n_between / total_edges), 1)

  # Get signs
  signs <- swiftelse(runif_xoshiro(n_between) < negative_proportion, -1, 1)

  # Set edges
  block_matrix[(!membership_matrix) & block_matrix] <- shuffle(remaining_weights[between_index]) * signs
  block_matrix[within_index] <- shuffle(
    c(sorted_edges[reserved_index], remaining_weights[remaining_index])
  )

  # Make symmetric
  network <- block_matrix + t(block_matrix)

  # Add names to network
  names(membership) <- row.names(network) <- colnames(network) <- paste0(
    "V", format_integer(seq_len(total_nodes), digits(total_nodes) - 1)
  )

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set lambda
  lambda <- NA

  # Check for positive definite
  if(anyNA(R) || !is_positive_definite(R)){

    # Try to condition network
    condition_output <- try(condition_network(network), silent = TRUE)

    # Check for positive definite again
    if(is(condition_output, "try-error")){
      stop("Ill-conditioned network")
    }

    # Store outputs
    R <- condition_output$R; lambda <- condition_output$lambda

    # Check for positive definite again
    if(anyNA(R) || !is_positive_definite(R)){
      stop("Not positive definite")
    }

    # Update network
    P <- cor2pcor(R)
    P[network == 0] <- 0
    network <- P

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

  # Return results
  return(
    list(
      R = R,
      network = network,
      membership = membership,
      params = attr(edges, "params"),
      lambda = lambda
    )
  )

}

