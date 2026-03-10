#' Simulates Stochastic Block Model Data
#'
#' Simulates data from a Gaussian Graphical Model (GGM) with a stochastic
#' block model (SBM) structure. Nodes are partitioned into communities, with
#' edge density controlled separately within and between communities. Edge
#' weights are derived empirically and the resulting network is used to
#' generate multivariate normal data. A \code{mixing} parameter controls the
#' proportion of the strongest edges that are allowed to appear between
#' communities rather than being reserved for within-community positions.
#' Parameters do not have default values (except \code{negative_proportion},
#' \code{mixing}, \code{target_condition}, and \code{max_iterations}) and
#' must each be set. See examples to get started
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
#' @param skew Numeric (length = 1 or \code{nodes * blocks}).
#' Skew to be included in categorical variables. It is randomly sampled from provided values.
#' Can be a single value or as many values as there are (total) variables.
#' Current skew implementation is between -2 and 2 in increments of 0.05.
#' Skews that are not in this sequence will be converted to their nearest
#' value in the sequence
#'
#' @param skew_range Numeric (length = 2).
#' Randomly selects skews within in the range.
#' Somewhat redundant with \code{skew} but more flexible
#'
#' @param mixing Numeric (length = 1).
#' Proportion of the strongest edges that are randomly assigned to
#' between-block positions rather than being reserved for within-block
#' positions. At \code{mixing = 0} (default), all of the strongest edges are
#' placed within communities, producing well-separated blocks. Higher values
#' redistribute more top-weighted edges to between-block positions, reducing
#' the weight separation between communities.
#' Must be between 0 and 1
#'
#' @param mixing_range Numeric (length = 2).
#' If provided, the mixing proportion is drawn uniformly from this range on
#' each call, overriding \code{mixing}. Useful for introducing variability
#' across simulation replications. For example,
#' \code{mixing_range = c(0.05, 0.20)} samples a mixing value between 5\%
#' and 20\% on each draw.
#' Both values must be between 0 and 1
#'
#' @param target_condition Numeric (length = 1).
#' Target condition number (using \code{\link{kappa}}) used
#' when ridge regularization is applied to an ill-conditioned precision matrix.
#' A lower value produces a better-conditioned (more stable) matrix.
#' Defaults to \code{30}.
#' For looser constraints, up to \code{100} is accepted but not recommended
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
#' \item{parameters}{
#' A list containing all input, derived, and estimated parameters:
#'
#' \itemize{
#'
#' \item \code{nodes} --- Number of nodes per community block (as supplied)
#'
#' \item \code{blocks} --- Number of community blocks
#'
#' \item \code{sample_size} --- Number of simulated cases
#'
#' \item \code{skew} --- Named numeric vector of per-variable skew values
#' actually applied during data generation
#'
#' \item \code{within_density} --- Within-block edge density (as supplied)
#'
#' \item \code{between_density} --- Between-block edge density (as supplied)
#'
#' \item \code{negative_proportion} --- Proportion of negative edges used
#'
#' \item \code{weibull} --- MLE parameters of the Weibull distribution fit
#' to the absolute edge weights, containing \code{shape} and \code{scale}
#'
#' \item \code{mixing} --- Numeric vector of length 2 giving the within-block
#' reservation range used during weight assignment. When \code{mixing_range}
#' is provided, this equals \code{1 - mixing_range}; otherwise both values
#' equal \code{1 - mixing}
#'
#' \item \code{Q} --- Modularity of the population network with respect to
#' the block membership structure, computed via \code{igraph::modularity}
#'
#' }
#'
#' }
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
#' \item \code{lambda} --- Ridge regularization parameter used to condition
#' the precision matrix when it was not initially positive definite;
#' \code{NA} if no conditioning was required
#'
#' \item \code{condition} --- Condition number (ratio of largest to smallest
#' eigenvalue) of the population correlation matrix \code{R}
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
#' # Control community separation via mixing
#' mixed_blocks <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block
#'   blocks = 3, # 3 community blocks
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = 0.90, # 90% edge probability within blocks
#'   between_density = 0.20, # 20% edge probability between blocks
#'   mixing = 0.10 # 10% of strongest edges can appear between blocks
#' )
#'
#' # Variable mixing across simulation replications
#' variable_mixing <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block
#'   blocks = 3, # 3 community blocks
#'   sample_size = 1000, # number of cases = 1000
#'   within_density = 0.90, # 90% edge probability within blocks
#'   between_density = 0.20, # 20% edge probability between blocks
#'   mixing_range = c(0.05, 0.20) # mixing sampled uniformly from 5-20%
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
# Updated 10.03.2026
simulate_sbm <- function(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size,
    skew = 0, skew_range = NULL,
    mixing = 0, mixing_range = NULL,
    target_condition = 30,
    max_iterations = 100
)
{

  # Check for missing negative proportion
  if(missing(negative_proportion)){

    # Compute proportion of negative edges based on empirical values
    negative_proportion <- pmin(
      pmax(
        0.3511293 + rnorm_ziggurat(1) * 0.08295474, # empirical mean +/- 1 SD
        0.08333333 # empirical minimum
      ),
      0.54945055 # empirical maximum
    )

  }

  # Check for input errors
  simulate_sbm_errors(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size, skew, mixing,
    target_condition, max_iterations
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

  # Collect errors
  errors <- character(length = max_iterations + 1)

  # Generate edges for the network
  while(bad_weights){

    # Try to get good weights
    output <- try(
      sbm_weights(
        block_matrix, membership, membership_matrix, sample_size,
        total_nodes, negative_proportion,
        mixing, mixing_range,
        target_condition, lower_triangle
      ), silent = TRUE
    )

    # Set good weights
    bad_weights <- is(output, "try-error")

    # Increase count
    weights_iter <- weights_iter + 1

    # If an error, collect it
    if(bad_weights){
      errors[weights_iter] <- trimws(strsplit(output[1], split = "\n")[[1]][2])
    }else if(output$condition > (target_condition + 10)){

      # Update bad weights and create error
      bad_weights <- TRUE
      errors[weights_iter] <- "Condition greater than target"

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
  data_output <- simulate_data(n = sample_size, R = output$R, skew = skew, skew_range = skew_range)

  # Return parameters
  return(
    list(
      data = data_output$data,
      parameters = list(
        nodes = nodes,
        blocks = blocks,
        sample_size = sample_size,
        skew = data_output$skew,
        within_density = within_density,
        between_density = between_density,
        negative_proportion = negative_proportion,
        weibull = output$params,
        mixing = output$mixing,
        Q = igraph::modularity(convert2igraph(abs(output$network)), output$membership)
      ),
      population = list(
        R = output$R,
        Omega = output$network,
        membership = output$membership
      ),
      convergence = list(
        connected = connected_iter,
        weights = weights_iter,
        lambda = output$lambda,
        condition = output$condition
      )
    )
  )

}

# Bug checking ----
# nodes = c(4, 6, 8); blocks = 3; sample_size = 1000
# within_density = 0.90; between_density = 0.20
# negative_proportion = 0.35; skew = 0; skew_range = NULL
# target_condition = 30; max_iterations = 100

#' @noRd
# Errors ----
# Updated 10.03.2026
simulate_sbm_errors <- function(
    nodes, blocks, within_density, between_density,
    negative_proportion, sample_size, skew, mixing,
    target_condition, max_iterations
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

  # Errors for 'skew'
  typeof_error(skew, "numeric")
  length_error(skew, c(1, swiftelse(length(nodes) == 1, nodes * blocks, sum(nodes))))
  range_error(skew, c(-2, 2))

  # Errors for 'mixing'
  typeof_error(mixing, "numeric")
  length_error(mixing, 1)
  range_error(mixing, c(0, 1))

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric")
  length_error(target_condition, 1)
  range_error(target_condition, c(1, 100))

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric")
  length_error(max_iterations, 1)
  range_error(max_iterations, c(1, Inf))

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
# Updated 10.03.2026
sbm_weights <- function(
    block_matrix, membership, membership_matrix, sample_size,
    total_nodes, negative_proportion,
    mixing, mixing_range,
    target_condition, lower_triangle
)
{

  # Check for mixing range
  if(!is.null(mixing_range)){
    typeof_error(mixing_range, "numeric") # object type error
    length_error(mixing_range, 2) # object length error
    range_error(mixing_range, c(0, 1)) # object range error
    mixing <- range(mixing_range)
  }else{
    mixing <- c(mixing, mixing)
  }

  # Reverse mixing
  mixing <- 1 - mixing

  # Generate edges
  total_edges <- sum(block_matrix)
  edges <- generate_edges(nonzero = total_edges, n = sample_size, p = total_nodes)

  # Sort edges
  sorted_edges <- sort(abs(edges), decreasing = TRUE)

  # Obtain edges
  within_index <- membership_matrix & block_matrix
  total_within <- sum(within_index)

  # Set mixing parameter
  within_reserved <- round(total_within * runif_xoshiro(1, min = mixing[2], max = mixing[1]))
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

  # Set condition and lambda
  lambda <- NA

  # Check for positive definite
  if(anyNA(R) || !is_positive_definite(R)){

    # Try to condition network
    output <- try(condition_network(network, target_condition), silent = TRUE)

    # Check for positive definite again
    if(is(output, "try-error")){
      stop("Ill-conditioned network")
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

  # Return results
  return(
    list(
      R = R,
      network = network,
      membership = membership,
      params = attr(edges, "params"),
      mixing = mixing,
      lambda = lambda,
      condition = kappa(R)
    )
  )

}

