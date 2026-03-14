#' Simulates Small-World GGM Data
#'
#' @description
#' Simulates data from a Gaussian Graphical Model (GGM) with a small-world
#' network structure using the Watts-Strogatz model. Nodes are arranged in a
#' ring lattice where each node is connected to its \code{k} nearest neighbors,
#' and edges are then randomly rewired with probability \code{rewire}. The
#' number of nearest neighbors \code{k} is derived from \code{nodes} and
#' \code{density}. Edge weights are assigned by structural priority: edges
#' closer to their original lattice positions receive larger partial correlation
#' weights, grounded in the empirical observation that shorter-distance
#' connections tend to carry stronger weights in psychometric networks.
#' The resulting network is used to generate multivariate normal data.
#' Parameters do not have default values (except \code{negative_proportion},
#' \code{target_condition}, and \code{max_iterations}) and must each be set.
#' See Details and Examples to get started.
#'
#' @param nodes Numeric (length = 1).
#' Number of nodes in the network.
#' Minimum of three nodes.
#'
#' @param density Numeric (length = 1).
#' Controls the initial connectivity of the ring lattice by determining the
#' number of nearest neighbors \code{k} each node is connected to before
#' rewiring. Specifically, \code{k} is derived as
#' \eqn{k = \mathrm{round}((\text{nodes} \times \frac{\text{nodes}-1}{2} \times
#' \text{density}) / \text{nodes})}, subject to a minimum of 2 and a maximum
#' of \eqn{\lfloor (\text{nodes}-1)/2 \rfloor}.
#' Must be between 0 and 1.
#'
#' @param rewire Numeric (length = 1).
#' Probability of rewiring each edge in the Watts-Strogatz model.
#' Values near 0 preserve the regular lattice structure; values near 1
#' produce approximately random networks. The small-world regime typically
#' occurs at intermediate values (roughly 0.01 to 0.30).
#' Must be between 0 and 1.
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
#' Proportion of edges that are negative (inhibitory).
#' Must be between 0 and 0.50. The upper bound of 0.50 is a mathematical
#' constraint of the sign-flipping procedure used to assign negative edges
#' (see Details).
#' If not provided, a value is sampled from an empirical distribution with
#' mean = 0.34 and SD = 0.086, bounded between 0.083 and 0.50.
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases to generate from the population multivariate normal
#' distribution.
#'
#' @param skew Numeric (length = 1 or \code{nodes}).
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
#' Maximum number of attempts to find (1) a connected network structure,
#' (2) a network within empirical small-world bounds, and (3) a valid set
#' of edge weights. Defaults to \code{100}.
#'
#' @details
#' \strong{Edge weight assignment}
#'
#' Edge weights are assigned by ranking edges according to their distance
#' in the original ring lattice before rewiring. Edges that were retained
#' from the lattice receive their original neighbor distance (1 = nearest
#' neighbor, 2 = second nearest, etc.); rewired edges are assigned a
#' distance of \code{k + 1}, placing them at the bottom of the priority
#' order. Pre-generated Weibull weights (sorted descending) are then
#' mapped onto this ranking so that the largest weights go to the
#' shortest-distance edges. This assignment is a structural heuristic
#' grounded in the Watts-Strogatz data generating process and is
#' consistent with empirical observations that shorter-distance local
#' connections tend to carry larger partial correlation weights in
#' psychometric networks.
#'
#' \strong{Negative edge assignment}
#'
#' Negative edges are introduced by flipping the signs of a subset of
#' nodes — multiplying all edges incident to those nodes by -1. The number
#' of nodes to flip is derived from \code{negative_proportion} via the
#' inversion formula \eqn{(1 - \sqrt{1 - 2p}) / 2}, where \eqn{p} is the
#' target proportion of negative edges. This formula assumes that an edge
#' is negative if and only if exactly one of its endpoints is flipped,
#' giving an expected negative proportion of \eqn{2 \cdot (k/n) \cdot
#' (1 - k/n)} for \eqn{k} flipped nodes out of \eqn{n} total. The formula
#' requires \eqn{p \leq 0.50}, which is why \code{negative_proportion} is
#' bounded at 0.50.
#'
#' \strong{Small-worldness screening}
#'
#' Each generated adjacency structure is screened using the omega statistic
#' (Telesford et al., 2011): \eqn{\omega = (L_r / L) - (C / C_l)}, where
#' \eqn{L} is the average shortest path length of the network, \eqn{L_r}
#' is the mean ASPL of degree-matched random networks, \eqn{C} is the
#' average clustering coefficient, and \eqn{C_l} is the clustering
#' coefficient of a lattice with the same \code{k}. Networks with
#' \eqn{|\omega| > 0.80} are rejected. Note that the lattice reference
#' uses the Watts-Strogatz construction (p = 0) rather than a maximally
#' regular lattice, which is internally consistent with the data generating
#' process.
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data matrix (\code{sample_size x nodes}) drawn
#' from the population GGM.}
#'
#' \item{parameters}{
#' A list of input, derived, and estimated parameters:
#' \itemize{
#'   \item \code{nodes} --- Number of nodes (as supplied)
#'   \item \code{density} --- Initial lattice density (as supplied)
#'   \item \code{neighbors} --- Number of nearest neighbors \code{k} in the
#'   initial ring lattice, derived from \code{nodes} and \code{density}
#'   \item \code{rewire} --- Edge rewiring probability (as supplied)
#'   \item \code{negative_proportion} --- Proportion of negative edges used
#'   \item \code{sample_size} --- Number of simulated cases
#'   \item \code{skew} --- Named numeric vector of per-variable skew values
#'   actually applied
#'   \item \code{weibull} --- Weibull \code{shape} and \code{scale} parameters
#'   of the absolute edge weight distribution
#'   \item \code{omega} --- Smallworldness omega statistic of the generated
#'   network; values near zero indicate small-world structure, negative values
#'   indicate lattice-like structure, and positive values indicate random-like
#'   structure
#' }
#' }
#'
#' \item{population}{
#' Population-level network parameters:
#' \itemize{
#'   \item \code{R} --- Population correlation matrix derived from the GGM
#'   \item \code{Omega} --- Population partial correlation matrix (GGM edge
#'   weights)
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
#' # Basic small-world network (moderate density, moderate rewiring)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.20,
#'   sample_size = 500
#' )
#'
#' # Lattice-like structure (low rewiring preserves local connectivity)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.01,
#'   sample_size = 500
#' )
#'
#' # Random-like structure (high rewiring destroys lattice regularity)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.80,
#'   sample_size = 500
#' )
#'
#' # Fix the proportion of negative edges
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.20,
#'   sample_size = 500,
#'   negative_proportion = 0.20
#' )
#'
#' # Larger network with higher density
#' result <- simulate_smallworld(
#'   nodes = 40,
#'   density = 0.50,
#'   rewire = 0.10,
#'   sample_size = 1000
#' )
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{Seminal introduction to the Watts-Strogatz small-world model} \cr
#' Watts, D. J., & Strogatz, S. H. (1998).
#' Collective dynamics of 'small-world' networks.
#' \emph{Nature}, \emph{393}(6684), 440--442.
#'
#' \strong{Omega statistic for smallworldness} \cr
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., &
#' Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367--375.
#'
#' @export
#'
# Simulate Small-world GGM data ----
# Updated 12.03.2026
simulate_smallworld <- function(
    nodes, density, rewire, snr = 1, negative_proportion,
    sample_size, skew = 0, skew_range = NULL,
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
      ), 0.50 # empirical maximum (0.54945055)
      # needs to be capped at 0.50 due to flipping variables
    )

  }

  # Check for input errors
  simulate_smallworld_errors(
    nodes, density, rewire, snr, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
  )

  # Set neighbors
  neighbors <- get_neighbors(nodes, density)

  # Create lattice network
  lattice <- smallworld_generate(nodes, neighbors, 0)

  # Obtain lower triangle
  lower_triangle <- lower.tri(lattice)

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
          "Reached maximum iterations: Could not find weights that for the small-world structure ",
          "Resulting rejections were due to:\n\n",
          paste0(names(rejection_table), " = ", rejection_table, collapse = "\n")
        )
      )

    }

    # Small-world adjacency
    A <- smallworld_generate(nodes, neighbors, rewire)

    # Check adjacency
    if(!igraph::is_connected(convert2igraph(A))){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where every node was connected."

      # Move to next iteration
      next

    }

    # Estimate small-worldness (|0.50| based on Telesford et al., 2011)
    omega <- smallworldness(A, nodes, neighbors)[["omega"]]

    # Check for smallworldness (based on empirical range of -0.580 to 0.70)
    if(abs(omega) > 0.80){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where network was small-world."

      # Move to next iteration
      next

    }

    # Try to get good weights
    output <- try(
      smallworld_weights(
        A, lattice, nodes, neighbors, sample_size,
        snr, negative_proportion, target_condition,
        max_correlation, omega, lower_triangle
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
    n = sample_size, p = nodes, R = output$R, skew = skew, skew_range = skew_range
  )

  # Return parameters
  return(
    list(
      data = data_output$data,
      parameters = list(
        nodes = nodes,
        density = density,
        neighbors = neighbors,
        rewire = rewire,
        negative_proportion = negative_proportion,
        sample_size = sample_size,
        skew = data_output$skew,
        weibull = output$params,
        omega = output$omega
      ),
      population = list(R = output$R, Omega = output$network),
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
# nodes = 20; density = 0.50; rewire = 0.10
# negative_proportion = 0.35; sample_size = 1000
# skew = 0; skew_range = NULL; snr = 1
# target_condition = 30;
# max_correlation = 0.80; max_iterations = 100

#' @noRd
# Errors ----
# Updated 14.03.2026
simulate_smallworld_errors <- function(
    nodes, density, rewire, snr, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
)
{

  # Errors for 'nodes'
  typeof_error(nodes, "numeric", "simulate_smallworld")
  length_error(nodes, 1, "simulate_smallworld")
  range_error(nodes, c(3, Inf), "simulate_smallworld")

  # Errors for 'density'
  typeof_error(density, "numeric", "simulate_smallworld")
  length_error(density, 1, "simulate_smallworld")
  range_error(density, c(0, 1), "simulate_smallworld")

  # Errors for 'rewire'
  typeof_error(rewire, "numeric", "simulate_smallworld")
  length_error(rewire, 1, "simulate_smallworld")
  range_error(rewire, c(0, 1), "simulate_smallworld")

  # Error for 'snr'
  typeof_error(snr, "numeric", "simulate_smallworld")
  length_error(snr, 1, "simulate_smallworld")
  range_error(snr, c(0.01, Inf), "simulate_smallworld")

  # Send warning if beyond bounds
  if((snr < 0.64) | (snr > 1.72)){
    warning(paste0(
      "The signal-to-ratio specified (", round(snr, 2), ") is beyond ",
      "the bounds of cataloged empirical data (0.64-1.72). "
    ))
  }

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric", "simulate_smallworld")
  length_error(negative_proportion, 1, "simulate_smallworld")
  range_error(negative_proportion, c(0, 0.50), "simulate_smallworld")

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric", "simulate_smallworld")
  length_error(sample_size, 1, "simulate_smallworld")
  range_error(sample_size, c(1, Inf), "simulate_smallworld")

  # Errors for 'skew'
  typeof_error(skew, "numeric", "simulate_smallworld")
  length_error(skew, c(1, nodes), "simulate_smallworld")
  range_error(skew, c(-2, 2), "simulate_smallworld")

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric", "simulate_smallworld")
  length_error(target_condition, 1, "simulate_smallworld")
  range_error(target_condition, c(1, 100), "simulate_smallworld")

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric", "simulate_smallworld")
  length_error(max_correlation, 1, "simulate_smallworld")
  range_error(max_correlation, c(0, 1), "simulate_smallworld")

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric", "simulate_smallworld")
  length_error(max_iterations, 1, "simulate_smallworld")
  range_error(max_iterations, c(1, Inf), "simulate_smallworld")

}

#' @noRd
# Get number of neighbors based on density ----
# Updated 08.03.2026
get_neighbors <- function(nodes, density)
{

  # Degrees of freedom
  half_df <- (nodes - 1) / 2

  # Total edges with density
  edges <- nodes * half_df * density

  # Return with minimum neighbor check
  return(
    max(
      min(
        max(round(edges / nodes), 1), # maximum neighbors
        floor(half_df) # minimum degrees of freedom
      ), 2 # at least two neighbors
    )
  )

}

#' @noRd
# Generate smallworld ----
# Updated 08.03.2026
smallworld_generate <- function(nodes, neighbors, rewire)
{

  return(
    silent_call(
      as.matrix(
        igraph::as_adjacency_matrix(
          igraph::sample_smallworld(
            dim = 1, size = nodes,
            nei = neighbors, p = rewire
          )
        )
      )
    )
  )

}

#' @noRd
# Compute smallworldness ----
# Updated 08.03.2026
smallworldness <- function (A, nodes, neighbors, iter = 100)
{

  # Obtain edges
  degree <- colSums(A, na.rm = TRUE)

  # Convert network to {igraph}
  I <- convert2igraph(A)

  # Obtain empirical values
  ASPL <- igraph::mean_distance(I)
  CC <- igraph::transitivity(I, type = "average")

  # Create lattice network
  lattice <- igraph::sample_smallworld(
    dim = 1, size = nodes,
    nei = neighbors, p = 0
  )

  # Obtain lattice values
  lattice_CC <- igraph::transitivity(lattice, type = "average")
  # Not technically correct based on Telesford et al. (2011);
  # however, matches the data generating process and therefore
  # is the correct choice for determine smallworldness

  # Collect random networks
  random_ASPL <- mean(
    nvapply(seq_len(iter), function(i){
      igraph::mean_distance(igraph::sample_degseq(out.deg = degree, method = "vl"))
    }), na.rm = TRUE
  )

  # Compute omega
  omega <- (random_ASPL / ASPL) - (CC / lattice_CC)

  # Return values
  return(
    c(
      "omega" = omega, "aspl" = ASPL, "cc" = CC,
      "random_aspl" = random_ASPL, "lattice_cc" = lattice_CC
    )
  )

}

#' @noRd
# Smallworld weight generation ----
# Updated 14.03.2026
smallworld_weights <- function(
    A, lattice, nodes, neighbors, sample_size,
    snr, negative_proportion, target_condition,
    max_correlation, omega, lower_triangle
)
{

  # Create sparse lattice network
  sparse_lattice <- sparse_network(lattice)
  sparse_lattice <- sparse_lattice[sparse_lattice$weight == 1,]
  distance <- abs(sparse_lattice$row - sparse_lattice$col)
  sparse_lattice$weight <- pmin(distance, nodes - distance)

  # Get not retained edges
  retained <- (lattice == 1) & (A == 1)

  # Get sparse retained edges
  sparse_retained <- sparse_network(retained)
  sparse_retained <- sparse_retained[sparse_retained$weight,]

  # Create character sequence
  sparse_lattice_edges <- paste(sparse_lattice$row, sparse_lattice$col)
  sparse_retained_edges <- paste(sparse_retained$row, sparse_retained$col)

  # Update sparse lattice
  sparse_lattice <- sparse_lattice[sparse_lattice_edges %in% sparse_retained_edges,]

  # Reset lattice network
  lattice[] <- neighbors + 1

  # Rebuild lattice
  for(i in 1:nrow(sparse_lattice)){
    lattice[sparse_lattice$row[i],sparse_lattice$col[i]] <-
      lattice[sparse_lattice$col[i],sparse_lattice$row[i]] <-
      sparse_lattice$weight[i]
  }

  # Obtain nonzero edges
  nonzero <- A[lower_triangle] != 0

  # Total edges
  total_edges <- sum(nonzero)

  # Initialize network
  network <- matrix(0, nrow = nodes, ncol = nodes)

  # Generate edges (returned edges are sorted)
  edges <- generate_edges(nonzero = total_edges, n = sample_size, p = nodes, snr = snr)

  # Set weights order
  weight_order <- rank(lattice[lower_triangle][nonzero], ties.method = "random")
  network[lower_triangle][nonzero] <- edges[weight_order]
  network <- network + t(network) # make symmetric

  # Get signs
  network <- determine_signs(network, nodes, negative_proportion)

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set lambda
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
      omega = omega,
      params = attr(edges, "params"),
      lambda = lambda,
      condition = condition
    )
  )

}

#' @noRd
# Determine balance of signs ----
# Updated 11.03.2026
determine_signs <- function(network, nodes, negative_proportion)
{

  # Check for zero negative weights
  if(negative_proportion == 0){
    return(network)
  }

  # Convert negative_proportion (edge-level) to node flip proportion
  flip_proportion <- (1 - sqrt(1 - 2 * negative_proportion)) / 2
  flip_number <- round(nodes * flip_proportion)

  # Flip at least one
  flip_number <- max(1, min(flip_number, floor(nodes / 2)))

  # Determine node indices
  index <- shuffle(seq_len(nodes), flip_number)

  # Flip nodes
  network[index,] <- -network[index,]
  network[,index] <- -network[,index]

  # Return the network
  return(network)

}