#' Simulates Small-World GGM Data
#'
#' Simulates data from a Gaussian Graphical Model (GGM) with a small-world
#' network structure using the Watts-Strogatz model. Nodes are arranged in a
#' ring lattice and edges are randomly rewired with probability \code{rewire}.
#' Edge weights are derived empirically and the resulting network is used to
#' generate multivariate normal data. Parameters do not have default values
#' (except \code{negative_proportion}, \code{target_condition}, and \code{max_iterations}) and must
#' each be set. See examples to get started
#'
#' @param nodes Numeric (length = 1).
#' Number of nodes in the network.
#' Minimum three nodes
#'
#' @param density Numeric (length = 1).
#' Controls the initial connectivity of the ring lattice by determining the
#' number of nearest neighbors each node is connected to before rewiring.
#' Must be between 0 and 1
#'
#' @param rewire Numeric (length = 1).
#' Probability of rewiring each edge in the Watts-Strogatz model.
#' Higher values produce more random networks; lower values preserve
#' the lattice structure.
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
#' @param skew Numeric (length = 1 or \code{nodes}).
#' Skew to be included in categorical variables. It is randomly sampled from provided values.
#' Can be a single value or as many values as there are (total) variables.
#' Current skew implementation is between -2 and 2 in increments of 0.05.
#' Skews that are not in this sequence will be converted to their nearest
#' value in the sequence. Not recommended to use with \code{variables_range}
#'
#' @param skew_range Numeric (length = 2).
#' Randomly selects skews within in the range.
#' Somewhat redundant with \code{skew} but more flexible
#'
#' @param target_condition Numeric (length = 1).
#' Target condition number (using \code{\link{kappa}}) used
#' when ridge regularization is applied to an ill-conditioned precision matrix.
#' A lower value produces a better-conditioned (more stable) matrix.
#' Defaults to \code{10}.
#' For looser constraints, up to \code{30} is accepted but not recommended
#'
#' @param max_correlation Numeric (length = 1).
#' Maximum allowed absolute correlation between any pair of nodes in the
#' population correlation matrix \code{R}.
#' If the largest off-diagonal entry of \code{abs(R)} exceeds this threshold,
#' the current weight draw is rejected and a new one is attempted.
#' Must be between 0 and 1.
#' Defaults to \code{0.80}
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of attempts to find (1) a connected network structure,
#' (2) a network within empirical small-world bounds, and (3) a valid set
#' of edge weights.
#' Defaults to \code{100}
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified small-world network model}
#'
#' \item{parameters}{
#' A list containing all input, derived, and estimated parameters:
#'
#' \itemize{
#'
#' \item \code{nodes} --- Number of nodes in the network (as supplied)
#'
#' \item \code{density} --- Initial lattice density (as supplied)
#'
#' \item \code{neighbors} --- Number of nearest neighbors per node in the
#' initial ring lattice, derived from \code{nodes} and \code{density}
#'
#' \item \code{rewire} --- Edge rewiring probability (as supplied)
#'
#' \item \code{negative_proportion} --- Proportion of negative edges used
#'
#' \item \code{sample_size} --- Number of simulated cases
#'
#' \item \code{skew} --- Named numeric vector of per-variable skew values
#' actually applied during data generation
#'
#' \item \code{weibull} --- MLE parameters of the Weibull distribution fit
#' to the absolute edge weights, containing \code{shape} and \code{scale}
#'
#' \item \code{omega} --- Smallworldness omega statistic of the generated
#' network (Telesford et al., 2011); values near zero indicate small-world
#' structure
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
#' }
#'
#' }
#'
#' \item{convergence}{
#' A list containing iteration counts and conditioning information:
#'
#' \itemize{
#'
#' \item \code{iterations} --- Number of iterations needed
#'
#' \item \code{rejections} --- Reasons each pass was rejected (informative for challenging structures)
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
#' # Generate small-world data with default settings
#' basic_smallworld <- simulate_smallworld(
#'   nodes = 20, # 20 nodes in the network
#'   density = 0.30, # moderate initial lattice connectivity
#'   rewire = 0.20, # 20% rewiring probability
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Generate a more lattice-like network (low rewiring)
#' lattice_like <- simulate_smallworld(
#'   nodes = 20, # 20 nodes in the network
#'   density = 0.50, # moderate initial lattice connectivity
#'   rewire = 0.05, # 5% rewiring probability (closer to regular lattice)
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Generate a more random network (high rewiring)
#' random_like <- simulate_smallworld(
#'   nodes = 20, # 20 nodes in the network
#'   density = 0.30, # moderate initial lattice connectivity
#'   rewire = 0.80, # 80% rewiring probability (closer to random graph)
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Control the proportion of negative edges
#' with_negatives <- simulate_smallworld(
#'   nodes = 20, # 20 nodes in the network
#'   density = 0.30, # moderate initial lattice connectivity
#'   rewire = 0.20, # 20% rewiring probability
#'   sample_size = 1000, # number of cases = 1000
#'   negative_proportion = 0.20 # 20% of edges are negative
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
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367--375.
#'
#' @export
#'
# Simulate Small-world GGM data ----
# Updated 10.03.2026
simulate_smallworld <- function(
    nodes, density, rewire, negative_proportion,
    sample_size, skew = 0, skew_range = NULL,
    target_condition = 10, max_correlation = 0.80,
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
  simulate_smallworld_errors(
    nodes, density, rewire, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
  )

  # Set neighbors
  neighbors <- get_neighbors(nodes, density)

  # Create lattice network
  lattice <- smallworld_generate(nodes, neighbors, 0)

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

    # Check for smallworldness
    if((omega < -0.30) | (omega > 0.70)){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where network was small-world."

      # Move to next iteration
      next

    }

    # Try to get good weights
    output <- try(
      smallworld_weights(
        A, lattice, nodes, sample_size, neighbors,
        negative_proportion, target_condition,
        max_correlation, omega
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
  data_output <- simulate_data(n = sample_size, R = output$R, skew = skew, skew_range = skew_range)

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
# skew = 0; skew_range = NULL
# target_condition = 10;
# max_correlation = 0.80; max_iterations = 100

#' @noRd
# Errors ----
# Updated 10.03.2026
simulate_smallworld_errors <- function(
    nodes, density, rewire, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
)
{

  # Errors for 'nodes'
  typeof_error(nodes, "numeric")
  length_error(nodes, c(1, Inf))
  range_error(nodes, c(3, Inf))

  # Errors for 'density'
  typeof_error(density, "numeric")
  length_error(density, 1)
  range_error(density, c(0, 1))

  # Errors for 'rewire'
  typeof_error(rewire, "numeric")
  length_error(rewire, 1)
  range_error(rewire, c(0, 1))

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
  length_error(skew, c(1, nodes))
  range_error(skew, c(-2, 2))

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric")
  length_error(target_condition, 1)
  range_error(target_condition, c(1, 30))

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric")
  length_error(max_correlation, 1)
  range_error(max_correlation, c(0, 1))

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric")
  length_error(max_iterations, 1)
  range_error(max_iterations, c(1, Inf))

}

#' @noRd
# Get number of neighbors based on density
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
# Updated 10.03.2026
smallworld_weights <- function(
    A, lattice, nodes, sample_size, neighbors,
    negative_proportion, target_condition,
    max_correlation, omega
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

  # Obtain lower triangle
  lower_triangle <- lower.tri(A)

  # Obtain nonzero edges
  nonzero <- A[lower_triangle] != 0

  # Total edges
  total_edges <- sum(nonzero)

  # Initialize network
  network <- matrix(0, nrow = nodes, ncol = nodes)

  # Get signs
  signs <- swiftelse(runif_xoshiro(total_edges) < negative_proportion, -1, 1)

  # Generate edges
  edges <- generate_edges(nonzero = total_edges, n = sample_size, p = nodes)
  sorted_edges <- sort(edges, decreasing = TRUE)

  # Set up for smallworld Schur complement
  ## ASPL component
  distance_normed <- lattice[lower_triangle][nonzero] / (neighbors + 1)
  ## Set weights order
  weight_order <- rank((1 - distance_normed), ties.method = "random")
  network[lower_triangle][nonzero] <- sorted_edges[weight_order] * signs
  network <- network + t(network) # make symmetric

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
  condition <- kappa(R)
  if(condition > (target_condition + 5)){
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

