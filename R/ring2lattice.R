#' @title Construct a Near-Degree-Preserving Ring Lattice
#'
#' @description Converts a partial correlation network matrix into a ring
#' lattice whose per-node degree sequence approximates the original as closely
#' as possible. An adjacency matrix is derived automatically from
#' \code{network} (non-zero entries become edges), so the function accepts any
#' weighted or binary network directly. Nodes are internally reordered by
#' degree before the lattice is built. Starting from a ring lattice dense
#' enough to cover the maximum target degree, the algorithm runs two
#' independent pruning passes — one visiting nodes from highest to lowest
#' surplus and one from lowest to highest — and retains whichever pass
#' produces the smaller total residual surplus (ties broken by the higher
#' average clustering coefficient). As a final fallback, if the empirical
#' network itself has a higher clustering coefficient than either pruned
#' lattice, the empirical adjacency is returned instead. The function is
#' designed for use as the lattice baseline in small-worldness calculations.
#'
#' @param network Matrix.
#' A square, symmetric numeric matrix representing a network (e.g., partial
#' correlations). Non-zero off-diagonal entries are treated as edges; the
#' binary adjacency is derived internally as \code{network != 0}.
#' Isolated nodes (degree zero) are supported and are disconnected from the
#' lattice before pruning begins.
#'
#' @return A square symmetric binary adjacency matrix of the same dimension as
#' \code{network} representing the resulting ring lattice. The average
#' clustering coefficient of the returned lattice is attached as the attribute
#' \code{"CC"} and can be retrieved with \code{attr(result, "CC")}.
#'
#' The degree of each node in the returned lattice is less than or equal to
#' the corresponding degree in \code{network}; exact equality is achieved
#' wherever the pruning algorithm can satisfy it within \code{nodes * 2}
#' iterations per pass.
#'
#' @details
#' ## Algorithm
#'
#' \strong{Initialisation.} The binary adjacency
#' \eqn{A = [\mathbf{1}(\text{network}_{ij} \neq 0)]} is constructed
#' internally. Nodes are reordered in decreasing degree before all subsequent
#' steps. The target degree of node \eqn{i} (after reordering) is
#' \eqn{d_i = \sum_j A_{ij}}. A ring lattice is seeded by connecting node
#' \eqn{i} to node \eqn{j} whenever the ring distance between them does not
#' exceed \eqn{\lceil (d_i + d_j) / 2 \rceil}, so every node begins with at
#' least as many edges as its target degree. Nodes with \eqn{d_i = 0} are
#' immediately disconnected.
#'
#' \strong{Circular distance.} A distance matrix recording the ring distance
#' between every pair of positions is precomputed and maintained throughout:
#' \deqn{\delta_{ij} = \min\!\bigl(|i - j|,\; p - |i - j|\bigr)}
#'
#' \strong{Surplus degree.} The surplus of node \eqn{i} is
#' \eqn{e_i = \hat{d}_i - d_i \geq 0}, where \eqn{\hat{d}_i} is the current
#' ring degree. Pruning continues while \eqn{\sum_i e_i > 0} or until
#' \eqn{p \times 2} iterations have elapsed.
#'
#' \strong{Dual-pass pruning.} Two independent pruning passes are run on the
#' same seed lattice:
#' \itemize{
#'   \item \emph{Highest-first pass} — nodes are visited in decreasing order
#'     of surplus.
#'   \item \emph{Lowest-first pass} — nodes are visited in increasing order
#'     of surplus.
#' }
#' Within each pass, for an over-degree node \eqn{i}, candidate removal
#' partners are current neighbours that also carry surplus (\eqn{e_j > 0}).
#' The partner to remove is selected by the following lexicographic priority:
#' \enumerate{
#'   \item \strong{Clustering constraint} — prefer edges whose removal reduces
#'     local clustering least, scored as \eqn{1 / (c_{ij} + 1)}, where
#'     \eqn{c_{ij}} is the number of common neighbours of \eqn{i} and \eqn{j}.
#'   \item \strong{Surplus} — among equally ranked candidates, prefer the
#'     neighbour with the largest surplus \eqn{e_j}.
#'   \item \strong{Ring distance} — as a final tiebreaker, prefer the
#'     neighbour furthest away on the ring (\eqn{\delta_{ij}}), keeping the
#'     lattice as local as possible.
#' }
#'
#' \strong{Pass selection.} After both passes complete, the one with the
#' smaller total residual surplus \eqn{\sum_i e_i} is retained. If the two
#' passes are tied on total surplus, the one with the higher average clustering
#' coefficient is chosen.
#'
#' \strong{Empirical fallback.} If the average clustering coefficient of the
#' winning lattice is still lower than that of the original empirical network,
#' the empirical adjacency (reordered by degree) is returned as the lattice
#' baseline instead.
#'
#' \strong{Termination.} Each pass exits when all surpluses reach zero or
#' \eqn{p \times 2} outer iterations have been completed without further
#' change.
#'
#' @examples
#' # Get network
#' network <- network_estimation(basic_smallworld)
#'
#' # Construct ring lattice
#' L <- ring2lattice(network)
#'
#' # Retrieve the attached clustering coefficient
#' attr(L, "CC")
#'
#' # Degree sequences should be identical (or very close)
#' cbind(target = colSums(network != 0), achieved = colSums(L))
#'
#' @references
#' \strong{Ring lattice and small-world networks} \cr
#' Watts, D. J., & Strogatz, S. H. (1998).
#' Collective dynamics of \sQuote{small-world} networks.
#' \emph{Nature}, \emph{393}(6684), 440--442.
#' \doi{10.1038/30918}
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Near-degree preserving lattice ----
# Updated 21.03.2026
ring2lattice <- function(network)
{

  # Automatically construct adjacency
  A <- network != 0

  # Get nodes and density
  nodes <- dim(A)[2]
  degree <- colSums(A)

  # Set degree order
  degree_order <- order(degree, decreasing = TRUE)
  degree <- degree[degree_order]

  # Set up distance matrix
  node_sequence <- seq_len(nodes)
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  ring <- distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Convert ring matrix to zeros and ones
  ring <- ring <= (ceiling(outer(degree, degree, "+") / 2))

  # Zero out zero degree nodes
  zero_degree <- degree == 0
  if(any(zero_degree)){
    ring[zero_degree,] <- ring[,zero_degree] <- FALSE
  }

  # Compute ring degree
  ring_degree <- colSums(ring)

  # Extra per node
  extra <- ring_degree - degree

  # Perform pass from highest extra to lowest
  highest <- ring_pruning(
    degree = degree, ring_degree = ring_degree, extra = ring_degree - degree,
    ring = ring, distance_matrix = distance_matrix, max_iter = nodes * 2,
    decreasing = TRUE
  )

  # Perform pass from lowest extra to highest
  lowest <- ring_pruning(
    degree = degree, ring_degree = ring_degree, extra = ring_degree - degree,
    ring = ring, distance_matrix = distance_matrix, max_iter = nodes * 2,
    decreasing = FALSE
  )

  # Determine extras
  lowest_extra <- sum(abs(lowest$extra))
  highest_extra <- sum(abs(highest$extra))

  # Check for ties
  if(lowest_extra == highest_extra){

    # Compute clustering coefficients
    highest_CC <- igraph::transitivity(convert2igraph(highest$ring), type = "average")
    lowest_CC <- igraph::transitivity(convert2igraph(lowest$ring), type = "average")

    # Set ring
    ring <- swiftelse(highest_CC > lowest_CC, highest$ring, lowest$ring)

  }else if(lowest_extra < highest_extra){
    ring <- lowest$ring
  }else{
    ring <- highest$ring
  }

  # Compute clustering coefficients
  lattice_CC <- igraph::transitivity(convert2igraph(ring), type = "average")

  # Compute empirical clustering coefficient
  empirical_CC <- igraph::transitivity(convert2igraph(A), type = "average")

  # Fallback: If ringed solution cannot beat empirical,
  # then empirical may be most optimized lattice already
  empirical_flag <- empirical_CC > lattice_CC
  ring <- swiftelse(empirical_flag, A[degree_order, degree_order], ring)

  # Attach clustering coefficient
  attr(ring, "CC") <- swiftelse(empirical_flag, empirical_CC, lattice_CC)

  # Return ring
  return(ring)

}

#' @noRd
# Core algorithm ----
# Updated 21.03.2026
ring_pruning <- function(degree, ring_degree, extra, ring, distance_matrix, max_iter, decreasing)
{

  # Initialize order
  initial_order <- order(extra, decreasing = decreasing)

  # Initialize iterations
  iter <- 0

  # Loop while there are still some left over
  while(iter < max_iter){

    # Increase iterations
    iter <- iter + 1

    # Obtain previous extra
    previous_extra <- extra

    # Loop over extra
    for(i in initial_order){

      # Skip if already satisfied by a previous iteration
      if(extra[i] <= 0){
        next
      }

      # Current neighbors with excess edges
      neighbors <- which(ring[i,] & (extra > 0))

      # Skip if there isn't any available
      if(length(neighbors) == 0){
        next
      }

      # Clustering constraint
      common_neighbors <- nvapply(neighbors, function(j){sum(ring[i,] & ring[j,])})

      # Prioritize clustering, excess, and then distance
      index <- order(
        1 / (common_neighbors + 1), extra[neighbors],
        distance_matrix[i, neighbors], # distance constraint
        decreasing = TRUE
      )[1]
      target <- neighbors[index]

      # Loop over target to update distance and edges
      ring[target, i] <- ring[i, target] <- FALSE

      # Update state
      extra <- colSums(ring) - degree

    }

    # Check for no change
    if(all(extra == 0) || all(previous_extra == extra)){
      break
    }

  }

  # Return the output
  return(list(ring = ring, extra = extra))

}
