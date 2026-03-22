#' @title Construct a Near-Degree-Preserving Ring Lattice
#'
#' @description Converts a network matrix into a ring lattice whose per-node
#' degree sequence approximates the original as closely as possible, while
#' maximising the average clustering coefficient. An adjacency matrix is
#' derived automatically from \code{network} (non-zero entries become edges),
#' so the function accepts any weighted or binary network directly.
#'
#' Starting from an overcomplete ring lattice seeded from the input degree
#' sequence, the algorithm runs \code{shuffles} independent pruning passes to
#' reduce each node's degree to its target. Edge-removal decisions are guided
#' by a lexicographic priority (clustering preservation, surplus, structural
#' distance). The pass producing the smallest residual surplus — further
#' broken by highest average clustering coefficient — is returned in the
#' original node ordering of \code{network}. An empirical fallback ensures
#' the result is never worse than the input network's own clustering coefficient.
#'
#' @param network Matrix.
#' A square, symmetric numeric matrix representing a network (e.g., partial
#' correlations). Non-zero off-diagonal entries are treated as edges; the
#' binary adjacency is derived internally as \code{network != 0}.
#' Isolated nodes (degree zero) are supported and are disconnected from the
#' lattice before pruning begins.
#'
#' @param shuffles Numeric (length = 1).
#' Total number of independent pruning passes to run on the seed lattice.
#' Must be at least 2. The first two passes use fixed node orderings
#' (highest-to-lowest surplus and lowest-to-highest surplus); the remaining
#' \code{shuffles - 2} passes use independently randomized node orderings.
#' Increasing \code{shuffles} improves the chance of finding a higher-
#' clustering solution at the cost of proportionally more computation.
#' Defaults to \code{100}.
#'
#' @return A square symmetric binary adjacency matrix of the same dimension
#' as \code{network}, in the original node ordering of \code{network},
#' representing the resulting ring lattice. The average clustering coefficient
#' of the returned lattice is attached as the attribute \code{"CC"} and can
#' be retrieved with \code{attr(result, "CC")}.
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
#' steps; results are un-permuted before return so that output rows and
#' columns correspond to the same nodes as the input \code{network}. The
#' target degree of node \eqn{i} (after reordering) is
#' \eqn{d_i = \sum_j A_{ij}}. A ring lattice is seeded by connecting node
#' \eqn{i} to node \eqn{j} whenever the circular ring distance does not
#' exceed \eqn{\lceil (d_i + d_j) / 2 \rceil}, so every node begins with at
#' least as many edges as its target degree. The diagonal is set to
#' \code{FALSE} before pruning to prevent self-loops. Nodes with
#' \eqn{d_i = 0} are immediately disconnected.
#'
#' \strong{Distance matrix.} After the seed ring is constructed, the pruning
#' distance matrix is replaced with the Euclidean distances between rows of
#' the seed adjacency matrix:
#' \deqn{D_{ij} = \lVert \mathbf{r}_i - \mathbf{r}_j \rVert_2}
#' where \eqn{\mathbf{r}_i} is the \eqn{i}-th row of the seed ring. This
#' measures structural dissimilarity between neighborhoods in the seed.
#'
#' \strong{Surplus degree.} The surplus of node \eqn{i} is
#' \eqn{e_i = \hat{d}_i - d_i \geq 0}, where \eqn{\hat{d}_i} is the current
#' ring degree. Pruning continues while \eqn{\sum_i e_i > 0} or until
#' \eqn{p \times 2} outer iterations have elapsed without change.
#'
#' \strong{Multi-pass pruning.} \code{shuffles} independent pruning passes
#' are run on the same seed lattice via compiled C code:
#' \itemize{
#'   \item \emph{Highest-first pass} — nodes visited in decreasing order of
#'     surplus.
#'   \item \emph{Lowest-first pass} — nodes visited in increasing order of
#'     surplus.
#'   \item \emph{Random passes} — \code{shuffles - 2} additional passes each
#'     visit nodes in an independently randomised order.
#' }
#' Within each pass, for an over-degree node \eqn{i}, candidate removal
#' partners are current neighbors that also carry surplus (\eqn{e_j > 0}).
#' The partner to remove is selected by the following lexicographic priority:
#' \enumerate{
#'   \item \strong{Clustering constraint} — prefer the neighbour with the
#'     fewest common neighbors with \eqn{i} (minimise \eqn{c_{ij}}), so
#'     that removing the edge disrupts local clustering as little as possible.
#'   \item \strong{Surplus} — among equally ranked candidates, prefer the
#'     neighbor with the largest surplus \eqn{e_j}.
#'   \item \strong{Structural distance} — as a final tiebreaker, prefer the
#'     neighbor with the largest \eqn{D_{ij}}, removing structurally
#'     dissimilar connections first. On a complete tie the first-encountered
#'     candidate is retained, consistent with R's \code{order()} stability.
#' }
#' Surplus values are updated incrementally after each edge removal (O(1)
#' per removal) and convergence is checked against a running scalar total.
#'
#' \strong{Pass selection.} After all \code{shuffles} passes complete, the
#' subset with the smallest total residual surplus \eqn{\sum_i |e_i|} is
#' identified. If multiple passes share that minimum, the one with the
#' highest average clustering coefficient is chosen.
#'
#' \strong{Empirical fallback.} If the average clustering coefficient of the
#' winning lattice is lower than that of the original empirical network,
#' the empirical adjacency is returned as the lattice baseline directly in
#' its original node ordering.
#'
#' \strong{Node ordering.} All internal computation uses degree-descending
#' node ordering. The returned matrix is un-permuted via
#' \code{original_order <- order(degree_order)} so that output rows and
#' columns correspond to the same nodes as the input \code{network}.
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
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Near-degree preserving lattice ----
# Updated 22.03.2026
ring2lattice <- function(network, shuffles = 100)
{

  # Automatically construct adjacency
  A <- network != 0

  # Get nodes and density
  nodes <- dim(A)[2]
  degree <- colSums(A)

  # Set degree order
  degree_order <- order(degree, decreasing = TRUE)
  degree <- degree[degree_order]
  original_order <- order(degree_order)

  # Set up distance matrix
  node_sequence <- seq_len(nodes)
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  ring <- distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Convert ring matrix to zeros and ones
  ring <- ring <= (ceiling(outer(degree, degree, "+") / 2))

  # Set diagonal to FALSE
  diag(ring) <- FALSE

  # Update distance matrix using Euclidean distance
  distance_matrix <- as.matrix(dist(ring))

  # Zero out zero degree nodes
  zero_degree <- degree == 0
  if(any(zero_degree)){
    ring[zero_degree,] <- ring[,zero_degree] <- FALSE
  }

  # Compute ring degree
  ring_degree <- colSums(ring)

  # Extra per node
  extra <- ring_degree - degree
  total_extra <- sum(extra)

  # Set up maximum iterations
  max_iter <- nodes * 2

  # Create orders
  node_orders <- c(
    list(
      # highest-to-lowest
      ring_pruning(
        nodes = nodes, ring_degree = ring_degree, extra = extra,
        total_extra = total_extra, ring = ring, distance_matrix = distance_matrix,
        max_iter = max_iter, initial_order = order(extra, decreasing = TRUE)
      ),
      # lowest-to-highest
      ring_pruning(
        nodes = nodes, ring_degree = ring_degree, extra = extra,
        total_extra = total_extra, ring = ring, distance_matrix = distance_matrix,
        max_iter = max_iter, initial_order = order(extra, decreasing = FALSE)
      )
    ),
    lapply(seq_len(shuffles - 2), function(i){

      ring_pruning(
        nodes = nodes, ring_degree = ring_degree, extra = extra,
        total_extra = total_extra, ring = ring, distance_matrix = distance_matrix,
        max_iter = max_iter, initial_order = shuffle(node_sequence)
      )

    })
  )

  # Identify best passes
  all_extras <- nvapply(node_orders, function(x){sum(x$extra)})
  minimum <- which.min(abs(all_extras))
  extra_minimum <- all_extras[minimum]
  best_passes <- node_orders[all_extras == extra_minimum]

  # Compute CCs
  best_ccs <- nvapply(best_passes, function(x){
    igraph::transitivity(convert2igraph(x$ring), type = "average")
  })

  # Collect maximum
  ring <- best_passes[[which.max(best_ccs)]]$ring

  # Warn if complete degree-preservation was not found
  if(extra_minimum != 0){

    warning(paste0(
      "Exact degree-preservation was not found. There's a difference of ",
      extra_minimum, " from the ring's degree minus the network's degree."
    ))

  }

  # Compute clustering coefficients
  lattice_CC <- igraph::transitivity(convert2igraph(ring), type = "average")

  # Compute empirical clustering coefficient
  empirical_CC <- igraph::transitivity(convert2igraph(A), type = "average")

  # Fallback: If ringed solution cannot beat empirical,
  # then empirical may be most optimized lattice already
  empirical_flag <- empirical_CC > lattice_CC
  ring <- swiftelse(
    empirical_flag, A, # return empirical
    ring[original_order, original_order] # return ring
  )

  # Warn on fallback
  if(empirical_flag){

    warning(paste0(
      "The lattice solution did not produce a better CC (", round(lattice_CC, 3),
      ") than the empirical (", round(empirical_CC, 3), ").\n",
      "Falling back to empirical solution..."
    ))

  }

  # Ensure named matrix
  dimnames(ring) <- dimnames(network)

  # Attach clustering coefficient
  attr(ring, "CC") <- swiftelse(empirical_flag, empirical_CC, lattice_CC)

  # Return ring
  return(ring)

}

#' @noRd
# Core algorithm ----
# Updated 21.03.2026
ring_pruning <- function(nodes, ring_degree, extra, total_extra, ring, distance_matrix, max_iter, initial_order)
{

  # # Initialize iterations
  # iter <- 0
  #
  # # Loop while there are still some left over
  # while(iter < max_iter){
  #
  #   # Increase iterations
  #   iter <- iter + 1
  #
  #   # Obtain previous total
  #   previous_total <- total_extra
  #
  #   # Loop over extra
  #   for(i in initial_order){
  #
  #     # Skip if already satisfied by a previous iteration
  #     if(extra[i] <= 0){
  #       next
  #     }
  #
  #     # Current neighbors with excess edges
  #     neighbors <- which(ring[i,] & (extra > 0))
  #
  #     # Skip if there isn't any available
  #     if(length(neighbors) == 0){
  #       next
  #     }
  #
  #     # Clustering constraint
  #     common_neighbors <- nvapply(neighbors, function(j){sum(ring[i,] & ring[j,])})
  #
  #     # Prioritize clustering, excess, and then distance
  #     index <- order(
  #       1 / (common_neighbors + 1), extra[neighbors],
  #       distance_matrix[i, neighbors], # distance constraint
  #       decreasing = TRUE
  #     )[1]
  #     target <- neighbors[index]
  #
  #     # Loop over target to update distance and edges
  #     ring[target, i] <- ring[i, target] <- FALSE
  #
  #     # Update state
  #     update <- c(i, target)
  #     extra[update] <- extra[update] - 1
  #     total_extra <- total_extra - 2
  #
  #   }
  #
  #   # Check for no change
  #   if((total_extra <= 0) || (previous_total == total_extra)){
  #     break
  #   }
  #
  # }

  # Obtain result
  result <- .Call(
    "ring_pruning_c",
    as.integer(extra),
    as.integer(total_extra),
    as.integer(ring),
    as.numeric(distance_matrix),
    as.integer(max_iter),
    as.integer(initial_order),
    PACKAGE = "L0ggm"
  )

  # Update ring
  result$ring <- matrix(as.logical(result$ring), nrow = nodes, ncol = nodes)

  # Return the output
  return(result)

}
