#' @title Computes Various Small-Worldness Metrics
#'
#' @description Computes the small-worldness of a network using one of four
#' methods. All simulation-based methods generate degree-preserving random
#' graphs as the null baseline. The lattice reference network (where required)
#' is constructed using \code{\link[L0ggm]{ring2lattice}}, which produces a
#' degree-preserving ring lattice that maximizes the clustering coefficient
#'
#' @param A Matrix or data frame.
#' An adjacency matrix representing the network.
#' Must be square, symmetric, and contain only non-negative values
#'
#' @param lattice Matrix (optional).
#' A pre-computed lattice adjacency matrix to use as the regular-network
#' reference for the \code{"omega"} and \code{"SWP"} methods.
#' If not provided, a lattice is generated automatically via
#' \code{\link[L0ggm]{ring2lattice}}
#'
#' @param method Character (length = 1).
#' The method used to compute small-worldness.
#' Defaults to \code{"SWP"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"analytical"} --- Computes the Humphries & Gurney (2008)
#' \eqn{S} metric using closed-form approximations for the Erdős–Rényi
#' random graph baseline:
#' \deqn{S = \frac{C / C_{\text{rand}}}{L / L_{\text{rand}}}}
#' where \eqn{C_{\text{rand}} \approx \langle k \rangle / n} and
#' \eqn{L_{\text{rand}} \approx \ln(n) / \ln(\langle k \rangle)}.
#' No random graphs are generated. Values greater than 1 indicate
#' small-world structure
#'
#' \item \code{"S"} --- Computes the Humphries & Gurney (2008) \eqn{S}
#' metric by simulation. Random graphs are generated under the
#' Erdős–Rényi model (\code{\link[igraph]{sample_gnm}}) with the same
#' number of nodes and edges as \code{A}. The clustering coefficient
#' uses global transitivity (\eqn{C^\Delta}). Values greater than 1
#' indicate small-world structure
#'
#' \item \code{"omega"} --- Computes the \eqn{\omega} metric of
#' Telesford et al. (2011):
#' \deqn{\omega = \frac{L_{\text{rand}}}{L} - \frac{C}{C_{\text{latt}}}}
#' Random graphs preserve the degree sequence via edge-switching
#' (\code{\link[igraph]{sample_degseq}}, \code{method = "edge.switching.simple"}).
#' Values near zero indicate small-world structure; negative values indicate
#' lattice-like structure; positive values indicate random-like structure.
#' Bounded in \eqn{[-1, 1]}
#'
#' \item \code{"SWP"} --- Computes the Small-World Propensity of
#' Muldoon et al. (2016):
#' \deqn{\phi = 1 - \sqrt{\frac{\Delta_C^2 + \Delta_L^2}{2}}}
#' where \eqn{\Delta_C = (C_{\text{latt}} - C) / (C_{\text{latt}} - C_{\text{rand}})}
#' and \eqn{\Delta_L = (L - L_{\text{rand}}) / (L_{\text{latt}} - L_{\text{rand}})}.
#' Both \eqn{\Delta_C} and \eqn{\Delta_L} are bounded to \eqn{[0, 1]},
#' guaranteeing \eqn{\phi \in [0, 1]}. Random graphs preserve the degree
#' sequence via edge-switching. Values close to 1 indicate strong
#' small-world structure; a pragmatic threshold of \eqn{\phi_T = 0.60}
#' has been suggested
#'
#' }
#'
#' @param iter Numeric (length = 1).
#' Number of random graphs to generate when estimating the null baseline.
#' Defaults to \code{100}.
#' Ignored when \code{method = "analytical"}.
#' Higher values produce more stable estimates at the cost of computation time
#'
#' @return Numeric (length = 1).
#' The small-worldness value computed using the chosen \code{method}:
#'
#' \itemize{
#'
#' \item \code{"analytical"} and \code{"S"} --- \eqn{S > 1} indicates small-world structure
#'
#' \item \code{"omega"} --- Values near \eqn{|\omega| < 0.50} indicate small-world structure
#' (ranges between -1 and 1)
#'
#' \item \code{"SWP"} --- Values \eqn{\phi > 0.60} indicate strong small-world structure
#' (ranges between 0 and 1)
#'
#' }
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{omega} \cr
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367--375.
#'
#' \strong{S and analytical} \cr
#' Humphries, M. D., & Gurney, K. (2008).
#' Network 'small-world-ness': A quantitative method for determining canonical network equivalence.
#' \emph{PLoS ONE}, \emph{3}(4), e0002051.
#'
#' \strong{SWP} \cr
#' Muldoon, S. F., Bridgeford, E. W., & Bassett, D. S. (2016).
#' Small-world propensity and weighted brain networks.
#' \emph{Scientific Reports}, \emph{6}(1), 22057.
#'
#' @examples
#' # Get network
#' network <- network_estimation(basic_smallworld)
#'
#' # Compute SWP (default)
#' swp <- smallworldness(A = network)
#'
#' # Compute omega
#' omega <- smallworldness(A = network, method = "omega")
#'
#' # Compute analytical S
#' S <- smallworldness(A = network, method = "analytical")
#'
#' # Compute simulated S with more iterations
#' S_sim <- smallworldness(A = network, method = "S", iter = 1000)
#'
#' @export
#'
# Compute smallworldness ----
# Updated 22.03.2026
smallworldness <- function(A, lattice, method = c("analytical", "omega", "S", "SWP"), iter = 100)
{

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  method <- set_default(method, "swp", smallworldness)

  # Ensure binary network
  A <- A != 0

  # Obtain edges
  degree <- colSums(A, na.rm = TRUE)

  # Convert network to {igraph}
  I <- convert2igraph(A)

  # Set up function switch
  smallworld_FUN <- switch(
    method,
    "analytical" = smallworldness_analytical,
    "omega" = smallworldness_omega,
    "s" = smallworldness_S,
    "swp" = smallworldness_SWP
  )

  # Obtain empirical values
  ASPL <- igraph::mean_distance(I)
  CC <- igraph::transitivity(I, type = swiftelse(method == "s", "global", "average"))
  # Switch CC based on method

  # Compute smallworldness
  return(
    smallworld_FUN(A = A, I = I, ASPL = ASPL, CC = CC, degree = degree, lattice = lattice, iter = iter)
  )

}

#' @noRd
# Analytical ----
# Updated 22.03.2026
smallworldness_analytical <- function(A, ASPL, CC, degree, ...)
{

  # Obtain nodes and average degree
  nodes <- dim(A)[2]
  mean_degree <- mean(degree, na.rm = TRUE)

  # Return analytical
  return(
    (CC / (mean_degree / nodes)) / (ASPL / (log(nodes) / log(mean_degree)))
  )

}

#' @noRd
# Telesford et al. (2011) ----
# Updated 22.03.2026
smallworldness_omega <- function(A, I, ASPL, CC, degree, lattice, iter, ...)
{

  # Collect random graphs
  random_graphs <- lapply(
    seq_len(iter), function(i){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }
  )

  # Compute ASPL
  random_ASPL <- mean(nvapply(random_graphs, igraph::mean_distance), na.rm = TRUE)

  # Check if lattice is missing
  if(missing(lattice)){
    lattice_CC <- attributes(ring2lattice(A))$CC
  }else{
    lattice_CC <- igraph::transitivity(convert2igraph(lattice), type = "average")
  }

  # Return omega
  return((random_ASPL / ASPL) - (CC / lattice_CC))

}

#' @noRd
# Humphries & Gurney (2008) ----
# Updated 22.03.2026
smallworldness_S <- function(A, I, ASPL, CC, iter, ...)
{

  # Obtain nodes and edges
  nodes <- dim(A)[2]
  edges <- sum(A[lower.tri(A)] != 0)

  # Obtain random graphs
  random_graphs <- lapply(seq_len(iter), function(i){
    igraph::sample_gnm(n = nodes, m = edges)
  })

  # Compute S
  S <- nvapply(
    random_graphs, function(x){
      (CC / igraph::transitivity(x, type = "global")) / (ASPL / igraph::mean_distance(x))
    }
  )

  # Return S
  return(mean(S, na.rm = TRUE))

}

#' @noRd
# Boundaries for SWP ----
# Updated 22.03.2026
boundary_values <- function(empirical, random, regular)
{
  any(is.infinite(c(empirical, random, regular))) | any(is.na(c(empirical, random, regular)))
}

#' @noRd
# Muldoon et al. (2016) ----
# Updated 22.03.2026
smallworldness_SWP <- function(A, I, ASPL, CC, degree, lattice, iter, ...)
{

  # Collect random graphs
  random_graphs <- lapply(
    seq_len(iter), function(i){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }
  )

  # Compute ASPL and CC
  random_ASPL <- mean(nvapply(random_graphs, igraph::mean_distance), na.rm = TRUE)
  random_CC <- mean(nvapply(random_graphs, igraph::transitivity, type = "average"), na.rm = TRUE)

  # Check if lattice is missing
  if(missing(lattice)){
    lattice <- ring2lattice(A)
  }

  # Ensure lattice is {igraph}
  lattice <- convert2igraph(lattice)

  # Compute ASPL and CC
  lattice_ASPL <- igraph::mean_distance(lattice)
  lattice_CC <- igraph::transitivity(lattice, type = "average")

  # Compute deltas
  delta_ASPL <- max(0, (ASPL - random_ASPL)) / (lattice_ASPL - random_ASPL)
  delta_CC <- max(0, (lattice_CC - CC)) / (lattice_CC - random_CC)

  # Boundary check
  delta_ASPL <- min(1, swiftelse(boundary_values(ASPL, random_ASPL, lattice_ASPL), 1, delta_ASPL))
  delta_CC <- min(1, swiftelse(boundary_values(CC, random_CC, lattice_CC), 1, delta_CC))

  # Return SWP
  return(1 - sqrt((delta_ASPL^2 + delta_CC^2) / 2))

}