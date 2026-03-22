
/*
 * ring_pruning.c
 * --------------
 * C implementation of the ring lattice pruning algorithm
 *
 * Matrices are passed from R as flat integer vectors in column-major order
 * (R's native layout). Because ring and distance_matrix are symmetric,
 * element [i, j] (0-indexed) is accessed as flat[i + j * n] throughout.
 *
 * Tie-breaking note
 * -----------------
 * When two candidate neighbours j1 < j2 are equal on all three keys
 * (common neighbors, surplus, ring distance), the comparison uses > for
 * the final distance check so that j2 (the FIRST found) wins.  This matches
 * R's `order` behavior.
 */
 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
 
SEXP ring_pruning_c(
  	SEXP extra_r, SEXP total_extra_r, SEXP ring_r, 
  	SEXP dist_matrix_r, SEXP max_iter_r, SEXP initial_order_r
){

    // Get number of nodes from length of extra vector
    int n = length(extra_r);
 
    // Get maximum number of iterations
    int max_iter = asInteger(max_iter_r);
 
    // Get total surplus across all nodes
    int total_extra = asInteger(total_extra_r);
 
    // Allocate output buffer for surplus vector
    SEXP extra_out = PROTECT(allocVector(INTSXP, n));
 
    // Allocate output buffer for adjacency matrix
    SEXP ring_out = PROTECT(allocVector(INTSXP, n * n));
 
    // Get pointer to output surplus buffer
    int *extra = INTEGER(extra_out);
 
    // Get pointer to output adjacency matrix buffer
    int *ring = INTEGER(ring_out);
 
    // Get pointer to input surplus vector
    int *ei = INTEGER(extra_r);
 
    // Copy input surplus values into output buffer
    for (int k = 0; k < n; k++) extra[k] = ei[k];
 
    // Get pointer to input adjacency matrix
    int *ri = INTEGER(ring_r);
 
    // Copy input adjacency matrix into output buffer
    for (int k = 0; k < n * n; k++) ring[k] = ri[k];
 
    // Get read-only pointer to distance matrix
    double *dist = REAL(dist_matrix_r);
 
    // Get read-only pointer to node visit order (1-indexed from R)
    int *ord = INTEGER(initial_order_r);
 
    // Iterate up to the maximum number of pruning passes
    for (int iter = 0; iter < max_iter; iter++) {
 
        // Record total surplus at the start of this pass for convergence check
        int previous_total = total_extra;
 
        // Visit each node in the prescribed order
        for (int oi = 0; oi < n; oi++) {
 
            // Convert 1-indexed R order to 0-indexed C index
            int i = ord[oi] - 1;
 
            // Skip nodes with no surplus
            if (extra[i] <= 0) continue;
 
            // Initialize best candidate index to root
            int best_j = -1;
 
            // Initialize best common-neighbor count to root
            int best_common = -1;
 
            // Initialize best surplus value to root
            int best_xtra = -1;
 
            // Initialize best ring distance to root
            double best_dist = -1.0;
 
            // Scan all nodes as potential removal targets
            for (int j = 0; j < n; j++) {
 
                // Skip self-loop
                if (j == i) continue;
 
                // Skip non-neighbours
                if (!ring[i + j * n]) continue;
 
                // Skip neighbours with no surplus
                if (extra[j] <= 0) continue;
 
                // Initialize common-neighbor counter
                int common = 0;
 
                // Count nodes adjacent to both i and j
                for (int k = 0; k < n; k++) {
                    if (ring[i + k * n] && ring[j + k * n]) common++;
                }
 
                // Initialize update flag
                int better = 0;
 
                // Always update when no candidate has been found yet
                if (best_j < 0) {
                    better = 1;
 
                // Prefer fewer common neighbors
                } else if (common < best_common) {
                    better = 1;
 
                } else if (common == best_common) {
 
                    // Prefer higher surplus
                    if (extra[j] > best_xtra) {
                        better = 1;
 
                    } else if (extra[j] == best_xtra) {
 
                        // Prefer larger ring distance; > so first-found wins on a full tie
                        if(dist[i + j * n] > best_dist) {
                            better = 1;
                        }
                    }
                }
 
                // Update best candidate if current node is better
                if (better) {
                    best_j = j;
                    best_common = common;
                    best_xtra = extra[j];
                    best_dist = dist[i + j * n];
                }
            }
 
            // Skip if no eligible neighbor was found
            if (best_j < 0) continue;
 
            // Remove edge from i to best_j
            ring[i + best_j * n] = 0;
 
            // Remove edge from best_j to i (symmetric)
            ring[best_j + i * n] = 0;
 
            // Decrement surplus for node i
            extra[i]--;
 
            // Decrement surplus for removed neighbor
            extra[best_j]--;
 
            // Update global surplus count
            total_extra -= 2;
        }
 
        // Stop if all surplus is eliminated
        if (total_extra <= 0) break;
 
        // Stop if no edges were removed this pass (convergence)
        if (total_extra == previous_total) break;
    }
 
    // Allocate named return list with two elements
    SEXP result = PROTECT(allocVector(VECSXP, 2));
 
    // Allocate character vector for list names
    SEXP names = PROTECT(allocVector(STRSXP, 2));
 
    // Set first list name to "ring"
    SET_STRING_ELT(names, 0, mkChar("ring"));
 
    // Set second list name to "extra"
    SET_STRING_ELT(names, 1, mkChar("extra"));
 
    // Store adjacency matrix in result list
    SET_VECTOR_ELT(result, 0, ring_out);
 
    // Store surplus vector in result list
    SET_VECTOR_ELT(result, 1, extra_out);
 
    // Attach names attribute to result list
    setAttrib(result, R_NamesSymbol, names);
 
    // Unprotect all four PROTECT calls: extra_out, ring_out, result, names
    UNPROTECT(4);
 
    // Return named list
    return result;
    
}
