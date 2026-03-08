#%%%%%%%%%%%%%%%%%%%%%%%#
#### {L0ggm} Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%#

# These helpers are internal functions that are repeatedly
# used across the package. They are implemented here once
# so that changes only need to be made here and not in
# multiple places across the package.
#
# The goal is to minimize redundant calls to improve reliability of code.
# Sections are separated by how they are associated with their utility
# in other functions.

#%%%%%%%%%%%%%%%%%%%%%%
# *APPLY FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

# These functions are merely wrappers over the `*apply`
# family to execute specific tasks that are often repeatedly
# done. Their purpose are to provide slightly cleaner code.

# These `*vapply` are specific versions of `vapply` that
# pre-specify the output. These versions are generally
# simpler (in function) and marginally faster than
# the more commonly used `sapply`. They are more strict
# in that they require knowing the output format (e.g., numeric)

#' @noRd
# Character vector/matrix output apply ----
# Updated 06.07.2023
cvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = character(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Logical vector/matrix output apply ----
# Updated 06.07.2023
lvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = logical(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Numeric vector/matrix output apply ----
# Updated 06.07.2023
nvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = numeric(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Unlist `lapply` ----
# Updated 14.07.2023
ulapply <- function(X, FUN, ..., recursive = TRUE)
{
  return(unlist(lapply(X, FUN, ...), recursive = recursive))
}

# The `row_apply` and `column_apply` functions are
# not necessary faster than `apply` but they are
# more strict in that they return a matrix with
# the same dimensions that was input into 'X'.
# Therefore, there is never a need to transpose
# the output because it will always have the same
# dimensions as what went in

#' @noRd
# Stricter `apply` by row ----
# Slightly faster than `apply` but not appreciably
# Updated 06.07.2023
row_apply <- function(X, FUN, ...)
{

  # Ensure matrix
  X <- as.matrix(X)

  # Get dimensions of data
  dimensions <- dim(X)

  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X), "\"",
        " is not available for `row_apply`. Only \"character\", \"logical\", and \"numeric\" are available"
      ), call. = FALSE
    )
  )

  # Return call
  return(
    matrix(
      vapply_FUN(
        split(X, seq_len(dimensions[1])),
        FUN = FUN, ..., LENGTH = dimensions[2],
        USE.NAMES = FALSE # names are supplied below
      ),
      dimensions, dimnames = dimnames(X), byrow = TRUE
    )
  )

}

#' @noRd
# Stricter `apply` by column ----
# Slightly faster than `apply` but not appreciably
# Updated 06.07.2023
column_apply <- function(X, FUN, ...)
{

  # Ensure matrix
  X <- as.matrix(X)

  # Get dimensions of data
  dimensions <- dim(X)

  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X), "\"",
        " is not available for `row_apply`. Only \"character\", \"logical\", and \"numeric\" are available"
      ), call. = FALSE
    )
  )

  # Return call
  return(
    matrix(
      vapply_FUN(
        split(t(X), seq_len(dimensions[2])),
        FUN = FUN, ..., LENGTH = dimensions[1],
        USE.NAMES = FALSE # names are supplied below
      ),
      dimensions, dimnames = dimnames(X), byrow = FALSE
    )
  )

}

#' @noRd
# 3D Array apply ----
# Apply's a function across a 3D array to obtain a 2D array
# About 2x faster than `apply(simplify2array(X), 1:2, FUN)`
# Updated 06.07.2023
symmetric_matrix_lapply <- function(X, FUN, ...){

  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X[[1]]),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X[[1]]), "\"",
        " is not available for `symmetric_matrix_lapply`. Only \"character\", \"logical\", and \"numeric\" are available"
      ), call. = FALSE
    )
  )

  # Get dimensions of single matrix in list
  dimensions <- dim(X[[1]])

  # Pre-obtain lower triangle indices
  # Equivalent to: (`lower.tri(matrix, diag = TRUE)`)
  lower_triangle_index <- .row(dimensions) >= .col(dimensions)

  # Get lower triangles
  lower_matrix <- do.call(rbind, lapply(X, function(x){x[lower_triangle_index]}))

  # Apply function
  values <- vapply_FUN(as.data.frame(lower_matrix), FUN, ...)

  # Initialize return matrix
  new_matrix <- matrix(
    nrow = dimensions[1], ncol = dimensions[2],
    dimnames = dimnames(X[[1]])
  )

  # Add values to lower triangle
  new_matrix[lower_triangle_index] <- values

  # Transpose
  new_matrix <- t(new_matrix)

  # Add values again to lower triagle
  new_matrix[lower_triangle_index] <- values

  # Return new matrix
  return(new_matrix)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# REPRODUCIBLE RANDOM GENERATION ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For these random generation functions, the goal
# is reproducibility, which runs counter to "random".
# However, for bootstrap results to be reproduced
# there is a need to generate "random" samples that
# can be reproduced so that the user arrives at the
# same results *without* influencing R's random number
# generation or changing a user-defined seed previous
# to, while, or after using {EGAnet} functions.
#
# Assistance in writing C and interfacing with R
# was provided by GPT-4

#' @noRd
# Generate integer seeds ----
# (non-randomly based on `seed`)
# `NULL` is reserved for actual pseudo-randomness
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 27.07.2023
reproducible_seeds <- function(n, seed = NULL)
{

  # Return call from C
  return(
    .Call(
      "r_xoshiro_seeds",
      as.integer(n), swiftelse(is.null(seed), 0, seed),
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Generate uniform data ----
# Allows adjustment of range
# Updated 26.10.2023
runif_xoshiro <- function(n, min = 0, max = 1, seed = NULL)
{

  # Get values
  values <- .Call(
    "r_xoshiro_uniform",
    as.integer(n),
    swiftelse(is.null(seed), 0, seed),
    PACKAGE = "L0ggm"
  )

  # Check for changes to minimum and maximum
  if(min != 0 || max != 1){ # transform
    values <- min + (max - min) * values
  }


  # Return call from C
  return(values)

}

#' @noRd
# Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 30.07.2023
shuffle <- function(x, size = length(x), seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle",
      as.integer(seq_along(x)),
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "L0ggm"
    )][seq_len(size)]
  )

}

#' @noRd
# Shuffle (with replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 30.07.2023
shuffle_replace <- function(x, size = length(x), seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle_replace",
      x, swiftelse(is.null(seed), 0, seed),
      PACKAGE = "L0ggm"
    )][seq_len(size)]
  )

}

#' @noRd
# Random normal generation with Ziggurat ----
# https://people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat/ziggurat.html
# Updated 27.07.2023
rnorm_ziggurat <- function(n, seed = NULL)
{

  # Return call from C
  return(
    .Call(
      "r_ziggurat",
      as.integer(n),
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "L0ggm"
    )
  )

}

#' @noRd
# Generate multivariate normal data ----
# Removes MASS package dependency (from version 7.3.60)
# Original code here for legacy and bug checking
# Updated 03.07.2023
MASS_mvrnorm <- function(
    n = 1, mu, Sigma,
    tol = 1e-06, empirical = FALSE, EISPACK = FALSE
)
{

  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)

}

#' @noRd
# Generate multivariate normal data (quicker) ----
# A very aggressive optimization
# Removes 'tol', 'empirical' and 'EISPACK' arguments
# Removes safeguards like positive definite (should be from `auto_correlate`)
# 'empirical' is ALWAYS = FALSE
# Pre-computes `np` (n * p)
# Avoids repeated calls to `eigen` in `MASS_mvrnorm`
# Pre-computes `eS$vectors %*% diag(sqrt(pmax(ev, 0)), p)`
# Updated 27.07.2023
mvrnorm_precompute <- function(cases, Sigma)
{

  # Get p (number of variables)
  p <- dim(Sigma)[2]

  # Compute eigenvalues
  eS <- eigen(Sigma, symmetric = TRUE)

  # Get eigenvalues
  eigenvalues <- eS$values

  # Get non-zero eigenvalues
  non_zero <- eigenvalues > 0

  # Set negative eigenvalues to zero
  eigenvalues[!non_zero] <- 0

  # Get square root of non-zero eigenvalues
  eigenvalues[non_zero] <- sqrt(eigenvalues[non_zero])

  # Return pre-computed values
  return(
    list(
      p = p,
      np = cases * p,
      coV = eS$vectors %*% diag(eigenvalues, p)
    )
  )

}

#' @noRd
# Generate multivariate normal data (quick) ----
# Updated 27.07.2023
MASS_mvrnorm_quick <- function(seed = NULL, p, np, coV)
{
  return(t(tcrossprod(coV, matrix(rnorm_ziggurat(np, seed), ncol = p))))
}

#%%%%%%%%%%%%%%%%%%%%%
# PARALLELIZATION ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Clear memory ----
# Updated 07.11.2023
clear_memory <- function()
{
  sink <- capture.output(
    gc(verbose = FALSE, reset = TRUE, full = TRUE)
  )
}

#' @noRd
# Get available memory ----
# Updated 18.10.2023
available_memory <- function()
{

  # Get operating system
  OS <- tolower(Sys.info()["sysname"])

  # Branch based on OS
  if(OS == "windows"){ # Windows

    # Alternative (outputs memory in kB)
    bytes <- as.numeric(
      trimws(system("wmic OS get FreePhysicalMemory", intern = TRUE))[2]
    ) * 1e+03

  }else if(OS == "linux"){ # Linux

    # Split system information
    info_split <- strsplit(system("free", intern = TRUE), split = " ")

    # Remove "Mem:" and "Swap:"
    info_split <- lapply(info_split, function(x){gsub("Mem:", "", x)})
    info_split <- lapply(info_split, function(x){gsub("Swap:", "", x)})

    # Get actual values
    info_split <- lapply(info_split, function(x){x[x != ""]})

    # Bind values
    info_split <- do.call(rbind, info_split[1:2])

    # Get free values (Linux reports in *kilo*bytes -- thanks, Aleksandar Tomasevic)
    bytes <- as.numeric(info_split[2, info_split[1,] == "available"]) * 1e+03

  }else{ # Mac

    # System information
    system_info <- system("top -l 1 -s 0 | grep PhysMem", intern = TRUE)

    # Get everything after comma
    unused <- gsub(" .*,", "", system_info)

    # Get values only
    value <- gsub(" unused.", "", gsub("PhysMem: ", "", unused))

    # Check for bytes
    if(grepl("M", value)){
      bytes <- as.numeric(gsub("M", "", value)) * 1e+06
    }else if(grepl("G", value)){
      bytes <- as.numeric(gsub("G", "", value)) * 1e+09
    }else if(grepl("K", value)){
      bytes <- as.numeric(gsub("K", "", value)) * 1e+03
    }else if(grepl("B", value)){ # edge case
      bytes <- as.numeric(gsub("B", "", value)) * 1
    }else if(grepl("T", value)){ # edge case
      bytes <- as.numeric(gsub("T", "", value)) * 1e+12
    }

  }

  # Return bytes
  return(bytes)

}

#' @noRd
# Store random state ----
# Updated 04.08.2023
store_state <- function()
{
  if(exists(".Random.seed", envir = parent.frame(2))){
    assign("saved_state", .Random.seed, envir = parent.frame(2))
  }
}

#' @noRd
# Restore and remove random state ----
# Updated 04.08.2023
restore_state <- function()
{
  if(exists("saved_state", envir = parent.frame(2))){
    saved_state <- get("saved_state", envir = parent.frame(2))
    assign(".Random.seed", saved_state, envir = parent.frame(2))
    rm("saved_state", envir = parent.frame(2))
  }
}

#%%%%%%%%%%%%%%%%%%%%%%
# FASTER SEQUENCES ----
#%%%%%%%%%%%%%%%%%%%%%%

# These functions aren't often used because
# dimensions of the data are usually used
# elsewhere in {EGAnet}; however, these
# functions use `seq_len` and `dim` which
# have minimal speed advantages over their
# respective `1:nrow` and `1:ncol`
#
# The primary benefit of using this functions
# is `seq_len` which throws an error with
# `seq_len(0)` where `1:0` will not
# (at least not at the top of the call)

#' @noRd
# Faster row sequences ----
# Updated 06.07.2023
nrow_sequence <- function(data)
{
  return(seq_len(dim(data)[1]))
}


#' @noRd
# Faster column sequences ----
# Updated 06.07.2023
ncol_sequence <- function(data)
{
  return(seq_len(dim(data)[2]))
}

#' @noRd
# Faster `ifelse` ----
# For single value replacements:
# 1.5x faster with 1 value
# 2.5x faster with 10 values
# >= 18x faster with >= 100 values
# Updated 26.02.2026
swiftelse <- function(condition, true, false)
{
  # Get condition length
  condition_length <- length(condition)

  # Check for single value
  if(condition_length == 1){

    # If TRUE
    if(condition){
      return(true)
    }else{ # Otherwise, FALSE
      return(false)
    }

  }

  # Initialize result
  result <- vector(typeof(false), condition_length)

  # Set TRUE condition
  if(length(true) == 1){
    result[condition] <- true
  }else{
    result[condition] <- true[condition]
  }

  # Set FALSE condition
  if(length(false) == 1){
    result[!condition] <- false
  }else{

    # Get opposite condition (slightly faster than repeated calls)
    opposite_condition <- !condition
    result[opposite_condition] <- false[opposite_condition]
  }


  # Return result
  return(result)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# FORMATTING FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine decimal places before non-zero ----
# From StackOverflow: https://stackoverflow.com/questions/35553244/count-leading-zeros-between-the-decimal-point-and-first-nonzero-digit
# Updated 10.11.2023
leading_zero <- function(number)
{
  floor(-log10(.Machine$double.eps + abs(number) - floor(abs(number))))
}

#' @noRd
# Determine number of digits in a number ----
# Updated 24.07.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(abs(number))) + 1)

}

#' @noRd
# Determine format of bytes ----
# Updated 22.07.2023
byte_digits <- function(bytes)
{

  # Get number of digits
  ndigits <- digits(bytes)

  # Loop over values
  if(ndigits > 11){
    value <- paste0(round(bytes / 1e+12, 2), " TB")
  }else if(ndigits > 8){
    value <- paste0(round(bytes / 1e+09, 2), " GB")
  }else if(ndigits > 5){
    value <- paste0(round(bytes / 1e+06, 2), " MB")
  }else if(ndigits > 2){
    value <- paste0(round(bytes / 1e+03, 2), " KB")
  }else{
    value <- paste0(round(bytes, 2), " B")
  }

  # Return configured value
  return(value)

}

#' @noRd
# Format number with certain decimals ----
# Mainly for naming and printing
# Updated 14.06.2024
format_decimal <- function(numbers, places)
{

  return(
    formatC(
      x = numbers,
      digits = places,
      format = "f", flag = "0"
    )
  )

}

#' @noRd
# Format number with certain integer ----
# Mainly for naming and printing
# Updated 14.06.2024
format_integer <- function(numbers, places)
{

  return(
    formatC(
      x = numbers,
      digits = places,
      format = "d", flag = "0"
    )
  )

}

#%%%%%%%%%%%%%%%%%%%%%%
# OBJECT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Faster data frame initialization ----
# Initializes a matrix and converts to a data frame
# Data frames are initialized as individual vectors
# carrying extra overhead
# Updated 06.07.2023
fast.data.frame <- function(
    data = NA, nrow = 1, ncol = 1,
    rownames = NULL, colnames = NULL,
    ...
)
{

  return(
    as.data.frame(
      matrix(
        data = data,
        nrow = nrow, ncol = ncol,
        dimnames = list(rownames, colnames)
      ), make.names = FALSE,
      ...
    )
  )

}

#' @noRd
# Determine object type ----
# Updated 05.07.2023
get_object_type <- function(object)
{

  # Remove attributes and class from object
  unclassed_object <- remove_attributes(object)

  # Add object type to object's attributes
  ## In order of precedence
  if(is(object, "tbl")){ # needs precedence over 'data.frame'
    return("tibble")
  }else if(is.data.frame(object)){ # needs precedence over 'list'
    return("data.frame")
  }else if(is.list(unclassed_object)){ # needs precedence over 'vector'
    return("list")
  }else if(is.factor(object)){
    return("factor")
  }else if(is.vector(unclassed_object)){
    return("vector")
  }else if(is.matrix(unclassed_object)){
    return("matrix")
  }else if(is.array(unclassed_object)){
    return("array")
  }

}

#' @noRd
# Force vector ----
# For 1D matrix or data frames only
# Updated 05.07.2023
force_vector <- function(object)
{

  # Get object type
  object_type <- get_object_type(object)

  # Branch across types
  if(object_type == "vector"){
    return(object) # already a vector
  }else if(object_type %in% c("matrix", "data.frame", "tibble", "factor")){

    # Ensure matrix
    new_matrix <- as.matrix(object)

    # Get dimensions
    dimensions <- dim(new_matrix)

    # Check for more than one dimension
    if(all(dimensions) > 1){
      stop("Cannot sufficiently force object into 'vector'.", call. = FALSE)
    }

    # Set as vector
    new_vector <- as.vector(new_matrix)

    # Set names based on dimension
    names(new_vector) <- unlist(dimnames(object)[dimensions != 1])

    # Return new vector
    return(new_vector)

  }

}

#' @noRd
# Force numeric ----
# For factors
# Updated 05.07.2023
force_numeric <- function(object)
{

  # Check for whether object is already numeric
  if(is.numeric(object)){
    return(object) # just return
  }else if(is.factor(object)){ # factors are trickier...

    # Switch based on type in the levels
    if(is.numeric(levels(object))){

      # Get original values rather than factor values
      return(as.numeric(as.character(object)))

    }else{

      # Just used factor assigned values
      return(as.numeric(object))

    }

  }

}

#' @noRd
# Convert matrix to {igraph} network
# Updated 09.08.2023
convert2igraph <- function (A, diagonal = 0)
{

  # Convert to matrix
  A <- as.matrix(A)

  # Change diagonal (to zero)
  diag(A) <- diagonal

  # Return {igraph} network
  return(
    silent_call(
      igraph::graph_from_adjacency_matrix(
        A, weighted = TRUE, mode = "undirected",
        add.colnames = FALSE
      )
    )
  )

}

#%%%%%%%%%%%%%%%%%%%%%
# COUNT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Unique Length ----
# Counts length of unique variables (without `NA`)
# Updated 03.07.2023
unique_length <- function(x)
{
  return(length(unique(x)) - anyNA(x))
}

#' @noRd
# Edge count ----
# Counts the number of edges
# (assumes network is _symmetric_ matrix)
# Updated 11.10.2023
edge_count <- function(network, nodes, diagonal = FALSE)
{
  return((sum(network != 0, na.rm = TRUE) - swiftelse(diagonal, nodes, 0)) * 0.50)
}

#' @noRd
# Faster table ----
# Skips many checks in `table`
# Updated 22.07.2023
fast_table <- function(values)
{

  # Get factor values
  factor_values <- swiftelse(
    is.factor(values), values, factor(values)
  )

  # Get factor levels
  factor_levels <- levels(factor_values)

  # Return with names
  return(
    setNames(
      tabulate(factor_values, length(factor_levels)),
      factor_levels
    )
  )

}

#' @noRd
# Converts entire vectors to a single value ----
# Same logic as {plyr}'s `id_var` and `ninteraction`
# Avoids {plyr} dependency for single function
# Updated 06.07.2023
vector2factor <- function(data)
{

  # Use "paste" method
  stringed_values <- do.call(
    # Get total unique values each variable takes
    paste, rev(lapply( # `ninteraction`

      # Ensure data frame
      as.data.frame(data), function(x){

        # Get matched ID (`id_var`)
        matched_id <- match(x, sort(unique(x), na.last = TRUE))

        # Attached number of unique as attribute
        attr(matched_id, "n") <- max(matched_id)

        # Return result
        return(matched_id)

      }))
  )

  # Return unique ID
  return(factor(match(stringed_values, unique(stringed_values))))

}

#' @noRd
# Count table ----
# Provides counts of repeated rows in a data frame
# Same logic as {plyr}'s `count`
# Avoids {plyr} dependency for single function
# Updated 22.07.2023
count_table <- function(data, proportion = FALSE)
{

  # Ensure data frame
  data <- as.data.frame(data)

  # Get unique ID
  unique_id <- vector2factor(data)

  # Get unique solutions
  unique_solutions <- data[!duplicated(unique_id),,drop = FALSE]

  # Tabulate data
  Value <- fast_table(unique_id) # tabulate(unique_id, length(levels(unique_id)))

  # Check for proportion
  if(proportion){
    Value <- Value / dim(data)[1]
  }

  # Return data frame
  return(
    as.data.frame(
      cbind(unique_solutions[order(unique(unique_id)),], Value)
    )
  )

}

#%%%%%%%%%%%%%%%%%%%%
# DATA FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine number of categories in data ----
# Updated 06.07.2023
data_categories <- function(data)
{
  return(nvapply(as.data.frame(data), unique_length))
}

#' @noRd
# All-purpose symmetric checker ----
# Faster but less robust than `isSymmetric`
# Defaults to floating point tolerance
# Updated 01.07.2023
is_symmetric <- function(data, tolerance = 1e-06){

  # Get dimensions
  dimensions <- dim(data)

  # Check for whether rows equal columns
  if(dimensions[1] != dimensions[2]){ # Not a (square) matrix
    return(FALSE)
  }else{

    # Convert to matrix
    data_matrix <- as.matrix(data)

    # Check that all are equal
    return(all(abs(data_matrix - t(data_matrix)) < tolerance))

  }

}

#' @noRd
# Ensure data has dimension names ----
# Updated 24.07.2023
ensure_dimension_names <- function(data)
{

  # Determine object type
  object_type <- get_object_type(data)

  # Branch for vector vs. matrix and data frame
  if(object_type == "vector"){

    # Get length of vector
    vector_length <- length(data)

    # Check for names
    if(is.null(names(data))){
      names(data) <- paste0(
        "V", format_integer(
          seq_len(vector_length),
          digits(vector_length) - 1
        )
      )
    }

  }else if(object_type %in% c("matrix", "data.frame")){

    # Get dimensions
    dimensions <- dim(data)

    # Get dimension names
    dimension_names <- dimnames(data)

    # Check for column names
    if(is.null(dimension_names[[2]])){

      # Standardize names
      dimnames(data)[[2]] <- paste0(
        "V", format_integer(
          seq_len(dimensions[2]),
          digits(dimensions[2]) - 1
        )
      )

    }

    # Check for symmetric matrix
    if(is_symmetric(data)){

      # Check for row names
      if(is.null(dimnames(data)[[1]]) | any(dimension_names[[1]] != dimension_names[[2]])){

        # Assign column names to row names
        dimnames(data)[[1]] <- dimnames(data)[[2]]

      }

    }

  }

  # Return named data
  return(data)

}

#' @noRd
# Transfer names from data to output ----
# Usually for matrix and data frame
# Updated 25.06.2023
transfer_names <- function(data, output)
{

  # Obtain column names
  column_names <- dimnames(data)[[2]]

  # Check for variable names
  if(!is.null(column_names)){

    # Add names to rows and columns
    dimnames(output) <- list(column_names, column_names)

  }

  # Return named output
  return(output)

}

#' @noRd
# Function to check for usable variables ----
# Updated 24.07.2023
usable_data <- function(data, verbose)
{

  # Turn data into a data matrix
  data_matrix <- data.matrix(data)

  # All missing after coercions
  remove_columns <- lvapply(as.data.frame(is.na(data_matrix)), all)

  # Send warning if there are any removals
  if(any(remove_columns)){

    # Send warning
    warning(
      paste0(
        "Some variables could not to be coerced to numeric values. These variables have been removed from the analysis:\n",
        paste0("'", dimnames(data_matrix)[[2]][remove_columns], "'", collapse = ", "),
        "\n\nIf these variables were not intended to be removed, then try converting them to numeric values before inputting the data into the function"
      ), call. = FALSE
    )

    # Remove these variables from `data` and `data_matrix`
    keep_columns <- !remove_columns
    data <- data[,keep_columns]
    data_matrix <- data_matrix[,keep_columns]

  }

  # Determine whether each variable was coerced
  coercions <- lvapply(as.data.frame(data_matrix != data), any, na.rm = TRUE)

  # Send warning for any coercions
  if(any(coercions) & verbose){

    # Send warning
    warning(
      paste0(
        "Several variables were coerced to numeric values. These variables were changed to numeric values:\n",
        paste0("'", dimnames(data_matrix)[[2]][coercions], "'", collapse = ", ")
      ), call. = FALSE
    )

  }

  # Send back usable data
  return(data_matrix)

}

#' @noRd
# Whether usable data is needed ----
# Updated 07.09.2023
needs_usable <- function(ellipse)
{
  return(!"needs_usable" %in% names(ellipse) || ellipse$needs_usable)
}

#' @noRd
# Obtain data, sample size, correlation matrix ----
# Generic function to get the usual needed inputs
# Updated 29.11.2025
obtain_sample_correlations <- function(data, n, corr, na_data, verbose, ...)
{

  # Check if data is a correlation matrix
  if(is_symmetric(data)){

    # Check for sample size
    if(is.null(n)){
      stop(
        "A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.",
        call. = FALSE
      )
    }

    # If symmetric and sample size is provided, then
    # data to correlation matrix
    correlation_matrix <- data

  }else{

    # Check for usable data
    if(needs_usable(list(...))){
      data <- usable_data(data, verbose)
    }

    # Obtain sample size
    n <- dim(data)[1]

    # Check for automatic correlations
    if(corr == "auto"){
      correlation_matrix <- auto_correlate(
        data = data, corr = "pearson", na_data = na_data,
        verbose = verbose, ...
      )
    }else{
      correlation_matrix <- cor(
        data, use = swiftelse(
          na_data == "pairwise",
          "pairwise.complete.obs",
          "complete.obs"
        ), method = corr
      )
    }

  }

  # Return results
  return(
    list(
      data = data, n = n,
      correlation_matrix = correlation_matrix
    )
  )

}

# Compute thresholds ----
#' @noRd
# Updated 22.07.2023
obtain_thresholds <- function(categorical_variable)
{

  # Obtain cumulative sums from frequency table
  cumulative_sum <- cumsum(fast_table(categorical_variable))

  # Obtain cumulative length
  cumsum_length <- length(cumulative_sum)

  # Obtain thresholds
  return(
    qnorm(
      cumulative_sum[-cumsum_length] /
        cumulative_sum[cumsum_length]
    )
  )

}

#%%%%%%%%%%%%%%%%%%%%%%%%%
# QUIET CALLS & LOADS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

# Lots of R packages and functions print
# messages that add to the noise of {EGAnet}
# outputs; these functions make other
# packages' messages go away
#
# NOTE: In some cases, some package messages
# are important and should not be ignored
# (don't use these functions in those cases)

#' @noRd
# General function to silently obtain output ----
# Updated 01.07.2023
silent_call <- function(...)
{

  # Make call
  sink <- capture.output(
    result <- suppressWarnings(
      suppressMessages(
        ...
      )
    )
  )

  # Return result
  return(result)

}

#' @noRd
# General function to silently load package ----
# Updated 01.07.2023
silent_load <- function(...)
{
  return(suppressPackageStartupMessages(...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FUNCTIONS & ARGUMENTS FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Strip classes and attributes
# Updated 09.07.2023
remove_attributes <- function(object)
{

  # Remove all except special attributes (see `?structure`)
  attributes(object) <- attributes(object)[
    names(attributes(object)) %in% c(
      "dim", "dimnames", "names", "tsp", "levels"
    )
  ]

  # Return unclassed object
  return(unclass(object))

}

#' @noRd
# Set default argument (cleaner missing function) ----
# Updated 28.07.2023
set_default <- function(argument, default, FUN, several.ok = FALSE)
{

  # Check for type error
  typeof_error(argument, c(typeof(default), "closure"), "set_default")

  # Return argument if it's a function
  if(is.function(argument)){
    return(argument)
  }

  # Return default if NULL
  if(is.null(argument)){
    return(default)
  }

  # Return default if length > 1 and several are NOT okay
  if(length(argument) > 1 & !several.ok){
    return(default)
  }

  # Get argument name
  argument_name <- as.character(substitute(argument))

  # Get choices for the function
  if(!is.function(FUN)){ # choices are provided in 'FUN'
    choices <- FUN
  }else{ # get formal arguments from the function itself

    # Get formal argument choices from default call
    choices <- formals(FUN)[[argument_name]]

    # Check if choices is a call
    if(is.call(choices)){

      # Get choices
      choices <- as.character(choices)

      # Remove "c" from call
      choices <- choices[-which(choices == "c")]

    }

  }

  # Make argument lowercase (if necessary)
  argument <- tolower(argument); choices <- tolower(choices)

  # Get arguments not in choices
  if(any(!argument %in% choices)){ # Repeated calls OK since returning error

    # Get arguments not in choices
    not_in_choices <- !argument %in% choices

    # Get number not in choices
    number_not <- sum(not_in_choices)

    # Send error for bad argument
    if(number_not > 0){

      # Only one not in choices
      if(number_not == 1){
        argument_text <- paste0(" = \"", argument[not_in_choices], "\"")
      }else{ # Multiple not in choices
        argument_text <- paste0(
          " = c(",
          paste0("\"", argument[not_in_choices], "\"", collapse = ", "),
          ")"
        )
      }

      # Throw error
      stop(
        paste0(
          paste0("Invalid argument: ", argument_name, argument_text),
          "\n\nPlease choose from: ",
          paste0("\"", choices, "\"", collapse = ", ")
        ), call. = FALSE
      )

    }

  }

  # Return argument
  return(argument)

}

#' @noRd
# Overwrite Arguments ----
# Updated 24.06.2023
overwrite_arguments <- function(main, ARGS)
{

  # Obtain replacement arguments
  target_args <- ARGS[names(ARGS) %in% names(main)]

  # Replace the arguments with the same name
  main[names(target_args)] <- target_args

  # Return main arguments
  return(main)

}

#' @noRd
# Function to obtain arguments ----
# Updated 06.07.2023
obtain_arguments <- function(FUN, FUN.args)
{

  # Overwrite arguments
  FUN.formals <- overwrite_arguments(formals(FUN), FUN.args)

  # Remove ellipses
  FUN.formals <- FUN.formals[names(FUN.formals) != "..."]

  # Remove call arguments (assume they are supplied elsewhere)
  return(FUN.formals[!lvapply(FUN.formals, is, "call")])

}

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Experimental warning ----
# Updated 03.08.2023
experimental <- function(function_name)
{
  warning(
    paste0(
      "This implementation of `", function_name,
      "` is ", styletext("experimental", "italics"), ". \n\nThe underlying function ",
      "and/or output may change until the results have been appropriately vetted and validated."
    ), call. = FALSE
  )
}

#' @noRd
# Error for class ----
# Updated 04.09.2023
class_error <- function(input, expected_class, function_name){

  # Check for object types
  if(!is(input, expected_class)){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is an object with ", paste0("'", class(input), "'", collapse = ", "), " class(es)",
        ". Input is expected to be ", paste0("'", expected_class, "'", collapse = ", "), " class(es)",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#class-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for object type ----
# Updated 04.09.2023
object_error <- function(input, expected_type, function_name){

  # Get input type
  input_type <- get_object_type(input)

  # Check for object types
  if(!input_type %in% expected_type){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is a ", paste0("'", input_type, "'", collapse = ", "), " object",
        ". Input is expected to be ", paste0("'", expected_type, "'", collapse = ", "), " object",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#object-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `typeof` ----
# Updated 04.09.2023
typeof_error <- function(input, expected_value, function_name){

  # Switch out "closure" with "function"
  if("closure" %in% expected_value){
    expected_value[expected_value == "closure"] <- "function"
  }

  # Get type of input
  typeof_input <- typeof(input)

  # Convert "closure" to "function"
  if("closure" %in% typeof_input){
    typeof_input[typeof_input == "closure"] <- "function"
  }

  # Convert "integer" and "double" to "numeric"
  ## Input
  typeof_input <- swiftelse(
    typeof_input %in% c("integer", "double"),
    "numeric", typeof_input
  )
  ## Expected value
  expected_value <- swiftelse(
    any(expected_value %in% c("integer", "double")),
    "numeric", expected_value
  )

  # Check for value
  if(!typeof_input %in% expected_value){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is ", paste("'", typeof_input, "'", sep = "", collapse = ", "), " type",
        ". Input is expected to be ",
        paste0("'", expected_value, "'", collapse = " or "), " type",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#typeof-error"
        # can use "or" because `typeof` won't ever be more than two
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `length` ----
# Updated 04.09.2023
length_error <- function(input, expected_lengths, function_name){

  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste0("'", expected_lengths, "'", collapse = " or "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#length-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `range` ----
# Updated 04.09.2023
range_error <- function(input, expected_ranges, function_name){

  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges, na.rm = TRUE)
  expected_minimum <- min(expected_ranges, na.rm = TRUE)

  # Obtain maximum and minimum values
  actual_maximum <- round(max(input, na.rm = TRUE), 3)
  actual_minimum <- round(min(input, na.rm = TRUE), 3)

  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Maximum value of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ),
      call = function_name
    )
  }

  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Minimum value of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ),
      call = function_name
    )
  }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COVARIANCE CONVERSION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Zero-order correlations to partial correlations ----
# Updated 29.06.2023
cor2pcor <- function(correlation_matrix)
{

  # Convert to inverse correlations to partial correlations
  partial_correlations <- -cov2cor(solve(correlation_matrix))

  # Set diagonal to zero
  diag(partial_correlations) <- 0

  # Return
  return(partial_correlations)

}

#' @noRd
# Inverse covariances to partial correlations ----
# Updated 07.03.2026
inv2pcor <- function(precision)
{

  # Convert to inverse correlations to partial correlations
  partial_correlations <- -cov2cor(precision)

  # Set diagonal to zero
  diag(partial_correlations) <- 0

  # Return
  return(partial_correlations)

}

#' @noRd
# Partial correlations to zero-order correlations ----
# Updated 29.06.2023
pcor2cor <- function(partial_correlations)
{

  # Set diagonal to negative 1
  diag(partial_correlations) <- -1

  # Return partial correlations as zero-order correlations
  return(cov2cor(solve(-partial_correlations)))

}

#' @noRd
# Partial correlations to inverse covariances ----
# Updated 25.09.2024
pcor2inv <- function(partial_correlations)
{

  # Set diagonal to negative 1
  diag(partial_correlations) <- -1

  # Return inverse covariance matrix
  return(solve(-partial_correlations))

}

#%%%%%%%%%%%%%%%%%%%%%%%%
## NETWORK FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Create sparse network ----
# Updated 08.11.2023
sparse_network <- function(network)
{

  # Get number of nodes
  nodes <- dim(network)[2]

  # Get node sequence
  node_sequence <- seq_len(nodes)

  # Create data frame
  sparse_df <- data.frame(
    row = rep(node_sequence, each = nodes),
    col = rep(node_sequence, times = nodes),
    weight = as.vector(network)
  )

  # Return lower triangle
  return(sparse_df[sparse_df$row < sparse_df$col,])

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MATH & STATS FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# There are some redundancies in the math functions because
# it's often faster to call some functions directly
# rather than nesting non-base R functions in functions

#' @noRd
# Positive definite matrix ----
# Logical for whether a matrix is positive definite
# Updated 05.11.2024
is_positive_definite <- function(data)
{
  return(all(eigen(data, symmetric = TRUE, only.values = TRUE)$values > .Machine$double.eps))
}

#' @noRd
# Wrapper for `eigen` ----
# Extracts eigenvalues only
# Updated 29.06.2023
matrix_eigenvalues <- function(data)
{
  return(eigen(data, symmetric = TRUE, only.values = TRUE)$values)
}

#' @noRd
# Trace of matrix ----
# Updated 25.06.2023
trace <- function(object)
{
  return(sum(diag(object), na.rm = TRUE))
}

#' @noRd
# r to Fisher's z ----
# Updated 27.03.2025
r2z <- function(r)
{
  return(0.5 * log((1 + r) / (1 - r)))
}

#' @noRd
# Fisher's z to r ----
# Updated 27.03.2025
z2r <- function(z)
{
  return((exp(2 * z) - 1) / (1 + exp(2 * z)))
}

#' @noRd
# Partial correlation p-values ----
# Updated 27.03.2025
p_partial <- function(z, variables, sample_size)
{

  # Get test statistic
  test <- z / sqrt(1 / (sample_size - variables - 5))

  return(pnorm(test, lower.tail = FALSE) * 2)

}

#%%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# OS and System Check ----
# Updated 28.11.2023
system.check <- function()
{

  # Get OS usage
  OS <- unname(tolower(Sys.info()["sysname"]))

  # Get RStudio usage
  RSTUDIO <- swiftelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)

  # Return list
  return(
    list(
      OS = OS,
      R = paste0(R.version$major, ".", R.version$minor),
      RSTUDIO = RSTUDIO,
      TEXT = swiftelse(!RSTUDIO & OS != "linux", FALSE, TRUE)
    )
  )

}

#' @noRd
# Colorize text ----
# Updated 08.09.2020
colortext <- function(text, number = NULL, defaults = NULL)
{
  # Check system
  sys.check <- system.check()

  if(sys.check$TEXT)
  {
    # Defaults for number (white text)
    if(is.null(number) || number < 0 || number > 231)
    {number <- 15}

    # Check for default color
    if(!is.null(defaults))
    {
      # Adjust highlight color based on background color
      if(defaults == "highlight")
      {
        number <- 208
      }else{

        number <- switch(defaults,
                         message = 204,
                         red = 9,
                         orange = 208,
                         yellow = 11,
                         "light green" = 10,
                         green = 34,
                         cyan = 14,
                         blue = 12,
                         magenta = 13,
                         pink = 211,
        )

      }

    }

    return(paste("\033[38;5;", number, "m", text, "\033[0m", sep = ""))

  }else{return(text)}
}

#' @noRd
# Style text ----
# Updated 26.07.2023
styletext <- function(
    text, defaults = c(
      "bold", "italics", "highlight",
      "underline", "strikethrough"
    )
)
{

  if(system.check()$TEXT)
  {

    if(missing(defaults)){
      number <- 0
    }else{

      number <- switch(
        defaults,
        bold = 1, italics = 3,
        underline = 4, highlight = 7,
        strikethrough = 9
      )

    }

    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))

  }else{
    return(text)
  }

}

#' @noRd
# Symbols ----
# Updated 26.07.2024
textsymbol <- function(
    symbol = c(
      "alpha", "beta", "chi", "delta",
      "eta", "gamma", "lambda", "omega",
      "phi", "pi", "rho", "sigma", "tau",
      "theta", "square root", "infinity",
      "check mark", "x", "bullet"
    )
)
{
  # Return code
  return(
    switch(
      symbol,
      alpha = "\u03B1", beta = "\u03B2", chi = "\u03C7",
      delta = "\u03B4", eta = "\u03B7", gamma = "\u03B3",
      lambda = "\u03BB,", omega = "\u03C9", phi = "\u03C6",
      pi = "\u03C0", rho = "\u03C1", sigma = "\u03C3", tau = "\u03C4",
      theta = "\u03B8", "square root" = "\u221A", infinity = "\u221E",
      "check mark" = "\u2713", x = "\u2717", bullet = "\u2022"
    )
  )

}

#' @noRd
# Title Case ----
# Not truly title case -- just capitializes
# the first letter of each word
# Uses `tools::toTitleCase` for actual title case
# Updated 25.06.2023
totitle <- function(string)
{

  # Split by spaces
  words <- unlist(strsplit(string, split = " "))

  # Set first letters to uppercase
  titleCased <- cvapply(words, function(x){

    # Stitch together letters
    return(
      paste0(
        toupper(substr(x, 1, 1)),
        tolower(substr(x, 2, nchar(x)))
      )
    )

  })

  # Paste words back together
  return(paste(titleCased, collapse = " "))

}

#' @noRd
# General function to check for packages ----
# Updated 13.07.2023
check_package <- function(packages)
{

  # Performs what original `installed.packages()` does
  # but without additional fluff
  installed <- packages %in% ulapply(.libPaths(), list.files)

  # Determine which packages are not installed
  not_installed <- packages[!installed]

  # Print error with missing packages
  if(length(not_installed) != 0){

    # Organize packages error output
    if(length(not_installed) > 1){

      # Get missing packages
      missing_packages <- paste0("{", packages , "}", collapse = ", ")
      packages <- paste0("\"", packages, "\"", collapse = ", ")

      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages,
          " are not installed but are required for this function. ",
          "Please run \n\n",
          "`install.packages(c(", packages, "))`",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        ), call. = FALSE
      )

    }else{

      # Get missing packages
      missing_packages <- paste0("{", packages, "}")
      packages <- paste0("\"", packages, "\"")

      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages,
          " is not installed but is required for this function. ",
          "Please run \n\n",
          "`install.packages(", packages, ")`",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        ), call. = FALSE
      )

    }

  }

}
