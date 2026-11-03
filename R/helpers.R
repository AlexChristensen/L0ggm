#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### {L0ggm} Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

# These helpers are internal functions that are repeatedly
# used across the package. They are implemented here once
# so that changes only need to be made here and not in
# multiple places across the package.

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

#' @noRd
# Determine number of digits in a number ----
# Updated 24.07.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(abs(number))) + 1)

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
# Positive definite matrix ----
# Logical for whether a matrix is positive definite
# Updated 05.11.2024
is_positive_definite <- function(data)
{
  return(all(eigen(data, symmetric = TRUE, only.values = TRUE)$values > .Machine$double.eps))
}

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

