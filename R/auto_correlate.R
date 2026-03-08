#' @title Automatic correlations
#'
#' @description Automatically computes the proper correlations between
#' continuous and categorical variables. \code{NA} values are not treated
#' as categories
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param corr Character (length = 1).
#' The standard correlation method to be used.
#' Defaults to \code{"pearson"}.
#' Using \code{"pearson"} will compute polychoric, tetrachoric, polyserial,
#' and biserial correlations for categorical and categorical/continuous correlations
#' by default. To obtain \code{"pearson"} correlations regardless, use \code{\link{cor}}.
#' Other options of \code{"kendall"} and \code{"spearman"} are provided for
#' completeness and use \code{\link{cor}}
#'
#' @param ordinal_categories Numeric (length = 1).
#' \emph{Up to} the number of categories \emph{before} a variable is considered continuous.
#' Defaults to \code{7} categories before \code{8} is considered continuous
#'
#' @param forcePD Boolean (length = 1).
#' Whether positive definite matrix should be enforced.
#' Defaults to \code{TRUE}
#'
#' @param na_data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"pairwise"} --- Computes correlation for all available
#' cases between two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete
#' cases in the dataset
#'
#' }
#'
#' @param empty_method Character (length = 1).
#' Method for empty cell correction in \code{\link[L0ggm]{polychoric_matrix}}.
#' Defaults to \code{"none"}
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"none"} --- Adds no value (\code{empty_value = "none"})
#' to the empirical joint frequency table between two variables
#'
#' \item \code{"zero"} --- Adds \code{empty_value} to the cells with
#' zero in the joint frequency table between two variables
#'
#' \item \code{"all"} --- Adds \code{empty_value} to all
#' in the joint frequency table between two variables
#'
#' }
#'
#' @param empty_value Character (length = 1).
#' Value to add to the joint frequency table cells in \code{\link[L0ggm]{polychoric_matrix}}.
#' Defaults to \code{"none"}.
#' Accepts numeric values between 0 and 1 or specific methods:
#'
#' \itemize{
#'
#' \item \code{"none"} --- Adds no value (\code{0}) to the empirical
#' joint  frequency table between two variables
#'
#' \item \code{"point_five"} --- Adds \code{0.5} to the cells
#' defined by \code{empty_method}
#'
#' \item \code{"one_over"} --- Adds \code{1 / n} where \code{n} equals the
#' number of cells based on \code{empty_method}. For
#' \code{empty_method = "zero"}, \code{n} equals the number of zero cells
#'
#' }
#'
#' @param forceReturn Boolean (length = 1).
#' Whether correlation matrix should be forced to return.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to receive the correlation matrix "as is"
#'
#' @param verbose Boolean (length = 1).
#' Whether messages should be printed.
#' Defaults to \code{FALSE}
#'
#' @param ...
#' Not actually used but makes it easier for general functionality
#' in the package
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Simulate data
#' wmt <- wmt2[,7:24]
#'
#' # Obtain correlations
#' wmt_corr <- auto_correlate(wmt)
#'
#' @export
#'
# Automatic correlations ----
# Updated 07.03.2026
auto_correlate <- function(
    data, # Matrix or data frame
    corr = c("kendall", "pearson", "spearman"), # allow changes to standard correlations
    ordinal_categories = 7, # consider ordinal up to 7 categories
    forcePD = TRUE, # ensure result is positive definite
    na_data = c("pairwise", "listwise"), # use available or complete values
    empty_method = c("none", "zero", "all"), # zero frequencies in categorical correlations
    empty_value = c("none", "point_five", "one_over"), # value to use in zero cells
    forceReturn = FALSE, # return even if bad correlation matrix
    verbose = FALSE, # don't print messages
    ... # not actually used
)
{

  # Argument errors (return data in case of tibble)
  data <- auto_correlate_errors(data, ordinal_categories, forcePD, forceReturn, verbose, ...)

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "pearson", auto_correlate)
  na_data <- set_default(na_data, "pairwise", auto_correlate)
  empty_method <- set_default(empty_method, "none", auto_correlate)
  empty_value <- set_default(empty_value, "none", auto_correlate)

  # Ensure matrix
  data <- as.matrix(data)

  # Ensure variable names
  data <- ensure_dimension_names(data)

  # Get variable names
  variable_names <- dimnames(data)[[2]]

  # Set dimensions
  dimensions <- dim(data)

  # Assumes some preprocessing has taken place
  # to only select for appropriate variables

  # Determine whether categorical correlations are necessary
  if(corr != "pearson"){

    # Obtain correlation matrix
    correlation_matrix <- cor(
      x = data, use = swiftelse(
        na_data == "pairwise",
        "pairwise.complete.obs",
        "complete.obs"
      ), method = corr
    )

  }else{ # Proceed with determination of categorical correlations

    # Obtain the number of categories for each variables
    categories <- data_categories(data)

    # Determine categorical variables
    categorical_variables <- categories <= ordinal_categories

    # Determine number of categorical variables
    categorical_number <- sum(categorical_variables)

    # Determine whether there are any categorical variables
    if(categorical_number != 0){

      # Determine continuous variables
      continuous_variables <- !categorical_variables

      # Set up correlation matrix
      correlation_matrix <- matrix(
        nrow = dimensions[2], ncol = dimensions[2],
        dimnames = list(variable_names, variable_names)
      )

      # Determine whether there are more than one categorical variables
      if(categorical_number > 1){

        # Compute categorical correlations (only correlation for categorical data)
        # Add correlations to correlation matrix
        correlation_matrix[
          categorical_variables, categorical_variables # ensure proper indexing
        ] <- polychoric_matrix(
          data = data[,categorical_variables], na_data = na_data,
          empty_method = empty_method, empty_value = empty_value,
          needs_usable = FALSE # skip usable data check
        )

      }

      # Determine whether there are more than one continuous variables
      if(sum(continuous_variables) > 1){

        # Compute continuous correlations (only correlation for continuous data)
        # Add correlations to correlation matrix
        correlation_matrix[
          continuous_variables, continuous_variables # ensure proper indexing
        ] <- cor(
          x = data[,continuous_variables],
          use = swiftelse(
            na_data == "pairwise",
            "pairwise.complete.obs",
            "complete.obs"
          ), method = corr
        )

      }

      # Determine whether there are mixed variables
      if(categorical_number != dimensions[2]){ # Check for mixed variables

        # Obtain continuous data (keep as matrix)
        continuous_data <- data[,continuous_variables, drop = FALSE]

        # Loop over categorical indices
        for(i in which(categorical_variables)){

          # Fill matrix
          correlation_matrix[continuous_variables, i] <-
            correlation_matrix[i, continuous_variables] <-
            polyserial_vector( # computes polyserial vector
              categorical_variable = data[,i], # drops to vector
              continuous_variables = continuous_data,
              na_data = na_data
            )

        }

      }

    }else{

      # Compute Pearson's correlations
      correlation_matrix <- cor(
        x = data, use = swiftelse(
          na_data == "pairwise",
          "pairwise.complete.obs",
          "complete.obs"
        ), method = corr
      )

    }

  }

  # Set diagonal as one
  diag(correlation_matrix) <- 1

  # Try to test for positive definite
  PD <- try(is_positive_definite(correlation_matrix), silent = TRUE)
  PD_error <- is(PD, "try-error")

  # Check for error
  if(PD_error && forceReturn){
    return(correlation_matrix)
  }

  # Determine whether matrix is positive definite
  if((PD_error || !PD) && forcePD){

    # Send warning to user (if `verbose`)
    if(verbose){
      warning(
        "Correlation matrix is not positive definite. Finding nearest positive definite matrix using `Matrix::nearPD`",
        call. = FALSE
      )
    }

    # Regardless, make matrix positive definite
    correlation_matrix <- as.matrix(
      Matrix::nearPD(
        correlation_matrix, corr = TRUE,
        ensureSymmetry = TRUE, keepDiag = TRUE
      )$mat
    )

  }

  # Return correlation matrix
  return(correlation_matrix)

}

# Bug checking ----
# ## Different categories
# set.seed(1234)
# data = latentFactoR::simulate_factors(
#   factors = 4, variables = 4,
#   loadings = 0.55, cross_loadings = 0.05,
#   correlations = 0.30, sample_size = 1000,
#   variable_categories = c(
#     rep(2, 4), rep(5, 4),
#     rep(7, 4), rep(Inf, 4)
#   )
# )$data
# ordinal_categories = 7
# corr = "pearson"; forcePD = TRUE
# na_data = "pairwise"; empty_method = "none"
# empty_value = "none"; verbose = FALSE
#
# # Compare against {qgraph}'s `cor_auto`
# qgraph_correlations <- qgraph::cor_auto(data)
# EGAnet_correlations <- auto_correlate(data)
#
# # Difference
# max(abs(EGAnet_correlations - qgraph_correlations))
# # Biggest difference is between polyserial (7 categories with continuous)
#
# ## Add missing data
# data[sample(1:length(data), 1000)] <- NA
# # Compare against {qgraph}'s `cor_auto`
# qgraph_correlations <- qgraph::cor_auto(
#   data,
#   ordinalLevelMax = 8
#   # Needs to have 8 levels to account for missing data!!
# )
# EGAnet_correlations <- auto_correlate(data)
#
# # Difference
# max(abs(EGAnet_correlations - qgraph_correlations))
# # Biggest difference is between polyserial (7 categories with continuous)
#
## Zero cell counts
# data <- cbind(
#   c(0, 1, 2, 3, 4, 1, 2, 3, 1),
#   c(0, 1, 2, 2, 1, 2, 2, 4, 2)
# ); ordinal_categories = 7;
# corr = "pearson"; forcePD = TRUE;
# na_data = "pairwise"; empty_method = "none";
# empty_value = "none"; verbose = FALSE;
#
# The above bug checks for categorical data
# have been verified to match the output
# of Turbofuns::PolychoricRM to a maximum
# difference less than or equal to 1.0e-06
# (or one step beyond floating point)

#' @noRd
# Errors ----
# Updated 29.11.2025
auto_correlate_errors <- function(data, ordinal_categories, forcePD, forceReturn, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "auto_correlate")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'ordinal_categories' errors
  length_error(ordinal_categories, 1, "auto_correlate")
  typeof_error(ordinal_categories, "numeric", "auto_correlate")
  range_error(ordinal_categories, c(2, 11), "auto_correlate")

  # 'forcePD' errors
  length_error(forcePD, 1, "auto_correlate")
  typeof_error(forcePD, "logical", "auto_correlate")

  # 'verbose' errors
  length_error(verbose, 1, "auto_correlate")
  typeof_error(verbose, "logical", "auto_correlate")

  # 'forceReturn' errors
  length_error(forceReturn, 1, "auto_correlate")
  typeof_error(forceReturn, "logical", "auto_correlate")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return usable data (in case of tibble)
  return(data)

}

# Compute polyserial correlation ----
#' For a single categorical variable, compute correlations
#' with \emph{n} continuous variables
#'
#' Uses two-step approximation from {polycor}'s \code{polyserial}
#'
#' @noRd
# Updated 16.03.2024
polyserial_vector <- function(
    categorical_variable, continuous_variables,
    na_data = c("pairwise", "listwise")
)
{

  # Determine cases based on `na_data` argument
  if(na_data == "pairwise"){

    # Pairwise cases
    categorical_cases <- colSums(
      !is.na(categorical_variable) *
      !is.na(continuous_variables)
    )

  }else if(na_data == "listwise"){

    # Complete cases
    complete_cases <- complete.cases(cbind(categorical_variable, continuous_variables))

    # Repeat cases for number of continuous variables
    categorical_cases <- rep(sum(complete_cases), dim(continuous_variables)[2])

  }

  # Compute correlation
  return(
    sqrt((categorical_cases - 1) / categorical_cases) *
      sd(categorical_variable, na.rm = TRUE) *
      # Correlations with scaled continuous variables
      cor(
        categorical_variable, scale(continuous_variables),
        use = swiftelse(
          na_data == "pairwise",
          "pairwise.complete.obs",
          "complete.obs"
        )
      ) /
      # Compute sum of thresholds
      sum(dnorm(obtain_thresholds(categorical_variable)), na.rm = TRUE)
  )

}
