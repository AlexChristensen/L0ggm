#' @title Confusion Matrix Metrics for Edge Comparison and Recovery
#'
#' @description Computes many commonly used confusion matrix metrics
#'
#' @param base Matrix or data frame.
#' Network that will be treated as the "ground truth" such that
#' a false positive represents an edge that is present in \code{comparison}
#' but not in this network
#'
#' @param comparison Matrix or data frame.
#' Network that will be treated as the estimator such that a
#' false positive represents an edge that is present in this network
#' but not in \code{base}
#'
#' @param metric Character vector.
#' Defaults to \code{"all"} metrics.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"all"} --- All available metrics (default)
#'
#' \item \code{"sen"} --- Sensitivity (True Positive Rate):
#' \deqn{\frac{TP}{TP + FN}}
#'
#' \item \code{"spec"} --- Specificity (True Negative Rate):
#' \deqn{\frac{TN}{TN + FP}}
#'
#' \item \code{"ppv"} --- Positive Predictive Value (Precision):
#' \deqn{\frac{TP}{TP + FP}}
#'
#' \item \code{"npv"} --- Negative Predictive Value:
#' \deqn{\frac{TN}{TN + FN}}
#'
#' \item \code{"fdr"} --- False Discovery Rate:
#' \deqn{1 - PPV = \frac{FP}{TP + FP}}
#'
#' \item \code{"fom"} --- False Omission Rate:
#' \deqn{1 - NPV = \frac{FN}{TN + FN}}
#'
#' \item \code{"ba"} --- Balanced Accuracy:
#' \deqn{\frac{Sensitivity + Specificity}{2}}
#'
#' \item \code{"f1"} --- F1 Score (harmonic mean of PPV and Sensitivity):
#' \deqn{\frac{2TP}{2TP + FP + FN}}
#'
#' \item \code{"csi"} --- Critical Success Index (Jaccard / Threat Score):
#' \deqn{\frac{TP}{TP + FP + FN}}
#'
#' \item \code{"mcc"} --- Matthews Correlation Coefficient:
#' \deqn{\frac{TP \times TN - FP \times FN}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}}
#'
#' }
#'
#' @param full.names Boolean (length = 1).
#' Whether full or abbreviated names should be used.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for full names
#'
#' @return A named numeric vector of the requested confusion matrix metrics.
#' Values are generally in \eqn{[0, 1]}, with the exception of \code{"mcc"}
#' (Matthews Correlation Coefficient), which ranges from \eqn{-1} to \eqn{1}
#' where 1 indicates perfect agreement, 0 indicates chance-level performance,
#' and -1 indicates perfect disagreement. The vector contains only the elements
#' specified by \code{metric} (all ten metrics when \code{metric = "all"},
#' which is the default). Names are abbreviated (e.g., \code{"sen"}) when
#' \code{full.names = FALSE} (default), or expanded
#' (e.g., \code{"Sensitivity"}) when \code{full.names = TRUE}. Any metric
#' whose denominator is zero (e.g., sensitivity when no edges exist in
#' \code{base}) is returned as \code{NA}.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Set split
#' split <- sample(
#'   1:nrow(basic_smallworld),
#'   round(nrow(basic_smallworld) / 2)
#' )
#'
#' # Estimate networks
#' split1 <- network_estimation(basic_smallworld[split,])
#' split2 <- network_estimation(basic_smallworld[-split,])
#'
#' # Estimate metrics
#' edge_confusion(split1, split2)
#'
#' @export
#'
# Compute confusion matrix metrics ----
# Updated 26.02.2026
edge_confusion <- function(base, comparison, metric = c(
  "all", "sen", "spec",
  "ppv", "npv", "fdr", "fom",
  "ba", "f1", "csi", "mcc"
), full.names = FALSE)
{

  # Set names
  nice_names <- c(
    "sen" = "Sensitivity", "spec" = "Specificity", "ppv" = "Positive Predictive Value",
    "npv" = "Negative Predictive Value", "fdr" = "False Discovery Rate", "fom" = "False Omission Rate",
    "ba" = "Balanced Accuracy", "f1" = "F1 Score", "csi" = "Critical Success Index (Jaccard)",
    "mcc" = "Matthew's Correlation Coefficient"
  )

  # Set default
  metric <- set_default(metric, "all", edge_confusion, several.ok = TRUE)

  # Check for all
  if("all" %in% metric){
    metric <- c("sen", "spec", "ppv", "npv", "fdr", "fom", "ba", "f1", "csi", "mcc")
  }

  # Obtain lower triangle indices and networks
  lower_triangle <- lower.tri(base)
  base_lower <- base[lower_triangle]
  comp_lower <- comparison[lower_triangle]

  # Set up indices
  base_present <- base_lower != 0
  base_absent <- !base_present
  comp_present <- comp_lower != 0
  comp_absent <- !comp_present

  # Set up basic confusion matrix metrics
  TP <- sum(base_present & comp_present)
  FP <- sum(base_absent & comp_present)
  TN <- sum(base_absent & comp_absent)
  FN <- sum(base_present & comp_absent)

  # Totals
  P <- TP + FN
  N <- TN + FP
  TP2 <- 2 * TP

  # Check for zero denominators
  TPFP <- TP + FP
  TNFN <- TN + FN
  FPFN <- FP + FN
  TP2FPFN <- TP2 + FPFN
  TPFPFN <- TP + FPFN

  # Compute metrics
  sen <- swiftelse(P == 0, NA, TP / P)
  spec <- swiftelse(N == 0, NA, TN / N)
  ppv <- swiftelse(TPFP == 0, NA, TP / TPFP)
  npv <- swiftelse(TNFN == 0, NA, TN / TNFN)

  # Create metric vector
  metric_vector <- c(
    "sen" = sen, "spec" = spec, "ppv" = ppv, "npv" = npv,
    "fdr" = 1 - ppv, "fom" = 1 - npv, "ba" = (sen + spec) / 2,
    "f1" = swiftelse(TP2FPFN == 0, NA, TP2 / TP2FPFN),
    "csi" = swiftelse(TPFPFN == 0, NA, TP / TPFPFN),
    "mcc" = silent_call(cor(base_present, comp_present)) # equivalent
  )

  # Return metrics
  return(
    structure(
      metric_vector[metric],
      names = swiftelse(full.names, nice_names[metric], metric)
    )
  )

}
