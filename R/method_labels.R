#' Manuscript method labels
#'
#' Converts short calibration keys or legacy result tags to manuscript labels.
#' Unknown labels are retained. With no argument, returns the seven labels in
#' manuscript order. Calibration API keys Li and P1-P4 remain unchanged.
#' @param x Character vector of keys or labels, or NULL for all labels.
#' @return A character vector of display labels.
#' @export
method_labels <- function(x = NULL) {
  labels <- c("Internal-only", "Adaptive lasso", "Precision-weighted L1",
              "Integrated-gate", "Information-adaptive MCP",
              "LR-weighted L1", "Naive pooled")
  if (is.null(x)) return(labels)
  keys <- c("IO", "Li", "P1", "P2", "P3", "P4", "NP",
            "InternalOnly", "LiAdaptiveLasso", "P1_SEScaledL1", "P2_GatedL1",
            "P3_SEScaledMCP", "P4_LRWeightedL1", "NaivePooled", "AL")
  values <- c(labels, labels, "Adaptive lasso")
  x <- as.character(x)
  i <- match(x, keys)
  x[!is.na(i)] <- values[i[!is.na(i)]]
  x
}
