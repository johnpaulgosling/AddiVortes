#' @title K-Nearest Neighbours Index Search
#' @description Find the indices of k nearest neighbours for each query point.
#'   This is a replacement for FNN::knnx.index to remove external dependencies.
#'
#' @param data A numeric matrix of reference points (tessellation centres) where each
#'   row is a point.
#' @param query A numeric matrix of query points where each row is a point to find
#'   neighbors for.
#' @param k An integer specifying the number of nearest neighbours to find.
#' @param dim The index susbet of covariates considered.
#' @param metric Either "Euclidean" or "Spherical".
#'
#' @return A matrix of integers where each row corresponds to a query point and
#'   contains the row indices of the k nearest neighbours in the data matrix.
#'   For k=1, returns a vector of integers.
#'
#' @keywords internal
#' @noRd
knnx_index <- function(data, query, k = 1, dim, metric) {
  # Input validation
  if (!is.matrix(data)) data <- as.matrix(data)
  if (!is.matrix(query)) query <- as.matrix(query)
  if (ncol(data) != ncol(query)) {
    stop("Number of columns in data and query must match")
  }
  if (k <= 0 || k > nrow(data)) {
    stop("k must be positive and not greater than number of reference points")
  }
  
  # Call C++ implementation
  result <- .Call("knnx_index_cpp", data, query, as.integer(k), dim, metric)
  
  # For k=1, return as vector (to match FNN::knnx.index behavior)
  if (k == 1) {
    result <- as.vector(result)
  }
  
  return(result)
}