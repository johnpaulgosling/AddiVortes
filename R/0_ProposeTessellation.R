#' @title Propose a New Tessellation via a Stochastic Modification
#' @description Generates a proposal for a new tessellation by randomly selecting
#'   and applying one of six modification types to an existing tessellation:
#'   adding/removing/swapping a dimension or adding/removing/changing a centre.
#'
#' @details The choice of modification is determined stochastically based on a
#'   uniformly drawn random number. This function is typically used to generate a
#'   candidate state within a single iteration of a Metropolis-Hastings MCMC algorithm.
#'
#' @param tess_j The current tessellation matrix, where each row represents a
#'   centre and each column corresponds to a dimension.
#' @param dim_j An integer vector specifying the column indices of the covariates
#'   that define the dimensions of the current tessellation.
#' @param var The variance parameter for the normal distribution (`rnorm`) used to
#'   generate new coordinate values for centres and dimensions.
#' @param covariateIndices This parameter is accepted by the function but is
#'   **not used** in the current implementation.
#' @param numCovariates An integer giving the total number of covariates
#'   available for selection in the full dataset.
#'
#' @return A list containing three elements:
#'   \enumerate{
#'     \item A numeric matrix representing the proposed new tessellation (`tess_j_star`).
#'     \item An integer vector of the dimensions for the new tessellation (`dim_j_star`).
#'     \item A character string indicating the type of modification applied (e.g., "AD", "RC", "Swap").
#'   }
#'
#' @keywords internal
#' @noRd
proposeTessellation <- function(tess_j, dim_j, var, covariateIndices, numCovariates) {
  p <- runif(1, 0, 1)
  d_j_length <- length(dim_j)
  tess_j_rows <- nrow(tess_j)
  
  dim_j_star <- dim_j
  tess_j_star <- tess_j
  Modification <- "Change" # Default modification type
  
  if (p < 0.2 && d_j_length != numCovariates ||
      d_j_length == 1 && p < 0.4) {
    # Add Dimension: Use a more efficient sample-and-check loop
    Modification <- "AD"
    repeat {
      new_dim <- sample.int(numCovariates, 1)
      if (!any(dim_j == new_dim)) break
    }
    dim_j_star <- c(dim_j, new_dim)
    tess_j_star <- cbind(tess_j, rnorm(tess_j_rows, 0, var))
    
  } else if (p < 0.4 && d_j_length > 1) {
    # Remove Dimension
    Modification <- "RD"
    removed_dim_idx <- sample.int(d_j_length, 1)
    dim_j_star <- dim_j[-removed_dim_idx]
    tess_j_star <- tess_j[, -removed_dim_idx, drop = FALSE]
    
  } else if (p < 0.6 || (p < 0.8 && tess_j_rows == 1)) {
    # Add Centre
    Modification <- "AC"
    tess_j_star <- rbind(tess_j, rnorm(d_j_length, 0, var))
    
  } else if (p < 0.8 && tess_j_rows > 1) {
    # Remove Centre
    Modification <- "RC"
    center_removed_idx <- sample.int(tess_j_rows, 1)
    tess_j_star <- tess_j[-center_removed_idx, , drop = FALSE]
    
  } else if (p < 0.9 || d_j_length == numCovariates) {
    # Change Centre (This is the default modification)
    centre_to_change_idx <- sample.int(tess_j_rows, 1)
    tess_j_star[centre_to_change_idx, ] <- rnorm(d_j_length, 0, var)
    
  } else {
    # Swap Dimension
    Modification <- "Swap"
    # Select dimension to replace
    dim_to_change_idx <- sample.int(d_j_length, 1)
    # Select new dimension to replace it with
    repeat {
      new_dim <- sample.int(numCovariates, 1)
      if (!any(dim_j == new_dim)) break
    }
    dim_j_star[dim_to_change_idx] <- new_dim
    tess_j_star[, dim_to_change_idx] <- rnorm(tess_j_rows, 0, var)
  }
  
  # Ensure the result is always a matrix
  if (!is.matrix(tess_j_star)) {
    tess_j_star <- matrix(tess_j_star, ncol = length(dim_j_star))
  }
  
  return(list(tess_j_star, dim_j_star, Modification))
}