#' @title Propose Tessellation
#'
#' @description Propose a new tessellation by adding, removing, changing or
#' swapping dimensions or centres.
#'
#' @param x A matrix of covariates.
#' @param j The index of the tessellation to be modified.
#' @param Tess A list of matrices of the tessellations.
#' @param Dim A list of vectors of the dimensions of the tessellations.
#' @param var The variance of the normal distribution used to sample
#' new coordinates.
#' @param covariateIndices A sequence of integers representing the column indices.
#' @param numCovariates The total number of covariates.
#'
#' @return A list containing the new tessellation, the new dimension matrix
#' and the modification made.
#' @keywords internal
#' @noRd
#'
proposeTessellation <- function(tess_j, dim_j, var, covariateIndices, numCovariates) {
  # This function now operates on a single tessellation matrix (tess_j)
  # and dimension vector (dim_j), not the entire lists.
  
  p <- runif(1, 0, 1)
  
  dim_j_star <- dim_j
  tess_j_star <- tess_j
  
  if (p < 0.2 && length(dim_j) != numCovariates ||
      length(dim_j) == 1 && p < 0.4) {
    # Add Dimension
    available_covs <- covariateIndices[-dim_j]
    if (length(available_covs) == 1) {
      new_dim <- available_covs
    } else {
      new_dim <- sample(available_covs, 1)
    }
    dim_j_star <- c(dim_j, new_dim)
    tess_j_star <- cbind(tess_j, rnorm(nrow(tess_j), 0, var))
    Modification <- "AD"
    
  } else if (p < 0.4) {
    # Remove Dimension
    removed_dim_idx <- sample.int(length(dim_j), 1)
    dim_j_star <- dim_j[-removed_dim_idx]
    tess_j_star <- tess_j[, -removed_dim_idx, drop = FALSE]
    Modification <- "RD"
    
  } else if (p < 0.6 || (p < 0.8 && nrow(tess_j) == 1)) {
    # Add Centre
    tess_j_star <- rbind(tess_j, rnorm(length(dim_j), 0, var))
    Modification <- "AC"
    
  } else if (p < 0.8) {
    # Remove Centre
    if (nrow(tess_j) > 1) {
      center_removed_idx <- sample.int(nrow(tess_j), 1)
      tess_j_star <- tess_j[-center_removed_idx, , drop = FALSE]
    }
    # if only one centre, do nothing, effectively a "Change"
    Modification <- "RC"
    
  } else if (p < 0.9 || length(dim_j) == numCovariates) {
    # Change Centre
    centre_to_change_idx <- sample.int(nrow(tess_j), 1)
    tess_j_star[centre_to_change_idx, ] <- rnorm(length(dim_j), 0, var)
    Modification <- "Change"
    
  } else {
    # Swap Dimension
    available_covs <- covariateIndices[-dim_j]
    dim_to_change_idx <- sample.int(length(dim_j), 1)
    
    if (length(available_covs) == 1) {
      new_dim <- available_covs
    } else {
      new_dim <- sample(available_covs, 1)
    }
    
    dim_j_star[dim_to_change_idx] <- new_dim
    tess_j_star[, dim_to_change_idx] <- rnorm(nrow(tess_j), 0, var)
    Modification <- "Swap"
  }
  
  # Ensure the result is always a matrix, even with one row/column
  tess_j_star <- matrix(tess_j_star, ncol = length(dim_j_star))
  
  return(list(tess_j_star, dim_j_star, Modification))
}