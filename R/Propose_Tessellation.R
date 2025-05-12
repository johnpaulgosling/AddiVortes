#' @title Propose_Tessellation
#'
#' @description Propose a new tessellation by adding, removing, changing or swapping dimensions or centers.
#'
#' @param x A matrix of covariates.
#' @param j The index of the tessellation to be modified.
#' @param Tess A list of matrices of the tessellations.
#' @param Dim A list of vectors of the dimensions of the tessellations.
#' @param var The variance of the normal distribution used to sample new coordinates.
#'
#' @return A list containing the new tessellation, the new dimension matrix and the modification made.
#'
#' @export
Propose_Tessellation <- function(x, j, Tess, Dim, var) { # Propose a new tessellation

  p <- runif(1, 0, 1) # Randomly sample p to decide the proposed modification to tessellation.

  DimStar <- Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar <- Tess # Similar for the tessellation matrix.

  if (p < 0.2 & length(Dim[[j]]) != length(x[1, ]) | length(Dim[[j]]) == 1 & p < 0.4) {
    # Add a dimension if p is less then 0.2 or if p is less then 0.4 when
    # there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
    NumberOfCovariates <- 1:length(x[1, ]) # Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates <- NumberOfCovariates[-Dim[[j]]] # Remove all values that for the covariates that are already in the tessellation.
    DimStar[[j]] <- c(Dim[[j]], sample(NumberOfCovariates, 1)) # Uniformly sample a new covariate and add it to the dimension matrix.
    TessStar[[j]] <- cbind(Tess[[j]], rnorm(length(Tess[[j]][, 1]), 0, var)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
    Modification <- "AD"
  } else if (p < 0.4) {
    # Remove a dimension if p is less then 0.4.
    RemovedDim <- sample(1:length(Dim[[j]]), 1) # Uniformly sample the dimension to be removed.
    DimStar[[j]] <- DimStar[[j]][-RemovedDim] # Remove the dimension from the dimesion Matrix.
    TessStar[[j]] <- matrix(TessStar[[j]][, -RemovedDim], ncol = length(DimStar[[j]])) # Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
    Modification <- "RD"
  } else if (p < 0.6 || p < 0.8 & length(Tess[[j]][, 1]) == 1) {
    # Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
    TessStar[[j]] <- rbind(Tess[[j]], rnorm(length(Dim[[j]]), 0, var)) # Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
    Modification <- "AC"
  } else if (p < 0.8) {
    # Add a centre if p is less then 0.8.
    CenterRemoved <- sample(1:length(TessStar[[j]][, 1]), 1) # Sample a row.
    TessStar[[j]] <- matrix(TessStar[[j]][-CenterRemoved, ], ncol = length(Dim[[j]])) # Remove row sampled.
    Modification <- "RC"
  } else if (p < 0.9 || length(Dim[[j]]) == length(x[1, ])) {
    # Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
    TessStar[[j]][sample(1:length(TessStar[[j]][, 1]), 1), ] <- rnorm(length(Dim[[j]]), 0, var) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
    Modification <- "Change"
  } else {
    # Swap a dimension.
    NumberOfCovariates <- 1:length(x[1, ]) # Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates <- NumberOfCovariates[-Dim[[j]]] # Remove all values that for the covariates that are already in the tessellation.
    DimToChange <- sample(1:length(Dim[[j]]), 1) # Uniformly sample a dimension to change.
    DimStar[[j]][DimToChange] <- sample(NumberOfCovariates, 1) # Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
    TessStar[[j]][, DimToChange] <- rnorm(length(Tess[[j]][, 1]), 0, var) # Add new normally sampled coordinates new dimension added.
    Modification <- "Swap"
  }

  TessStar[[j]] <- matrix(TessStar[[j]], ncol = length(DimStar[[j]])) # Ensure the the Tessellation matrix is a "matrix" type.

  return(list(TessStar, DimStar, Modification)) # Return new proposed tessellation.
}
