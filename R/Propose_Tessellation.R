#' @title Propose_Tessellation
#'
#' @description Propose a new tessellation by adding, removing, changing or
#' swapping dimensions or centres.
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
Propose_Tessellation <- function(x, j,
                                 Tess, Dim, var) {

  p <- runif(1, 0, 1) # Randomly sample p to decide the proposed modification to tessellation.

  DimStar <- Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar <- Tess # Similar for the tessellation matrix.

  if (p < 0.2 & length(Dim[[j]]) != length(x[1, ]) | length(Dim[[j]]) == 1 & p < 0.4) {
    # Add a dimension if p is less then 0.2 or if p is less then 0.4 when
    # there is only one dimension in the Tessellation due to adjustments.

    # Let NumberOfCovariates be a vector from 1 to the number of covariates
    # considered.
    NumberOfCovariates <- 1:length(x[1, ])
    # Remove all values that for the covariates that are already in the
    # tessellation.
    NumberOfCovariates <- NumberOfCovariates[-Dim[[j]]]
    # Uniformly sample a new covariate and add it to the dimension matrix.
    DimStar[[j]] <- c(Dim[[j]], sample(NumberOfCovariates, 1))
    # Sample new coordinates from Normal distribution for the new dimension and
    # add it to the Tessellation matrix.
    TessStar[[j]] <- cbind(Tess[[j]], rnorm(length(Tess[[j]][, 1]), 0, var))
    Modification <- "AD"
  } else if (p < 0.4) {
    # Remove a dimension if p is less then 0.4.

    # Uniformly sample the dimension to be removed.
    RemovedDim <- sample(1:length(Dim[[j]]), 1)
    # Remove the dimension from the dimension matrix.
    DimStar[[j]] <- DimStar[[j]][-RemovedDim]
    # Remove the coordinates in the Tessellation matrix corresponding to the
    # dimension removed.
    TessStar[[j]] <- matrix(TessStar[[j]][, -RemovedDim],
                            ncol = length(DimStar[[j]]))
    Modification <- "RD"
  } else if (p < 0.6 || p < 0.8 & length(Tess[[j]][, 1]) == 1) {
    # Add a centre if p is less then 0.6 or if p is less then 0.4
    # when there is only one center in the Tessellation due to adjustments.

    # Add a new row of coordinates, sampled from a normal distribution,
    # to the Tessellation matrix to add a center.
    TessStar[[j]] <- rbind(Tess[[j]], rnorm(length(Dim[[j]]), 0, var))
    Modification <- "AC"
  } else if (p < 0.8) {
    # Add a centre if p is less then 0.8.

    # Sample a row.
    CenterRemoved <- sample(1:length(TessStar[[j]][, 1]), 1)
    # Remove row sampled.
    TessStar[[j]] <- matrix(TessStar[[j]][-CenterRemoved, ],
                            ncol = length(Dim[[j]]))
    Modification <- "RC"
  } else if (p < 0.9 || length(Dim[[j]]) == length(x[1, ])) {
    # Change a center if p is less then 0.9 or if the all the covariates are in
    # the tessellation.

    # Sample a row in the tessellaion matrix and change the coordinates of the
    # centre by sampling from a normal distribution.
    TessStar[[j]][sample(1:length(TessStar[[j]][, 1]),
                         1), ] <- rnorm(length(Dim[[j]]), 0, var)
    Modification <- "Change"
  } else {
    # Swap a dimension.

    # Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates <- 1:length(x[1, ])
    # Remove all values that for the covariates that are already in the tessellation.
    NumberOfCovariates <- NumberOfCovariates[-Dim[[j]]]
    # Uniformly sample a dimension to change.
    DimToChange <- sample(1:length(Dim[[j]]), 1)
    # Replace the Dimension to a new uniformly sampled covariate that is not
    # already in the tessellation.
    DimStar[[j]][DimToChange] <- sample(NumberOfCovariates, 1)
    # Add new normally sampled coordinates new dimension added.
    TessStar[[j]][, DimToChange] <- rnorm(length(Tess[[j]][, 1]), 0, var)
    Modification <- "Swap"
  }

  # Ensure the the Tessellation matrix is a "matrix" type.
  TessStar[[j]] <- matrix(TessStar[[j]], ncol = length(DimStar[[j]]))

  return(list(TessStar, DimStar, Modification))
}
