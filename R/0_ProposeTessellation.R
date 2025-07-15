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
proposeTessellation <- function(x, j, Tess, Dim, var, covariateIndices, numCovariates) { # <-- ADDED ARGUMENT
  p <- runif(1, 0, 1)
  
  DimStar <- Dim
  TessStar <- Tess
  if (p < 0.2 && length(Dim[[j]]) != numCovariates || # <-- MODIFIED LINE
      length(Dim[[j]]) == 1 && p < 0.4) {
    NumberOfCovariates <- covariateIndices[-Dim[[j]]]
    DimStar[[j]] <- c(Dim[[j]], sample(NumberOfCovariates, 1))
    TessStar[[j]] <- cbind(Tess[[j]], rnorm(length(Tess[[j]][, 1]), 0, var))
    Modification <- "AD"
  } else if (p < 0.4) {
    RemovedDim <- sample(seq_along(Dim[[j]]), 1)
    DimStar[[j]] <- DimStar[[j]][-RemovedDim]
    TessStar[[j]] <- matrix(TessStar[[j]][, -RemovedDim],
                            ncol = length(DimStar[[j]])
    )
    Modification <- "RD"
  } else if (p < 0.6 || p < 0.8 && length(Tess[[j]][, 1]) == 1) {
    TessStar[[j]] <- rbind(Tess[[j]], rnorm(length(Dim[[j]]), 0, var))
    Modification <- "AC"
  } else if (p < 0.8) {
    CenterRemoved <- sample(seq_along(TessStar[[j]][, 1]), 1)
    TessStar[[j]] <- matrix(TessStar[[j]][-CenterRemoved, ],
                            ncol = length(Dim[[j]])
    )
    Modification <- "RC"
  } else if (p < 0.9 || length(Dim[[j]]) == numCovariates) { # <-- MODIFIED LINE
    TessStar[[j]][sample(
      seq_along(TessStar[[j]][, 1]),
      1
    ), ] <- rnorm(length(Dim[[j]]), 0, var)
    Modification <- "Change"
  } else {
    NumberOfCovariates <- covariateIndices[-Dim[[j]]]
    DimToChange <- sample(seq_along(Dim[[j]]), 1)
    DimStar[[j]][DimToChange] <- sample(NumberOfCovariates, 1)
    TessStar[[j]][, DimToChange] <- rnorm(length(Tess[[j]][, 1]), 0, var)
    Modification <- "Swap"
  }
  
  TessStar[[j]] <- matrix(TessStar[[j]], ncol = length(DimStar[[j]]))
  
  return(list(TessStar, DimStar, Modification))
}