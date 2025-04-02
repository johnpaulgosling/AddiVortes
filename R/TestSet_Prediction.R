#' @title TestSet_Prediction
#'
#' @description A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
#'
#' @param x A matrix with the observations.
#' @param m The number of iterations in AddiVortes.
#' @param Tess A list with the tessellations.
#' @param Dim A list with the dimensions of the tessellations.
#' @param Pred A list with the predictions.
#'
#' @return A vector with the predictions.
#'
#' @export
TestSet_Prediction <- function(x, m, Tess, Dim, Pred) {
  Prediction <- rep(0, length(x[, 1]))
  for (j in 1:m) {
    NewTessIndexes <- Cell_Indexes(x, Tess[[j]], Dim[[j]])
    Prediction <- Prediction + Pred[[j]][NewTessIndexes]
  }
  
  return(Prediction)
}
