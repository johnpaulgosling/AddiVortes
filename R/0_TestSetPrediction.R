#' @title Test Set Prediction
#'
#' @description A function that derives a prediction for the output variable
#' give an observation, for one iteration in AddiVortes.
#'
#' @param x A matrix with the observations.
#' @param m The number of iterations in AddiVortes.
#' @param Tess A list with the tessellations.
#' @param Dim A list with the dimensions of the tessellations.
#' @param Pred A list with the predictions.
#'
#' @return A vector with the predictions.
#'
#' @keywords internal
#' @noRd
#'
testSetPrediction <- function(x, m,
                              Tess, Dim, Pred) {
  # Use lapply to iterate from 1 to m
  predictionList <- lapply(1:m, function(j) {
    NewTessIndexes <- cellIndices(x, Tess[[j]], Dim[[j]])
    Pred[[j]][NewTessIndexes]
  })

  # Sum up the predictions from each iteration
  # If prediction_list is empty (e.g. if m = 0), rowSums would error.
  # So, initialize Prediction and then sum if there's anything to sum.
  Prediction <- rep(0, nrow(x))
  if (length(predictionList) > 0) {
    Prediction <- rowSums(do.call(cbind, predictionList))
  }

  return(Prediction)
}
