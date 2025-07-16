#' @title Calculate Partial Residuals for a Tessellation
#' @description Calculates the vector of partial residuals by excluding the effect
#'  of the jth tessellation. It then aggregates these residuals and observation
#'  counts for both the current and the proposed new tessellation structures.
#'
#' @param y The numeric vector of the outcome or response variable.
#' @param j The index of the tessellation currently being updated.
#' @param SumOfAllTess A numeric vector representing the sum of predictions from
#'  all tessellations.
#' @param Pred A list where each element contains the predicted values for a
#'  tessellation.
#' @param lastTessPred The prediction from the jth tessellation in the previous
#'  iteration, used to correctly update the sum of all tessellation predictions.
#' @param indexes A numeric vector of integers mapping each observation to a
#'  cell in the current tessellation.
#' @param indexesStar A numeric vector of integers mapping each observation to a
#'  cell in the proposed new tessellation.
#' @param num_centres_new The number of centres (cells) in the proposed new
#'  tessellation.
#'
#' @return A list containing seven elements in the following order:
#'  \enumerate{
#'    \item A numeric vector of summed residuals for each cell in the current tessellation (`R_ijOld`).
#'    \item An integer vector of observation counts for each cell in the current tessellation (`n_ijOld`).
#'    \item A numeric vector of summed residuals for each cell in the proposed tessellation (`R_ijNew`).
#'    \item An integer vector of observation counts for each cell in the proposed tessellation (`n_ijNew`).
#'    \item The updated sum of all tessellation predictions (`SumOfAllTess`).
#'    \item The integer vector of cell indexes for the proposed tessellation (`indexesStar`).
#'    \item The integer vector of cell indexes for the current tessellation (`indexes`).
#'  }
#'
#' @keywords internal
#' @noRd
#' 
#' @useDynLib AddiVortes, .registration = TRUE
calculateResiduals <- function(y, j, SumOfAllTess, Pred, lastTessPred,
                               indexes, indexesStar, num_centres_new) {
  # This part of the logic remains in R
  if (j == 1) {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess - CurrentTessPred
  } else {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess + lastTessPred - CurrentTessPred
  }
  R_j <- y - SumOfAllTess
  
  # --- Call the C++ function via the .Call interface ---
  num_levels_old <- length(Pred[[j]])
  cpp_results <- .Call("calculate_residuals_cpp",
                       as.double(R_j),
                       as.integer(indexes),
                       as.integer(indexesStar),
                       as.integer(num_levels_old),
                       as.integer(num_centres_new))
  
  # --- Return the results in the original list format ---
  # The .Call function returns a named list, which we access using $
  return(list(
    cpp_results$R_ijOld, cpp_results$n_ijOld,
    cpp_results$R_ijNew, cpp_results$n_ijNew,
    SumOfAllTess, indexesStar, indexes
  ))
}