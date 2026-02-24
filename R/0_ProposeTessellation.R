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
#' @param sd The standard deviation parameter for the normal distribution (`rnorm`) 
#'   used to generate new coordinate values for centres and dimensions.
#' @param mu The mean parameter for the normal distribution used to generate new
#'   coordinate values for centres and dimensions.
#' @param covariateIndices This parameter is accepted by the function but is
#'   **not used** in the current implementation.
#' @param NumCovariates An integer giving the total number of covariates
#'   available for selection in the full dataset.
#' @param metric Either "Euclidean" or "Spherical".
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
#' 
#' @useDynLib AddiVortes, .registration = TRUE
proposeTessellation <- function(tess_j, dim_j, sd, mu, covariateIndices,
                                NumCovariates, metric) {
  # Call the C++ implementation via the .Call interface
  results <- .Call("propose_tessellation_cpp",
                   tess_j,
                   dim_j,
                   as.double(sd),
                   as.double(mu),
                   as.integer(NumCovariates),
                   as.integer(metric))
  
  return(results)
}