#' @title Create Covariate Structure
#'
#' @description
#' Converts the structure of the covariates to reflect geometric properties.
#' 
#' @details
#' A covariate data.frame can comprise coordinates taken from an outer product
#' of arbitrarily many Euclidean, Spherical, and Categorical spaces, which must
#' be checked and sanitised before use in tessellations. This function performs
#' the required steps to do so: checking the types of the covariates to ensure
#' compatibility, reordering them to combine Euclidean covariates and categorical
#' covariates, and grouping respective spherical coordinates in the required
#' form for later use.
#' 
#' If membership is not provided, it is inferred from the structure vector (often
#' passed as `metric` in AddiVortes). If the argument `one.hot` is `TRUE`, then
#' the categorical coordinates will be one-hot encoded and so the structure and
#' membership vectors are augmented for later use, with each categorical covariate
#' having (n-1) membership and structure elements, where n is the number of levels.
#' 
#' @param data The data.frame of covariates.
#' @param structure The structure vector, containing combinations of "E", "S", 
#' and "C" for Euclidean, Spherical, and Categorical respectively.
#' @param membership The specific membership groups, as integers; if `NULL`, the
#' membership is inferred from the structure vector.
#' @param one.hot Boolean: are the categorical covariates to be one-hot encoded?
#' 
#' @return A list consisting of:
#'   \item{data}{The (potentially reordered) covariate data.frame.}
#'   \item{structure}{The structure vector detailing the types of covariates.}
#'   item{membership}{The membership vector indicating which covariates belong 
#'   to which group.}
#' 
#' @keywords internal
#' @noRd
covariateStructure_internal <- function(data, structure, membership = NULL, one.hot = TRUE) {
  ## If membership is not provided, it's inferred from the structure
  if (is.null(membership)) {
    membership <- numeric(length(structure))
    membership[structure == "E" | structure == 0] <- 1
    membership[structure == "S" | structure == 1] <- 2
    membership[structure == "C" | structure == 2] <- 3
    structure <- c("E", "S", "C")[membership]
  }
  if (length(structure) != length(membership)) {
    if (length(structure) == 1) structure <- rep(structure, length(membership))
    else structure <- c(structure, rep(structure[1], length(membership)-length(structure)))
  }
  ## A fiddle to ensure that categorical variables live at the end of the list.
  if (any(structure == "C")) {
    new_structure <- character(0)
    new_member <- membership
    member_number <- 1
    which_euc <- which(structure == "E")
    for (idx in which_euc) {
      new_member[which(membership == idx)] <- member_number
      member_number <- member_number + 1
      new_structure <- c(new_structure, "E")
    }
    which_sph <- which(structure == "S")
    for (idx in which_sph) {
      new_member[which(membership == idx)] <- member_number
      member_number <- member_number + 1
      new_structure <- c(new_structure, "S")
    }
    which_cat <- which(structure == "C")
    for (idx in which_cat) {
      new_member[which(membership == idx)] <- member_number
      member_number <- member_number + 1
      new_structure <- c(new_structure, "C")
    }
    membership <- new_member
    structure <- new_structure
  }
  reduced_data <- data[,seq_along(membership),drop=FALSE]
  for (i in seq_along(structure)) {
    struct <- structure[i]
    if (struct == "E") {
      has_lev <- is.factor(reduced_data[,i]) || is.character(reduced_data[,i])
      if (has_lev) {
        struct <- structure[i] <- "C"
      }
    }
    if (struct == "C") {
      has_lev <- is.factor(reduced_data[,i])
      if (!has_lev) {
        is_chr <- is.character(reduced_data[,i])
        if (!is_chr) stop(paste("Covariate", names(reduced_data)[i],
                                "specified as categorical but has neither", 
                                "levels nor is a character."))
        reduced_data[,i] <- as.factor(reduced_data[,i])
      }
    }
  }
  if (any(structure == "S")) {
    members <- membership[structure == "S"]
    for (mem in unique(members)) {
      sphere_data <- reduced_data[,which(membership == mem)]
      has_lev <- sapply(seq_len(ncol(sphere_data)), function(idx) is.factor(sphere_data[,idx]))
      if(any(has_lev)) stop("Covariate(s) specified as spherical have factors.")
      param_extent <- apply(apply(sphere_data, 2, range), 2, diff)
      if (sum(param_extent > pi) > 1) 
        stop(paste("More than one spherical parameter in membership group",
                                                 mem, "has range >pi."))
      if (sum(param_extent > pi) == 1) {
        which_polar <- which(param_extent > pi)
        reduced_data[,membership == mem] <- sphere_data[,c(which(param_extent <= pi),
                                                     which(param_extent > pi))]
      }
    }
  }
  m_order <- order(membership)
  structure <- structure[order(membership)]
  membership <- membership[order(membership)]
  ## If categoricals are going to be one-hot-encoded, then want to pre-emptively
  # expand the structure and membership
  if (one.hot && any(structure == "C")) {
    n_onehot <- sum(sapply(seq_along(reduced_data)[structure == "C"], function(i) length(levels(reduced_data[,i]))-1))
    membership <- membership[structure != "C"]
    if (length(membership) == 0)
      membership <- rep(1, n_onehot)
    else
      membership <- c(membership, rep(max(membership)+1, n_onehot))
    structure <- c(structure[structure != "C"], rep("C", n_onehot))
    structure[structure == "C"] <- "E"
  }
  return(list(
    data = reduced_data,
    structure = structure,
    membership = membership
  ))
}
