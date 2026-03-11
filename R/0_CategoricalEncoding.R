#' @title Encode Categorical Covariates as Binary Variables
#'
#' @description
#' Converts categorical (character or factor) columns into d-1 binary indicator
#' variables (one-hot encoding with a reference category dropped). Binary columns
#' take values 0 or \code{catScaling} (default 1).
#'
#' @details
#' For a categorical variable with d unique levels, d-1 binary columns are created
#' by dropping the first level (sorted alphabetically for characters, or by factor
#' level order for factors) as the reference category. The resulting binary columns
#' are named \code{<original_col>_<level>}.
#'
#' When \code{encoding} is supplied (e.g., during prediction), the stored levels
#' from the training data are used. Any level in new data not seen during training
#' is treated as the reference level (all indicator columns equal 0).
#'
#' @param x A matrix or data frame potentially containing character or factor columns.
#' @param catScaling Numeric scalar. The value to assign to the "active" category
#'   in the binary indicator (default 1). Controls how much weight categorical
#'   differences receive in distance calculations.
#' @param encoding If \code{NULL} (default), the encoding is inferred from \code{x}.
#'   If provided (as returned previously from this function), that encoding is
#'   applied to \code{x} instead.
#'
#' @return A list with:
#'   \item{encoded}{A numeric matrix with categorical columns replaced by binary
#'     indicator columns.}
#'   \item{encoding}{A list of encoding metadata needed to re-apply the same
#'     transformation to new data, or \code{NULL} if no categorical columns were
#'     found.}
#'
#' @keywords internal
#' @noRd
encodeCategories_internal <- function(x, catScaling = 1, encoding = NULL) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  nCols <- ncol(x)
  colNames <- names(x)

  # Identify categorical columns (character or factor)
  isCat <- sapply(x, function(col) is.character(col) || is.factor(col))
  catColIndices <- unname(which(isCat))

  # If no categorical columns and no stored encoding, return as-is
  if (length(catColIndices) == 0 && is.null(encoding)) {
    encoded <- matrix(as.numeric(as.matrix(x)), nrow = nrow(x), ncol = nCols)
    colnames(encoded) <- colNames
    return(list(encoded = encoded, encoding = NULL))
  }

  isBuilding <- is.null(encoding)
  if (isBuilding) {
    # Build encoding from current data
    colEncodings <- vector("list", nCols)
    for (j in catColIndices) {
      colData <- x[[j]]
      levs <- if (is.factor(colData)) {
        levels(colData)
      } else {
        sort(unique(as.character(colData)))
      }
      colEncodings[[j]] <- list(levels = levs)
    }
  } else {
    catColIndices <- encoding$catColIndices
    colEncodings <- encoding$colEncodings
    catScaling <- encoding$catScaling
  }

  # Build encoded matrix column by column
  resultCols <- list()
  resultNames <- character(0)
  encodedBinaryCols <- integer(0)
  colIdx <- 0L

  for (j in seq_len(nCols)) {
    if (j %in% catColIndices) {
      levs <- colEncodings[[j]]$levels
      levs_non_ref <- levs[-1]  # drop first level (reference category)
      colData <- as.character(x[[j]])
      for (lev in levs_non_ref) {
        colIdx <- colIdx + 1L
        resultCols[[length(resultCols) + 1L]] <- as.numeric(colData == lev) * catScaling
        resultNames <- c(resultNames, paste0(colNames[j], "_", lev))
        encodedBinaryCols <- c(encodedBinaryCols, colIdx)
      }
    } else {
      colIdx <- colIdx + 1L
      resultCols[[length(resultCols) + 1L]] <- as.numeric(x[[j]])
      resultNames <- c(resultNames, colNames[j])
    }
  }

  encoded <- do.call(cbind, resultCols)
  colnames(encoded) <- resultNames

  if (isBuilding) {
    encoding <- list(
      catColIndices    = catColIndices,
      colEncodings     = colEncodings,
      catScaling       = catScaling,
      origNCols        = nCols,
      origColNames     = colNames,
      encodedBinaryCols = encodedBinaryCols
    )
  }

  return(list(encoded = encoded, encoding = encoding))
}
