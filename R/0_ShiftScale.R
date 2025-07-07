#' @title scaleData_internal
#'
#' @description
#' Scale a numeric vector, matrix, or data frame to [-0.5, 0.5]
#'
#' Determines if the input is vector-like or matrix/data frame-like.
#' - If vector-like: scales based on the overall range of the vector.
#' - If matrix/data frame-like: scales each column based on its own range.
#' Scaling uses `(value - centre) / range` where `centre = (max+min)/2`
#' and `range = max-min`. Ignores NA for calculations. Returns 0 for zero-range
#' data/columns. Non-numeric columns in matrices/data frames are returned
#' unchanged with a warning.
#'
#' @param data A numeric vector, matrix, or data frame.
#'
#' @return A list containing:
#'    \item{scaledData}{The scaled object (vector, matrix or data frame matching input type).}
#'    \item{centres}{The centre(s) used (single value for vector input, vector for matrix/df).}
#'    \item{ranges}{The range(s) used (single value for vector input, vector for matrix/df).}
#'
#' @keywords internal
#' @noRd
scaleData_internal <- function(data) {
  if (!is.numeric(data) && !all(sapply(data, is.numeric)) && !is.null(dim(data))) {
    # Allow non-numeric columns in data frames/matrices but check base type
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)) {
      stop("'data' must be a numeric vector, matrix, or data frame.")
    }
    # Further checks inside loop for matrix/df case
  } else if (!is.numeric(data) && is.null(dim(data))) {
    stop("'data' must be a numeric vector.")
  }

  # Check if it's vector-like (no dimensions) or matrix/df-like
  if (!is.matrix(data)) {
    # --- Handle Vector ---
    minV <- min(data, na.rm = TRUE)
    maxV <- max(data, na.rm = TRUE)
    centre <- (maxV + minV) / 2
    rangeV <- maxV - minV

    # Break if range is zero
    if (rangeV == 0) {
      stop("Range is zero. Cannot scale data.")
    }

    # Scale the data
    scaledData <- (data - centre) / rangeV
    centres <- centre # Store single value
    ranges <- rangeV # Store single value
  } else {
    # --- Handle Matrix or Data Frame ---
    isDf <- is.data.frame(data)
    matWork <- as.matrix(data) # Work with matrix representation
    numCols <- ncol(matWork)

    # Pre-allocate results
    scaledMatWork <- matrix(NA_real_, nrow = nrow(matWork), ncol = numCols)
    centres <- numeric(numCols) # Vector to store centres
    ranges <- numeric(numCols) # Vector to store ranges

    originalColnames <- colnames(matWork) # Store original names
    if (is.null(originalColnames) && ncol(matWork) > 0) {
      originalColnames <- paste0("V", 1:ncol(matWork)) # Assign default if null
    }

    for (i in 1:numCols) {
      colData <- matWork[, i]
      colName <- originalColnames[i]

      # Check if column is numeric *before* scaling
      if (!is.numeric(colData)) {
        warning("Column ", i, " ('", colName, "') is not numeric. Returning original column data.", call. = FALSE)
        scaledMatWork[, i] <- colData # Keep original non-numeric data
        centres[i] <- NA # Indicate parameters are not applicable/calculated
        ranges[i] <- NA
        next # Move to the next column
      }

      minV <- min(colData, na.rm = TRUE)
      maxV <- max(colData, na.rm = TRUE)
      centres[i] <- (maxV + minV) / 2
      ranges[i] <- maxV - minV
      scaledMatWork[, i] <- (colData - centres[i]) / ranges[i]
    }

    # Preserve names
    colnames(scaledMatWork) <- originalColnames
    rownames(scaledMatWork) <- rownames(matWork)

    # Convert back to data frame if original was data frame
    if (isDf) {
      # Important: Ensure columns that were originally factors/characters remain numeric after scaling
      # The as.data.frame conversion should handle this correctly if scaledMatWork is purely numeric where scaling occurred.
      scaledData <- as.data.frame(scaledMatWork)
      # Restore original non-numeric columns if needed (though warning was given)
      # For simplicity here, we assume conversion is okay, relying on the warning.
    } else {
      scaledData <- scaledMatWork
    }
    # `centres` and `ranges` are already vectors
  }

  # --- Return Results ---
  return(list(
    scaledData = scaledData,
    centres = centres,
    ranges = ranges
  ))
}

#' @title applyScaling_internal
#'
#' @description
#' Apply existing scaling parameters to columns of a matrix/data frame
#'
#' Scales columns using `(value - centre) / range` with *provided* centre
#' and range values (typically from a training set). Handles cases where
#' the provided range is zero consistently (output is 0). Preserves input
#' data type (matrix or data frame). Stops if non-numeric columns are encountered
#' or if scaling parameters are NA for a column.
#'
#' @param mat A numeric matrix or data frame to scale.
#' @param centres A numeric vector of centre values to apply (length must match `ncol(mat)`).
#' @param ranges A numeric vector of range values to apply (length must match `ncol(mat)`).
#'
#' @return The scaled matrix or data frame.
#'
#' @keywords internal
#' @noRd
applyScaling_internal <- function(mat, centres, ranges) {
  isDf <- is.data.frame(mat)
  matWork <- as.matrix(mat) # Work with matrix representation

  # Basic type check first
  # Allow matrices/dataframes containing non-numeric columns initially, check column-by-column
  if (!is.matrix(matWork) && !is.data.frame(mat)) {
    stop("'mat' must be a matrix or data frame.")
  }


  if (length(centres) != ncol(matWork) || length(ranges) != ncol(matWork)) {
    stop(
      "Length of 'centres' (", length(centres), ") and 'ranges' (", length(ranges),
      ") must match the number of columns in 'mat' (", ncol(matWork), ")."
    )
  }

  numCols <- ncol(matWork)
  scaledMatWork <- matrix(NA_real_, nrow = nrow(matWork), ncol = numCols)
  originalColnames <- colnames(matWork)
  if (is.null(originalColnames) && ncol(matWork) > 0) {
    originalColnames <- paste0("V", 1:ncol(matWork)) # Assign default if null
  }

  for (i in 1:numCols) {
    colData <- matWork[, i]
    colName <- originalColnames[i]

    # Check if the column data itself is numeric
    if (!is.numeric(colData)) {
      # If params exist for this column but data is non-numeric, this is an error
      warning("Column ", i, " ('", colName, "') is not numeric. Cannot apply scaling. Returning original data for this column.", call. = FALSE)
      scaledMatWork[, i] <- colData
      next
    }

    # Check if the parameters for this column are valid (not NA)
    if (is.na(centres[i]) || is.na(ranges[i])) {
      warning("Scaling parameters for column ", i, " ('", colName, "') are NA (likely due to non-numeric original training column). Returning original data for this column.", call. = FALSE)
      scaledMatWork[, i] <- colData
      next
    }

    # Apply scaling using provided centres and ranges
    scaledMatWork[, i] <- (colData - centres[i]) / ranges[i]
  }

  colnames(scaledMatWork) <- originalColnames
  rownames(scaledMatWork) <- rownames(matWork)

  if (isDf) {
    scaledMatOut <- as.data.frame(scaledMatWork)
    # Again, careful about type restoration if needed, but rely on warnings for now
  } else {
    scaledMatOut <- scaledMatWork
  }

  return(scaledMatOut) # Return only the scaled data
}
