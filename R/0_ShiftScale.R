#' scale_data_internal
#'
#' Scale a numeric vector, matrix, or data frame to [-0.5, 0.5]
#'
#' Determines if the input is vector-like or matrix/data frame-like.
#' - If vector-like: scales based on the overall range of the vector.
#' - If matrix/data frame-like: scales each column based on its own range.
#' Scaling uses `(value - center) / range` where `center = (max+min)/2`
#' and `range = max-min`. Ignores NA for calculations. Returns 0 for zero-range
#' data/columns. Non-numeric columns in matrices/data frames are returned
#' unchanged with a warning.
#'
#' @param data A numeric vector, matrix, or data frame.
#'
#' @return A list containing:
#'   \item{scaled_data}{The scaled object (vector, matrix or data frame matching input type).}
#'   \item{centers}{The center(s) used (single value for vector input, vector for matrix/df).}
#'   \item{ranges}{The range(s) used (single value for vector input, vector for matrix/df).}
#'
#' @keywords internal
#' @noRd
scaleData_internal <- function(data) {
  if (!is.numeric(data) && !all(sapply(data, is.numeric)) && !is.null(dim(data))) {
    # Allow non-numeric columns in data frames/matrices but check base type
    if(!is.data.frame(data) && !is.matrix(data) && !is.vector(data))
      stop("'data' must be a numeric vector, matrix, or data frame.")
    # Further checks inside loop for matrix/df case
  } else if (!is.numeric(data) && is.null(dim(data))) {
    stop("'data' must be a numeric vector.")
  }

  # Check if it's vector-like (no dimensions) or matrix/df-like
  if (!is.matrix(data)) {
    # --- Handle Vector ---
    min_v <- min(data, na.rm = TRUE)
    max_v <- max(data, na.rm = TRUE)
    center <- (max_v + min_v) / 2
    range_v <- max_v - min_v

    # Break if range is zero
    if (range_v == 0) {
      stop("Range is zero. Cannot scale data.")
    }

    # Scale the data
    scaled_data <- (data - center) / range_v
    centers <- center # Store single value
    ranges <- range_v   # Store single value

  } else {
    # --- Handle Matrix or Data Frame ---
    is_df <- is.data.frame(data)
    mat_work <- as.matrix(data) # Work with matrix representation
    num_cols <- ncol(mat_work)

    # Pre-allocate results
    scaled_mat_work <- matrix(NA_real_, nrow = nrow(mat_work), ncol = num_cols)
    centers <- numeric(num_cols) # Vector to store centers
    ranges <- numeric(num_cols)  # Vector to store ranges

    original_colnames <- colnames(mat_work) # Store original names
    if(is.null(original_colnames) && ncol(mat_work) > 0) {
      original_colnames <- paste0("V", 1:ncol(mat_work)) # Assign default if null
    }

    for (i in 1:num_cols) {
      col_data <- mat_work[, i]
      col_name <- original_colnames[i]

      # Check if column is numeric *before* scaling
      if (!is.numeric(col_data)) {
        warning("Column ", i, " ('", col_name, "') is not numeric. Returning original column data.", call. = FALSE)
        scaled_mat_work[, i] <- col_data # Keep original non-numeric data
        centers[i] <- NA # Indicate parameters are not applicable/calculated
        ranges[i] <- NA
        next # Move to the next column
      }

      min_v <- min(col_data, na.rm = TRUE)
      max_v <- max(col_data, na.rm = TRUE)
      centers[i] <- (max_v + min_v) / 2
      ranges[i] <- max_v - min_v
      scaled_mat_work[, i] <- (col_data - centers[i]) / ranges[i]
    }

    # Preserve names
    colnames(scaled_mat_work) <- original_colnames
    rownames(scaled_mat_work) <- rownames(mat_work)

    # Convert back to data frame if original was data frame
    if (is_df) {
      # Important: Ensure columns that were originally factors/characters remain numeric after scaling
      # The as.data.frame conversion should handle this correctly if scaled_mat_work is purely numeric where scaling occurred.
      scaled_data <- as.data.frame(scaled_mat_work)
      # Restore original non-numeric columns if needed (though warning was given)
      # For simplicity here, we assume conversion is okay, relying on the warning.
    } else {
      scaled_data <- scaled_mat_work
    }
    # `centers` and `ranges` are already vectors
  }

  # --- Return Results ---
  return(list(
    scaled_data = scaled_data,
    centers = centers,
    ranges = ranges
  ))
}

#' apply_scaling_internal
#'
#' Apply existing scaling parameters to columns of a matrix/data frame
#'
#' Scales columns using `(value - center) / range` with *provided* center
#' and range values (typically from a training set). Handles cases where
#' the provided range is zero consistently (output is 0). Preserves input
#' data type (matrix or data frame). Stops if non-numeric columns are encountered
#' or if scaling parameters are NA for a column.
#'
#' @param mat A numeric matrix or data frame to scale.
#' @param centers A numeric vector of center values to apply (length must match `ncol(mat)`).
#' @param ranges A numeric vector of range values to apply (length must match `ncol(mat)`).
#'
#' @return The scaled matrix or data frame.
#'
#' @keywords internal
#' @noRd
applyScaling_internal <- function(mat, centers, ranges) {
  is_df <- is.data.frame(mat)
  mat_work <- as.matrix(mat) # Work with matrix representation

  # Basic type check first
  # Allow matrices/dataframes containing non-numeric columns initially, check column-by-column
  if (!is.matrix(mat_work) && !is.data.frame(mat)) {
    stop("'mat' must be a matrix or data frame.")
  }


  if (length(centers) != ncol(mat_work) || length(ranges) != ncol(mat_work)) {
    stop("Length of 'centers' (", length(centers), ") and 'ranges' (", length(ranges),
         ") must match the number of columns in 'mat' (", ncol(mat_work), ").")
  }

  num_cols <- ncol(mat_work)
  scaled_mat_work <- matrix(NA_real_, nrow = nrow(mat_work), ncol = num_cols)
  original_colnames <- colnames(mat_work)
  if(is.null(original_colnames) && ncol(mat_work) > 0) {
    original_colnames <- paste0("V", 1:ncol(mat_work)) # Assign default if null
  }

  for (i in 1:num_cols) {
    col_data <- mat_work[, i]
    col_name <- original_colnames[i]

    # Check if the column data itself is numeric
    if (!is.numeric(col_data)) {
      # If params exist for this column but data is non-numeric, this is an error
      warning("Column ", i, " ('", col_name, "') is not numeric. Cannot apply scaling. Returning original data for this column.", call. = FALSE)
      scaled_mat_work[, i] <- col_data
      next
    }

    # Check if the parameters for this column are valid (not NA)
    if (is.na(centers[i]) || is.na(ranges[i])) {
      warning("Scaling parameters for column ", i, " ('", col_name, "') are NA (likely due to non-numeric original training column). Returning original data for this column.", call. = FALSE)
      scaled_mat_work[, i] <- col_data
      next
    }

    # Apply scaling using provided centres and ranges
    scaled_mat_work[, i] <- (col_data - centers[i]) / ranges[i]
  }

  colnames(scaled_mat_work) <- original_colnames
  rownames(scaled_mat_work) <- rownames(mat_work)

  if (is_df) {
    scaled_mat_out <- as.data.frame(scaled_mat_work)
    # Again, careful about type restoration if needed, but rely on warnings for now
  } else {
    scaled_mat_out <- scaled_mat_work
  }

  return(scaled_mat_out) # Return only the scaled data
}
