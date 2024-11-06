#' Add Two Numbers
#'
#' Adds two numbers together.
#'
#' @param x The first number.
#' @param y The second number.
#'
#' @return The sum of `x` and `y`.
#'
#' @examples
#' add_numbers(2, 3)
#'
#' @export
add_numbers <- function(x, y) {
  # Somewhere to store error messages
  errors <- character()

  # Check x is a number
  if (!is.numeric(x)) errors <- c(errors, '`x` must be a number')

  # Check y is a number
  if (!is.numeric(y)) errors <- c(errors, '`y` must be a number')

  # Return any errors
  if (length(errors) > 0) stop(paste(errors,
                                     '\n  '))

  # Return the sum of x and y
  return(x + y)
}
