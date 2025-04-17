#' Generate Random Deviates from the Inverse Gamma Distribution (Internal)
#'
#' This function generates random deviates from the inverse gamma distribution,
#' providing the same functionality as `invgamma::rinvgamma` without adding
#' an external dependency.
#'
#' The inverse gamma distribution with shape $a$ (alpha) and rate $b$ (beta)
#' has the probability density function:
#' $f(x) = \frac{b^a}{\Gamma(a)} x^{-a-1} \exp(-b/x)$
#' for $x > 0$, $a > 0$, $b > 0$.
#' This implementation uses the relationship that if $Y \sim Gamma(shape=a, rate=b)$,
#' then $1/Y \sim InvGamma(shape=a, rate=b)$.
#'
#' @param n Numeric, the number of observations required.
#' @param shape Numeric, the shape parameter ($a$). Must be positive.
#' @param rate Numeric, the rate parameter ($b$). Must be positive.
#'   Cannot be specified if `scale` is specified. Defaults to 1.
#' @param scale Numeric, the scale parameter, equivalent to `1/rate`. Must be
#'   positive. Cannot be specified if `rate` is specified. Defaults to `1/rate`.
#'
#' @return A numeric vector of length `n` containing random deviates from the
#'   Inverse Gamma(shape, rate) distribution.
#'
#' @references Based on the source code of `invgamma::rinvgamma` by Stefan Kloppenborg.
#'             Relies on `stats::rgamma`.
#'
#' @keywords internal
#' @noRd
#'
rinvgamma_internal <- function(n, shape, rate = 1, scale = 1/rate) {
  # Check n
  if (length(n) != 1 || !is.numeric(n) || n < 0 || floor(n) != n) {
    stop("'n' must be a non-negative integer.")
  }
  # Return numeric(0) if n=0, consistent with rgamma
  if (n == 0) return(numeric(0))

  # Check shape
  if (length(shape) != 1 || !is.numeric(shape) || shape <= 0) {
    stop("'shape' must be a single positive number.")
  }

  # Handle rate vs scale parameters, ensuring only one effective parameter is used
  # and it's positive.
  if (!missing(rate) && !missing(scale)) {
    if (!is.numeric(rate) || length(rate) != 1 || rate <= 0) {
      stop("'rate' must be a single positive number.")
    }
    if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
      stop("'scale' must be a single positive number.")
    }
  }

  # Determine the rate parameter to use, based on logic from invgamma source
  if (missing(rate) && !missing(scale)) {
    if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
      stop("'scale' must be a single positive number.")
    }
    actual_rate <- 1 / scale
  } else {
    # Use the value 'rate' has (either user-provided or default=1)
    if (!is.numeric(rate) || length(rate) != 1 || rate <= 0) {
      stop("'rate' must be a single positive number.")
    }
    actual_rate <- rate
  }

  gamma_sample <- stats::rgamma(n = n,
                                shape = shape,
                                rate = actual_rate)
  inv_gamma_sample <- 1 / gamma_sample

  return(inv_gamma_sample)
}
