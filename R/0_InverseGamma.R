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

#' Calculate Quantiles of the Inverse Gamma Distribution (Internal)
#'
#' This function computes quantiles (the inverse CDF) of the inverse gamma distribution,
#' providing similar functionality to potentially available functions in external packages
#' (like `invgamma::qinvgamma`) without adding an external dependency.
#'
#' The inverse gamma distribution with shape $a$ (alpha) and rate $b$ (beta)
#' has the probability density function:
#' $f(x) = \frac{b^a}{\Gamma(a)} x^{-a-1} \exp(-b/x)$
#' for $x > 0$, $a > 0$, $b > 0$.
#' This implementation uses the relationship that if $Y \sim Gamma(shape=a, rate=b)$,
#' and $Q_Y(p, lower.tail=FALSE)$ is the quantile function of $Y$ giving the value $y$
#' such that $P(Y \ge y) = p$, then the quantile $x$ for the inverse gamma
#' distribution such that $P(X \le x) = p$ is given by $x = 1 / Q_Y(p, lower.tail=FALSE)$.
#' It relies on `stats::qgamma`.
#'
#' @param p Numeric vector of probabilities. Values outside $[0, 1]$ (or the
#'  corresponding range if `log.p = TRUE`) will result in `NaN`.
#' @param shape Numeric, the shape parameter ($a$). Must be positive.
#' @param rate Numeric, the rate parameter ($b$). Must be positive.
#'  Cannot be specified if `scale` is specified. Defaults to 1.
#' @param scale Numeric, the scale parameter, equivalent to `1/rate`. Must be
#'  positive. Cannot be specified if `rate` is specified. Defaults to `1/rate`.
#' @param lower.tail Logical; if TRUE (default), probabilities are $P[X \le x]$,
#'  otherwise, $P[X > x]$.
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#'
#' @return A numeric vector of quantiles corresponding to the probabilities `p`.
#'
#' @references Based on the relationship between Gamma and Inverse Gamma distributions
#'  and relies on `stats::qgamma`.
#'
#' @keywords internal
#' @noRd
#'
qinvgamma_internal <- function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
  # Check shape
  if (length(shape) != 1 || !is.numeric(shape) || shape <= 0 || !is.finite(shape)) {
    stop("'shape' must be a single positive finite number.")
  }

  # Handle rate vs scale parameters, ensuring only one effective parameter is used
  # and it's positive and finite.
  rate_missing <- missing(rate)
  scale_missing <- missing(scale)

  if (!rate_missing && !scale_missing) {
    # Both provided, check consistency implicitly via default scale=1/rate
    # Need explicit check if defaults were different or if strict error desired
    # Using rate's value if both provided (matches rinvgamma default behaviour)
    if (!is.numeric(rate) || length(rate) != 1 || rate <= 0 || !is.finite(rate)) {
      stop("'rate' must be a single positive finite number.")
    }
    actual_rate <- rate
  } else if (rate_missing && !scale_missing) {
    # Only scale provided
    if (!is.numeric(scale) || length(scale) != 1 || scale <= 0 || !is.finite(scale)) {
      stop("'scale' must be a single positive finite number.")
    }
    actual_rate <- 1 / scale
  } else {
    # Use rate (either provided or default=1)
    if (!is.numeric(rate) || length(rate) != 1 || rate <= 0 || !is.finite(rate)) {
      stop("'rate' must be a single positive finite number.")
    }
    actual_rate <- rate
  }

  # Check p
  if (!is.numeric(p)) {
    stop("'p' must be numeric.")
  }

  # Vectorised calculation using stats::qgamma
  # The quantile x for InvGamma(a, b) such that P(X <= x) = p
  # corresponds to 1 / y where y is the quantile for Gamma(a, b)
  # such that P(Y >= y) = p.
  # This means y = qgamma(p, shape, rate, lower.tail = FALSE)
  # So, x = 1 / qgamma(p, shape, rate, lower.tail = FALSE)

  # Need to handle log.p and lower.tail for the input p *before* passing to qgamma
  # stats::qgamma handles log.p and lower.tail internally, but we need the
  # *final* target probability p (corresponding to P(X <= x)) to apply the formula.

  # Let qgamma handle the log.p and lower.tail transformations on the probability
  # Argument to qgamma needs to be the probability p for the *Gamma* distribution's
  # upper tail (P(Y >= y)) which is numerically the same as the probability p
  # for the *Inverse Gamma* distribution's lower tail (P(X <= x)).
  # So we pass p directly to qgamma, but specify lower.tail=FALSE

  # Handle edge cases for p before calling qgamma which might warn or error
  # qgamma handles p=0/p=1 (or log equivalents) correctly for lower.tail=FALSE
  # qgamma(0, ..., lower.tail=FALSE) = Inf => 1/Inf = 0
  # qgamma(1, ..., lower.tail=FALSE) = 0   => 1/0 = Inf

  # Catch issues where p might be outside valid range after log.p/lower.tail logic
  # This is implicitly handled by qgamma which returns NaN for invalid p

  # Note: We must pass log.p to qgamma as well!
  quantiles_gamma_upper <- stats::qgamma(p,
                                         shape = shape,
                                         rate = actual_rate,
                                         lower.tail = FALSE, # Key part: P(Y >= y) = p
                                         log.p = log.p)

  # Inverse relationship
  quantiles_invgamma <- 1 / quantiles_gamma_upper

  # Explicitly set InvGamma quantile for p=0 (lower tail) to 0
  # and p=1 (lower tail) to Inf, respecting log.p and lower.tail flags.

  # Determine the effective probability values corresponding to P(X <= x)
  if (log.p) {
    if (any(p > 0 & !is.na(p))) {
      warning("Negative values expected for p when log.p = TRUE")
      # qgamma will likely return NaN for p > 0 anyway
    }
    prob_eff <- exp(p)
  } else {
    if (any((p < 0 | p > 1) & !is.na(p) )) {
      warning("Probabilities p are expected to be in [0, 1]")
      # qgamma will likely return NaN anyway
    }
    prob_eff <- p
  }

  if (!lower.tail) {
    prob_eff <- 1 - prob_eff # Convert P(X > x) probability to P(X <= x) probability
  }

  # Apply edge case corrections based on effective probability P(X <= x)
  quantiles_invgamma[prob_eff == 0] <- 0
  quantiles_invgamma[prob_eff == 1] <- Inf

  return(quantiles_invgamma)
}
