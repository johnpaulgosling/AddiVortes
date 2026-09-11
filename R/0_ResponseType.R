#' @title Infer the modelling task from a response vector
#'
#' @description
#' Detects whether `y` should be treated as a regression, binary classification
#' or multinomial classification response.
#'
#' Factor, character and logical vectors are treated as classification (two
#' levels: binary; three or more: multinomial). Numeric `y` with exactly two
#' unique values in `{0, 1}` is binary classification. Any other numeric `y`
#' is regression.
#'
#' @param y A response vector.
#'
#' @return A list with `task`, `classLevels`, `yClass` (0-based integer codes,
#'   or `NULL` for regression) and `nLatents`.
#'
#' @keywords internal
#' @noRd
infer_response_type_internal <- function(y) {
  if (is.null(y)) {
    stop("'y' must be provided.", call. = FALSE)
  }
  if (length(y) < 1L) {
    stop("'y' must have length greater than 0.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("'y' must not contain missing values.", call. = FALSE)
  }

  if (is.factor(y) || is.character(y) || is.logical(y)) {
    if (is.character(y)) {
      y_fac <- factor(y)
    } else if (is.logical(y)) {
      y_fac <- factor(y, levels = c(FALSE, TRUE))
    } else {
      y_fac <- droplevels(y)
    }
    levels_y <- levels(y_fac)
    if (length(levels_y) < 2L) {
      stop("Classification responses must have at least two classes.",
           call. = FALSE)
    }
    codes <- as.integer(y_fac) - 1L
    task <- if (length(levels_y) == 2L) "binary" else "multinomial"
    return(list(
      task = task,
      classLevels = levels_y,
      yClass = codes,
      nLatents = length(levels_y) - 1L
    ))
  }

  if (is.numeric(y)) {
    unique_y <- unique(y)
    if (length(unique_y) == 2L && all(unique_y %in% c(0, 1))) {
      return(list(
        task = "binary",
        classLevels = c("0", "1"),
        yClass = as.integer(y == 1),
        nLatents = 1L
      ))
    }
    return(list(
      task = "regression",
      classLevels = NULL,
      yClass = NULL,
      nLatents = 1L
    ))
  }

  stop("'y' must be numeric, factor, character, or logical.", call. = FALSE)
}

#' @title Independent-probit class probabilities from latent means
#'
#' @description
#' For a list of `nLatents` matrices of latent means `G_d` (each `n` by
#' `nDraws`), estimate class probabilities under independent N(G, I) latents.
#' The reference class (index 1) is chosen when every latent is negative;
#' otherwise the class is `which.max(z) + 1`.
#'
#' @param G_list A list of numeric matrices, one per latent dimension.
#' @param nMC Number of Monte Carlo draws of `z ~ N(G, I)` per posterior draw.
#' @param perDraw If `TRUE`, also return an `n x K x nDraws` array of
#'   per-draw probabilities.
#'
#' @return A list with `mean` (`n x K` matrix) and optionally `draw`.
#'
#' @keywords internal
#' @noRd
multinomialProbabilities_internal <- function(G_list, nMC = 32L,
                                              perDraw = FALSE) {
  n <- nrow(G_list[[1]])
  n_draws <- ncol(G_list[[1]])
  n_latents <- length(G_list)
  k_classes <- n_latents + 1L
  counts <- matrix(0, n, k_classes)
  draw_probs <- if (perDraw) {
    array(0, dim = c(n, k_classes, n_draws))
  } else {
    NULL
  }

  row_idx <- seq_len(n)
  for (s in seq_len(n_draws)) {
    g <- matrix(0, n, n_latents)
    for (d in seq_len(n_latents)) {
      g[, d] <- G_list[[d]][, s]
    }
    class_counts <- matrix(0, n, k_classes)
    for (ignored in seq_len(nMC)) {
      z <- g + matrix(stats::rnorm(n * n_latents), n, n_latents)
      max_z <- z[, 1]
      if (n_latents > 1L) {
        for (d in 2:n_latents) {
          max_z <- pmax(max_z, z[, d])
        }
      }
      winner <- max.col(z, ties.method = "first")
      cls <- ifelse(max_z < 0, 1L, winner + 1L)
      class_counts[cbind(row_idx, cls)] <- class_counts[cbind(row_idx, cls)] + 1
    }
    if (perDraw) {
      draw_probs[, , s] <- class_counts / nMC
    }
    counts <- counts + class_counts
  }

  list(mean = counts / (n_draws * nMC), draw = draw_probs)
}

#' @title Posterior latent-link matrices for each ensemble
#'
#' @description
#' Evaluates `G_d(x)` for each latent ensemble by calling the compiled
#' predictor on the corresponding subset of tessellations.
#'
#' @keywords internal
#' @noRd
latentLinkMatrices_internal <- function(object, x_scaled, showProgress = FALSE) {
  n_latents <- object$nLatents
  if (is.null(n_latents) || n_latents < 1L) {
    n_latents <- 1L
  }
  m_per <- object$mPerLatent
  if (is.null(m_per) || is.na(m_per)) {
    n_tess <- if (length(object$posteriorTess) > 0) {
      length(object$posteriorTess[[1]])
    } else {
      0L
    }
    m_per <- as.integer(n_tess / n_latents)
  }

  G_list <- vector("list", n_latents)
  for (d in seq_len(n_latents)) {
    idx <- ((d - 1L) * m_per + 1L):(d * m_per)
    tess_d <- lapply(object$posteriorTess, function(draw) draw[idx])
    dim_d <- lapply(object$posteriorDim, function(draw) draw[idx])
    pred_d <- lapply(object$posteriorPred, function(draw) draw[idx])
    G_list[[d]] <- .Call(
      "addi_vortes_predict_cpp",
      x_scaled,
      tess_d,
      dim_d,
      pred_d,
      as.integer(object$metric_red),
      as.integer(object$member_red),
      as.logical(showProgress && d == n_latents)
    )
  }
  G_list
}

#' @title Whether an AddiVortes object is a classification fit
#'
#' @keywords internal
#' @noRd
isClassification_internal <- function(object) {
  !is.null(object$task) && object$task %in% c("binary", "multinomial")
}
