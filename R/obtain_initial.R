#' Obtain initial beta coefficients using MVMR.dIVW with ridge penalty.
#'
#' @param beta.exposure A data.frame or matrix of SNP-exposure association estimates.
#' @param se.exposure A data.frame or matrix of standard errors for beta.exposure.
#' @param beta.outcome A numeric vector of SNP-outcome association estimates.
#' @param se.outcome A numeric vector of standard errors for beta.outcome.
#' @param iv_strength_parameter Sample IV strength parameter.
#' @param P A K-by-K exposure correlation matrix when overlap = FALSE, or a (K+1)-by-(K+1) joint exposure-outcome correlation matrix when overlap = TRUE.
#' @param overlap Logical. Whether exposure and outcome GWAS summary statistics overlap. Default = FALSE.
#' @param CV_fold Number of data-thinning cross-validation folds. Default = 5.
#' @param n_times Number of repeated cross-validation runs. Default = 1.
#' @param epsilon Tolerance for ADMM projection. Default = 1e-4.
#' @param seed Optional random seed.
#' @param lambda.1se Whether to use the one-standard-error rule. Default = TRUE.
#' @param lambda.length Number of candidate ridge penalties. Default = 20.
#'
#' @return A list with elements beta, phi, phi.index, phi.grid,
#'   cv.loss.mean, cv.loss.se, cv.valid.n, and cv.loss.
#' @export
obtain_initial <- function(
    beta.exposure,
    se.exposure,
    beta.outcome,
    se.outcome,
    iv_strength_parameter,
    P,
    overlap = FALSE,
    CV_fold = 5,
    n_times = 1,
    epsilon = 1e-4,
    seed = NULL,
    lambda.1se = TRUE,
    lambda.length = 20
) {
  dat <- validate_initial_inputs(
    beta.exposure, se.exposure, beta.outcome, se.outcome,
    iv_strength_parameter, P, overlap, CV_fold, n_times, epsilon,
    lambda.length
  )
  beta.exposure <- dat$beta.exposure
  se.exposure <- dat$se.exposure
  beta.outcome <- dat$beta.outcome
  se.outcome <- dat$se.outcome
  P <- dat$P
  P.exposure <- dat$P.exposure
  rho.xy <- dat$rho.xy
  p <- dat$p
  K <- dat$K

  e <- 1 / CV_fold
  phi_cand <- make_phi_grid(
    iv_strength_parameter = iv_strength_parameter,
    p = p,
    length.out = lambda.length
  )

  cv_loss <- matrix(
    NA_real_, nrow = CV_fold * n_times, ncol = length(phi_cand)
  )

  if (isTRUE(overlap)) {
    Sig.joint <- make_joint_thinning_covariance(se.exposure, se.outcome, P)
    joint.summary <- cbind(beta.exposure, beta.outcome)
  } else {
    Sig <- make_thinning_covariance(se.exposure, P.exposure)
  }

  for (m in seq_len(n_times)) {
    if (!is.null(seed)) set.seed(seed + m)

    if (isTRUE(overlap)) {
      dt.joint <- datathin::datathin(
        joint.summary, family = "mvnormal", arg = Sig.joint, K = CV_fold
      )
      dt.exposure <- dt.joint[, seq_len(K), , drop = FALSE]
      dt.outcome <- dt.joint[, K + 1L, , drop = FALSE]
    } else {
      dt.exposure <- datathin::datathin(
        beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold
      )
      dt.outcome <- datathin::datathin(
        beta.outcome, family = "normal", arg = se.outcome^2, K = CV_fold
      )
    }

    for (l in seq_len(CV_fold)) {
      train_idx <- setdiff(seq_len(CV_fold), l)
      row_id <- l + (m - 1L) * CV_fold

      X.train <- apply(
        dt.exposure[, , train_idx, drop = FALSE], c(1, 2), sum
      )
      y.train <- apply(
        dt.outcome[, 1, train_idx, drop = FALSE], 1, sum
      )
      seX.train <- se.exposure * sqrt(1 - e)
      sey.train <- se.outcome * sqrt(1 - e)

      X.test <- dt.exposure[, , l, drop = FALSE][, , 1]
      y.test <- dt.outcome[, 1, l]
      seX.test <- se.exposure * sqrt(e)
      sey.test <- se.outcome * sqrt(e)

      train_cache <- tryCatch(
        prepare_mvmr_ridge_system(
          beta.exposure = X.train,
          se.exposure = seX.train,
          se.outcome = sey.train,
          P = P.exposure,
          epsilon = epsilon
        ),
        error = function(e) NULL
      )
      test_cache <- tryCatch(
        prepare_mvmr_ridge_system(
          beta.exposure = X.test,
          se.exposure = seX.test,
          se.outcome = sey.test,
          P = P.exposure,
          epsilon = epsilon
        ),
        error = function(e) NULL
      )

      if (is.null(train_cache) || is.null(test_cache)) next

      B.test <- compute_overlap_bias(seX.test, sey.test, rho.xy)
      B.train <- compute_overlap_bias(seX.train, sey.train, rho.xy)
      c.test <- as.numeric(test_cache$XtW %*% y.test) - B.test

      for (jphi in seq_along(phi_cand)) {
        fit <- tryCatch({
          ridge_system <- factorize_mvmr_ridge_system(
            train_cache$MV.plus, phi_cand[jphi]
          )
          solve_mvmr_ridge_system(
            ridge_system,
            as.numeric(train_cache$XtW %*% y.train) - B.train
          )
        }, error = function(e) NULL)

        if (is.null(fit) || any(!is.finite(fit))) next

        loss <- drop(crossprod(fit, test_cache$MV.plus %*% fit)) -
          2 * drop(crossprod(c.test, fit))
        if (is.finite(loss)) cv_loss[row_id, jphi] <- loss
      }
    }
  }

  cv_summary <- summarize_initial_cv(cv_loss, margin = 2)
  valid_full <- cv_summary$valid.n == CV_fold * n_times &
    is.finite(cv_summary$mean)

  if (!any(valid_full)) {
    stop("No ridge tuning value had valid fits across all CV folds.")
  }

  loss_valid <- cv_summary$mean
  loss_valid[!valid_full] <- Inf
  min_idx <- which.min(loss_valid)
  selected_idx <- min_idx

  if (isTRUE(lambda.1se) && is.finite(cv_summary$se[min_idx])) {
    cutoff <- cv_summary$mean[min_idx] + cv_summary$se[min_idx]
    eligible <- which(
      valid_full & is.finite(cv_summary$mean) &
        cv_summary$mean <= cutoff
    )
    if (length(eligible) > 0L) {
      selected_idx <- eligible[which.max(phi_cand[eligible])]
    }
  }

  final_cache <- prepare_mvmr_ridge_system(
    beta.exposure, se.exposure, se.outcome, P.exposure, epsilon
  )
  final_system <- factorize_mvmr_ridge_system(
    final_cache$MV.plus, phi_cand[selected_idx]
  )
  beta_hat <- solve_mvmr_ridge_system(
    final_system,
    as.numeric(final_cache$XtW %*% beta.outcome) -
      compute_overlap_bias(se.exposure, se.outcome, rho.xy)
  )

  list(
    beta = as.numeric(beta_hat),
    phi = phi_cand[selected_idx],
    phi.index = selected_idx,
    phi.grid = phi_cand,
    cv.loss.mean = cv_summary$mean,
    cv.loss.se = cv_summary$se,
    cv.valid.n = cv_summary$valid.n,
    cv.loss = cv_loss
  )
}

#' Obtain pleiotropy-robust initial beta coefficients.
#'
#' @inheritParams obtain_initial
#' @param lambda.alpha Optional candidate penalties for SNP-specific pleiotropy.
#' @param lambda.alpha.length Length of the lambda.alpha grid. Default = 20.
#' @param lambda.alpha.min.ratio Smallest-to-largest lambda.alpha ratio.
#' @param max.iter Maximum alternating updates. Default = 200.
#' @param err Joint absolute-relative convergence tolerance. Default = 1e-4.
#'
#' @return A list with beta, alpha, invalid.snp, selected tuning parameters,
#'   tuning grids, convergence information, and CV diagnostics.
#' @export
obtain_initial_pleiotropy <- function(
    beta.exposure,
    se.exposure,
    beta.outcome,
    se.outcome,
    iv_strength_parameter,
    P,
    overlap = FALSE,
    CV_fold = 5,
    n_times = 1,
    epsilon = 1e-4,
    seed = NULL,
    lambda.1se = TRUE,
    lambda.length = 20,
    lambda.alpha = NULL,
    lambda.alpha.length = 20,
    lambda.alpha.min.ratio = 0.05,
    max.iter = 200,
    err = 1e-4
) {
  dat <- validate_initial_inputs(
    beta.exposure, se.exposure, beta.outcome, se.outcome,
    iv_strength_parameter, P, overlap, CV_fold, n_times, epsilon,
    lambda.length
  )
  beta.exposure <- dat$beta.exposure
  se.exposure <- dat$se.exposure
  beta.outcome <- dat$beta.outcome
  se.outcome <- dat$se.outcome
  P <- dat$P
  P.exposure <- dat$P.exposure
  rho.xy <- dat$rho.xy
  p <- dat$p
  K <- dat$K

  if (length(lambda.alpha.length) != 1L || !is.finite(lambda.alpha.length) ||
      lambda.alpha.length < 1L) {
    stop("lambda.alpha.length must be a positive integer.")
  }
  if (length(lambda.alpha.min.ratio) != 1L ||
      !is.finite(lambda.alpha.min.ratio) ||
      lambda.alpha.min.ratio <= 0 || lambda.alpha.min.ratio > 1) {
    stop("lambda.alpha.min.ratio must be in (0, 1].")
  }
  if (length(max.iter) != 1L || !is.finite(max.iter) || max.iter < 1L) {
    stop("max.iter must be a positive integer.")
  }
  if (length(err) != 1L || !is.finite(err) || err <= 0) {
    stop("err must be positive.")
  }

  e <- 1 / CV_fold
  phi_cand <- make_phi_grid(iv_strength_parameter, p, lambda.length)

  pilot_cache <- prepare_mvmr_ridge_system(
    beta.exposure, se.exposure, se.outcome, P.exposure, epsilon
  )
  pilot_system <- factorize_mvmr_ridge_system(
    pilot_cache$MV.plus, max(phi_cand)
  )
  pilot_beta <- solve_mvmr_ridge_system(
    pilot_system,
    as.numeric(pilot_cache$XtW %*% beta.outcome) -
      compute_overlap_bias(se.exposure, se.outcome, rho.xy)
  )

  if (is.null(lambda.alpha)) {
    pilot_residual <- as.numeric(beta.outcome - beta.exposure %*% pilot_beta)
    lambda.alpha.max <- max(abs(pilot_residual) / se.outcome^2)
    if (!is.finite(lambda.alpha.max) || lambda.alpha.max <= 0) {
      stop("Cannot construct a valid lambda.alpha grid.")
    }
    lambda.alpha <- lambda.alpha.max * exp(seq(
      log(1), log(lambda.alpha.min.ratio),
      length.out = as.integer(lambda.alpha.length)
    ))
  } else {
    lambda.alpha <- sort(unique(as.numeric(lambda.alpha)), decreasing = TRUE)
    if (length(lambda.alpha) == 0L || any(!is.finite(lambda.alpha)) ||
        any(lambda.alpha < 0)) {
      stop("lambda.alpha must contain finite nonnegative values.")
    }
  }

  if (isTRUE(overlap)) {
    Sig.joint <- make_joint_thinning_covariance(
      se.exposure, se.outcome, P
    )
    joint.summary <- cbind(beta.exposure, beta.outcome)
  } else {
    Sig <- make_thinning_covariance(se.exposure, P.exposure)
  }

  cv_loss <- array(
    NA_real_,
    dim = c(CV_fold * n_times, length(phi_cand), length(lambda.alpha))
  )

  for (m in seq_len(n_times)) {
    if (!is.null(seed)) set.seed(seed + m)

    if (isTRUE(overlap)) {
      dt.joint <- datathin::datathin(
        joint.summary,
        family = "mvnormal",
        arg = Sig.joint,
        K = CV_fold
      )
      dt.exposure <- dt.joint[, seq_len(K), , drop = FALSE]
      dt.outcome <- dt.joint[, K + 1L, , drop = FALSE]
    } else {
      dt.exposure <- datathin::datathin(
        beta.exposure,
        family = "mvnormal",
        arg = Sig,
        K = CV_fold
      )
      dt.outcome <- datathin::datathin(
        beta.outcome,
        family = "normal",
        arg = se.outcome^2,
        K = CV_fold
      )
    }

    for (l in seq_len(CV_fold)) {
      train_idx <- setdiff(seq_len(CV_fold), l)
      row_id <- l + (m - 1L) * CV_fold

      X.train <- apply(
        dt.exposure[, , train_idx, drop = FALSE], c(1, 2), sum
      )
      y.train <- apply(
        dt.outcome[, 1, train_idx, drop = FALSE], 1, sum
      )
      seX.train <- se.exposure * sqrt(1 - e)
      sey.train <- se.outcome * sqrt(1 - e)

      X.test <- dt.exposure[, , l, drop = FALSE][, , 1]
      y.test <- dt.outcome[, 1, l]
      seX.test <- se.exposure * sqrt(e)
      sey.test <- se.outcome * sqrt(e)

      train_cache <- tryCatch(
        prepare_mvmr_ridge_system(X.train, seX.train, sey.train, P.exposure, epsilon),
        error = function(e) NULL
      )
      test_cache <- tryCatch(
        prepare_mvmr_ridge_system(X.test, seX.test, sey.test, P.exposure, epsilon),
        error = function(e) NULL
      )
      if (is.null(train_cache) || is.null(test_cache)) next

      B.train <- compute_overlap_bias(seX.train, sey.train, rho.xy)
      B.test <- compute_overlap_bias(seX.test, sey.test, rho.xy)
      wtest <- 1 / sey.test^2

      for (jphi in seq_along(phi_cand)) {
        ridge_system <- tryCatch(
          factorize_mvmr_ridge_system(train_cache$MV.plus, phi_cand[jphi]),
          error = function(e) NULL
        )
        if (is.null(ridge_system)) next

        beta.warm <- tryCatch(
          solve_mvmr_ridge_system(
            ridge_system,
            as.numeric(train_cache$XtW %*% y.train) - B.train
          ),
          error = function(e) NULL
        )
        if (is.null(beta.warm) || any(!is.finite(beta.warm))) next
        alpha.warm <- rep(0, p)

        for (ja in seq_along(lambda.alpha)) {
          fit <- tryCatch(
            mvmr_ridge_pleiotropy(
              beta.exposure = X.train,
              beta.outcome = y.train,
              se.outcome = sey.train,
              lambda.alpha = lambda.alpha[ja],
              max.iter = max.iter,
              err = err,
              beta.init = beta.warm,
              alpha.init = alpha.warm,
              XtW = train_cache$XtW,
              ridge.system = ridge_system,
              overlap.bias = B.train
            ),
            error = function(e) NULL
          )
          if (is.null(fit) || any(!is.finite(fit$beta)) ||
              any(!is.finite(fit$alpha))) next

          beta.warm <- fit$beta
          alpha.warm <- fit$alpha

          alpha.test <- fit$alpha * e / (1 - e)
          y.adjusted <- y.test - alpha.test
          c.test <- as.numeric(test_cache$XtW %*% y.adjusted) - B.test
          loss <- drop(crossprod(fit$beta, test_cache$MV.plus %*% fit$beta)) -
            2 * drop(crossprod(c.test, fit$beta)) +
            sum(wtest * y.adjusted^2)

          if (is.finite(loss)) cv_loss[row_id, jphi, ja] <- loss
        }
      }
    }
  }

  cv_summary <- summarize_initial_cv(cv_loss, margin = c(2, 3))
  valid_full <- cv_summary$valid.n == CV_fold * n_times &
    is.finite(cv_summary$mean)
  if (!any(valid_full)) {
    stop("No pleiotropy-aware tuning combination had valid fits across all CV folds.")
  }

  loss_valid <- cv_summary$mean
  loss_valid[!valid_full] <- Inf
  best_idx <- which(loss_valid == min(loss_valid), arr.ind = TRUE)[1, ]

  if (isTRUE(lambda.1se) &&
      is.finite(cv_summary$se[best_idx[1], best_idx[2]])) {
    cutoff <- cv_summary$mean[best_idx[1], best_idx[2]] +
      cv_summary$se[best_idx[1], best_idx[2]]
    eligible <- which(
      valid_full & is.finite(cv_summary$mean) &
        cv_summary$mean <= cutoff,
      arr.ind = TRUE
    )
    if (nrow(eligible) > 0L) {
      max_alpha <- max(lambda.alpha[eligible[, 2]])
      eligible <- eligible[
        lambda.alpha[eligible[, 2]] == max_alpha, , drop = FALSE
      ]
      max_phi <- max(phi_cand[eligible[, 1]])
      eligible <- eligible[
        phi_cand[eligible[, 1]] == max_phi, , drop = FALSE
      ]
      best_idx <- eligible[1, ]
    }
  }

  final_cache <- prepare_mvmr_ridge_system(
    beta.exposure, se.exposure, se.outcome, P.exposure, epsilon
  )
  final_system <- factorize_mvmr_ridge_system(
    final_cache$MV.plus, phi_cand[best_idx[1]]
  )
  final_beta_init <- solve_mvmr_ridge_system(
    final_system,
    as.numeric(final_cache$XtW %*% beta.outcome) -
      compute_overlap_bias(se.exposure, se.outcome, rho.xy)
  )
  final_fit <- mvmr_ridge_pleiotropy(
    beta.exposure = beta.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    lambda.alpha = lambda.alpha[best_idx[2]],
    max.iter = max.iter,
    err = err,
    beta.init = final_beta_init,
    alpha.init = rep(0, p),
    XtW = final_cache$XtW,
    ridge.system = final_system,
    overlap.bias = compute_overlap_bias(se.exposure, se.outcome, rho.xy)
  )

  invalid.snp <- which(abs(final_fit$alpha) > 1e-7)

  list(
    beta = as.numeric(final_fit$beta),
    alpha = as.numeric(final_fit$alpha),
    invalid.snp = invalid.snp,
    phi = phi_cand[best_idx[1]],
    phi.index = best_idx[1],
    phi.grid = phi_cand,
    lambda.alpha = lambda.alpha[best_idx[2]],
    lambda.alpha.index = best_idx[2],
    lambda.alpha.grid = lambda.alpha,
    iter = final_fit$iter,
    converged = final_fit$converged,
    cv.loss.mean = cv_summary$mean,
    cv.loss.se = cv_summary$se,
    cv.valid.n = cv_summary$valid.n,
    cv.loss = cv_loss
  )
}

validate_initial_inputs <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome,
    iv_strength_parameter, P, overlap, CV_fold, n_times, epsilon,
    lambda.length
) {
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  beta.outcome <- as.numeric(beta.outcome)
  se.outcome <- as.numeric(se.outcome)
  P <- as.matrix(P)

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)
  if (p < 1L || K < 1L) stop("beta.exposure must have positive dimensions.")
  if (!identical(dim(se.exposure), c(p, K))) {
    stop("se.exposure must have the same dimensions as beta.exposure.")
  }
  if (length(beta.outcome) != p || length(se.outcome) != p) {
    stop("beta.outcome and se.outcome must have one value per SNP.")
  }
  overlap.info <- parse_overlap_correlation(P, K, overlap)
  if (any(!is.finite(beta.exposure)) || any(!is.finite(se.exposure)) ||
      any(!is.finite(beta.outcome)) || any(!is.finite(se.outcome)) ||
      any(!is.finite(P))) {
    stop("All association estimates, standard errors, and P must be finite.")
  }
  if (any(se.exposure <= 0) || any(se.outcome <= 0)) {
    stop("All standard errors must be positive.")
  }
  if (!isTRUE(all.equal(P, t(P), tolerance = 1e-8))) {
    stop("P must be symmetric.")
  }
  if (min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("P must be positive semidefinite.")
  }
  if (length(iv_strength_parameter) != 1L ||
      !is.finite(iv_strength_parameter)) {
    stop("iv_strength_parameter must be a finite scalar.")
  }
  if (length(CV_fold) != 1L || !is.finite(CV_fold) || CV_fold < 2L) {
    stop("CV_fold must be at least 2.")
  }
  if (length(n_times) != 1L || !is.finite(n_times) || n_times < 1L) {
    stop("n_times must be at least 1.")
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be positive.")
  }
  if (length(lambda.length) != 1L || !is.finite(lambda.length) ||
      lambda.length < 2L) {
    stop("lambda.length must be at least 2.")
  }

  list(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    P = P,
    P.exposure = overlap.info$P.exposure,
    rho.xy = overlap.info$rho.xy,
    p = p,
    K = K
  )
}

make_phi_grid <- function(iv_strength_parameter, p, length.out) {
  upper <- (max(iv_strength_parameter, 0) * sqrt(p) + p)^(2 / 5)
  if (!is.finite(upper) || upper <= 0) {
    stop("Failed to construct a valid ridge tuning grid.")
  }
  exp(seq(
    log(upper * 1e2), log(upper * 1e-4),
    length.out = as.integer(length.out)
  ))
}

make_thinning_covariance <- function(se.exposure, P) {
  p <- nrow(se.exposure)
  K <- ncol(se.exposure)
  Sig <- array(0, dim = c(p, K, K))
  for (j in seq_len(p)) {
    Sig[j, , ] <- P * tcrossprod(se.exposure[j, ])
  }
  Sig
}

summarize_initial_cv <- function(cv_loss, margin) {
  loss_mean <- apply(cv_loss, margin, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  })
  loss_se <- apply(cv_loss, margin, function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) Inf else stats::sd(x) / sqrt(length(x))
  })
  valid_n <- apply(cv_loss, margin, function(x) sum(is.finite(x)))
  list(mean = loss_mean, se = loss_se, valid.n = valid_n)
}

prepare_mvmr_ridge_system <- function(
    beta.exposure, se.exposure, se.outcome, P, epsilon = 1e-4
) {
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  se.outcome <- as.numeric(se.outcome)

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)
  w <- 1 / se.outcome^2
  XtW <- t(beta.exposure) * rep(w, each = K)
  M <- XtW %*% beta.exposure

  weighted_se <- se.exposure * sqrt(w)
  V <- P * crossprod(weighted_se)
  MV <- M - V
  MV.plus <- BDcocolasso::ADMM_proj(
    MV / sqrt(p), epsilon = epsilon
  )$mat * sqrt(p)

  list(XtW = XtW, MV.plus = MV.plus)
}

factorize_mvmr_ridge_system <- function(MV.plus, phi) {
  if (length(phi) != 1L || !is.finite(phi) || phi < 0) {
    stop("phi must be a finite nonnegative scalar.")
  }
  K <- nrow(MV.plus)
  A <- as.matrix(MV.plus) + phi * diag(K)
  chol.A <- tryCatch(chol(A), error = function(e) NULL)
  list(A = A, chol = chol.A)
}

solve_mvmr_ridge_system <- function(ridge.system, rhs) {
  rhs <- as.numeric(rhs)
  if (!is.null(ridge.system$chol)) {
    return(as.numeric(backsolve(
      ridge.system$chol,
      forwardsolve(t(ridge.system$chol), rhs)
    )))
  }
  out <- tryCatch(
    solve(ridge.system$A, rhs),
    error = function(e) MASS::ginv(ridge.system$A) %*% rhs
  )
  as.numeric(out)
}

mvmr_ridge_pleiotropy <- function(
    beta.exposure,
    beta.outcome,
    se.outcome,
    lambda.alpha,
    max.iter = 200,
    err = 1e-4,
    beta.init = NULL,
    alpha.init = NULL,
    XtW,
    ridge.system,
    overlap.bias = NULL
) {
  beta.exposure <- as.matrix(beta.exposure)
  beta.outcome <- as.numeric(beta.outcome)
  se.outcome <- as.numeric(se.outcome)
  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)
  overlap.bias <- if (is.null(overlap.bias)) rep(0, K) else as.numeric(overlap.bias)

  beta.current <- if (is.null(beta.init)) {
    solve_mvmr_ridge_system(
      ridge.system, as.numeric(XtW %*% beta.outcome) - overlap.bias
    )
  } else {
    as.numeric(beta.init)
  }
  alpha.current <- if (is.null(alpha.init)) rep(0, p) else as.numeric(alpha.init)

  if (length(beta.current) != K || length(alpha.current) != p) {
    stop("Invalid beta.init or alpha.init dimension.")
  }

  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(max.iter)) {
    rhs <- as.numeric(XtW %*% (beta.outcome - alpha.current)) - overlap.bias
    beta.new <- solve_mvmr_ridge_system(ridge.system, rhs)

    residual <- as.numeric(beta.outcome - beta.exposure %*% beta.new)
    threshold <- lambda.alpha * se.outcome^2
    alpha.new <- sign(residual) * pmax(abs(residual) - threshold, 0)

    max.change <- max(
      max(abs(beta.new - beta.current)),
      max(abs(alpha.new - alpha.current))
    )
    parameter.scale <- max(
      1,
      max(abs(beta.current)),
      max(abs(alpha.current))
    )

    beta.current <- beta.new
    alpha.current <- alpha.new

    if (max.change <= err * parameter.scale) {
      converged <- TRUE
      break
    }
  }

  list(
    beta = as.numeric(beta.current),
    alpha = as.numeric(alpha.current),
    iter = iter,
    converged = converged
  )
}
