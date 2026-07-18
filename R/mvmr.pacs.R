#' Function to perform the MVMR-PACS estimator for signal-group identification and variable selection in two-sample multivariable Mendelian randomization with highly correlated exposures.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effects of a SNP on the K exposures, usually obtained from a GWAS.
#' @param se.exposure A data.frame or matrix of estimated standard errors corresponding to beta.exposure.
#' @param beta.outcome A numeric vector of the estimated marginal effects of SNPs on the outcome, usually obtained from a GWAS.
#' @param se.outcome A numeric vector of estimated standard errors corresponding to beta.outcome.
#' @param P A K-by-K matrix for the estimated shared correlation matrix among the SNP-exposure association estimates, where K is the number of exposures.
#' @param type Weighting scheme in penalization: 1 = weights independent of the correlation between beta.exposure; 2 = weights inversely related to the correlation between beta.exposure; 3 = weights nonzero only if the absolute correlation between beta.exposure is greater than rr_cut_off; 4 = weights inversely related to the correlation between beta.exposure only if the absolute correlation between beta.exposure is greater than rr_cut_off.
#' @param rr_cut_off Pre-specified correlation cutoff used in the weighting schemes specified by type. Default = 0.
#' @param fix.cor Logical. Whether to use the full-data correlation matrix of beta.exposure in cross-validation instead of recomputing the correlation matrix within each training fold. Default = FALSE.
#' @param betawt Numeric vector of initial beta coefficients. Default = NULL, in which case an initial estimator is obtained internally.
#' @param tau Numeric vector of candidate values for the adaptive weighting parameter. Default = c(0.5, 1, 2, 3).
#' @param CV_fold Number of folds for cross-validation based on multi-fold data thinning. Default = 5.
#' @param n_times Number of repeated cross-validation runs. Default = 1.
#' @param lambda.length Length of the candidate tuning parameter sequence for the PACS penalty. Default = 20.
#' @param eps Tolerance level for ADMM projection. Default = 1e-4.
#' @param digit Number of digits to keep when determining signal-groups. Default = 3, corresponding to decimal precision 0.001.
#' @param seed Seed for reproducibility. Default = NULL.
#' @param lambda.1se Logical. Whether to use the one-standard-error rule when selecting tuning parameters. Default = TRUE.
#' @param pleiotropy Logical. Whether to allow and select SNP-specific pleiotropic effects using an additional l1 penalty. Default = TRUE.
#' @param lambda.alpha Optional numeric vector of candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects. Default = NULL, in which case a data-dependent decreasing sequence is constructed.
#' @param lambda.alpha.length Length of the candidate tuning parameter sequence for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Default = 20.
#' @param lambda.alpha.min.ratio Ratio between the smallest and largest candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Default = 0.05.
#'
#' @return A list with the following elements:
#' \item{beta}{Estimated direct effects of each exposure on the outcome. Exposures with estimated effects of the same magnitude, up to the specified decimal precision, belong to the same signal-group.}
#' \item{alpha}{Estimated SNP-specific pleiotropic effects. If pleiotropy = FALSE, this is a vector of zeros.}
#' \item{invalid.snp}{Indices of SNPs with nonzero estimated pleiotropic effects.}
#' \item{signal_group}{Identified signal-groups.}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies instrument strength in the sample.}
#' \item{lambda.index}{Index of the selected PACS penalty tuning parameter.}
#' \item{lambda}{Selected PACS penalty tuning parameter.}
#' \item{tau.index}{Index of the selected adaptive weighting parameter.}
#' \item{tau}{Selected adaptive weighting parameter.}
#' \item{lambda.alpha.index}{Index of the selected tuning parameter for the l1 penalty on SNP-specific pleiotropic effects.}
#' \item{lambda.alpha}{Selected tuning parameter for the l1 penalty on SNP-specific pleiotropic effects.}
#' \item{lambda.alpha.grid}{Candidate tuning parameter sequence used for the l1 penalty on SNP-specific pleiotropic effects.}
#' \item{pleiotropy}{Logical indicator of whether pleiotropic effects were modeled.}
#' \item{beta.init}{Initial beta estimates used to construct adaptive weights.}
#' \item{cv.loss.mean}{Array of mean cross-validation losses over repeated data-thinning folds.}
#' \item{cv.loss.sd}{Array of estimated standard errors of the cross-validation losses.}
#'
#' @export
mvmr.pacs <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome, P,
    type = 2, rr_cut_off = 0, fix.cor = FALSE,
    betawt = NULL, tau = c(0.5, 1, 2, 3),
    CV_fold = 5, n_times = 1, lambda.length = 20,
    eps = 1e-4, digit = 3, seed = NULL, lambda.1se = TRUE,
    pleiotropy = TRUE,                   # NEW
    lambda.alpha = NULL,                 # NEW
    lambda.alpha.length = 20,            # NEW
    lambda.alpha.min.ratio = 0.05        # NEW
) {

  e <- 1 / CV_fold

  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  beta.outcome <- as.numeric(beta.outcome)
  se.outcome <- as.numeric(se.outcome)

  RR_obs <- stats::cor(beta.exposure)

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)

  iv_strength_parameter <- mr.divw::mvmr.ivw(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    gen_cor = P
  )$iv_strength_parameter

  initial.fit <- NULL
  if (is.null(betawt)) {
    if (pleiotropy) {
      initial.fit <- obtain_initial_pleiotropy(
        beta.exposure, se.exposure, beta.outcome, se.outcome,
        iv_strength_parameter, P,
        lambda.length = lambda.length,
        seed = seed,
        epsilon = eps,
        CV_fold = CV_fold,
        n_times = n_times,
        lambda.1se = lambda.1se,
        lambda.alpha = lambda.alpha,
        lambda.alpha.length = lambda.alpha.length,
        lambda.alpha.min.ratio = lambda.alpha.min.ratio
      )
    } else {
      initial.fit <- obtain_initial(
        beta.exposure, se.exposure, beta.outcome, se.outcome,
        iv_strength_parameter, P,
        lambda.length = lambda.length,
        seed = seed,
        epsilon = eps,
        CV_fold = CV_fold,
        n_times = n_times,
        lambda.1se = lambda.1se
      )
    }
    betawt <- initial.fit$beta
  } else {
    betawt <- as.numeric(betawt)
    if (length(betawt) != K || any(!is.finite(betawt))) {
      stop("betawt must be a finite numeric vector of length K.")
    }
  }

  lambda_cand <- generate_tuning_para(
    beta.exposure, se.exposure, se.outcome, K, P,
    lambda.length = lambda.length
  )

  # NEW: construct recommended lambda.alpha grid
  if (pleiotropy) {
    if (is.null(lambda.alpha)) {
      r0 <- as.numeric(beta.outcome - beta.exposure %*% betawt)
      lambda.alpha.max <- max(abs(r0) / se.outcome^2, na.rm = TRUE)

      if (!is.finite(lambda.alpha.max) || lambda.alpha.max <= 0) {
        stop("Cannot construct lambda.alpha grid: lambda.alpha.max is not finite or non-positive.")
      }

      lambda.alpha <- lambda.alpha.max *
        exp(seq(log(1), log(lambda.alpha.min.ratio), length.out = lambda.alpha.length))
    } else {
      lambda.alpha <- sort(unique(lambda.alpha), decreasing = TRUE)
    }
  } else {
    lambda.alpha <- 0
  }

  # UPDATED: grid includes lambda.alpha
  grid <- expand.grid(
    lambda = lambda_cand,
    tau = tau,
    lambda.alpha = lambda.alpha
  )

  obj_f_pacs <- array(
    NA_real_,
    dim = c(n_times, CV_fold, length(lambda_cand), length(tau), length(lambda.alpha))
  )

  numCores <- parallel::detectCores() - 1
  numCores <- max(1L, numCores)

  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)

  for (m in 1:n_times) {
    cat("CV times = ", m, "\n")

    Sig <- array(0, dim = c(p, K, K))
    for (j in 1:p) {
      Sig[j, , ] <- diag(se.exposure[j, ]) %*% P %*% diag(se.exposure[j, ])
    }

    if (!is.null(seed)) set.seed(seed + m)

    dt.exposure <- datathin::datathin(
      beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold
    )
    dt.outcome <- datathin::datathin(
      beta.outcome, family = "normal", arg = se.outcome^2, K = CV_fold
    )

    cv_results <- foreach::foreach(
      l = 1:CV_fold,
      .packages = c(
        "stats", "BDcocolasso", "mr.divw", "expm",
        "MASS", "datathin", "MVMR.PACS", "Matrix"
      )
    ) %dopar% {

      beta.exposure.train <- apply(dt.exposure[, , -l], c(1, 2), sum)
      se.exposure.train <- se.exposure * sqrt(1 - e)

      beta.outcome.train <- apply(dt.outcome[, 1, -l], 1, sum)
      se.outcome.train <- se.outcome * sqrt(1 - e)

      if (fix.cor) {
        rr.train <- RR_obs
      } else {
        rr.train <- stats::cor(beta.exposure.train)
      }

      beta.exposure.test <- dt.exposure[, , l]
      se.exposure.test <- se.exposure * sqrt(e)

      beta.outcome.test <- dt.outcome[, 1, l]
      se.outcome.test <- se.outcome * sqrt(e)

      V.test <- Reduce("+", lapply(1:nrow(beta.exposure.test), function(j) {
        diag(se.exposure.test[j, ]) %*% P %*% diag(se.exposure.test[j, ]) *
          se.outcome.test[j]^(-2)
      }))

      W.test <- diag(1 / se.outcome.test^2)

      MV.test <- t(beta.exposure.test) %*% W.test %*% beta.exposure.test - V.test
      MV_plus.test <- BDcocolasso::ADMM_proj(
        MV.test / sqrt(p), epsilon = eps
      )$mat * sqrt(p)

      losses_array <- array(
        NA_real_,
        dim = c(length(lambda_cand), length(tau), length(lambda.alpha))
      )

      for (i in 1:nrow(grid)) {
        lambda_val <- grid[i, "lambda"]
        tau_val <- grid[i, "tau"]
        lambda_alpha_val <- grid[i, "lambda.alpha"]

        fit <- tryCatch({
          dIVW_PACS_cluster(
            beta.exposure = beta.exposure.train,
            se.exposure = se.exposure.train,
            beta.outcome = beta.outcome.train,
            se.outcome = se.outcome.train,
            P = P,
            RR = rr.train,
            type = type,
            rr = rr_cut_off,
            betawt = betawt,
            eps = eps,
            lambda = lambda_val,
            tau = tau_val,
            digit = digit,
            pleiotropy = pleiotropy,
            lambda.alpha = lambda_alpha_val
          )
        }, error = function(e) {
          NULL
        })

        lambda_idx <- which(lambda_cand == lambda_val)
        tau_idx <- which(tau == tau_val)
        lambda_alpha_idx <- which(lambda.alpha == lambda_alpha_val)

        # NEW: safely skip failed or divergent fits
        if (is.null(fit) || is.null(fit$coefficients)) {
          losses_array[lambda_idx, tau_idx, lambda_alpha_idx] <- NA_real_
          next
        }

        beta_hat <- as.numeric(fit$coefficients)

        if (any(!is.finite(beta_hat))) {
          losses_array[lambda_idx, tau_idx, lambda_alpha_idx] <- NA_real_
          next
        }

        # NEW: validation loss with alpha if pleiotropy = TRUE
        if (pleiotropy) {
          alpha_hat_train <- as.numeric(fit$alpha)

          if (length(alpha_hat_train) != p || any(!is.finite(alpha_hat_train))) {
            losses_array[lambda_idx, tau_idx, lambda_alpha_idx] <- NA_real_
            next
          }

          # Data thinning scales the mean signal.
          # Training alpha estimates correspond to (1-e) * alpha_j.
          # Test fold has mean e * alpha_j.
          alpha_hat_test <- alpha_hat_train * e / (1 - e)

        } else {
          alpha_hat_test <- rep(0, p)
        }

        y_alpha_test <- beta.outcome.test - alpha_hat_test

        XW_y <- t(beta.exposure.test) %*% (W.test %*% y_alpha_test)

        loss <- t(beta_hat) %*% (MV_plus.test %*% beta_hat) -
          2 * crossprod(XW_y, beta_hat) +
          as.numeric(t(y_alpha_test) %*% W.test %*% y_alpha_test)

        if (!is.finite(loss)) {
          loss <- NA_real_
        }

        losses_array[lambda_idx, tau_idx, lambda_alpha_idx] <- as.numeric(loss)
      }

      losses_array
    }

    for (l in 1:CV_fold) {
      obj_f_pacs[m, l, , , ] <- cv_results[[l]]
    }
  }

  parallel::stopCluster(cl)

  valid_n <- apply(obj_f_pacs, c(3, 4, 5), function(x) sum(is.finite(x)))

  res_mean <- apply(obj_f_pacs, c(3, 4, 5), function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x)
  })

  res_sd <- apply(obj_f_pacs, c(3, 4, 5), function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(Inf)
    stats::sd(x) / sqrt(length(x))
  })

  if (all(!is.finite(res_mean))) {
    stop("All CV fits failed or diverged. Try increasing lambda.alpha.min.ratio or lambda.alpha.")
  }

  n_cv_total <- CV_fold * n_times
  valid_full <- valid_n == n_cv_total & is.finite(res_mean)

  if (!any(valid_full)) {
    stop(
      "No tuning parameter combination had fully valid CV fits across all ",
      n_cv_total,
      " data-thinning runs. Try increasing lambda.alpha.min.ratio, ",
      "using a larger lambda.alpha grid, or checking convergence failures."
    )
  }

  warning(
    "Selecting tuning parameters only among combinations with fully valid CV fits across all ",
    n_cv_total,
    " data-thinning runs, so that the CV standard error for the minimizer is well defined."
  )

  res_mean_valid <- res_mean

  res_mean_valid[!valid_full] <- Inf

  best_idx <- which(
    res_mean_valid == min(res_mean_valid, na.rm = TRUE),
    arr.ind = TRUE
  )[1, ]

  if (lambda.1se) {
    best_se <- res_sd[best_idx[1], best_idx[2], best_idx[3]]

    if (!is.finite(best_se)) {
      warning(
        "The selected model does not have a finite CV standard error; ",
        "using minimum-CV rule instead of 1SE rule."
      )
    } else {
      cutoff <- res_mean[best_idx[1], best_idx[2], best_idx[3]] + best_se

      true_indices <- which(
        res_mean <= cutoff &
          is.finite(res_mean) &
          valid_n == n_cv_total,
        arr.ind = TRUE
      )

      if (nrow(true_indices) > 0) {
        lambda_alpha_values_ok <- lambda.alpha[true_indices[, 3]]

        max_lambda_alpha <- max(lambda_alpha_values_ok, na.rm = TRUE)
        keep1 <- true_indices[lambda_alpha_values_ok == max_lambda_alpha, , drop = FALSE]

        lambda_values_ok <- lambda_cand[keep1[, 1]]

        max_lambda <- max(lambda_values_ok, na.rm = TRUE)
        keep2 <- keep1[lambda_values_ok == max_lambda, , drop = FALSE]

        tau_values_ok <- tau[keep2[, 2]]

        max_tau <- max(tau_values_ok, na.rm = TRUE)
        keep3 <- keep2[tau_values_ok == max_tau, , drop = FALSE]

        best_idx <- keep3[1, ]
      } else {
        warning(
          "No fully valid tuning parameter combination was within the 1SE cutoff; ",
          "using minimum-CV rule instead."
        )
      }
    }
  }

  pacs.model <- dIVW_PACS_cluster(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    P = P,
    RR = RR_obs,
    type = type,
    rr = rr_cut_off,
    lambda = lambda_cand[best_idx[1]],
    betawt = betawt,
    tau = tau[best_idx[2]],
    eps = eps,
    digit = digit,
    pleiotropy = pleiotropy,
    lambda.alpha = lambda.alpha[best_idx[3]]
  )

  beta.pacs <- pacs.model$coefficients
  cluster.pacs <- pacs.model$cls

  res <- list()
  res$beta <- beta.pacs
  res$alpha <- pacs.model$alpha
  res$invalid.snp <- pacs.model$invalid.snp
  res$signal_group <- cluster.pacs
  res$iv_strength_parameter <- iv_strength_parameter

  res$lambda.index <- best_idx[1]
  res$lambda <- lambda_cand[best_idx[1]]

  res$tau.index <- best_idx[2]
  res$tau <- tau[best_idx[2]]

  res$lambda.alpha.index <- best_idx[3]
  res$lambda.alpha <- lambda.alpha[best_idx[3]]
  res$lambda.alpha.grid <- lambda.alpha

  res$pleiotropy <- pleiotropy
  res$beta.init <- betawt
  res$initial.fit <- initial.fit
  res$initial.invalid.snp <- if (!is.null(initial.fit$invalid.snp)) {
    initial.fit$invalid.snp
  } else {
    integer(0)
  }
  res$cv.loss.mean <- res_mean
  res$cv.loss.sd <- res_sd

  set.seed(NULL)

  return(res)
}
