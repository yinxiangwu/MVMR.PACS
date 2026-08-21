#' Function to perform the MVMR-PACS estimator for signal-group identification and variable selection in two-sample multivariable Mendelian randomization with highly correlated exposures.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effects of a SNP on the K exposures, usually obtained from a GWAS.
#' @param se.exposure A data.frame or matrix of estimated standard errors corresponding to beta.exposure.
#' @param beta.outcome A numeric vector of the estimated marginal effects of SNPs on the outcome, usually obtained from a GWAS.
#' @param se.outcome A numeric vector of estimated standard errors corresponding to beta.outcome.
#' @param P A K-by-K exposure correlation matrix when overlap = FALSE, or a (K+1)-by-(K+1) joint exposure-outcome correlation matrix when overlap = TRUE.
#' @param overlap Logical. Whether exposure and outcome GWAS summary statistics overlap. Default = FALSE.
#' @param initial.overlap Logical. Whether the internally estimated initial coefficients should account for exposure-outcome sample overlap. By default, inherits overlap. It may be set to FALSE when overlap = TRUE so that the initial estimator ignores overlap while MVMR-PACS still accounts for it. Setting initial.overlap = TRUE when overlap = FALSE is not supported because P does not contain exposure-outcome correlations in that setting. This argument is only used when betawt = NULL.
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
#' @param initial.pleiotropy Logical. Whether the internally estimated initial coefficients should use the pleiotropy-robust initial estimator. By default, inherits pleiotropy. This argument is only used when betawt = NULL.
#' @param lambda.alpha Optional numeric vector of candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects. This grid is used by each pleiotropy-enabled stage, including the initial estimator when initial.pleiotropy = TRUE. Default = NULL, in which case a data-dependent decreasing sequence is constructed.
#' @param lambda.alpha.length Length of the candidate tuning parameter sequence for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Applies to each pleiotropy-enabled stage. Default = 20.
#' @param lambda.alpha.min.ratio Ratio between the smallest and largest candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Applies to each pleiotropy-enabled stage. Default = 0.05.
#' @param lambda.alpha.max.multiplier Multiplier applied to the data-dependent upper end of the lambda.alpha grid for each pleiotropy-enabled stage. Default = 4.
#' @param max.pleiotropy.prop Prespecified hard upper bound on the mean proportion of nonzero alpha estimates across CV folds for each pleiotropy-enabled stage. The admissible tuning set is restricted to combinations satisfying this bound before minimum-CV and 1SE selection. Default = 0.10.
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
#' \item{initial.overlap}{Requested overlap setting for the internal initial estimator; operative when betawt = NULL.}
#' \item{initial.pleiotropy}{Requested pleiotropy setting for the internal initial estimator; operative when betawt = NULL.}
#' \item{beta.init}{Initial beta estimates used to construct adaptive weights.}
#' \item{cv.loss.mean}{Array of mean cross-validation losses over repeated data-thinning folds.}
#' \item{cv.loss.sd}{Array of estimated standard errors of the cross-validation losses.}
#' \item{cv.pleiotropy.prop.mean}{Array of mean proportions of nonzero alpha estimates across CV runs for each tuning combination.}
#' \item{selected.cv.pleiotropy.prop}{Mean CV pleiotropic proportion for the selected tuning combination; when pleiotropy = TRUE this is constrained to be no larger than max.pleiotropy.prop.}
#' \item{final.pleiotropy.prop}{Proportion of nonzero alpha estimates after refitting the selected tuning combination to the full data.}
#' \item{final.pleiotropy.within.cap}{Logical indicator of whether the full-data refit also satisfies max.pleiotropy.prop.}
#' \item{sparse_pleiotropy_supported}{Logical indicator that an admissible CV tuning combination satisfying max.pleiotropy.prop was available. The function stops if none is available.}
#' \item{constrained.cv.loss}{CV loss of the minimum-loss model within the prespecified pleiotropy-constrained model class.}
#' \item{unconstrained.cv.loss}{CV loss of the unrestricted minimum-loss model, retained for diagnostics.}
#'
#' @export
mvmr.pacs <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome, P,
    overlap = FALSE,
    type = 2, rr_cut_off = 0, fix.cor = FALSE,
    betawt = NULL, tau = c(0.5, 1, 2, 3),
    CV_fold = 5, n_times = 1, lambda.length = 20,
    eps = 1e-4, digit = 3, seed = NULL, lambda.1se = TRUE,
    pleiotropy = TRUE,                   # NEW
    lambda.alpha = NULL,                 # NEW
    lambda.alpha.length = 20,            # NEW
    lambda.alpha.min.ratio = 0.05,       # NEW
    lambda.alpha.max.multiplier = 4,     # NEW
    max.pleiotropy.prop = 0.10,          # NEW
    initial.overlap = overlap,
    initial.pleiotropy = pleiotropy
) {

  e <- 1 / CV_fold

  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  beta.outcome <- as.numeric(beta.outcome)
  se.outcome <- as.numeric(se.outcome)

  validate_initial_options(
    overlap = overlap,
    pleiotropy = pleiotropy,
    initial.overlap = initial.overlap,
    initial.pleiotropy = initial.pleiotropy,
    initial.required = is.null(betawt)
  )

  if (isTRUE(pleiotropy) || (is.null(betawt) && isTRUE(initial.pleiotropy))) {
    if (length(lambda.alpha.max.multiplier) != 1L ||
        !is.finite(lambda.alpha.max.multiplier) ||
        lambda.alpha.max.multiplier < 1) {
      stop("lambda.alpha.max.multiplier must be at least 1.")
    }
    if (length(max.pleiotropy.prop) != 1L ||
        !is.finite(max.pleiotropy.prop) ||
        max.pleiotropy.prop <= 0 || max.pleiotropy.prop >= 1) {
      stop("max.pleiotropy.prop must be in (0, 1).")
    }
  }

  RR_obs <- stats::cor(beta.exposure)

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)
  overlap.info <- parse_overlap_correlation(P, K, overlap)
  P.exposure <- overlap.info$P.exposure
  rho.xy <- overlap.info$rho.xy
  P.initial <- if (isTRUE(initial.overlap)) P else P.exposure

  iv_strength_parameter <- mr.divw::mvmr.ivw(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    gen_cor = P.exposure
  )$iv_strength_parameter

  initial.fit <- NULL
  if (is.null(betawt)) {
    if (isTRUE(initial.pleiotropy)) {
      initial.fit <- obtain_initial_pleiotropy(
        beta.exposure, se.exposure, beta.outcome, se.outcome,
        iv_strength_parameter, P.initial,
        overlap = initial.overlap,
        lambda.length = lambda.length,
        seed = seed,
        epsilon = eps,
        CV_fold = CV_fold,
        n_times = n_times,
        lambda.1se = lambda.1se,
        lambda.alpha = lambda.alpha,
        lambda.alpha.length = lambda.alpha.length,
        lambda.alpha.min.ratio = lambda.alpha.min.ratio,
        lambda.alpha.max.multiplier = lambda.alpha.max.multiplier,
        max.pleiotropy.prop = max.pleiotropy.prop
      )
    } else {
      initial.fit <- obtain_initial(
        beta.exposure, se.exposure, beta.outcome, se.outcome,
        iv_strength_parameter, P.initial,
        overlap = initial.overlap,
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
    beta.exposure, se.exposure, se.outcome, K, P.exposure,
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
        exp(seq(
          log(lambda.alpha.max.multiplier),
          log(lambda.alpha.min.ratio),
          length.out = lambda.alpha.length
        ))
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
  prop_alpha_pacs <- array(
    NA_real_,
    dim = c(n_times, CV_fold, length(lambda_cand), length(tau), length(lambda.alpha))
  )

  numCores <- parallel::detectCores() - 1
  numCores <- max(1L, numCores)

  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)

  for (m in 1:n_times) {
    cat("CV times = ", m, "\n")

    if (isTRUE(overlap)) {
      Sig.joint <- make_joint_thinning_covariance(se.exposure, se.outcome, P)
      joint.summary <- cbind(beta.exposure, beta.outcome)
    } else {
      Sig <- make_thinning_covariance(se.exposure, P.exposure)
    }

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
        diag(se.exposure.test[j, ]) %*% P.exposure %*% diag(se.exposure.test[j, ]) *
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
      prop_alpha_array <- array(
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
            P = P.exposure,
            RR = rr.train,
            type = type,
            rr = rr_cut_off,
            betawt = betawt,
            eps = eps,
            lambda = lambda_val,
            tau = tau_val,
            digit = digit,
            pleiotropy = pleiotropy,
            lambda.alpha = lambda_alpha_val,
            overlap.bias = compute_overlap_bias(
              se.exposure.train, se.outcome.train, rho.xy
            )
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

          prop_alpha_array[
            lambda_idx, tau_idx, lambda_alpha_idx
          ] <- mean(abs(alpha_hat_train) > 1e-7)

          # Data thinning scales the mean signal.
          # Training alpha estimates correspond to (1-e) * alpha_j.
          # Test fold has mean e * alpha_j.
          alpha_hat_test <- alpha_hat_train * e / (1 - e)

        } else {
          alpha_hat_test <- rep(0, p)
          prop_alpha_array[
            lambda_idx, tau_idx, lambda_alpha_idx
          ] <- 0
        }

        y_alpha_test <- beta.outcome.test - alpha_hat_test

        XW_y <- t(beta.exposure.test) %*% (W.test %*% y_alpha_test) -
          compute_overlap_bias(se.exposure.test, se.outcome.test, rho.xy)

        loss <- t(beta_hat) %*% (MV_plus.test %*% beta_hat) -
          2 * crossprod(XW_y, beta_hat) +
          as.numeric(t(y_alpha_test) %*% W.test %*% y_alpha_test)

        if (!is.finite(loss)) {
          loss <- NA_real_
        }

        losses_array[lambda_idx, tau_idx, lambda_alpha_idx] <- as.numeric(loss)
      }

      list(loss = losses_array, prop_alpha = prop_alpha_array)
    }

    for (l in 1:CV_fold) {
      obj_f_pacs[m, l, , , ] <- cv_results[[l]]$loss
      prop_alpha_pacs[m, l, , , ] <- cv_results[[l]]$prop_alpha
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

  prop_alpha_mean <- apply(prop_alpha_pacs, c(3, 4, 5), function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  })

  prop_alpha_max <- apply(prop_alpha_pacs, c(3, 4, 5), function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x)
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

  # Keep the unconstrained minimum only as a diagnostic.
  res_mean_valid <- res_mean
  res_mean_valid[!valid_full] <- Inf

  unconstrained_idx <- which(
    res_mean_valid == min(res_mean_valid, na.rm = TRUE),
    arr.ind = TRUE
  )[1, ]

  # HARD SPARSITY CONSTRAINT:
  # Restrict the admissible tuning set before finding the CV minimum.
  # Thus max.pleiotropy.prop is honored by the selected CV model rather
  # than used only as a preference inside the unconstrained 1SE region.
  if (pleiotropy) {
    sparse_valid <- valid_full &
      is.finite(prop_alpha_mean) &
      prop_alpha_mean <= max.pleiotropy.prop

    if (!any(sparse_valid)) {
      min_prop <- suppressWarnings(min(
        prop_alpha_mean[valid_full & is.finite(prop_alpha_mean)],
        na.rm = TRUE
      ))
      if (!is.finite(min_prop)) min_prop <- NA_real_
      stop(
        "No fully valid MVMR-PACS tuning combination satisfies ",
        "max.pleiotropy.prop = ", max.pleiotropy.prop, ". ",
        "The smallest mean CV pleiotropic proportion among fully valid ",
        "candidates was ", signif(min_prop, 4), ". Increase the upper ",
        "range of lambda.alpha (for example via ",
        "lambda.alpha.max.multiplier) or reconsider the prespecified cap."
      )
    }
  } else {
    sparse_valid <- valid_full
  }

  # Minimum CV loss within the admissible (constrained) model class.
  res_mean_constrained <- res_mean
  res_mean_constrained[!sparse_valid] <- Inf
  constrained_idx <- which(
    res_mean_constrained == min(res_mean_constrained, na.rm = TRUE),
    arr.ind = TRUE
  )[1, ]
  best_idx <- constrained_idx

  # If requested, apply the 1SE rule within the constrained model class.
  # Preference order follows the scientific objective:
  #   1) largest PACS lambda (stronger beta grouping/selection),
  #   2) largest lambda.alpha (sparser SNP-specific pleiotropy),
  #   3) largest tau as the remaining tie-breaker.
  if (lambda.1se) {
    best_se <- res_sd[
      constrained_idx[1], constrained_idx[2], constrained_idx[3]
    ]

    if (!is.finite(best_se)) {
      warning(
        "The constrained minimum-CV model does not have a finite CV ",
        "standard error; using the constrained minimum-CV model instead ",
        "of the 1SE rule."
      )
    } else {
      cutoff <- res_mean[
        constrained_idx[1], constrained_idx[2], constrained_idx[3]
      ] + best_se

      eligible <- which(
        sparse_valid & is.finite(res_mean) & res_mean <= cutoff,
        arr.ind = TRUE
      )

      if (nrow(eligible) > 0L) {
        # 1. Stronger beta shrinkage/grouping first.
        lambda_values_ok <- lambda_cand[eligible[, 1]]
        max_lambda <- max(lambda_values_ok, na.rm = TRUE)
        keep1 <- eligible[
          lambda_values_ok == max_lambda, , drop = FALSE
        ]

        # 2. Then stronger alpha shrinkage.
        if (pleiotropy) {
          lambda_alpha_values_ok <- lambda.alpha[keep1[, 3]]
          max_lambda_alpha <- max(lambda_alpha_values_ok, na.rm = TRUE)
          keep2 <- keep1[
            lambda_alpha_values_ok == max_lambda_alpha, , drop = FALSE
          ]
        } else {
          keep2 <- keep1
        }

        # 3. Remaining tie-breaker: larger tau.
        tau_values_ok <- tau[keep2[, 2]]
        max_tau <- max(tau_values_ok, na.rm = TRUE)
        keep3 <- keep2[
          tau_values_ok == max_tau, , drop = FALSE
        ]

        best_idx <- keep3[1, ]
      }
    }
  }

  # By construction, the selected CV model satisfies the prespecified cap.
  sparse_pleiotropy_supported <- TRUE

  pacs.model <- dIVW_PACS_cluster(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    P = P.exposure,
    RR = RR_obs,
    type = type,
    rr = rr_cut_off,
    lambda = lambda_cand[best_idx[1]],
    betawt = betawt,
    tau = tau[best_idx[2]],
    eps = eps,
    digit = digit,
    pleiotropy = pleiotropy,
    lambda.alpha = lambda.alpha[best_idx[3]],
    overlap.bias = compute_overlap_bias(se.exposure, se.outcome, rho.xy)
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

  res$max.pleiotropy.prop <- max.pleiotropy.prop
  res$cv.pleiotropy.prop.mean <- prop_alpha_mean
  res$cv.pleiotropy.prop.max <- prop_alpha_max
  res$selected.cv.pleiotropy.prop <- prop_alpha_mean[
    best_idx[1], best_idx[2], best_idx[3]
  ]
  res$final.pleiotropy.prop <- if (pleiotropy) {
    mean(abs(pacs.model$alpha) > 1e-7)
  } else {
    0
  }
  res$final.pleiotropy.within.cap <-
    !pleiotropy || res$final.pleiotropy.prop <= max.pleiotropy.prop
  res$sparse_pleiotropy_supported <- sparse_pleiotropy_supported

  # Constrained-minimum diagnostics.
  res$constrained.lambda <- lambda_cand[constrained_idx[1]]
  res$constrained.tau <- tau[constrained_idx[2]]
  res$constrained.lambda.alpha <- lambda.alpha[constrained_idx[3]]
  res$constrained.cv.loss <- res_mean[
    constrained_idx[1], constrained_idx[2], constrained_idx[3]
  ]
  res$constrained.cv.pleiotropy.prop <- prop_alpha_mean[
    constrained_idx[1], constrained_idx[2], constrained_idx[3]
  ]

  # Unconstrained optimum retained only for sensitivity/diagnostic reporting.
  res$unconstrained.lambda <- lambda_cand[unconstrained_idx[1]]
  res$unconstrained.tau <- tau[unconstrained_idx[2]]
  res$unconstrained.lambda.alpha <- lambda.alpha[unconstrained_idx[3]]
  res$unconstrained.cv.loss <- res_mean[
    unconstrained_idx[1], unconstrained_idx[2], unconstrained_idx[3]
  ]
  res$unconstrained.cv.pleiotropy.prop <- prop_alpha_mean[
    unconstrained_idx[1], unconstrained_idx[2], unconstrained_idx[3]
  ]

  if (pleiotropy && !res$final.pleiotropy.within.cap) {
    warning(
      "The selected tuning combination satisfies the mean CV pleiotropy ",
      "cap (", signif(res$selected.cv.pleiotropy.prop, 4), " <= ",
      max.pleiotropy.prop, "), but the full-data refit selected a ",
      "pleiotropic proportion of ", signif(res$final.pleiotropy.prop, 4),
      ". The tuning parameters were not changed post hoc; inspect ",
      "final.pleiotropy.within.cap as a stability diagnostic."
    )
  }

  res$pleiotropy <- pleiotropy
  res$overlap <- overlap
  res$initial.pleiotropy <- initial.pleiotropy
  res$initial.overlap <- initial.overlap
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
