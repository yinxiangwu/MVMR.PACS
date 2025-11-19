#' Function to perform MVMR-PACS estimator for signal-group identification and variable selection in two-sample multivariable Mendeliean randomization with highly-correlated exposures.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param P A K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure.
#' @param type Weighting scheme in penalization: 1 = weights independent of correlation between beta.exposure; 2 = weights inversely related to correlation between beta.exposure; 3 = weights non-zero only if absolute correlation between beta.exposure greater than rr_cut_off; 4 = weights inversely related to correlation between beta.exposure only if absolute correlation between beta.exposure greater than rr_cut_off.
#' @param rr_cut_off Pre-specified cut-off used in different weighting schemes, specified by type argument.
#' @param fix.cor Whether or not recompute correlation of beta.exposure in cross validation. Default is FALSE.
#' @param betawt Initial beta coefficients. Default = NULL, meaning MVMR.dIVW with ridge penalty is used.
#' @param tau Tuning parameter range. Default = c(0.5, 1, 2, 3).
#' @param CV_fold The number of folds for cross validation based on multifold data thinning. Default = 5.
#' @param n_times The number of repeated cross validation. Default = 1.
#' @param lambda.length Length of candidate penalization parameter. Default = 20.
#' @param eps Tolerance level for ADMM.proj. Default = 1e-4.
#' @param digit Number of digits to keep when determining signal-groups. Default = 3, i.e., decimal precision to 0.001.
#' @param seed Seed for reproducibility. Default = NULL.
#' @param lambda.1se Whether or not to perform 1SE rules when selecting tuning parameters. Default = TRUE.
#'
#' @returns A list with elements
#' \item{beta}{Estimated direct effects of each exposure on the outcome. Exposures with estimated effects of same magnitude (up to a prespecified decimal precision) belong to the same signal-group.}
#' \item{signal_group}{Identified signal-groups}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' \item{lambda.index}{Index for the best penalization parameter lambda}
#' \item{lambda}{Value for the best penalization parameter lambda}
#' \item{tau.index}{Index for the best penalization parameter tau}
#' \item{tau}{Value for the best penalization parameter tau}
#' \item{beta.init}{Initial beta estimates}
#' @export
#'
#' @examples
mvmr.pacs <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P, type = 2, rr_cut_off = 0, fix.cor = FALSE,
                                  betawt = NULL, tau = c(0.5,1,2,3), CV_fold = 5, n_times = 1, lambda.length = 20,
                                  eps = 1e-4, digit = 3, seed = NULL, lambda.1se = TRUE) {
  # Define the fraction for CV splits
  e <- 1/CV_fold

  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)

  RR_obs <- stats::cor(beta.exposure)

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)

  iv_strength_parameter <- mr.divw::mvmr.ivw(beta.exposure = beta.exposure,
                                    se.exposure = se.exposure,
                                    beta.outcome = beta.outcome,
                                    se.outcome = se.outcome,
                                    gen_cor = P)$iv_strength_parameter

  if (is.null(betawt)) {
    # obtain initial beta estimates for subsequent estimators
    betawt <- obtain_initial(beta.exposure, se.exposure, beta.outcome, se.outcome, iv_strength_parameter, P)
  }

  # Fit dIVW PACS: create weight matrix and grid of tuning parameters
  W <- diag(1 / se.outcome^2)
  lambda_cand <- generate_tuning_para(beta.exposure, se.exposure, se.outcome, K, P, lambda.length = lambda.length)
  grid <- expand.grid(lambda = lambda_cand, tau = tau)

  # Prepare array to store CV losses
  obj_f_pacs <- array(0, dim = c(n_times, CV_fold, length(lambda_cand), length(tau)))

  ## ---- set up parallel backend (using available cores minus one) ----
  numCores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(numCores)

  ## the helper functions now come from the MVMR.PACS namespace on each worker

  doParallel::registerDoParallel(cl)

  for (m in 1:n_times) {
    cat('CV times = ', m, "\n")
    # Compute Sig once per m (independent of folds)
    Sig <- array(0, dim = c(p, K, K))
    for (j in 1:p) {
      Sig[j, , ] <- diag(se.exposure[j, ]) %*% P %*% diag(se.exposure[j, ])
    }
    if (!is.null(seed)) set.seed(seed + m)
    # Generate the CV data splits
    dt.exposure <- datathin::datathin(beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold)
    dt.outcome  <- datathin::datathin(beta.outcome, family = "normal", arg = se.outcome^2, K = CV_fold)

    # Parallelize over CV folds
    cv_results <- foreach::foreach(
      l = 1:CV_fold,
      .packages = c("stats", "BDcocolasso", "mr.divw", "expm",
                    "MASS", "datathin", "MVMR.PACS")  # <-- add your package here
    ) %dopar% {
      cat('Processing CV fold = ', l, "\n")
      # Training data
      beta.exposure.train <- apply(dt.exposure[, , -l], c(1, 2), sum)
      se.exposure.train   <- se.exposure * sqrt(1 - e)
      beta.outcome.train  <- apply(dt.outcome[, 1, -l], 1, sum)
      se.outcome.train    <- se.outcome * sqrt(1 - e)

      if (fix.cor) {
        rr.train <- RR_obs
      } else {
        rr.train <- stats::cor(beta.exposure.train)
      }

      # Testing data
      beta.exposure.test <- dt.exposure[, , l]
      se.exposure.test   <- se.exposure * sqrt(e)
      beta.outcome.test  <- dt.outcome[, 1, l]
      se.outcome.test    <- se.outcome * sqrt(e)

      # Compute test variance component and adjusted MV matrix
      V.test <- Reduce("+", lapply(1:nrow(beta.exposure.test), function(j) {
        diag(se.exposure.test[j, ]) %*% P %*% diag(se.exposure.test[j, ]) * (se.outcome.test[j]^(-2))
      }))
      W.test <- diag(1 / se.outcome.test^2)
      MV.test <- t(beta.exposure.test) %*% W.test %*% beta.exposure.test - V.test
      MV_plus.test <- BDcocolasso::ADMM_proj(MV.test / sqrt(p), epsilon = eps)$mat * sqrt(p)

      # Loop over grid of tuning parameters to compute the loss
      losses_matrix <- matrix(0, nrow = length(lambda_cand), ncol = length(tau))
      for (i in 1:nrow(grid)) {
        lambda_val <- grid[i, "lambda"]
        tau_val    <- grid[i, "tau"]
        beta_hat   <- tryCatch({
          dIVW_PACS_cluster(
            beta.exposure = beta.exposure.train,
            se.exposure   = se.exposure.train,
            beta.outcome  = beta.outcome.train,
            se.outcome    = se.outcome.train,
            P             = P,
            RR            = rr.train,
            type          = type,
            rr            = rr_cut_off,
            betawt        = betawt,
            eps           = eps,
            lambda        = lambda_val,
            tau           = tau_val,
            digit         = digit
          )$coefficients
        }, error = function(e) { return(rep(NA, K)) })

        XW_y   <- t(beta.exposure.test) %*% (W.test %*% beta.outcome.test)
        loss   <- t(beta_hat) %*% (MV_plus.test %*% beta_hat) - 2 * crossprod(XW_y, beta_hat)
        # Determine the corresponding indices for lambda and tau
        lambda_idx <- which(lambda_cand == lambda_val)
        tau_idx    <- which(tau == tau_val)
        losses_matrix[lambda_idx, tau_idx] <- as.numeric(loss)
      }
      return(losses_matrix)
    } # end foreach

    # Combine the results from each CV fold
    for (l in 1:CV_fold) {
      obj_f_pacs[m, l, , ] <- cv_results[[l]]
    }
  } # end for n_times

  # Shut down the cluster
  parallel::stopCluster(cl)

  # Average the loss over CV folds and repetitions
  res_mean <- apply(obj_f_pacs, c(3, 4), mean, na.rm = TRUE)
  res_sd   <- apply(obj_f_pacs, c(3, 4), stats::sd, na.rm = TRUE) / sqrt(CV_fold * n_times)

  # Choose the tuning parameters with the minimum CV loss
  best_idx <- which(res_mean == min(res_mean, na.rm = TRUE), arr.ind = TRUE)[1, ]

  if (lambda.1se) {
    # Find indices of TRUE entries
    true_indices <- which(res_mean <= res_mean[best_idx[1],best_idx[2]] + res_sd[best_idx[1],best_idx[2]], arr.ind = TRUE)

    # Smallest lambda index
    min_row <- min(true_indices[, "row"])

    # Largest tau index
    max_col_in_min_row <- max(true_indices[true_indices[, "row"] == min_row, "col"])

    # tuning parameters with 1se rule
    best_idx <- c(min_row, max_col_in_min_row)
  }

  # Re-fit the PACS model on the entire dataset with the selected tuning parameters
  pacs.model <- dIVW_PACS_cluster(beta.exposure = beta.exposure,
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
                                  digit = digit)

  beta.pacs    <- pacs.model$coefficients
  cluster.pacs <- pacs.model$cls

  # Package the results
  res <- list()
  res$beta                <- beta.pacs
  res$signal_group             <- cluster.pacs
  res$iv_strength_parameter <- iv_strength_parameter
  res$lambda.index        <- best_idx[1]
  res$lambda              <- lambda_cand[best_idx[1]]
  res$tau.index           <- best_idx[2]
  res$tau                 <- tau[best_idx[2]]
  res$beta.init           <- betawt

  set.seed(NULL)

  return(res)
}
