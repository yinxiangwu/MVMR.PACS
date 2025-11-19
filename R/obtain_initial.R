#' Obtain initial beta coefficients using MVMR.dIVW with ridge penalty.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param iv_strength_parameter Sample IV strength parameter. Output from mvmr.srivw or mvmr.ivw estimators.
#' @param P A K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure.
#' @param CV_fold The number of folds for cross validation based on multifold data thinning. Default = 5.
#' @param n_times The number of repeated cross validation. Default = 1.
#' @param epsilon Tolerance level for ADMM.proj. Default = 1e-4.
#' @param seed Seed for reproducibility. Default = NULL.
#' @param lambda.1se Whether or not to perform 1SE rules when selecting tuning parameters. Default = TRUE.
#' @param lambda.length Length of candidate penalization parameter. Default = 20.
#'
#' @returns
#' \item{betawt}{Estimated direct effects of each exposure on the outcome.}
#' @export
#'
#' @examples
obtain_initial <- function(
    beta.exposure,    # p x K
    se.exposure,      # p x K (SEs, >0)
    beta.outcome,     # p x 1 (or vector)
    se.outcome,       # p x 1 (or vector, >0)
    iv_strength_parameter,
    P,                # K x K (symmetric PSD)
    CV_fold = 5,
    n_times = 1,
    epsilon = 1e-4,
    seed = NULL,
    lambda.1se = TRUE,
    lambda.length = 20
){
  # --- basic checks ---
  p <- nrow(beta.exposure); K <- ncol(beta.exposure)
  stopifnot(nrow(se.exposure)==p, ncol(se.exposure)==K)
  stopifnot(length(se.outcome)==p, length(beta.outcome)==p)
  stopifnot(is.matrix(P), nrow(P)==K, ncol(P)==K)

  e <- 1/CV_fold

  # tuning grid
  upper <- (max(iv_strength_parameter, 0)*sqrt(p) + p)^(2/5)
  phi_cand <- exp(seq(log(upper*1e2), log(upper*1e-4), length.out = lambda.length))

  obj_f_ridge <- matrix(NA_real_, nrow = CV_fold*n_times, ncol = length(phi_cand))

  # precompute Sig array for mvnormal thinning
  Sig <- array(0, dim = c(p, K, K))
  for (j in 1:p) {
    D <- diag(se.exposure[j,], nrow=K, ncol=K)
    Sig[j,,] <- D %*% P %*% D
  }

  for (m in 1:n_times) {
    if (!is.null(seed)) set.seed(seed + m)

    dt.exposure <- datathin::datathin(beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold)
    dt.outcome  <- datathin::datathin(beta.outcome,  family = "normal",   arg = se.outcome^2, K = CV_fold)

    for (l in 1:CV_fold) {
      # train (sum over all folds except l)
      beta.exposure.train <- apply(dt.exposure[,,-l, drop=FALSE], c(1,2), sum)
      se.exposure.train   <- se.exposure * sqrt(1 - e)
      beta.outcome.train  <- apply(dt.outcome[,1,-l, drop=FALSE], 1, sum)
      se.outcome.train    <- se.outcome * sqrt(1 - e)

      # test (the l-th slice)
      beta.exposure.test <- dt.exposure[,,l]
      se.exposure.test   <- se.exposure * sqrt(e)
      beta.outcome.test  <- dt.outcome[,1,l]
      se.outcome.test    <- se.outcome * sqrt(e)

      # W.test as diagonal: prefer weights vector to avoid forming big diag
      wtest <- 1/(se.outcome.test^2)

      # MV.test = X' W X - sum_j Sigma_xj P Sigma_xj / se_yj^2
      # Vectorized V.test:
      # sum_j diag(se_exj) %*% P %*% diag(se_exj) * wtest[j]
      V.test <- matrix(0.0, K, K)
      for (j in 1:p) {
        Dj <- diag(se.exposure.test[j,], K, K)
        V.test <- V.test + (wtest[j]) * (Dj %*% P %*% Dj)
      }

      XtW  <- t(beta.exposure.test) * rep(wtest, each = K)  # K x p (each column j scaled by wtest[j])
      MV.test <- XtW %*% beta.exposure.test - V.test
      MV_plus.test <- BDcocolasso::ADMM_proj(MV.test/sqrt(p), epsilon = epsilon)$mat * sqrt(p)

      # c.test = X' W y
      c.test <- as.numeric(t(beta.exposure.test) %*% (wtest * beta.outcome.test))

      for (jphi in seq_along(phi_cand)) {
        b <- mvmr_ridge(
          beta.exposure = beta.exposure.train,
          se.exposure   = se.exposure.train,
          beta.outcome  = beta.outcome.train,
          se.outcome    = se.outcome.train,
          P             = P,
          phi           = phi_cand[jphi]
        )
        b <- as.numeric(b)  # ensure vector
        obj_f_ridge[l + (m-1)*CV_fold, jphi] <- drop(Matrix::crossprod(b, MV_plus.test %*% b)) - 2 * sum(c.test * b)
      }
    }
  }

  # CV summary (guard against NA)
  col_means <- colMeans(obj_f_ridge, na.rm = TRUE)
  col_se    <- apply(obj_f_ridge, 2, function(x) sd(x, na.rm = TRUE) / sqrt(CV_fold*n_times))

  # min + 1se rule
  s_min <- which.min(col_means)
  thresh <- col_means[s_min] + col_se[s_min]
  s_1se <- min(which(col_means <= thresh))

  phi_final <- if (isTRUE(lambda.1se)) phi_cand[s_1se] else phi_cand[s_min]

  # final fit on all data
  betawt <- mvmr_ridge(
    beta.exposure = beta.exposure,
    se.exposure   = se.exposure,
    beta.outcome  = beta.outcome,
    se.outcome    = se.outcome,
    P             = P,
    phi           = phi_final
  )
  betawt <- as.numeric(betawt)

  return(betawt)
}
