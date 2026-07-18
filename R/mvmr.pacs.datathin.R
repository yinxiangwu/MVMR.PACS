#' Function to perform MVMR-PACS estimator for signal-group identification and variable selection, along with two-fold data thinning for post-selection inference.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effects of a SNP on the K exposures, usually obtained from a GWAS.
#' @param se.exposure A data.frame or matrix of estimated standard errors corresponding to beta.exposure.
#' @param beta.outcome A numeric vector of the estimated marginal effects of SNPs on the outcome, usually obtained from a GWAS.
#' @param se.outcome A numeric vector of estimated standard errors corresponding to beta.outcome.
#' @param P A K-by-K matrix for the estimated shared correlation matrix among the SNP-exposure association estimates, where K is the number of exposures.
#' @param type Weighting scheme in penalization: 1 = weights independent of the correlation between beta.exposure; 2 = weights inversely related to the correlation between beta.exposure; 3 = weights nonzero only if the absolute correlation between beta.exposure is greater than rr_cut_off; 4 = weights inversely related to the correlation between beta.exposure only if the absolute correlation between beta.exposure is greater than rr_cut_off.
#' @param rr_cut_off Pre-specified correlation cutoff used in the weighting schemes specified by type. Default = 0.
#' @param tau Numeric vector of candidate values for the adaptive weighting parameter. Default = c(0.5, 1, 2, 3).
#' @param fold Number of folds used for cross-validation based on multi-fold data thinning within the selection dataset. Default = 5.
#' @param over.dispersion Logical. Whether SRIVW should account for balanced horizontal pleiotropy when performing inference. Default = FALSE.
#' @param re Number of repeated two-fold data-thinning runs. Default = 1. Increase this value to assess robustness of selected signal-groups and post-selection inference.
#' @param n_times Number of repeated cross-validation runs used when selecting tuning parameters within the selection dataset. Default = 1.
#' @param epsilon Tolerance level for ADMM projection. Default = 1e-4.
#' @param digit Number of digits to keep when determining signal-groups. Default = 3, corresponding to decimal precision 0.001.
#' @param include_zero_group Logical. Whether to treat the non-causal group identified by MVMR-PACS from the selection dataset as a group for post-selection inference. Default = FALSE.
#' @param lambda.1se Logical. Whether to use the one-standard-error rule when selecting tuning parameters. Default = TRUE.
#' @param lambda.length Length of the candidate tuning parameter sequence for the PACS penalty. Default = 20.
#' @param pleiotropy Logical. Whether to allow and select SNP-specific pleiotropic effects using an additional l1 penalty in the selection dataset. Default = TRUE.
#' @param lambda.alpha Optional numeric vector of candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects. Default = NULL, in which case a data-dependent decreasing sequence is constructed.
#' @param lambda.alpha.length Length of the candidate tuning parameter sequence for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Default = 20.
#' @param lambda.alpha.min.ratio Ratio between the smallest and largest candidate tuning parameter values for the l1 penalty on SNP-specific pleiotropic effects when lambda.alpha = NULL. Default = 0.05.
#'
#' @return A list with the following elements:
#' \item{grouping.hist}{Identified signal-groups across repeated data-thinning runs.}
#' \item{beta.hist}{Matrix of estimated direct effects of each signal-group on the outcome from SRIVW across repeated data-thinning runs.}
#' \item{se.hist}{Matrix of estimated standard errors of direct effects of each signal-group on the outcome from SRIVW across repeated data-thinning runs.}
#' \item{iv.hist.dt1}{Sample IV strength parameter in the selection dataset across repeated data-thinning runs, after excluding selected pleiotropic SNPs if pleiotropy = TRUE.}
#' \item{iv.hist.dt2}{Sample IV strength parameter in the inferential dataset based on identified signal-groups across repeated data-thinning runs, after excluding selected pleiotropic SNPs if pleiotropy = TRUE.}
#' \item{condF.dt1}{Conditional F-statistics in the selection dataset across repeated data-thinning runs, after excluding selected pleiotropic SNPs if pleiotropy = TRUE.}
#' \item{condF.dt2}{Conditional F-statistics in the inferential dataset based on identified signal-groups across repeated data-thinning runs, after excluding selected pleiotropic SNPs if pleiotropy = TRUE.}
#' \item{alpha.hist}{Matrix of estimated SNP-specific pleiotropic effects from the selection dataset across repeated data-thinning runs. If pleiotropy = FALSE, entries are zero.}
#' \item{invalid.snp.hist}{A list containing indices of SNPs identified as pleiotropic in the selection dataset for each repeated data-thinning run.}
#' \item{n.invalid.snp}{Number of SNPs identified as pleiotropic in the selection dataset for each repeated data-thinning run.}
#' \item{pleiotropy}{Logical indicator of whether SNP-specific pleiotropic effects were modeled.}
#'
#' @export
mvmr.pacs.datathin <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome, P,
    type = 2, rr_cut_off = 0, tau = c(0.5, 1, 2, 3),
    fold = 5, over.dispersion = FALSE,
    re = 1, n_times = 1, epsilon = 1e-4, digit = 3,
    include_zero_group = FALSE, lambda.length = 20, lambda.1se = TRUE,
    pleiotropy = TRUE,
    lambda.alpha = NULL,
    lambda.alpha.length = 20,
    lambda.alpha.min.ratio = 0.05
) {

  e1 <- 1 / 2

  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)

  grouping_hist <- rep(NA, re)

  beta_est <- matrix(NA, nrow = re, ncol = K)
  beta_se  <- matrix(NA, nrow = re, ncol = K)

  iv.hist.dt1 <- rep(NA, re)
  iv.hist.dt2 <- rep(NA, re)

  condF.dt1 <- matrix(NA, nrow = re, ncol = K)
  condF.dt2 <- matrix(NA, nrow = re, ncol = K)

  alpha_hist <- matrix(NA, nrow = re, ncol = p)       # NEW
  invalid.snp.hist <- vector("list", re)              # NEW
  n.invalid.snp <- rep(NA_integer_, re)               # NEW
  initial.invalid.snp.hist <- vector("list", re)

  for (r in 1:re) {

    cat("start analysis run:", r, "\n")

    Sig <- array(0, dim = c(p, K, K))
    for (j in 1:p) {
      Sig[j, , ] <- diag(se.exposure[j, ]) %*% P %*% diag(se.exposure[j, ])
    }

    dt.exposure <- datathin::datathin(
      beta.exposure, family = "mvnormal", arg = Sig, K = 2
    )

    dt.outcome <- datathin::datathin(
      beta.outcome, family = "normal", arg = se.outcome^2, K = 2
    )

    beta.exposure.train <- apply(dt.exposure[, , -2], c(1, 2), sum)
    se.exposure.train   <- se.exposure * sqrt(1 - e1)

    beta.outcome.train <- as.vector(apply(dt.outcome[, 1, -2, drop = FALSE], 1, sum))
    se.outcome.train   <- se.outcome * sqrt(1 - e1)

    beta.exposure.test <- dt.exposure[, , 2]
    se.exposure.test   <- se.exposure * sqrt(e1)

    beta.outcome.test <- dt.outcome[, 1, 2]
    se.outcome.test   <- se.outcome * sqrt(e1)

    iv_str_train <- mr.divw::mvmr.ivw(
      beta.exposure = beta.exposure.train,
      se.exposure = se.exposure.train,
      beta.outcome = beta.outcome.train,
      se.outcome = se.outcome.train,
      gen_cor = P
    )$iv_strength_parameter

    if (pleiotropy) {
      initial.fit <- obtain_initial_pleiotropy(
        beta.exposure.train, se.exposure.train,
        beta.outcome.train, se.outcome.train,
        iv_strength_parameter = iv_str_train,
        n_times = 5, P = P,
        epsilon = epsilon,
        lambda.length = lambda.length,
        lambda.1se = lambda.1se,
        lambda.alpha = lambda.alpha,
        lambda.alpha.length = lambda.alpha.length,
        lambda.alpha.min.ratio = lambda.alpha.min.ratio
      )
      initial.invalid.snp.hist[[r]] <- initial.fit$invalid.snp
    } else {
      initial.fit <- obtain_initial(
        beta.exposure.train, se.exposure.train,
        beta.outcome.train, se.outcome.train,
        iv_strength_parameter = iv_str_train,
        n_times = 5, P = P,
        epsilon = epsilon,
        lambda.length = lambda.length,
        lambda.1se = lambda.1se
      )
      initial.invalid.snp.hist[[r]] <- integer(0)
    }
    betawt <- initial.fit$beta

    cat("Perform PACS on the training set to obtain grouping \n")

    pacs_res <- mvmr.pacs(
      beta.exposure = beta.exposure.train,
      se.exposure = se.exposure.train,
      beta.outcome = beta.outcome.train,
      se.outcome = se.outcome.train,
      P = P,
      type = type,
      rr_cut_off = rr_cut_off,
      betawt = betawt,
      CV_fold = fold,
      tau = tau,
      lambda.length = lambda.length,
      eps = epsilon,
      digit = digit,
      n_times = n_times,
      lambda.1se = lambda.1se,
      pleiotropy = pleiotropy,
      lambda.alpha = lambda.alpha,
      lambda.alpha.length = lambda.alpha.length,
      lambda.alpha.min.ratio = lambda.alpha.min.ratio
    )

    if (all(round(pacs_res$beta, digit) == 0)) {
      pacs_res$beta <- rep(10^(-digit), K)
    }

    grouping_hist[r] <- paste(pacs_res$signal_group, collapse = "-")

    # NEW: store pleiotropy results from selection dataset
    if (pleiotropy) {
      alpha_hist[r, ] <- pacs_res$alpha
      invalid.snp <- pacs_res$invalid.snp
    } else {
      alpha_hist[r, ] <- rep(0, p)
      invalid.snp <- integer(0)
    }

    invalid.snp.hist[[r]] <- invalid.snp
    n.invalid.snp[r] <- length(invalid.snp)

    # NEW: exclude selected pleiotropic SNPs before SRIVW inference
    keep.snp <- setdiff(seq_len(p), invalid.snp)

    if (length(keep.snp) <= K) {
      warning(
        "Too few SNPs remain after excluding pleiotropic SNPs in run ",
        r, ". SRIVW inference skipped for this run."
      )
      next
    }

    beta.exposure.train.keep <- beta.exposure.train[keep.snp, , drop = FALSE]
    se.exposure.train.keep   <- se.exposure.train[keep.snp, , drop = FALSE]
    beta.outcome.train.keep  <- beta.outcome.train[keep.snp]
    se.outcome.train.keep    <- se.outcome.train[keep.snp]

    beta.exposure.test.keep <- beta.exposure.test[keep.snp, , drop = FALSE]
    se.exposure.test.keep   <- se.exposure.test[keep.snp, , drop = FALSE]
    beta.outcome.test.keep  <- beta.outcome.test[keep.snp]
    se.outcome.test.keep    <- se.outcome.test[keep.snp]

    res_dt1 <- SRIVW_group_est(
      beta.exposure = beta.exposure.train.keep,
      se.exposure = se.exposure.train.keep,
      beta.outcome = beta.outcome.train.keep,
      se.outcome = se.outcome.train.keep,
      P = P,
      est_beta = pacs_res$beta,
      digit = digit,
      include_zero_group = include_zero_group,
      over.dispersion = over.dispersion
    )

    iv.hist.dt1[r] <- res_dt1$iv_strength
    condF.dt1[r, ] <- res_dt1$condF

    cat("Using SRIVW to obtain causal effect estimates of grouped risk factors \n")

    res_dt2 <- SRIVW_group_est(
      beta.exposure = beta.exposure.test.keep,
      se.exposure = se.exposure.test.keep,
      beta.outcome = beta.outcome.test.keep,
      se.outcome = se.outcome.test.keep,
      P = P,
      est_beta = pacs_res$beta,
      digit = digit,
      include_zero_group = include_zero_group,
      over.dispersion = over.dispersion
    )

    beta_est[r, ] <- res_dt2$beta_est
    beta_se[r, ]  <- res_dt2$beta_se

    iv.hist.dt2[r] <- res_dt2$iv_strength
    condF.dt2[r, ] <- res_dt2$condF
  }

  res <- list()
  res$grouping.hist <- grouping_hist
  res$beta.hist <- beta_est
  res$se.hist <- beta_se
  res$iv.hist.dt1 <- iv.hist.dt1
  res$iv.hist.dt2 <- iv.hist.dt2
  res$condF.dt1 <- condF.dt1
  res$condF.dt2 <- condF.dt2

  res$alpha.hist <- alpha_hist
  res$invalid.snp.hist <- invalid.snp.hist
  res$n.invalid.snp <- n.invalid.snp
  res$initial.invalid.snp.hist <- initial.invalid.snp.hist
  res$pleiotropy <- pleiotropy

  return(res)
}
