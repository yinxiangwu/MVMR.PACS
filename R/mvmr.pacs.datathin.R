#' Function to perform MVMR-PACS estimator for signal-group identification and variable selection, along with two-fold data thinning for post-selection inference.
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param P A K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure.
#' @param type Weighting scheme in penalization: 1 = weights independent of correlation between beta.exposure; 2 = weights inversely related to correlation between beta.exposure; 3 = weights non-zero only if absolute correlation between beta.exposure greater than rr_cut_off; 4 = weights inversely related to correlation between beta.exposure only if absolute correlation between beta.exposure greater than rr_cut_off.
#' @param rr_cut_off Pre-specified cut-off used in different weighting schemes, specified by type argument.
#' @param tau Tuning parameter range. Default = c(0.5, 1, 2, 3).
#' @param fold The number of folds for cross validation based on multifold data thinning. Default = 5.
#' @param over.dispersion Should SRIVW consider balanced horizontal pleiotropy when performing inference? Default is FALSE.
#' @param re The number of repeated data thinning procedure. Default = 1. Increase this number of assess robustness of results.
#' @param n_times The number of cross validation runs when selecting tuning parameters within the selection dataset. Default = 1.
#' @param epsilon Tolerance level for ADMM.proj. Default = 1e-4.
#' @param digit Number of digits to keep when determining signal-groups. Default = 3, i.e., decimal precision to 0.001.
#' @param include_zero_group Whether or not consider non-causal-group identified by MVMR-PACS from the selection dataset as a group for post-selection inference. Default = FALSE.
#' @param lambda.1se Whether or not to perform 1SE rules when selecting tuning parameters. Default = TRUE.
#'
#' @returns A list with elements
#' \item{grouping.hist}{Identified signal-groups across repeated data thinning runs.}
#' \item{beta.hist}{Matrix of estimated direct effect of each signal-group on the outcome from SRIVW across repeated data thinning runs.}
#' \item{se.hist}{Matrix of estimated standard error of direct effect of each signal-group on the outcome from SRIVW across repeated data thinning runs.}
#' \item{iv.hist.dt1}{Sample IV strength parameter in the selection dataset across repeated data thinning runs.}
#' \item{iv.hist.dt1}{Sample IV strength parameter in the inferential dataset based on signal-groups across repeated data thinning runs.}
#' \item{condF.dt1}{Conditional F-statistics in the selection dataset across repeated data thinning runs.}
#' \item{condF.dt2}{Conditional F-statistics in the inferential dataset based on signal-groups across repeated data thinning runs.}
#' @export
#'
#' @examples
mvmr.pacs.datathin <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P, type = 2, rr_cut_off = 0, tau = c(0.5,1,2,3),fold = 5,over.dispersion = FALSE,
                                  re = 1, n_times = 1,  epsilon = 1e-4, digit = 3, include_zero_group = FALSE, lambda.1se = TRUE) {

  # default is splitting the data into 10 folds, using half to pick grouping, the other half to perform estimation
  e1 <- 1/2
  e2 <- 1/fold

  # number of SNPs
  p <- nrow(beta.exposure)

  # number of risk factors
  K <- ncol(beta.exposure)

  # history of clustering
  grouping_hist <- rep(NA, re)

  # SRIVW estimates on the testing dataset
  beta_est <- matrix(NA,nrow = re, ncol = K)
  beta_se <- matrix(NA,nrow = re, ncol = K)

  # IV strength on the training dataset
  iv.hist.dt1 <- rep(NA, re)

  # IV strength on the testing dataset
  iv.hist.dt2 <- rep(NA, re)

  # conditional F-statistics on the training dataset
  condF.dt1 <- matrix(NA, nrow = re, ncol = K)

  # conditional F-statistics on the testing dataset
  condF.dt2 <- matrix(NA, nrow = re, ncol = K)

  # repeated process

  for (r in 1:re) {

    cat("start analysis run:", r, '\n')

    # data thinning

    Sig <- array(0, dim = c(p, K, K))

    for (j in 1:p) {
      Sig[j,,] = diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,])
    }

    dt.exposure <- datathin::datathin(beta.exposure, family="mvnormal", arg = Sig, K = 2)

    dt.outcome <- datathin::datathin(beta.outcome, family = 'normal', arg = se.outcome^2, K = 2)

    # using the first 5 folds to perform grouping selection
    # data based on the first half of the data
    beta.exposure.train <- apply(dt.exposure[, , -2], c(1, 2), sum)
    se.exposure.train <- se.exposure * sqrt(1-e1)
    beta.outcome.train <- as.vector(apply(dt.outcome[, 1, -2, drop = FALSE], 1, sum))
    se.outcome.train <- se.outcome * sqrt(1-e1)
    iv_str_train <- mr.divw::mvmr.ivw(beta.exposure = beta.exposure.train,se.exposure = se.exposure.train,beta.outcome = beta.outcome.train,se.outcome = se.outcome.train,gen_cor = P)$iv_strength_parameter

    # data based on the second half of the data for inference
    beta.exposure.test <- dt.exposure[,,2]
    se.exposure.test <- se.exposure * sqrt(e1)
    beta.outcome.test <- dt.outcome[,1,2]
    se.outcome.test <- se.outcome * sqrt(e1)

    # initial estimates
    betawt <- obtain_initial(beta.exposure.train,se.exposure.train,beta.outcome.train,se.outcome.train,iv_strength_parameter = iv_str_train, n_times = 5, P = P, lambda.1se = lambda.1se)

    # perform PACS on the training set
    cat("Perform PACS on the training set to obtain grouping \n")


    pacs_res <- mvmr.pacs(beta.exposure = beta.exposure.train,
                                      se.exposure = se.exposure.train,
                                      beta.outcome = beta.outcome.train,
                                      se.outcome = se.outcome.train,
                                      P = P,
                                      type = type,
                                      rr_cut_off = rr_cut_off,
                                      betawt = betawt,
                                      CV_fold = fold,
                                      tau = tau,
                                      eps = epsilon,
                                      digit = digit,
                                      n_times = n_times,
                                      lambda.1se = lambda.1se)

    if (all(round(pacs_res$beta,digit) == 0)) {pacs_res$beta <- rep(10^(-digit),K)} # this consider all risk factors as a group even if they jointly have no effect

    grouping_hist[r] <- paste(pacs_res$signal_group, collapse = "-")

    res_dt1 <- SRIVW_group_est(beta.exposure = beta.exposure.train,se.exposure = se.exposure.train,
                               beta.outcome = beta.outcome.train,se.outcome = se.outcome.train,
                               P = P, est_beta = pacs_res$beta, digit = digit, include_zero_group = include_zero_group, over.dispersion = over.dispersion)

    iv.hist.dt1[r] <- res_dt1$iv_strength

    condF.dt1[r,] <- res_dt1$condF

    # estimation on the other data
    cat("Using SRIVW to obtain causal effect estimates of grouped risk factors \n")

    res_dt2 <- SRIVW_group_est(beta.exposure = beta.exposure.test,se.exposure = se.exposure.test,
                               beta.outcome = beta.outcome.test,se.outcome = se.outcome.test,
                               P = P, est_beta = pacs_res$beta, digit = digit, include_zero_group = include_zero_group, over.dispersion = over.dispersion)

    beta_est[r,] <- res_dt2$beta_est

    beta_se[r,] <- res_dt2$beta_se

    iv.hist.dt2[r] <- res_dt2$iv_strength

    condF.dt2[r,] <- res_dt2$condF

  }

  # output
  res <- NULL
  res$grouping.hist <- grouping_hist
  res$beta.hist <- beta_est
  res$se.hist <- beta_se
  res$iv.hist.dt1 <- iv.hist.dt1
  res$iv.hist.dt2 <- iv.hist.dt2
  res$condF.dt1 <- condF.dt1
  res$condF.dt2 <- condF.dt2
  return(res)
}
