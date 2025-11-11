strength_mvmr2 <- function (r_input, gencov) 
{
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  if (!("mvmr_format" %in% class(r_input))) {
    stop("The class of the data object must be \"mvmr_format\", please resave the object with the output of format_mvmr().")
  }
  if (missing(gencov)) {
    gencov <- as.numeric(0)
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  Wj <- 1/r_input[, 3]^2
  exp.number <- length(names(r_input)[-c(1, 2, 3)])/2
  A <- summary(lm(as.formula(paste("betaYG~ -1 +", paste(names(r_input)[seq(4, 
                                                                            3 + exp.number, by = 1)], collapse = "+"))), data = r_input))$coef
  for (i in 1:exp.number) {
    dimnames(A)[[1]][i] <- paste0("exposure", i, collapse = "")
  }
  delta_mat <- matrix(0, ncol = exp.number, nrow = exp.number - 
                        1)
  for (i in 1:exp.number) {
    regressand <- names(r_input[3 + i])
    regressors <- names(r_input)[-c(1, 2, 3, 4 + exp.number:length(names(r_input)))]
    C <- paste(regressand, "~", "-1 +", paste(regressors[-i], 
                                              collapse = "+"))
    D.reg <- lm(C, data = r_input)
    delta_mat[, i] <- D.reg$coefficients
  }
  sigma2xj_dat <- matrix(ncol = exp.number, nrow = length(r_input[, 
                                                                  1]), 0)
  if (length(gencov) < 2) {
    sebetas <- r_input[, (exp.number + 4):length(r_input)]
    for (i in 1:exp.number) {
      se.temp <- as.matrix(sebetas[, -i])
      for (j in 1:(exp.number - 1)) {
        sigma2xj_dat[, i] <- sigma2xj_dat[, i] + (se.temp[, 
                                                          j]^2 * delta_mat[j, i]^2)
      }
      sigma2xj_dat[, i] <- sigma2xj_dat[, i] + sebetas[, 
                                                       i]^2
    }
  }
  if (length(gencov) > 2) {
    sigma2xj_dat <- matrix(ncol = exp.number, nrow = length(r_input[, 
                                                                    1]), 0)
    delta.temp <- matrix(0, ncol = exp.number, nrow = exp.number)
    for (i in 1:exp.number) {
      if (i == 1) {
        delta.temp[, i] <- c(-1, delta_mat[, i])
      }
      if (i > 1 & i < exp.number) {
        delta.temp[, i] <- c(delta_mat[1:(i - 1), i], 
                             -1, delta_mat[i:(exp.number - 1), i])
      }
      if (i == exp.number) {
        delta.temp[, i] <- c(delta_mat[, i], -1)
      }
      for (l in 1:length(r_input[, 1])) {
        sigma2xj_dat[l, i] <- sigma2xj_dat[l, i] + t(delta.temp[, 
                                                                i]) %*% gencov[[l]] %*% delta.temp[, i]
      }
    }
  }
  Q_strength <- matrix(ncol = exp.number, nrow = 1, 0)
  for (i in 1:exp.number) {
    betas <- r_input[, c(4:(3 + exp.number))]
    betas <- data.frame(betas[, -i])
    temp.sub <- 0
    for (j in 1:(exp.number - 1)) {
      temp.sub <- temp.sub + (delta_mat[j, i] * betas[, 
                                                      j])
    }
    Q_strength[i] <- sum((1/sigma2xj_dat[, i]) * ((r_input[, 
                                                           3 + i] - temp.sub)^2))
    Q_strength[i] <- Q_strength[i]/nrow(r_input)
  }
  Q_strength <- data.frame(Q_strength)
  names(Q_strength) <- dimnames(A)[[1]]
  rownames(Q_strength) <- "F-statistic"
  #cat("\n")
  #cat("Conditional F-statistics for instrument strength\n")
  #cat("\n")
  #print(Q_strength)
  return(Q_strength)
}