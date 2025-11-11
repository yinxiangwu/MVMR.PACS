#####################################################################
#                                                                   #
# Function to apply MVMR-PACS approach                              #
# Code adapted from  Sharma et al. Consistent group identification  #
# and variable selection in regression with correlated predictors.  #
# Pubmed: 23772171                                                  #
#####################################################################

library(doParallel)
library(foreach)

dIVW_PACS_cluster=function(beta.exposure,se.exposure,beta.outcome,se.outcome,lambda,betawt,P,RR,tau = 1,cluster = NULL,type=1,rr=0,digit=3,eps=1e-4,max.iter = 1000,err=1e-4)
{
  require(mvtnorm)
  require(MASS)
  require(Matrix)
  if(lambda <=0) {return(cat("ERROR: Lambda must be > 0. \n"))}
  if(rr <0) {return(cat("ERROR: RR must be >=0. \n"))}
  if(rr >1) {return(cat("ERROR: RR must be <=1. \n"))}
  if(eps <=0) {return(cat("ERROR: Eps must be > 0. \n"))}
  if(!(type %in% 1:4)){return(cat("ERROR: Type must be 1, 2, 3 or 4. \n"))}
  # data preparation
  p <- length(beta.outcome)
  K <- ncol(beta.exposure)
  littleeps=10^-7
  # require p>1
  qp=K*(K-1)/2
  vc=0
  row.vec1=0
  col.vec1=0
  for(w in 1:(K-1))
  {
    row.vec1=c(row.vec1,c((vc+1):(vc+(K-w))))
    col.vec1=c(col.vec1,rep(w,length(c((vc+1):(vc+(K-w))))))
    vc=vc+(K-w)
  }
  c1=1
  c2=K-1
  row.vec2=0
  col.vec2=0
  for(w in 1:(K-1))
  {
    row.vec2=c(row.vec2,c(c1:c2))
    col.vec2=c(col.vec2,c((w+1):K))
    c1=c2+1
    c2=c2+K-w-1
  }
  dm=sparseMatrix(i=c(row.vec1[-1],row.vec2[-1]),j=c(col.vec1[-1],col.vec2[-1]),x=c(rep(1,qp),rep(-1,qp)))
  dp=sparseMatrix(i=c(row.vec1[-1],row.vec2[-1]),j=c(col.vec1[-1],col.vec2[-1]),x=c(rep(1,qp),rep(1,qp)))
  rm(c1,c2,w,vc,row.vec1,col.vec1,row.vec2,col.vec2)
  # — Apply clustering mask if requested —
  if (!is.null(cluster)) {
    tmp_x <- 1:K
    check_dat <- data.frame(x = tmp_x, y = cluster)
    # For each row in dp/dm, if the two corresponding clusters differ, zero that row.
    for (i in 1:nrow(dp)) {
      cls_pair <- check_dat[cols[which(as.numeric(dm[i,]) != 0)], "y"]
      if (cls_pair[1] != cls_pair[2]) {
        dp[i,] <- rep(0, K)
      }
    }
    for (i in 1:nrow(dm)) {
      cls_pair <- check_dat[cols[which(as.numeric(dm[i,]) != 0)], "y"]
      if (cls_pair[1] != cls_pair[2]) {
        dm[i,] <- rep(0, K)
      }
    }
  }
  
  # — Construct penalty weight vector ascvec based on type —
  abs_dm_b <- as.vector(abs(dm %*% betawt))
  abs_dp_b <- as.vector(abs(dp %*% betawt))
  
  if (type == 1) {
    ascvec <- c(1/abs(betawt)^tau, 1/abs_dm_b^tau, 1/abs_dp_b^tau)
  } else if (type == 2) {
    # For symmetric matrices, lower.tri extracts the unique pairs.
    crm <- 1/(1 - RR[lower.tri(RR)])^tau
    crp <- 1/(1 + RR[lower.tri(RR)])^tau
    ascvec <- c(1/abs(betawt)^tau, crm/abs_dm_b^tau, crp/abs_dp_b^tau)
  } else if (type == 3) {
    corp <- ifelse(RR[lower.tri(RR)] > rr, 1, 0)
    corm <- ifelse(RR[lower.tri(RR)] < -rr, 1, 0)
    ascvec <- c(1/abs(betawt)^tau, corm/abs_dm_b^tau, corp/abs_dp_b^tau)
  } else if (type == 4) {
    corp <- ifelse(RR[lower.tri(RR)] > rr, 1, 0)
    crm  <- corp/(1 - RR[lower.tri(RR)])^tau
    corm <- ifelse(RR[lower.tri(RR)] < -rr, 1, 0)
    crp  <- corm/(1 + RR[lower.tri(RR)])^tau
    ascvec <- c(1/abs(betawt)^tau, crm/abs_dm_b^tau, crp/abs_dp_b^tau)
  }
  
  # — Precompute weight matrices and related quantities —
  W   <- Matrix::Diagonal(x = 1/se.outcome^2)
  M   <- crossprod(beta.exposure, W %*% beta.exposure)
  A_weighted <- se.exposure / matrix(se.outcome, nrow = p, ncol = K, byrow = FALSE)
  S <- t(A_weighted) %*% A_weighted
  V <- P * S
  MV  <- M - V
  MVplus <- ADMM_proj(as.matrix(MV)/sqrt(p), epsilon = eps)$mat * sqrt(p)
  XWy <- crossprod(beta.exposure, W %*% beta.outcome)
  
  # — Iterative update —
  betal <- betawt
  for (iter in seq_len(max.iter)) {
    old <- betal
    
    # Create diagonal penalty matrices using sparse Diagonals
    D1 <- Matrix::Diagonal(x = ascvec[1:K] / (abs(old) + littleeps))
    D2 <- Matrix::Diagonal(x = ascvec[(K + 1):(K + qp)] / (as.numeric(abs(dm %*% old)) + littleeps))
    D3 <- Matrix::Diagonal(x = ascvec[(K + qp + 1):(K + 2 * qp)] / (as.numeric(abs(dp %*% old)) + littleeps))
    
    A <- MVplus + lambda * (D1 + t(dm) %*% D2 %*% dm + t(dp) %*% D3 %*% dp)
    
    # Use solve() if A is well-conditioned; otherwise use MASS::ginv()
    if (rcond(as.matrix(A)) > .Machine$double.eps) {
      betal_new <- as.numeric(solve(A, XWy))
    } else {
      betal_new <- as.numeric(MASS::ginv(as.matrix(A)) %*% XWy)
    }
    
    if (max(abs((betal_new - old) / (abs(old) + littleeps))) < err) {
      betal <- betal_new
      break
    }
    betal <- betal_new
  }
  
  cls <- groupAssignment(betal, digit = digit)
  
  list(
    coefficients = betal,
    cls          = cls,
    cls_chr      = paste(cls, collapse = "-"),
    lambda       = lambda,
    rr           = abs(rr),
    init         = betawt,
    type         = type,
    eps          = eps,
    littleeps    = littleeps,
    iter         = iter
  )
}

SRIVW <- function (beta.exposure, se.exposure, beta.outcome, se.outcome, G, gen_cor = NULL, phi_cand = 0) 
{
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {
    stop("this function is developed for multivariable MR")
  }
  K <- ncol(beta.exposure)
  if (nrow(beta.exposure) != length(beta.outcome)) {
    stop("The number of SNPs in beta.exposure and beta.outcome is different")
  }
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  p <- nrow(beta.exposure)
  W <- diag(se.outcome^(-2))
  Vj <- lapply(1:p, function(j) diag(se.exposure[j, ]) %*% 
                 gen_cor[[j]] %*% diag(se.exposure[j, ]))
  Vj_root_inv <- lapply(1:p, function(j) {
    P_eigen <- eigen(gen_cor[[j]])
    P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% 
      t(P_eigen$vectors)
    P_root_inv %*% diag(1/se.exposure[j, ])
  })
  IV_strength_matrix <- Reduce("+", lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j, ]
    beta.exposure.V %*% t(beta.exposure.V)
  })) - p * diag(K)
  iv_strength_parameter <- min(eigen(IV_strength_matrix/sqrt(p))$values)
  V <- Reduce("+", lapply(1:p, function(j) {
    Vj[[j]] * (se.outcome[j]^(-2))
  }))
  M <- t(beta.exposure) %*% W %*% beta.exposure
  MV <- M - V
  MV_eigvalues <- eigen(MV)$values
  MV_eigen <- eigen(MV)
  if (is.null(phi_cand)) {
    phi_cand <- c(0, exp(seq(0, 17, by = 0.5) - min(iv_strength_parameter)))
  }
  phi_length <- length(phi_cand)
  MV.l.inv.long <- Reduce(rbind, lapply(1:phi_length, function(l) {
    MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_cand[l]/MV_eigvalues)) %*% 
      t(MV_eigen$vectors)
  }))
  beta.est <- MV.l.inv.long %*% t(beta.exposure) %*% W %*% 
    (beta.outcome)
  prof.lik <- sapply(1:phi_length, function(l) {
    beta.hat <- beta.est[(1 + (l - 1) * K):(l * K)]
    bvb <- sapply(1:p, function(j) t(beta.hat) %*% Vj[[j]] %*% 
                    beta.hat)
    tau2 <- 0
    S <- diag(1/(se.outcome^2 + bvb + tau2))
    1/p * t(beta.outcome - beta.exposure %*% beta.hat) %*% 
      S %*% (beta.outcome - beta.exposure %*% beta.hat)
  })
  phi_selected <- phi_cand[which.min(prof.lik)]
  MV.l.inv <- MV_eigen$vectors %*% diag(1/(MV_eigvalues + 
                                             phi_selected/MV_eigvalues)) %*% t(MV_eigen$vectors)
  mvmr.adIVW <- MV.l.inv %*% t(beta.exposure) %*% W %*% 
    beta.outcome
  tau2_adivw <- 0
  adIVW_Vt <- Reduce("+", lapply(1:p, function(j) {
    m <- beta.exposure[j, ] %*% t(beta.exposure[j, ]) * 
      (se.outcome[j]^(-2))
    v <- Vj[[j]] * (se.outcome[j]^(-2))
    bvb <- as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
    vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
    m * (1 + bvb + tau2_adivw * se.outcome[j]^(-2)) + 
      vbbv
  }))
  mvmr.adIVW.se <- sqrt(diag(MV.l.inv %*% adIVW_Vt %*% 
                               MV.l.inv))
  return(list(beta.hat = as.vector(mvmr.adIVW), beta.se = as.vector(mvmr.adIVW.se), 
              iv_strength_parameter = iv_strength_parameter, phi_selected = phi_selected, 
              tau.square = tau2_adivw))
}

groupAssignment <- function(beta, digit = 3) {
  beta <- round(beta, digit)
  tol <- 10^(-digit)
  K <- length(beta)
  groups <- rep(NA, K)
  group_map <- list()  # to store mapping: key (as character) -> group number
  group_counter <- 0
  
  for (k in seq_len(K)) {
    # Use a key based on absolute value. For zero, use "0".
    key <- if (abs(beta[k]) < tol) {
      "0"
    } else {
      # Use a rounded absolute value as a key.
      format(abs(beta[k]), nsmall = digit)
    }
    
    # If this key has been seen before, assign the corresponding group;
    # otherwise, create a new group.
    if (key %in% names(group_map)) {
      groups[k] <- group_map[[key]]
    } else {
      group_counter <- group_counter + 1
      group_map[[key]] <- group_counter
      groups[k] <- group_counter
    }
  }
  
  return(groups)
}

constructG <- function(x, y, include_zero_group = FALSE) {
  if (length(x) != length(y)) {
    stop("`x` and `y` must be the same length.")
  }
  if (all(y==0) & include_zero_group == FALSE) {
    stop("all betas are zero")
  }
  
  # 1. the groups with any non-zero signs
  all_groups     <- unique(x)
  nonzero_groups <- all_groups[sapply(all_groups, function(g) any(y[x == g] != 0))]
  
  K <- length(x)
  L <- length(nonzero_groups)
  
  # 2. allocate L × K zero matrix
  G <- matrix(0L, nrow = L, ncol = K,
              dimnames = list(nonzero_groups, NULL))
  
  # 3. fill in each nonzero-sign group
  for (i in seq_along(nonzero_groups)) {
    g    <- nonzero_groups[i]
    idx  <- which(x == g)
    p_l  <- y[idx[1]]            # sign of the first member
    G[i, idx] <- p_l * y[idx]
  }
  
  # 4. optionally append the "zero-sign" group as row L+1
  if (include_zero_group) {
    zero_idx <- which(y == 0)
    zero_row <- integer(K)
    zero_row[zero_idx] <- 1L
    G <- rbind(G, zero_row)
    # give it a name so you can see it’s the zero-group
    rownames(G)[nrow(G)] <- "zero_signs"
  }
  
  return(G)
}

computeGroupedSE <- function(se.exposure, P, G) {
  # se.exposure : p x K matrix of SNP–exposure SEs
  # P    : K x K exposure–exposure correlation matrix
  # G    : L x K grouping/transformation matrix
  
  p <- nrow(se.exposure)
  J <- nrow(G)
  se_grouped <- matrix(NA_real_, nrow = p, ncol = J)
  Corr_group <- vector(mode = 'list', length = p)
  
  for(i in seq_len(p)) {
    Sigma_i    <- diag(se.exposure[i, ]) %*% P %*% diag(se.exposure[i, ])     # K×K covariance for SNP i
    cov_group  <- G %*% Sigma_i %*% t(G)                        # J×J covariance of grouped exposures
    se_grouped[i, ] <- sqrt(diag(cov_group))
    Corr_group[[i]] <- diag(1 / se_grouped[i,],nrow = J) %*% cov_group %*% diag(1 / se_grouped[i,],nrow = J)
  }
  
  colnames(se_grouped) <- paste0("group", seq_len(J))
  
  return(list(se_grouped = se_grouped,
              Corr_group = Corr_group))
}

computeGroupedData <- function(beta.exposure, se.exposure, P, G) {
  beta.exposure <- beta.exposure %*% t(G)
  out <- computeGroupedSE(se.exposure,P,G)
  se.exposure <- out$se_grouped
  P_grouped <- out$Corr_group
  return(list(beta.exposure = beta.exposure,
              se.exposure = se.exposure,
              P_grouped = P_grouped))
}

SRIVW_group_est <- function(beta.exposure,se.exposure,beta.outcome,se.outcome,P,est_beta,include_zero_group = FALSE, digit = 3) {
  K <- ncol(beta.exposure)
  p <- nrow(beta.exposure)
  # construct transformation matrix for beta-exposure associations
  group <- groupAssignment(est_beta, digit = digit)
  G_collapse <- constructG(x = group, y = sign(round(est_beta,digit)), include_zero_group = include_zero_group)
  L <- nrow(G_collapse)
  stopifnot(L>0)
  
  # Use another data to do inference
  group_dat <- computeGroupedData(beta.exposure,se.exposure,P,G = G_collapse)
  
  if (L == 1) {
    final_fit <- mr.divw(beta.exposure = group_dat$beta.exposure,beta.outcome = beta.outcome,
                         se.exposure = group_dat$se.exposure,se.outcome = se.outcome)
    beta_est <- final_fit$beta.hat %*% G_collapse
    beta_se <- final_fit$beta.se %*% G_collapse
    iv.strength.dt <- final_fit$condition
  } else {
    final_fit <- SRIVW(beta.exposure = group_dat$beta.exposure, beta.outcome = beta.outcome,
                       se.exposure = group_dat$se.exposure, se.outcome = se.outcome, 
                       phi_cand = NULL, gen_cor = group_dat$P_grouped)
    beta_est <- final_fit$beta.hat %*% G_collapse
    beta_se <- abs(final_fit$beta.se %*% G_collapse)
    iv.strength.dt <- final_fit$iv_strength_parameter
  }
  
  # calculate minimum conditional F-stats
  if (L > 1) {
    F.data <- format_mvmr(BXGs = group_dat$beta.exposure,
                          BYG = beta.outcome,
                          seBXGs = group_dat$se.exposure,
                          seBYG = se.outcome,
                          RSID = 1:p)
    
    invisible(capture.output(fres <- strength_mvmr(r_input = F.data, gencov = lapply(1:p, function(j) 
    {diag(group_dat$se.exposure[j,]) %*% group_dat$P_grouped[[j]] %*% diag(group_dat$se.exposure[j,])}))
    ))
    
    condF <- abs(as.numeric(fres) %*% G_collapse)
  } else {
    condF <- NA
  }
  
  res <- NULL
  res$beta_est <- beta_est
  res$beta_se <- abs(beta_se)
  res$iv_strength <- iv.strength.dt
  res$condF <- condF
  return(res)
}

generate_tuning_para <- function(beta.exposure, se.exposure, P, lambda.length) {
  
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  p <- nrow(beta.exposure)
  W <- diag(se.outcome^(-2))
  Vj <- lapply(1:p, function(j) diag(se.exposure[j, ]) %*% 
                 P %*% diag(se.exposure[j, ]))
  Vj_root_inv <- lapply(1:p, function(j) {
    P_eigen <- eigen(P)
    P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% 
      t(P_eigen$vectors)
    P_root_inv %*% diag(1/se.exposure[j, ])
  })
  IV_strength_matrix <- Reduce("+", lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j, ]
    beta.exposure.V %*% t(beta.exposure.V)
  })) - p * diag(K)
  lambda_min <- min(eigen(IV_strength_matrix)$values)
  C <- ifelse(lambda_min > p, (lambda_min/sqrt(lambda_min + p))^(2/3), (p/2)^(1/3))
  exp(seq(log(C), log(C*1e-4), 
          length.out = lambda.length))
  
}

mvmr_pacs_cv_parallel <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P, RR_obs, type = 2, rr_cut_off = 0, fix.cor = FALSE,
                                  betawt = NULL, tau = c(0.5,1,2,3), CV_fold = 5, n_times = 1, lambda.length = 20, 
                                  eps = 1e-4, digit = 3, seed = NULL, lambda.1se = TRUE) {
  # Define the fraction for CV splits
  e <- 1/CV_fold
  
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  
  p <- nrow(beta.exposure)
  K <- ncol(beta.exposure)
  
  iv_strength_parameter <- mvmr.ivw(beta.exposure = beta.exposure,
                                    se.exposure = se.exposure,
                                    beta.outcome = beta.outcome,
                                    se.outcome = se.outcome,
                                    gen_cor = P)$iv_strength_parameter
  
  if (is.null(betawt)) {
    # obtain initial beta estimates for subsequent estimators
    betawt <- obtain_initial(beta.exposure, se.exposure, beta.outcome, se.outcome,
                             iv_strength_parameter = iv_strength_parameter)
  }
  
  # Fit dIVW PACS: create weight matrix and grid of tuning parameters
  W <- diag(1 / se.outcome^2)
  lambda_cand <- generate_tuning_para(beta.exposure, se.exposure, P, lambda.length = lambda.length)
  grid <- expand.grid(lambda = lambda_cand, tau = tau)
  
  # Prepare array to store CV losses
  obj_f_pacs <- array(0, dim = c(n_times, CV_fold, length(lambda_cand), length(tau)))
  
  # Set up parallel backend (using available cores minus one)
  numCores <- parallel::detectCores() - 1
  cl <- makeCluster(numCores)
  
  # Source your custom functions on each worker.
  clusterEvalQ(cl, source("PACS_funs_091125.R"))
  
  registerDoParallel(cl)
  
  for (m in 1:n_times) {
    cat('CV times = ', m, "\n")
    # Compute Sig once per m (independent of folds)
    Sig <- array(0, dim = c(p, K, K))
    for (j in 1:p) {
      Sig[j, , ] <- diag(se.exposure[j, ]) %*% P %*% diag(se.exposure[j, ])
    }
    if (!is.null(seed)) set.seed(seed + m)
    # Generate the CV data splits
    dt.exposure <- datathin(beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold)
    dt.outcome  <- datathin(beta.outcome, family = "normal", arg = se.outcome^2, K = CV_fold)
    
    # Parallelize over CV folds
    cv_results <- foreach(l = 1:CV_fold, .packages = c("stats","BDcocolasso","mr.divw","expm",
                                                       "MASS","datathin")) %dopar% {
                                                         cat('Processing CV fold = ', l, "\n")
                                                         # Training data
                                                         beta.exposure.train <- apply(dt.exposure[, , -l], c(1, 2), sum)
                                                         se.exposure.train   <- se.exposure * sqrt(1 - e)
                                                         beta.outcome.train  <- apply(dt.outcome[, 1, -l], 1, sum)
                                                         se.outcome.train    <- se.outcome * sqrt(1 - e)
                                                         
                                                         if (fix.cor) {
                                                           rr.train <- RR_obs
                                                         } else {
                                                           rr.train <- cor(beta.exposure.train)
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
                                                         MV_plus.test <- ADMM_proj(MV.test / sqrt(p), epsilon = eps)$mat * sqrt(p)
                                  
                                                         # Loop over grid of tuning parameters to compute the loss
                                                         losses_matrix <- matrix(0, nrow = length(lambda_cand), ncol = length(tau))
                                                         for (i in 1:nrow(grid)) {
                                                           lambda_val <- grid[i, "lambda"]
                                                           tau_val    <- grid[i, "tau"]
                                                           beta_hat   <- tryCatch({dIVW_PACS_cluster(beta.exposure = beta.exposure.train,
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
                                                                                                     digit         = digit)$coefficients
                                                           },error = function(e) {return(rep(NA,K))})
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
  stopCluster(cl)
  
  # Average the loss over CV folds and repetitions
  res_mean <- apply(obj_f_pacs, c(3, 4), mean, na.rm = TRUE)
  res_sd   <- apply(obj_f_pacs, c(3, 4), sd, na.rm = TRUE) / sqrt(CV_fold * n_times)
  
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
  res$cluster             <- cluster.pacs
  res$iv_strength_parameter <- iv_strength_parameter
  res$lambda.index        <- best_idx[1]
  res$lambda              <- lambda_cand[best_idx[1]]
  res$tau.index           <- best_idx[2]
  res$tau                 <- tau[best_idx[2]]
  res$beta.init           <- betawt
  res$best_idx <- best_idx
  
  set.seed(NULL)
  
  return(res)
}


DT_PACS_dIVW_parallel <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P, RR_obs, cluster = NULL, tau = c(0.5,1,2,3), fold1 = 2, fold2 = 5, 
                                  re = 10, n_times = 1, rr_cut_off = 0, type = 2, epsilon = 1e-4, digit = 3, include_zero_group = FALSE, iv_strength_cut = NULL, lambda.1se = TRUE) {
  
  # default is splitting the data into 10 folds, using half to pick grouping, the other half to perform estimation
  e1 <- 1/fold1
  e2 <- 1/fold2
  
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
    
    dt.exposure <- datathin(beta.exposure, family="mvnormal", arg = Sig, K = fold1)
    
    dt.outcome <- datathin(beta.outcome, family = 'normal', arg = se.outcome^2, K = fold1)
    
    # using the first 5 folds to perform grouping selection
    # data based on the first half of the data
    beta.exposure.train <- apply(dt.exposure[, , -fold1], c(1, 2), sum)
    se.exposure.train <- se.exposure * sqrt(1-e1)
    beta.outcome.train <- apply(dt.outcome[, 1, -fold1, drop = FALSE], 1, sum) %>% as.vector(.)
    se.outcome.train <- se.outcome * sqrt(1-e1)
    iv_str_train <- mvmr.ivw(beta.exposure = beta.exposure.train,se.exposure = se.exposure.train,beta.outcome = beta.outcome.train,se.outcome = se.outcome.train,gen_cor = P)$iv_strength_parameter
    
    # data based on the second half of the data for inference
    beta.exposure.test <- dt.exposure[,,fold1]
    se.exposure.test <- se.exposure * sqrt(e1)
    beta.outcome.test <- dt.outcome[,1,fold1]
    se.outcome.test <- se.outcome * sqrt(e1)
    
    # initial estimates
    betawt <- obtain_initial(beta.exposure.train,se.exposure.train,beta.outcome.train,se.outcome.train,iv_strength_parameter = iv_str_train, n_times = 5, P = P, lambda.1se = lambda.1se)
    
    # perform PACS on the training set
    cat("Perform PACS on the training set to obtain grouping \n")
    
      
    pacs_res <- mvmr_pacs_cv_parallel(beta.exposure = beta.exposure.train, 
                                        se.exposure = se.exposure.train, 
                                        beta.outcome = beta.outcome.train, 
                                        se.outcome = se.outcome.train, 
                                        P = P, 
                                        type = type, 
                                        rr_cut_off = rr_cut_off, 
                                        RR_obs = RR_obs, 
                                        betawt = betawt, 
                                        CV_fold = fold2, 
                                        tau = tau, 
                                        eps = epsilon,
                                        digit = digit, 
                                        n_times = n_times,
                                        lambda.1se = lambda.1se) 
    
    if (all(round(pacs_res$beta,digit) == 0)) {pacs_res$beta <- rep(10^(-digit),K)} # this consider all risk factors as a group even if they jointly have no effect
    
    grouping_hist[r] <- paste(pacs_res$cluster, collapse = "-")
    
    res_dt1 <- SRIVW_group_est(beta.exposure = beta.exposure.train,se.exposure = se.exposure.train,
                               beta.outcome = beta.outcome.train,se.outcome = se.outcome.train,
                               P = P, est_beta = pacs_res$beta, digit = digit, include_zero_group = include_zero_group)
    
    iv.hist.dt1[r] <- res_dt1$iv_strength
    
    condF.dt1[r,] <- res_dt1$condF
    
    # estimation on the other data
    cat("Using SRIVW to obtain causal effect estimates of grouped risk factors \n")
    
    res_dt2 <- SRIVW_group_est(beta.exposure = beta.exposure.test,se.exposure = se.exposure.test,
                               beta.outcome = beta.outcome.test,se.outcome = se.outcome.test,
                               P = P, est_beta = pacs_res$beta, digit = digit, include_zero_group = include_zero_group)
    
    beta_est[r,] <- res_dt2$beta_est
    
    beta_se[r,] <- res_dt2$beta_se
    
    iv.hist.dt2[r] <- res_dt2$iv_strength
    
    condF.dt2[r,] <- res_dt2$condF
    
  }
  
  # grouping counts
  group_counts <- table(grouping_hist)
  # IV strength in the training set
  avg_iv_dt1 <- tapply(iv.hist.dt1, grouping_hist, mean, na.rm = TRUE)[names(group_counts)]
  # final grouping
  if (is.null(iv_strength_cut)) {
    final_cls <- names(which.max(group_counts))
  } else {
    final_cls <- names(which.max(group_counts[avg_iv_dt1 > iv_strength_cut]))
  }
  
  # output
  res <- NULL
  res$major_group <- final_cls
  res$beta <- apply(matrix(beta_est[grouping_hist == final_cls,],ncol = K), 2, median, na.rm = TRUE)
  res$se<- apply(matrix(beta_se[grouping_hist == final_cls,],ncol = K), 2, median, na.rm = TRUE)
  res$iv.hist.dt1 <- avg_iv_dt1
  res$iv.hist.dt2 <- iv.hist.dt2
  res$condF.dt1 <- apply(matrix(condF.dt1[grouping_hist == final_cls,],ncol = K), 2, median, na.rm = TRUE)
  res$condF.dt2 <- apply(matrix(condF.dt2[grouping_hist == final_cls,],ncol = K), 2, median, na.rm = TRUE)
  res$grouping_hist <- grouping_hist
  res$beta.hist <- beta_est
  res$se.hist <- beta_se
  return(res)
}

mvmr_ridge <- function(beta.exposure, se.exposure,
                       beta.outcome,  se.outcome,
                       P, phi, epsilon = 1e-4) {
  # Dimensions and basic checks
  p <- nrow(beta.exposure); K <- ncol(beta.exposure)
  stopifnot(nrow(se.exposure) == p, ncol(se.exposure) == K)
  stopifnot(length(beta.outcome) == p, length(se.outcome) == p)
  stopifnot(is.matrix(P), nrow(P) == K, ncol(P) == K)
  if (any(se.outcome <= 0) || any(se.exposure <= 0)) {
    stop("Non-positive standard errors detected.")
  }
  if (phi < 0) stop("phi must be nonnegative.")
  
  # Coerce shapes and symmetrize P
  X <- beta.exposure
  y <- as.numeric(beta.outcome)
  
  # Outcome weights
  w <- 1 / (se.outcome^2)                    # length p
  
  XtW <- t(X) * rep(w, each = K)             # K x p
  M   <- XtW %*% X                           # K x K
  
  # V = sum_j diag(se_xj) P diag(se_xj) * w_j
  V <- matrix(0.0, K, K)
  for (j in 1:p) {
    Dj <- diag(se.exposure[j,], K, K)
    V  <- V + w[j] * (Dj %*% P %*% Dj)
  }
  
  # MV = M - V (symmetrize for safety)
  MV <- M - V
  
  # PSD projection with scaling for stability
  MVs      <- MV / sqrt(p)
  MV_plusS <- ADMM_proj(MVs, epsilon = epsilon)$mat
  MV_plus  <- MV_plusS * sqrt(p)   # re-symmetrize + scale back
  
  # Ridge system: A b = rhs
  A   <- MV_plus + (phi * diag(K))
  rhs <- as.numeric(XtW %*% y)               # X' W y
  
  b <- solve(A) %*% rhs
  
  # Ensure vector output
  as.numeric(b)
}


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
    
    dt.exposure <- datathin(beta.exposure, family = "mvnormal", arg = Sig, K = CV_fold)
    dt.outcome  <- datathin(beta.outcome,  family = "normal",   arg = se.outcome^2, K = CV_fold)
    
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
      MV_plus.test <- ADMM_proj(MV.test/sqrt(p), epsilon = epsilon)$mat * sqrt(p)
      
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
        obj_f_ridge[l + (m-1)*CV_fold, jphi] <- drop(crossprod(b, MV_plus.test %*% b)) - 2 * sum(c.test * b)
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