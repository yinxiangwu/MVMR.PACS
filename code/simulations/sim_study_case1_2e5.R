task_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# simulation study for MVMR-PACS
# Create date: 090924
# Number of risk factors = 10, two clusters

# load required packages and functions
library(mvtnorm)
library(CVXR)
library(BDcocolasso)
library(magrittr)
library(glmnet)
library(mr.divw)
library(expm)
library(MVMR)
library(foreach)
library(MASS)
library(datathin)
library(doParallel)
library(MVMRcMLSuSiE)
library(MRcML)
library(MVMRcML)
source("~/MR PACS/my functions/summary_mvMR_SSS.R")
source("~/MR PACS/my functions/summary_mvMR_BF.R")
source("~/MR PACS/my functions/data_gen_individual_test.R")
source("~/MR PACS/my functions/PACS_funs_091125.R")
source("~/MR PACS/my functions/F_stats_calculator.R")
source("~/MR PACS/my functions/mvMR_glmnet.R")

# simulation study case 1 - vary true beta0

# factors to consider

# number of snps

p = 500

# number of risk factors

K = 10

# true beta

beta0 = c(1,1,1,0,0,0,0,0,0.5,0)

# sample size of individual-level data

N = 2e5

# number of simulations

n_rep = 1

# sigma of gamma

sigma_gamma = 0.001

# sigma of confounder u 

sigma_u = 2

# sigma of residual error

sigma_e = 1

epsilon = 1e-4

# vectors to store the results

res_ivw = matrix(NA, n_rep, ncol = K)
res_ivw_inclusion = matrix(NA, n_rep, ncol = K)
res_srivw = matrix(NA, n_rep, ncol = K)
res_srivw_inclusion = matrix(NA, n_rep, ncol = K)
res_ivw_lasso = matrix(NA, n_rep, ncol = K)
res_dlasso = matrix(NA, n_rep, ncol = K)
res_pacs_0.8 = matrix(NA, n_rep, ncol = K)
res_pacs_v2 = matrix(NA, n_rep, ncol = K)
res_pacs_v2_1se = matrix(NA, n_rep, ncol = K)
res_pacs_v2_gen = matrix(NA, n_rep, ncol = K)
res_mrbma = matrix(NA, n_rep, ncol = K)
mrbma_mip = matrix(NA, n_rep, ncol = K)
susie_num_cluster = rep(NA, n_rep)
res_susie = matrix(NA, n_rep, ncol = K)
res_iv_strength = rep(NA, n_rep)
res_cond_F = matrix(NA, n_rep, ncol = K)
cluster_pacs_0.8 = matrix(NA, n_rep, ncol = K)
cluster_pacs_v2 = matrix(NA, n_rep, ncol = K)
cluster_pacs_v2_1se = matrix(NA, n_rep, ncol = K)
cluster_pacs_v2_gen = matrix(NA, n_rep, ncol = K)

for (i in 1:n_rep) {
  tryCatch({
    
    # generate data
    
    data = data_gen_individual(beta0 = beta0, sigma_gamma = sigma_gamma, sigma_u = sigma_u, sigma_e = sigma_e, n = N)
    
    full_df <- data$full_df
    
    beta.exposure <- as.matrix(full_df[,1:10])
    
    se.exposure <- as.matrix(full_df[,11:20])
    
    beta.outcome <- as.numeric(full_df$beta.outcome)
    
    se.outcome <- as.numeric(full_df$se.outcome)
    
    P <- as.matrix(data$P)
    
    gen_cor <- cor(data$gamma_coef)
    
    # IVW estimator
    
    ivw.obj <- mvmr.ivw(beta.exposure,se.exposure = se.exposure,beta.outcome = beta.outcome,se.outcome = se.outcome,gen_cor = P)
    
    res_ivw[i,] <- ivw.obj$beta.hat
    
    res_ivw_inclusion[i,] <- (2 * pnorm(abs(ivw.obj$beta.hat/ivw.obj$beta.se),lower.tail = FALSE)) < 0.05/K
    
    res_iv_strength[i] <- ivw.obj$iv_strength_parameter
    
    # SRIVW estimator
    
    srivw.obj <- mvmr.srivw(beta.exposure,se.exposure = se.exposure,beta.outcome = beta.outcome,se.outcome = se.outcome,gen_cor = P,phi_cand = NULL)
    
    res_srivw[i,] <- srivw.obj$beta.hat
    
    res_srivw_inclusion[i,] <- (2 * pnorm(abs(srivw.obj$beta.hat/srivw.obj$beta.se),lower.tail = FALSE)) < 0.05/K
    
    # IVW + LASSO penalty, use existing code from the author of MRBMA https://github.com/verena-zuber/demo_AMD
    
    ivw_lasso_input=new("mvMRInput", betaX = as.matrix(beta.exposure/se.outcome), betaY = as.matrix(beta.outcome/se.outcome), snps=as.character(1:p), exposure=as.character(1:K))
    
    lasso=glmnet_l1_mr(ivw_lasso_input, cv=TRUE, cv.param="lambda.1se")			
    
    res_ivw_lasso[i,] <- lasso@Estimate
    
    # obtain initial beta estimates for subsequent estimators
    
    betawt <- obtain_initial(beta.exposure,se.exposure,beta.outcome,se.outcome,iv_strength_parameter = ivw.obj$iv_strength_parameter, epsilon = epsilon, P=P, n_times = 3, lambda.1se = TRUE)
    
    # dIVW + LASSO penalty
    
    dlasso <- mvmr_pacs_cv_parallel(beta.exposure = beta.exposure, se.exposure = se.exposure, beta.outcome = beta.outcome, se.outcome = se.outcome, 
                                    P = P, RR_obs = cor(beta.exposure), tau = c(1,2,3), type = 3, rr_cut_off = 1, n_times = 3,
                                    betawt = betawt, eps = epsilon, lambda.1se = TRUE)
    
    res_dlasso[i,] = dlasso$beta
    
    # PACS
    
    pacs_0.8 <- mvmr_pacs_cv_parallel(beta.exposure = beta.exposure, se.exposure = se.exposure, beta.outcome = beta.outcome, se.outcome = se.outcome, 
                                      P = P, RR_obs = cor(beta.exposure), tau = c(1,2,3), type = 4, rr_cut_off = 0.8, n_times = 3,
                                      betawt = betawt, eps = epsilon, lambda.1se = TRUE)
    
    res_pacs_0.8[i,] <- pacs_0.8$beta
    
    cluster_pacs_0.8[i,] <- pacs_0.8$cluster
    
    # PACS
    
    pacs_v2 <- mvmr_pacs_cv_parallel(beta.exposure = beta.exposure, se.exposure = se.exposure, beta.outcome = beta.outcome, se.outcome = se.outcome, 
                                     P = P, RR_obs = cor(beta.exposure), tau = c(1,2,3), type = 2, n_times = 3,
                                     betawt = betawt, eps = epsilon, lambda.1se = TRUE)
    
    res_pacs_v2[i,] <- pacs_v2$beta
    
    cluster_pacs_v2[i,] <- pacs_v2$cluster
    
    # MR-BMA
    
    betaX_ivw = beta.exposure / se.outcome
    
    betaY_ivw = beta.outcome / se.outcome
    
    mrbma_input=new("mvMRInput", betaX = as.matrix(betaX_ivw), betaY = as.matrix(betaY_ivw), snps=as.character(1:p), exposure=as.character(1:K), outcome = "amd")
    
    BMA_output=summarymvMR_SSS(mrbma_input,kmin=1,kmax=10, prior_prob=0.5, max_iter=100000)
    
    mr.bma.out = sss.report.mr.bma(BMA_output, top = 10, write.out = FALSE)
    
    res_mrbma[i,] <- as.numeric(mr.bma.out[order(as.numeric(mr.bma.out[,1])),3])
    
    mrbma_mip[i,] <- as.numeric(mr.bma.out[order(as.numeric(mr.bma.out[,1])),2])
    
    # SuSiE
    
    pval.exp <- 2 * pnorm(abs(beta.exposure/se.exposure), lower.tail = FALSE)
    
    tryCatch({
      
      # step 1 screening using univariable MR-cML
      
      p.vec <- rep(NA, K)
      
      for (j in 1:K) {
        
        ind <- pval.exp[,j] < 0.05/p
        
        p.vec[j] <- mr_cML(as.vector(beta.exposure[ind,j]),beta.outcome[ind],se_exp = as.vector(se.exposure[ind,j]),se_out = se.outcome[ind],n = N)$MA_BIC_p
        
      }
      
      subset.idx <- which(p.vec < 0.05 / K)
      
      K_sub <- length(subset.idx)
      
      ind <- apply(pval.exp[,subset.idx], 1, function(x) any(x < 0.05/p/K_sub))
      
      beta.exposure.mat <- beta.exposure[ind,subset.idx]
      
      se.exposure.mat <- se.exposure[ind,subset.idx]
      
      beta.outcome.vec <- beta.outcome[ind]
      
      se.outcome.vec <- se.outcome[ind]
      
      pval.exposure.mat <- 2*pnorm(abs(beta.exposure[ind,subset.idx]/se.exposure[ind,subset.idx]),lower.tail = FALSE)
      
      # step 2 obtain initial estimates using MVMR-cML
      
      step2.res <- mvmr.cml.susie.step2(sample.sizes.subset = rep(N, K_sub+1), beta.exposure.mat = beta.exposure.mat, se.exposure.mat = se.exposure.mat, beta.outcome.vec = beta.outcome.vec, se.outcome.vec = se.outcome.vec, pval.exposure.mat = pval.exposure.mat, use.openGWAS = FALSE)
      
      rho.mat <- matrix(0, K+1, K+1)
      
      rho.mat[1:K,1:K] <- P
      
      # step 3 susie
      
      susie.res <- mvmr.cml.susie.step3(step2.res$mvdat, step2.res$invalid.idx, step2.res$theta.vec, rho.mat[c(subset.idx,K+1),c(subset.idx,K+1)])$alpha
      
      susie_num_cluster[i] <- sum(apply(susie.res,1,function(x) any(x>1/K_sub+.Machine$double.eps)))
      
      tmp.vec <- subset.idx
      
      susie.include <- rep(0,K)
      
      susie.include[unique(tmp.vec[Reduce(c,apply(susie.res,1,function(x) which(x>1/K_sub+.Machine$double.eps)))])] <- 1
      
      res_susie[i,] <- susie.include
      
    }, error = function(e) {})
    
    # Conditional F-stats
    F.data <- format_mvmr(BXGs = beta.exposure,
                          BYG = beta.outcome,
                          seBXGs = se.exposure,
                          seBYG = se.outcome,
                          RSID = 1:p)
    
    fres <- strength_mvmr(r_input = F.data, gencov = lapply(1:p, function(j) 
    {diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]) }))
    
    res_cond_F[i, ] <- as.numeric(fres)
  },error = function(e) {print(e)})
}
write.csv(res_ivw,paste0("Res2e5/res_ivw_",task_id,".csv"),row.names = FALSE)
write.csv(res_ivw_inclusion,paste0("Res2e5/res_ivw_inclusion_",task_id,".csv"),row.names = FALSE)
write.csv(res_srivw,paste0("Res2e5/res_srivw_",task_id,".csv"),row.names = FALSE)
write.csv(res_srivw_inclusion,paste0("Res2e5/res_srivw_inclusion_",task_id,".csv"),row.names = FALSE)
write.csv(res_ivw_lasso,paste0("Res2e5/res_ivw_lasso_",task_id,".csv"),row.names = FALSE)
write.csv(res_dlasso,paste0("Res2e5/res_dlasso_",task_id,".csv"),row.names = FALSE)
write.csv(res_pacs_0.8,paste0("Res2e5/res_pacs_0.8_",task_id,".csv"),row.names = FALSE)
write.csv(res_pacs_v2,paste0("Res2e5/res_pacs_v2_",task_id,".csv"),row.names = FALSE)

write.csv(res_mrbma,paste0("Res2e5/res_mrbma_",task_id,".csv"),row.names = FALSE)
write.csv(mrbma_mip,paste0("Res2e5/mrbma_mip_",task_id,".csv"),row.names = FALSE)
write.csv(res_susie,paste0("Res2e5/res_susie_",task_id,".csv"),row.names = FALSE)
write.csv(susie_num_cluster,paste0("Res2e5/susie_num_cluster_",task_id,".csv"),row.names = FALSE)
write.csv(res_iv_strength,paste0("Res2e5/res_iv_strength_",task_id,".csv"),row.names = FALSE)
write.csv(res_cond_F,paste0("Res2e5/res_cond_F_",task_id,".csv"),row.names = FALSE)

write.csv(cluster_pacs_0.8,paste0("Res2e5/cluster_pacs_0.8_",task_id,".csv"),row.names = FALSE)
write.csv(cluster_pacs_v2,paste0("Res2e5/cluster_pacs_v2_",task_id,".csv"),row.names = FALSE)





