library(magrittr)
library(mr.divw)
library(corrplot)
library(BDcocolasso)
library(glmnet)
library(MVMR)
library(datathin)
library(MVMRcMLSuSiE)
library(MRcML)
library(MVMRcML)
library(doParallel)
library(foreach)
library(MVMR.PACS)
setwd('~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/code/my functions/')
link_to_dat <- '~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/lipids_processed_data/'
source('mvMR_glmnet.R')
source("summary_mvMR_BF.R")
source("summary_mvMR_SSS.R")

dat <- read.csv(paste0(link_to_dat,'lipids_total32_5e-08.csv'))
P <- read.csv(paste0(link_to_dat,'lipids_total32_5e-08_cor_mat.csv'))[1:32,1:32] %>% as.matrix(.)
traits <- c('XS-VLDL-PL','XS-VLDL-L','LDL-D','L-LDL-PL','L-LDL-P',
            'L-LDL-L','L-LDL-FC','L-LDL-C','L-LDL-CE','M-LDL-P',
            'M-LDL-L','M-LDL-CE','S-LDL-P','S-LDL-L','S-LDL-C',
            'IDL-L','IDL-FC','IDL-C','HDL-D','XL-HDL-TG',
            'L-HDL-PL','M-HDL-PL','M-HDL-P','M-HDL-L','M-HDL-CE',
            'S-HDL-P','S-HDL-L','ApoA1','ApoB','HDL-C','LDL-C','TG')

beta.exposure <- dat[,paste0('gamma_exp',1:32)] %>% as.matrix(.)
se.exposure <- dat[,paste0('se_exp',1:32)] %>% as.matrix(.)
beta.outcome <- dat$gamma_out1
se.outcome <- dat$se_out1
p <- nrow(beta.exposure)
K <- 32

# calculate IV strength

F.data <- format_mvmr(BXGs = beta.exposure,
                      BYG = beta.outcome,
                      seBXGs = se.exposure,
                      seBYG = se.outcome,
                      RSID = dat$SNP)
fres <- strength_mvmr(r_input = F.data, gencov = lapply(1:p, function(j)
{diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]) }))

# correlation matrix of SNP-exposure associations

cor_mat <- cor(beta.exposure)
rownames(cor_mat) <- traits
colnames(cor_mat) <- traits
jpeg("~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/code/table figure/cor mat/lipids32_cor_mat.jpeg",height = 4000, width = 4000, res = 600)
corrplot(cor_mat,order = 'hclust')
dev.off()

# Analysis

# MVMR-IVW

ivw.obj <- mvmr.ivw(beta.exposure,se.exposure = se.exposure,beta.outcome = beta.outcome,se.outcome = se.outcome,gen_cor = P)
ivw.res <- ivw.obj$beta.hat

# SRIVW

srivw.obj <- mvmr.srivw(beta.exposure,se.exposure = se.exposure,beta.outcome = beta.outcome,se.outcome = se.outcome,gen_cor = P,phi_cand = NULL)
srivw.res <- srivw.obj$beta.hat

# IVW-LASSO

set.seed(1234)

ivw_lasso_input=new("mvMRInput", betaX = as.matrix(beta.exposure/se.outcome), betaY = as.matrix(beta.outcome/se.outcome), snps=as.character(1:p), exposure=as.character(1:K))

lasso=glmnet_l1_mr(ivw_lasso_input, cv=TRUE)

beta.lasso.res <- lasso@Estimate

lasso@Estimate

# Initial estimates

eps <- 1e-4

# MVMR-dLASSO

dlasso <- mvmr.pacs(beta.exposure = beta.exposure, se.exposure = se.exposure, beta.outcome = beta.outcome, se.outcome = se.outcome, P = P, type = 3, rr_cut_off = 1, tau = c(0.5, 1, 2, 3), eps = 1e-4, n_times = 5,lambda.1se = TRUE, lambda.length = 30, seed = 1234)

beta.dlasso.res <- dlasso$beta

pacs <- mvmr.pacs(beta.exposure = beta.exposure, se.exposure = se.exposure, beta.outcome = beta.outcome, se.outcome = se.outcome, P = P, type = 2, tau = c(0.5, 1, 2, 3), eps = 1e-4, n_times = 5, seed = 1234, lambda.length = 30, lambda.1se = TRUE)

beta.pacs.res <- pacs$beta

# MRMBA
betaX_ivw = beta.exposure / se.outcome
betaY_ivw = beta.outcome / se.outcome

mrbma_input=new("mvMRInput", betaX = as.matrix(betaX_ivw), betaY = as.matrix(betaY_ivw), snps=as.character(1:p), exposure=as.character(1:K), outcome = "amd")
BMA_output=summarymvMR_SSS(mrbma_input,prior_prob=0.5, max_iter=10000)
mr.bma.out = sss.report.mr.bma(BMA_output, top = 32, write.out = FALSE)
beta.mrbma.res <- as.numeric(mr.bma.out[order(as.numeric(mr.bma.out[,1])),3])
# marginal inclusion probability
mr.bma.out = sss.report.mr.bma(BMA_output, top = 32, write.out = FALSE)
mr.bma.out[,1] <- traits[as.numeric(mr.bma.out[,1])]

best.model.out = sss.report.best.model(BMA_output, prior_sigma=0.5, top = 10, write.out = FALSE, csv.file.name="amd_best_10models_n145")
print(best.model.out[,-1])
write.csv(print(best.model.out[,-1]),'Res/fig 5e-08/mr.bma.top10.csv')

# summarizing results in a Table
tbl_res <- rbind(ivw.res[,1], srivw.res[,1], beta.lasso.res, beta.dlasso.res,
                 beta.pacs.res,beta.mrbma.res)

rownames(tbl_res) <- c('IVW','adIVW','IVW-LASSO','MVMR-dLASSO','MVMR-PACS','MR-BMA')

colnames(tbl_res) <- traits

print(round(t(tbl_res),3))

write.csv(round(t(tbl_res),3),'~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/code/real data application/main_analysis_090925.csv')

### MVMR-PACS with data thinning for post-selection inference

res_dt <- mvmr.pacs.datathin(beta.exposure = beta.exposure,
                             se.exposure = se.exposure,
                             beta.outcome = beta.outcome,
                             se.outcome = se.outcome,
                             P = P,
                             tau = c(0.5,1,2,3),
                             n_times = 1,
                             re = 100,
                             rr_cut_off = 0,
                             type = 2,
                             epsilon = 1e-4,
                             lambda.1se = TRUE)

