### Code to perform MVMR-cML-SuSiE on lipids data
.libPaths(c("/home/CHSCC/jeco/R/x86_64-redhat-linux-gnu-library/4.1","/chru/analysis/yinxiangwu/R/x86_64-redhat-linux-gnu-library/4.3"))

library(MVMRcMLSuSiE)
library(TwoSampleMR)
library(MRcML)
library(MVMRcML)
library(data.table)
library(genetics.binaRies)
library(unixtools)
library(GRAPPLE)
set.tempdir('/chru/analysis/yinxiangwu/tmp')

traits <- c('XS-VLDL-PL','XS-VLDL-L','LDL-D','L-LDL-PL','L-LDL-P',
            'L-LDL-L','L-LDL-FC','L-LDL-C','L-LDL-CE','M-LDL-P',
            'M-LDL-L','M-LDL-CE','S-LDL-P','S-LDL-L','S-LDL-C',
            'IDL-L','IDL-FC','IDL-C','HDL-D','XL-HDL-TG',
            'L-HDL-PL','M-HDL-PL','M-HDL-P','M-HDL-L','M-HDL-CE',
            'S-HDL-P','S-HDL-L','ApoA1','ApoB','HDL_C','LDL_C','TG')
traits <- gsub("-", "_", traits)

### Step 1

sample.size.vec <- c(rep(115082,27),rep(439214,2),rep(1.32*1e6,3), 185000)

### prepare data for univariable screening

beta.exposure.ls <- vector('list',32)
se.exposure.ls <- vector('list',32)
beta.outcome.ls <- vector('list',32)
se.outcome.ls <- vector('list',32)

sel.files <- c(paste0('Richardson/Processed/',traits[1:27],'_exp.csv'),
                           'IEU/ApoA1_exp.csv',
                           'IEU/ApoB_exp.csv',
                           'GLGC/HDL_C_exp.csv',
                           'GLGC/LDL_C_exp.csv',
                           'GLGC/TG_exp.csv')

sel.files <- paste0("/chru/analysis/yinxiangwu/lipids/data/", sel.files)

tmp_out <- fread('/chru/analysis/yinxiangwu/lipids/data/CAD/cardiogramplusc4d_ukbb_cad.csv')

tmp_out <- format_data(as.data.frame(tmp_out), type = 'outcome')

bfile = '/chru/analysis/pipeline/MR/MR_forward_ARICexposure_IEUoutcome/reference_files/1000G/EUR'  ## https://github.com/MRCIEU/gwasvcf

for (i in 1:32) {
  
  tmp_exp <- fread(sel.files[i])

  tmp_exp <- format_data(as.data.frame(tmp_exp))
    
  clump_exp <- clump_data(tmp_exp,clump_p1 = 5e-8,clump_kb = 10000,clump_r2 = 0.001,bfile = bfile,plink_bin = genetics.binaRies::get_plink_binary())

  tmp_exp_clumped <- tmp_exp[tmp_exp$SNP %in% clump_exp$SNP, ]
  
  harmonized <- harmonise_data(
    tmp_exp_clumped, 
    tmp_out, 
    action = 2)
  
  harmonized <- harmonized[harmonized$mr_keep == TRUE,]
  
  fwrite(harmonized, paste0('/chru_subs/analysis_wuyx/lipids/For MVMR-cML-SuSiE/',traits[i],'_univariable_MR.csv'), row.names= FALSE)
  
  beta.exposure.ls[[i]] <- harmonized$beta.exposure
  se.exposure.ls[[i]] <- harmonized$se.exposure
  beta.outcome.ls[[i]] <- harmonized$beta.outcome
  se.outcome.ls[[i]] <- harmonized$se.outcome
  
  rm(tmp_exp)
    
}

step1.res <- mvmr.cml.susie.step1(sample.sizes = sample.size.vec, beta.exposure.ls = beta.exposure.ls, se.exposure.ls = se.exposure.ls, beta.outcome.ls = beta.outcome.ls, se.outcome.ls = se.outcome.ls, use.openGWAS = FALSE)

### Step 2

subset.idx <- which(step1.res < 0.05 / 32) # remove 22, 26, 27
sample.sizes.subset <- sample.size.vec[subset.idx]
sample.sizes.subset <- c(sample.sizes.subset, 185000)

######## MVMR all covariates ########
sel.files2 = exp.files2 <- c(paste0('Richardson/Processed/',traits[1:27],'_exp.csv'),
                           'IEU/ApoA1_exp.csv',
                           'IEU/ApoB_exp.csv',
                           'GLGC/HDL_C_exp.csv',
                           'GLGC/LDL_C_exp.csv',
                           'GLGC/TG_exp.csv')
sel.files2 <- paste0("/chru/analysis/yinxiangwu/lipids/data/", sel.files2)[subset.idx]
exp.files2 <- paste0("/chru/analysis/yinxiangwu/lipids/data/", exp.files2)[subset.idx]
out.file <- '/chru/analysis/yinxiangwu/lipids/data/CAD/cardiogramplusc4d_ukbb_cad.csv'

# p-value threshold 5e-8 # 592 independent SNPs
tmp <- getInput(sel.files2, exp.files2, out.file,
                plink_exe = genetics.binaRies::get_plink_binary(), plink_refdat = bfile, 
                max.p.thres = 5e-8, cal.cor = TRUE, get.marker.candidates = FALSE)
write.csv(tmp$data,file = "/chru/analysis/yinxiangwu/lipids/For MVMR-cML-SuSiE/lipids_total29_5e-08.csv",row.names = FALSE)
rownames(tmp$cor.mat) <- c(traits[subset.idx],'CAD')
colnames(tmp$cor.mat) <- c(traits[subset.idx],'CAD')
write.csv(tmp$cor.mat,file = "/chru/analysis/yinxiangwu/lipids/processed_data/lipids_total29_5e-08_cor_mat.csv",row.names = FALSE)


beta.exposure.mat <- tmp$data[,paste0('gamma_exp',1:29)] %>% as.matrix(.)
se.exposure.mat <- tmp$data[,paste0('se_exp',1:29)] %>% as.matrix(.)
beta.outcome.vec <- tmp$data[,paste0('gamma_out',1)]
se.outcome.vec <- tmp$data[,paste0('se_out',1)]
pval.exposure.mat <- 2 * pnorm(abs(beta.exposure.mat/se.exposure.mat), lower.tail = FALSE)

step2.res <- mvmr.cml.susie.step2(sample.sizes.subset = sample.sizes.subset, beta.exposure.mat = beta.exposure.mat, se.exposure.mat = se.exposure.mat, beta.outcome.vec = beta.outcome.vec, se.outcome.vec = se.outcome.vec, pval.exposure.mat = pval.exposure.mat, use.openGWAS = FALSE)

rho.mat <- tmp$cor.mat

### Step 3

step3.res <- mvmr.cml.susie.step3(step2.res$mvdat, step2.res$invalid.idx, step2.res$theta.vec, rho.mat)

short.traits <- traits[-c(22,26,27)]

short.traits[which(step3.res$alpha[1,] > 1/29)] # LDL-C
short.traits[which(step3.res$alpha[2,] > 1/29)] # ApoA1, HDL-C
short.traits[which(step3.res$alpha[3,] > 1/29)] # L-LDL-FC, IDL-L, IDL-FC, IDL-C

### Model enumeration

models <- list(c(28,25,7),c(28,25,16),c(28,25,17),c(28,25,18),c(28,27,7),c(28,27,16),c(28,27,17),c(28,27,18))

names(models) <- c('LDL-C+ApoA1+L_LDL-FC','LDL-C+ApoA1+IDL-L','LDL-C+ApoA1+IDL-FC','LDL-C+ApoA1+IDL-C',
                   'LDL-C+HDL-C+L_LDL-FC','LDL-C+HDL-C+IDL-L','LDL-C+HDL-C+IDL-FC','LDL-C+HDL-C+IDL-C')

res <- data.frame(comb = names(models),
                  exp1 = NA,
                  exp2 = NA,
                  exp3 = NA)

for (i in 1:8) {
  
  tmp.sel.files <- sel.files2[models[[i]]]
  tmp.exp.files <- tmp.sel.files
  
  # p-value threshold 5e-8 # 592 independent SNPs
  tmp <- getInput(tmp.sel.files, tmp.exp.files, out.file,
                  plink_exe = genetics.binaRies::get_plink_binary(), plink_refdat = bfile, 
                  max.p.thres = 5e-8, cal.cor = TRUE, get.marker.candidates = FALSE)
  write.csv(tmp$data,file = paste0("/chru/analysis/yinxiangwu/lipids/For MVMR-cML-SuSiE/model.comb/",names(models)[i],".csv"),row.names = FALSE)
  rownames(tmp$cor.mat) <- c(traits[models[[i]]],'CAD')
  colnames(tmp$cor.mat) <- c(traits[models[[i]]],'CAD')
  write.csv(tmp$cor.mat,file = paste0("/chru/analysis/yinxiangwu/lipids/For MVMR-cML-SuSiE/model.comb/",names(models)[i],"_cor.csv"),row.names = FALSE)
  
  tmp.beta.exposure <- tmp$data[,paste0('gamma_exp',1:3)] %>% as.matrix(.)
  tmp.se.exposure <- tmp$data[,paste0('se_exp',1:3)] %>% as.matrix(.)
  tmp.beta.outcome <- tmp$data[,paste0('gamma_out',1)] %>% as.matrix(.)
  tmp.se.outcome <- tmp$data[,paste0('se_out',1)] %>% as.matrix(.)
  
  tmp.rho.mat <- tmp$cor.mat
  
  Sig.inv <- invcov_mvmr(se_bx = tmp.se.exposure,
                                se_by = tmp.se.outcome,
                                rho_mat = tmp.rho.mat)
  
  n <- min(sample.sizes.subset[models[[i]]])
  
  MVcML.res.sub1 <- MVmr_cML(b_exp = tmp.beta.exposure,
                                b_out = tmp.beta.outcome,
                                se_bx = tmp.se.exposure,
                                Sig_inv_l = Sig.inv, n = n,
                                K_vec = 0:(nrow(tmp.beta.exposure)-3))
  
  res[i,2:4] <- MVcML.res.sub1$BIC_theta
  
}

write.csv(res, '/chru/analysis/yinxiangwu/lipids/For MVMR-cML-SuSiE/res.csv', row.names = FALSE)
