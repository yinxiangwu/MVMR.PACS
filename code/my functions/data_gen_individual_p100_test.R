data_gen_individual <- function(beta0, K=10, p = 100, n = 1e4, sigma_gamma = 0.01, sigma_u = 0.8, sigma_e = 0.8, seed = 123) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # SNPs associated with the first cluster
  n_p1 = n_p2 = 45
  sigma1 <- matrix(0.5, 6, 6)
  sigma1[1:3,1:3] <- 0.995
  sigma1[4:6,4:6] <- 0.9
  diag(sigma1) <- 1
  sigma1 <- sigma1 * sigma_gamma
  
  gamma_coef1 <- cbind(rmvnorm(n_p1, mean = rep(0, 6), sigma = sigma1), matrix(0, nrow = n_p1, ncol = 4))
  
  sigma2 <- matrix(0.3, 4, 4)
  diag(sigma2) <- 1
  sigma2 <- sigma2 * sigma_gamma
  
  gamma_coef2 <- cbind(matrix(0, nrow = n_p2, ncol = 6), rmvnorm(n_p2, mean = rep(0, 4), sigma = sigma2))
  
  # SNPs associated with all risk factors
  n_p3 <- 10
  sigma3 <- matrix(0, K, K)
  sigma3[1:6, 1:6] <- sigma1/sigma_gamma
  sigma3[7:10,7:10] <- 0.3
  sigma3[1:6,7:10] <- 0.3
  sigma3[7:10,1:6] <- 0.3
  diag(sigma3) <- 1
  sigma3 <- sigma3 * sigma_gamma
  
  gamma_coef3 <- rmvnorm(n_p3, mean = rep(0, K), sigma = sigma3)
  
  gamma_coef <- rbind(gamma_coef1, gamma_coef2, gamma_coef3, matrix(0,nrow = p,ncol = K)) # fixed across simulations
  
  # stop setting seed
  set.seed(NULL)
  
  # MAF
  maf <- runif(2*p, 0.01, 0.5)
  
  data_gen_onesample <- function() {
    z <- matrix(nrow = n, ncol = 2*p)
    for (j in 1:(2*p)) {
      z[, j] <- matrix(rbinom(n, size = 2, prob = maf[j]), ncol = 1) # z's are independent
    }
    u <- rnorm(n, 0, sigma_u)
    x <- z %*% gamma_coef + u + rmvnorm(n, mean = rep(0, K), sigma = diag(K) * sigma_e)
    x <- scale(x)
    y <- 10 + x %*% beta0 + u + rnorm(n)
    return(list(z = z, x = x, y = y))
  }
  
  tmp1 <- data_gen_onesample() # exposure dataset
  tmp2 <- data_gen_onesample() # outcome dataset
  
  full_df <- data.frame(beta.exposure1 = numeric(2*p))
  full_df[, c(paste0('beta.exposure', 1:K), paste0('se.exposure', 1:K))] <- sapply(1:(2*p), function(i) {
    m <- lm(tmp1$x ~ tmp1$z[, i])
    m_summary <- summary(m)
    return(c(m$coef[2, ],
             m_summary$`Response Y1`$coefficients[2, 2],
             m_summary$`Response Y2`$coefficients[2, 2],
             m_summary$`Response Y3`$coefficients[2, 2],
             m_summary$`Response Y4`$coefficients[2, 2],
             m_summary$`Response Y5`$coefficients[2, 2],
             m_summary$`Response Y6`$coefficients[2, 2],
             m_summary$`Response Y7`$coefficients[2, 2],
             m_summary$`Response Y8`$coefficients[2, 2],
             m_summary$`Response Y9`$coefficients[2, 2],
             m_summary$`Response Y10`$coefficients[2, 2]))
  }) %>% t(.)
  
  full_df[, c('beta.outcome', 'se.outcome')] <- sapply(1:(2*p), function(i) {
    m <- lm(tmp2$y ~ tmp2$z[, i])
    m_summary <- summary(m)
    return(c(m$coefficients[-1], m_summary$coefficients[2, 2]))
  }) %>% t(.)
  
  #####################
  for (i in 1:K) {
    full_df[, paste0('pval.exposure', i)] <- 2 * pnorm(abs(full_df[, paste0('beta.exposure', i)]) / full_df[, paste0('se.exposure', i)], lower.tail = FALSE)
  }
  
  null_snps <- (p+1):(2*p)
  
  # Estimate the Z-score matrix
  z_mat <- full_df[null_snps, paste0('beta.exposure', 1:K)] /
    full_df[null_snps, paste0('se.exposure', 1:K)]
  
  P <- cor(z_mat)
  
  return(list(full_df = full_df[-null_snps,],
              P = P,
              gamma_coef = gamma_coef))
}
