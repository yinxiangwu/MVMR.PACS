#
# 11th April 2019
# Calculating the lasso estimate using the multivariable input data 


# cv = True:  Estimate the regularisation parameter bY cross-validation
# cv = False: Parameter can be specified manually

#Note: Intercept is always true





library(glmnet)


#
# class mv-mr input
#

setClass("mvMRInput",
         representation(betaX = "matrix",
                        betaY = "matrix",
                        betaXse = "matrix",
                        betaYse = "matrix",
                        exposure = "character",
                        outcome = "character",
                        snps = "character",
                        effect_allele = "character",
                        other_allele  = "character",
                        eaf           = "numeric",                        
			correlation = "matrix")
)


#
# output class
#



setClass("MRl1",
         representation(Exposure = "character",
                        Outcome = "character",
                        Estimate = "numeric",
			Lambda = "numeric")
)


glmnet_l1_mr = function(object, cv=TRUE, lambda=0.1, nfold = 5, cv.param="lambda.1se"){


	bX = object@betaX
	bY = object@betaY

	if(cv==TRUE){
		cv = cv.glmnet(bX, bY, family = "gaussian", nfold = nfold, type.measure = "mse", intercept=FALSE, alpha = 1,  standardize = FALSE)
		if(cv.param=="lambda.1se"){bestlambda=cv$lambda.1se}
		if(cv.param=="lambda.min"){bestlambda=cv$lambda.min}
		g.out =  glmnet(bX, bY, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = 1,  standardize = FALSE)
		l1.coeff = coef(g.out)[2:(ncol(bX)+1)]
	}

	else{
		bestlambda = lambda
		g.out =  glmnet(bX, bY, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = 1,  standardize = FALSE)
		l1.coeff = coef(g.out)[2:(ncol(bX)+1)]
	}

	return(new("MRl1", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		Estimate =  l1.coeff,
		Lambda = bestlambda
	))

}



glmnet_debiased_l1_mr = function(beta.exposure, se.exposure, beta.outcome, se.outcome, betawt, P, epsilon = 1e-2, lambda=0.1){
  
  K <- ncol(beta.exposure)
  
  if (is.null(betawt)) {w <- rep(1,K)} else {w <- 1/abs(betawt)}
  
  p <- nrow(beta.exposure)
  
  # create pseudo data
  
  W <- diag(1/se.outcome^2)
  
  M <- t(beta.exposure) %*% W %*% beta.exposure
  
  V <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,])/se.outcome[j]^2) %>% Reduce("+",.)
  
  MV <- M-V
  
  M_plus <- ADMM_proj(MV/sqrt(p), epsilon = epsilon)$mat * sqrt(p)
  
  Xtilde <- chol(M_plus)
  
  ytilde <- solve(t(Xtilde), t(beta.exposure) %*% W %*% beta.outcome)
  
  bestlambda = lambda
  
  g.out =  glmnet(Xtilde, ytilde, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = 1,  standardize = FALSE, penalty.factor = w)
  
  l1.coeff = coef(g.out)[2:(ncol(Xtilde)+1)]
  
  return(list(Estimate =  l1.coeff,
              Lambda = bestlambda
  ))
  
}


glmnet_debiased_enet_mr = function(beta.exposure, se.exposure, beta.outcome, se.outcome, betawt, P, epsilon = 1e-2,  lambda=0.1, alpha=0.1){
  
  K <- ncol(beta.exposure)
  
  if (is.null(betawt)) {w <- rep(1,K)} else {w <- 1/abs(betawt)}
  
  p <- nrow(beta.exposure)
  
  # create pseudo data
  
  W <- diag(1/se.outcome^2)
  
  M <- t(beta.exposure) %*% W %*% beta.exposure
  
  V <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,])/se.outcome[j]^2) %>% Reduce("+",.)
  
  MV <- M-V
  
  M_plus <- ADMM_proj(MV/sqrt(p), epsilon = epsilon)$mat * sqrt(p)
  
  Xtilde <- chol(M_plus)
  
  ytilde <- solve(t(Xtilde), t(beta.exposure) %*% W %*% beta.outcome)
  
  bestlambda = lambda
  
  bestalpha = alpha
  
  enet.out =  glmnet(Xtilde, ytilde, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = bestalpha,  standardize = FALSE, penalty.factor = w)
  
  enet.coeff = coef(enet.out)[2:(ncol(Xtilde)+1)]
  
  return(list(Estimate =  enet.coeff,
              Lambda = bestlambda,
              Alpha = bestalpha
  ))
  
  
}






