#' @importFrom stats cor median sd
#' @importFrom utils capture.output
#' @importFrom Matrix sparseMatrix
#' @importFrom foreach %dopar%
#' @importFrom datathin datathin
#' @importFrom magrittr %>%
NULL

# Split the shared correlation matrix into the exposure block and the
# exposure-outcome correlations. When overlap = FALSE, P is K-by-K.
parse_overlap_correlation <- function(P, K, overlap = FALSE) {
  P <- as.matrix(P)
  expected <- if (isTRUE(overlap)) K + 1L else K
  if (!identical(dim(P), c(expected, expected))) {
    stop(
      "P must be ", expected, "-by-", expected,
      if (isTRUE(overlap)) " when overlap = TRUE." else " when overlap = FALSE."
    )
  }
  if (any(!is.finite(P)) || !isTRUE(all.equal(P, t(P), tolerance = 1e-8))) {
    stop("P must be a finite symmetric correlation matrix.")
  }
  if (min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("P must be positive semidefinite.")
  }
  list(
    P.exposure = P[seq_len(K), seq_len(K), drop = FALSE],
    rho.xy = if (isTRUE(overlap)) P[seq_len(K), K + 1L] else rep(0, K)
  )
}

compute_overlap_bias <- function(se.exposure, se.outcome, rho.xy) {
  se.exposure <- as.matrix(se.exposure)
  se.outcome <- as.numeric(se.outcome)
  as.numeric(colSums(se.exposure * (rho.xy / se.outcome)))
}

make_joint_thinning_covariance <- function(se.exposure, se.outcome, P) {
  se.exposure <- as.matrix(se.exposure)
  se.outcome <- as.numeric(se.outcome)
  p <- nrow(se.exposure)
  K <- ncol(se.exposure)
  Sig <- array(0, dim = c(p, K + 1L, K + 1L))
  for (j in seq_len(p)) {
    sj <- c(se.exposure[j, ], se.outcome[j])
    Sig[j, , ] <- P * tcrossprod(sj)
  }
  Sig
}


#####################################################################
#                                                                   #
# Code adapted from  Sharma et al. Consistent group identification  #
# and variable selection in regression with correlated predictors.  #
# Pubmed: 23772171                                                  #
#####################################################################
dIVW_PACS_cluster <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome,
    lambda, betawt, P, RR,
    tau = 1, type = 1, rr = 0, digit = 3,
    eps = 1e-4, max.iter = 1000, err = 1e-4,
    pleiotropy = FALSE,              # NEW
    lambda.alpha = NULL,             # NEW
    alpha.init = NULL,                # NEW
    overlap.bias = NULL
) {
  if (lambda <= 0) stop("lambda must be > 0.")
  if (rr < 0) stop("rr must be >= 0.")
  if (rr > 1) stop("rr must be <= 1.")
  if (eps <= 0) stop("eps must be > 0.")
  if (!(type %in% 1:4)) stop("type must be 1, 2, 3 or 4.")

  # NEW: check pleiotropy tuning parameter
  if (pleiotropy) {
    if (is.null(lambda.alpha)) {
      stop("ERROR: lambda.alpha must be provided when pleiotropy = TRUE. \n")
    }
    if (lambda.alpha < 0) {
      stop("ERROR: lambda.alpha must be >= 0. \n")
    }
  }

  # data preparation
  p <- length(beta.outcome)
  K <- ncol(beta.exposure)
  overlap.bias <- if (is.null(overlap.bias)) rep(0, K) else as.numeric(overlap.bias)
  if (length(overlap.bias) != K || any(!is.finite(overlap.bias))) {
    stop("overlap.bias must be a finite numeric vector of length K.")
  }
  littleeps <- 10^-7

  qp <- K * (K - 1) / 2
  vc <- 0
  row.vec1 <- 0
  col.vec1 <- 0
  for(w in 1:(K - 1)) {
    row.vec1 <- c(row.vec1, c((vc + 1):(vc + (K - w))))
    col.vec1 <- c(col.vec1, rep(w, length(c((vc + 1):(vc + (K - w))))))
    vc <- vc + (K - w)
  }

  c1 <- 1
  c2 <- K - 1
  row.vec2 <- 0
  col.vec2 <- 0
  for(w in 1:(K - 1)) {
    row.vec2 <- c(row.vec2, c(c1:c2))
    col.vec2 <- c(col.vec2, c((w + 1):K))
    c1 <- c2 + 1
    c2 <- c2 + K - w - 1
  }

  dm <- Matrix::sparseMatrix(
    i = c(row.vec1[-1], row.vec2[-1]),
    j = c(col.vec1[-1], col.vec2[-1]),
    x = c(rep(1, qp), rep(-1, qp))
  )

  dp <- Matrix::sparseMatrix(
    i = c(row.vec1[-1], row.vec2[-1]),
    j = c(col.vec1[-1], col.vec2[-1]),
    x = c(rep(1, qp), rep(1, qp))
  )

  rm(c1, c2, w, vc, row.vec1, col.vec1, row.vec2, col.vec2)

  # Construct penalty weight vector ascvec based on type
  abs_b <- pmax(abs(betawt), littleeps)
  abs_dm_b <- pmax(as.vector(abs(dm %*% betawt)), littleeps)
  abs_dp_b <- pmax(as.vector(abs(dp %*% betawt)), littleeps)

  if (type == 1) {
    ascvec <- c(1 / abs_b^tau, 1 / abs_dm_b^tau, 1 / abs_dp_b^tau)
  } else if (type == 2) {
    crm <- 1 / (1 - RR[lower.tri(RR)])^tau
    crp <- 1 / (1 + RR[lower.tri(RR)])^tau
    ascvec <- c(1 / abs_b^tau, crm / abs_dm_b^tau, crp / abs_dp_b^tau)
  } else if (type == 3) {
    corp <- ifelse(RR[lower.tri(RR)] > rr, 1, 0)
    corm <- ifelse(RR[lower.tri(RR)] < -rr, 1, 0)
    ascvec <- c(1 / abs_b^tau, corm / abs_dm_b^tau, corp / abs_dp_b^tau)
  } else if (type == 4) {
    corp <- ifelse(RR[lower.tri(RR)] > rr, 1, 0)
    crm  <- corp / (1 - RR[lower.tri(RR)])^tau
    corm <- ifelse(RR[lower.tri(RR)] < -rr, 1, 0)
    crp  <- corm / (1 + RR[lower.tri(RR)])^tau
    ascvec <- c(1 / abs_b^tau, crm / abs_dm_b^tau, crp / abs_dp_b^tau)
  }

  # Precompute weight matrices and debiased covariance
  W <- Matrix::Diagonal(x = 1 / se.outcome^2)

  M <- Matrix::crossprod(beta.exposure, W %*% beta.exposure)

  A_weighted <- se.exposure / matrix(se.outcome, nrow = p, ncol = K, byrow = FALSE)
  S <- Matrix::t(A_weighted) %*% A_weighted
  V <- P * S

  MV <- M - V
  MVplus <- BDcocolasso::ADMM_proj(as.matrix(MV) / sqrt(p), epsilon = eps)$mat * sqrt(p)

  # NEW: initialize alpha
  if (pleiotropy) {
    if (is.null(alpha.init)) {
      alphal <- rep(0, p)
    } else {
      if (length(alpha.init) != p) {
        return(cat("ERROR: alpha.init must have length equal to length(beta.outcome). \n"))
      }
      alphal <- as.numeric(alpha.init)
    }
  } else {
    alphal <- rep(0, p)
    lambda.alpha <- 0
  }

  # Iterative update
  betal <- betawt

  for (iter in seq_len(max.iter)) {
    old.beta <- betal
    old.alpha <- alphal

    # beta update given alpha
    # NEW: replace beta.outcome by beta.outcome - alpha
    XWy <- Matrix::crossprod(beta.exposure, W %*% (beta.outcome - alphal)) -
      overlap.bias

    D1 <- Matrix::Diagonal(x = ascvec[1:K] / (abs(old.beta) + littleeps))
    D2 <- Matrix::Diagonal(
      x = ascvec[(K + 1):(K + qp)] /
        (as.numeric(abs(dm %*% old.beta)) + littleeps)
    )
    D3 <- Matrix::Diagonal(
      x = ascvec[(K + qp + 1):(K + 2 * qp)] /
        (as.numeric(abs(dp %*% old.beta)) + littleeps)
    )

    A <- MVplus +
      lambda * (
        D1 +
          Matrix::t(dm) %*% D2 %*% dm +
          Matrix::t(dp) %*% D3 %*% dp
      )

    if (rcond(as.matrix(A)) > .Machine$double.eps) {
      betal_new <- as.numeric(Matrix::solve(A, XWy))
    } else {
      betal_new <- as.numeric(MASS::ginv(as.matrix(A)) %*% XWy)
    }

    # NEW: alpha update given beta via soft-thresholding
    if (pleiotropy) {
      r_alpha <- as.numeric(beta.outcome - beta.exposure %*% betal_new)

      thresh <- lambda.alpha * se.outcome^2

      alphal_new <- sign(r_alpha) * pmax(abs(r_alpha) - thresh, 0)
    } else {
      alphal_new <- rep(0, p)
    }

    beta.diff <- max(abs(betal_new - old.beta)) / (max(abs(old.beta)) + littleeps)
    alpha.diff <- max(abs(alphal_new - old.alpha)) / (max(abs(old.alpha)) + littleeps)

    betal <- betal_new
    alphal <- alphal_new

    if (max(beta.diff, alpha.diff) < err) {
      break
    }
  }

  cls <- groupAssignment(betal, digit = digit)

  list(
    coefficients = betal,
    alpha        = alphal,
    invalid.snp  = which(abs(alphal) > littleeps),
    cls          = cls,
    cls_chr      = paste(cls, collapse = "-"),
    lambda       = lambda,
    lambda.alpha = lambda.alpha,
    pleiotropy   = pleiotropy,
    overlap.bias = overlap.bias,
    rr           = abs(rr),
    init         = betawt,
    alpha.init   = alpha.init,
    type         = type,
    eps          = eps,
    littleeps    = littleeps,
    iter         = iter
  )
}

SRIVW <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome,
    gen_cor = NULL, phi_cand = 0,
    over.dispersion = FALSE, overlap = FALSE
) {
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  beta.outcome <- as.numeric(beta.outcome)
  se.outcome <- as.numeric(se.outcome)

  K <- ncol(beta.exposure)
  p <- nrow(beta.exposure)
  if (K < 1L || ncol(se.exposure) != K) {
    stop("beta.exposure and se.exposure must have the same positive number of columns.")
  }
  if (nrow(se.exposure) != p || length(beta.outcome) != p ||
      length(se.outcome) != p) {
    stop("Exposure and outcome inputs must contain the same number of SNPs.")
  }
  if (isTRUE(overlap) && isTRUE(over.dispersion)) {
    stop(
      "Simultaneous over.dispersion = TRUE and overlap = TRUE is not ",
      "currently supported by the SRIVW variance formula."
    )
  }

  if (is.null(gen_cor)) {
    gen_cor <- replicate(
      p, diag(if (isTRUE(overlap)) K + 1L else K), simplify = FALSE
    )
  } else if (is.matrix(gen_cor)) {
    gen_cor <- replicate(p, as.matrix(gen_cor), simplify = FALSE)
  }
  if (!is.list(gen_cor) || length(gen_cor) != p) {
    stop("gen_cor must be a correlation matrix or a list with one matrix per SNP.")
  }

  expected_dim <- if (isTRUE(overlap)) K + 1L else K
  if (any(vapply(gen_cor, function(x) {
    !identical(dim(as.matrix(x)), c(expected_dim, expected_dim))
  }, logical(1)))) {
    stop(
      "Each gen_cor matrix must be ", expected_dim, "-by-", expected_dim,
      if (isTRUE(overlap)) " when overlap = TRUE." else "."
    )
  }

  W <- diag(se.outcome^(-2), p, p)
  Sigma.joint <- lapply(seq_len(p), function(j) {
    if (isTRUE(overlap)) {
      sj <- c(se.exposure[j, ], se.outcome[j])
    } else {
      sj <- se.exposure[j, ]
    }
    diag(sj) %*% as.matrix(gen_cor[[j]]) %*% diag(sj)
  })
  Sigma.X <- lapply(Sigma.joint, function(x) x[seq_len(K), seq_len(K), drop = FALSE])

  Vj_root_inv <- lapply(seq_len(p), function(j) {
    P.X <- as.matrix(gen_cor[[j]])[seq_len(K), seq_len(K), drop = FALSE]
    ee <- eigen(P.X, symmetric = TRUE)
    if (min(ee$values) <= 0) stop("Exposure correlation matrix must be positive definite.")
    P_root_inv <- ee$vectors %*% diag(1 / sqrt(ee$values), K) %*% t(ee$vectors)
    P_root_inv %*% diag(1 / se.exposure[j, ], K)
  })

  IV_strength_matrix <- Reduce("+", lapply(seq_len(p), function(j) {
    bx.std <- Vj_root_inv[[j]] %*% beta.exposure[j, ]
    bx.std %*% t(bx.std)
  })) - p * diag(K)
  iv_strength_parameter <- min(
    eigen(IV_strength_matrix / sqrt(p), symmetric = TRUE, only.values = TRUE)$values
  )

  V <- Reduce("+", lapply(seq_len(p), function(j) {
    Sigma.X[[j]] * se.outcome[j]^(-2)
  }))
  M <- crossprod(beta.exposure, W %*% beta.exposure)
  MV <- M - V
  MV_eigen <- eigen(MV, symmetric = TRUE)
  MV_eigvalues <- MV_eigen$values

  if (is.null(phi_cand)) {
    phi_cand <- c(0, exp(seq(0, 17, by = 0.5) - iv_strength_parameter))
  }
  phi_cand <- as.numeric(phi_cand)
  phi_length <- length(phi_cand)

  regularized_inverse <- function(phi) {
    denom <- MV_eigvalues + phi / MV_eigvalues
    MV_eigen$vectors %*% diag(1 / denom, K) %*% t(MV_eigen$vectors)
  }
  MV.l.inv.list <- lapply(phi_cand, regularized_inverse)

  score <- as.numeric(crossprod(beta.exposure, W %*% beta.outcome))
  if (isTRUE(overlap)) {
    Adj_term <- Reduce("+", lapply(seq_len(p), function(j) {
      Sigma.joint[[j]][seq_len(K), K + 1L, drop = FALSE] *
        se.outcome[j]^(-2)
    }))
    score <- score - as.numeric(Adj_term)
  }
  beta.est <- lapply(MV.l.inv.list, function(Ainv) as.numeric(Ainv %*% score))

  prof.lik <- vapply(beta.est, function(beta.hat) {
    bvb <- vapply(seq_len(p), function(j) {
      drop(t(beta.hat) %*% Sigma.X[[j]] %*% beta.hat)
    }, numeric(1))
    if (isTRUE(overlap)) {
      bvxy <- vapply(seq_len(p), function(j) {
        sigmaxy <- Sigma.joint[[j]][seq_len(K), K + 1L, drop = FALSE]
        drop(t(beta.hat) %*% sigmaxy)
      }, numeric(1))
      residual_var <- se.outcome^2 + bvb - 2 * bvxy
    } else {
      residual_var <- se.outcome^2 + bvb
    }
    if (any(!is.finite(residual_var)) || any(residual_var <= 0)) return(Inf)
    residual <- beta.outcome - as.numeric(beta.exposure %*% beta.hat)
    mean(residual^2 / residual_var)
  }, numeric(1))

  selected_idx <- which.min(prof.lik)
  phi_selected <- phi_cand[selected_idx]
  MV.l.inv <- MV.l.inv.list[[selected_idx]]
  mvmr.adIVW <- beta.est[[selected_idx]]

  if (!isTRUE(overlap)) {
    if (isTRUE(over.dispersion)) {
      tau2_adivw <- Reduce("+", lapply(seq_len(p), function(j) {
        v <- Sigma.X[[j]] * se.outcome[j]^(-2)
        residual <- beta.outcome[j] - drop(beta.exposure[j, ] %*% mvmr.adIVW)
        residual^2 * se.outcome[j]^(-2) - 1 -
          drop(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
      })) / sum(diag(W))
      tau2_adivw <- max(as.numeric(tau2_adivw), 0)
    } else {
      tau2_adivw <- 0
    }

    adIVW_Vt <- Reduce("+", lapply(seq_len(p), function(j) {
      bx <- matrix(beta.exposure[j, ], ncol = 1)
      m <- bx %*% t(bx) * se.outcome[j]^(-2)
      v <- Sigma.X[[j]] * se.outcome[j]^(-2)
      bvb <- drop(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
      vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
      m * (1 + bvb + tau2_adivw * se.outcome[j]^(-2)) + vbbv
    }))
    meat <- adIVW_Vt
  } else {
    tau2_adivw <- NULL
    meat <- Reduce("+", lapply(seq_len(p), function(j) {
      bx <- matrix(beta.exposure[j, ], ncol = 1)
      bhat <- matrix(mvmr.adIVW, ncol = 1)
      m <- bx %*% t(bx) * se.outcome[j]^(-2)
      v <- Sigma.X[[j]] * se.outcome[j]^(-2)
      bvb <- drop(t(bhat) %*% v %*% bhat)
      vbbv <- v %*% bhat %*% t(bhat) %*% v
      sigmaxy <- Sigma.joint[[j]][seq_len(K), K + 1L, drop = FALSE]
      sy4 <- se.outcome[j]^(-4)
      A1 <- sigmaxy %*% t(sigmaxy) * sy4
      A2 <- Sigma.X[[j]] %*% bhat %*% t(sigmaxy) * sy4
      A3 <- sigmaxy %*% t(bhat) %*% Sigma.X[[j]] * sy4
      beta_sigmaxy_w <- drop(t(bhat) %*% sigmaxy) * se.outcome[j]^(-2)
      A4 <- beta_sigmaxy_w * m
      A5 <- beta_sigmaxy_w * sigmaxy %*% t(sigmaxy) * sy4
      m * (1 + bvb) + vbbv + A1 - A2 - A3 - 2 * A4 + 4 * A5
    }))
  }

  covariance <- MV.l.inv %*% meat %*% MV.l.inv
  beta.se <- sqrt(pmax(diag(covariance), 0))

  list(
    beta.hat = as.numeric(mvmr.adIVW),
    beta.se = as.numeric(beta.se),
    iv_strength_parameter = iv_strength_parameter,
    phi_selected = phi_selected,
    tau.square = tau2_adivw
  )
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

computeGroupedData <- function(
    beta.exposure, se.exposure, se.outcome, P, G, overlap = FALSE
) {
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  se.outcome <- as.numeric(se.outcome)
  K <- ncol(beta.exposure)
  p <- nrow(beta.exposure)

  overlap.info <- parse_overlap_correlation(P, K, overlap)
  P.exposure <- overlap.info$P.exposure
  beta.grouped <- beta.exposure %*% t(G)
  out <- computeGroupedSE(se.exposure, P.exposure, G)

  if (isTRUE(overlap)) {
    joint_cor <- vector("list", p)
    for (j in seq_len(p)) {
      Sigma.X <- diag(se.exposure[j, ], K) %*% P.exposure %*%
        diag(se.exposure[j, ], K)
      Sigma.XY <- se.exposure[j, ] * se.outcome[j] * overlap.info$rho.xy
      Sigma.joint <- rbind(
        cbind(Sigma.X, Sigma.XY),
        c(Sigma.XY, se.outcome[j]^2)
      )
      T.group <- matrix(0, nrow = nrow(G) + 1L, ncol = K + 1L)
      T.group[seq_len(nrow(G)), seq_len(K)] <- G
      T.group[nrow(G) + 1L, K + 1L] <- 1
      Sigma.group <- T.group %*% Sigma.joint %*% t(T.group)
      sd.group <- sqrt(diag(Sigma.group))
      joint_cor[[j]] <- Sigma.group / tcrossprod(sd.group)
    }
    P.grouped <- joint_cor
  } else {
    P.grouped <- out$Corr_group
  }

  list(
    beta.exposure = beta.grouped,
    se.exposure = out$se_grouped,
    P_grouped = P.grouped
  )
}

SRIVW_group_est <- function(
    beta.exposure, se.exposure, beta.outcome, se.outcome, P, est_beta,
    include_zero_group = FALSE, digit = 3,
    over.dispersion = FALSE, overlap = FALSE
) {
  K <- ncol(beta.exposure)
  p <- nrow(beta.exposure)
  group <- groupAssignment(est_beta, digit = digit)
  G_collapse <- constructG(
    x = group,
    y = sign(round(est_beta, digit)),
    include_zero_group = include_zero_group
  )
  L <- nrow(G_collapse)
  if (L < 1L) stop("No signal group is available for SRIVW inference.")

  group_dat <- computeGroupedData(
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    se.outcome = se.outcome,
    P = P,
    G = G_collapse,
    overlap = overlap
  )

  final_fit <- SRIVW(
    beta.exposure = group_dat$beta.exposure,
    se.exposure = group_dat$se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    phi_cand = if (L == 1L) 0 else NULL,
    gen_cor = group_dat$P_grouped,
    over.dispersion = over.dispersion,
    overlap = overlap
  )
  beta_est <- as.numeric(final_fit$beta.hat %*% G_collapse)
  beta_se <- abs(as.numeric(final_fit$beta.se %*% G_collapse))
  iv.strength.dt <- final_fit$iv_strength_parameter

  if (L > 1L) {
    P.exposure.grouped <- lapply(group_dat$P_grouped, function(x) {
      as.matrix(x)[seq_len(L), seq_len(L), drop = FALSE]
    })
    F.data <- MVMR::format_mvmr(
      BXGs = group_dat$beta.exposure,
      BYG = beta.outcome,
      seBXGs = group_dat$se.exposure,
      seBYG = se.outcome,
      RSID = seq_len(p)
    )
    invisible(capture.output(
      fres <- MVMR::strength_mvmr(
        r_input = F.data,
        gencov = lapply(seq_len(p), function(j) {
          diag(group_dat$se.exposure[j, ], L) %*%
            P.exposure.grouped[[j]] %*%
            diag(group_dat$se.exposure[j, ], L)
        })
      )
    ))
    condF <- abs(as.numeric(fres) %*% G_collapse)
  } else {
    condF <- rep(NA_real_, K)
  }

  list(
    beta_est = beta_est,
    beta_se = beta_se,
    iv_strength = iv.strength.dt,
    condF = condF
  )
}

generate_tuning_para <- function(beta.exposure, se.exposure, se.outcome, K, P, lambda.length) {

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
  MV_plusS <- BDcocolasso::ADMM_proj(MVs, epsilon = epsilon)$mat
  MV_plus  <- MV_plusS * sqrt(p)   # re-symmetrize + scale back

  # Ridge system: A b = rhs
  A   <- MV_plus + (phi * diag(K))
  rhs <- as.numeric(XtW %*% y)               # X' W y

  b <- solve(A) %*% rhs

  # Ensure vector output
  as.numeric(b)
}
