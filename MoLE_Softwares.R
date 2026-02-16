library(mixtools)
if(!require("nnet")){install.packages("nnet")}else{library(nnet)}

###Initialization function####
initialise.random=function(x,y,K){
  # Initialization
  count=0
  R=50
  counter=0
  res=list()
  while(count<R){
    fit=NULL
    try({fit=GMLRs.fit(y,x,k=K)},silent=F)
    if(!is.null(fit)) {count=count+1;res[[count]]=fit;}
    counter=counter+1
    if(counter==1e2) count=R
  }
  return(res)
}

##
invlogit <- function(z) 1/(1+exp(-z))

####GMLRs####
GMLRs.fit=function(y,x,k=2,init.pi=NULL,init.beta=NULL,init.sigma=NULL){
  n=length(y)
  p=ncol(as.matrix(x))+1
  ##Initialization
  if(!is.null(init.beta)){init.beta=init.beta}else{init.beta=matrix(rnorm(k*p),p,k)}  # intercept, slope
  if(!is.null(init.sigma)){init.sigma=init.sigma}else{init.sigma=rgamma(k,1,1)}
  if(!is.null(init.pi)){init.pi=init.pi}else{init.pi=rep(1/k,k)}
  
  fit_mlr <- regmixEM(y, x, k = k,beta=init.beta,sigma=init.sigma,lambda = init.pi)

  df=k*(3+nrow(fit_mlr$beta))+(k-1)
  loglik=fit_mlr$loglik
  # Information criteria
  AIC<- -2 * loglik + 2 * df
  BIC<- -2 * loglik + log(n) * df
  return(list(fit=fit_mlr,AIC=AIC,BIC=BIC))
}

####CGMLRs####
CGMLRs.fit <- function(x, y, K = 2, max_iter = 500, tol = 1e-6, init.pi=NULL,init.beta=NULL,init.sigma=NULL,init.eta=NULL,init.alpha=NULL) {
  n <- length(y)
  X <- cbind(1, x)
  p=ncol(X)
  
  # Initialize
  if(!is.null(init.pi)){pi_old=init.pi}else{pi_old=rep(1/K,K)}
  if(!is.null(init.beta)){beta_old=init.beta}else{beta_old=matrix(rnorm(K*p),p,K)}  # intercept, slope
  if(!is.null(init.sigma)){sigma_old=init.sigma}else{sigma_old=rgamma(K,1,1)}
  if(!is.null(init.alpha)){alpha_old=init.alpha}else{alpha_old <- rep(0.85, K)}  # proportion "good"
  if(!is.null(init.eta)){eta_old=init.eta}else{eta_old <- rep(2.0,  K)}  # variance inflation for "bad"
  
  z <- matrix(0, n, K)
  v <- matrix(0, n, K)
  
  for (iter in 1:max_iter) {
    # --- E-step ---
    for (k in 1:K) {
      mean_k <- X%*%beta_old[,k]
      z[, k] <- pi_old[k] * (
        alpha_old[k] * dnorm(y, mean_k, sigma_old[k]) +
          (1 - alpha_old[k]) * dnorm(y, mean_k, sqrt(eta_old[k]) * sigma_old[k])
      )
    }
    z <- z / rowSums(z)
    
    # --- v-step ---
    for (k in 1:K) {
      mean_k <- X%*%beta_old[, k]
      denom <- alpha_old[k] * dnorm(y, mean_k, sigma_old[k]) +
        (1 - alpha_old[k]) * dnorm(y, mean_k, sqrt(eta_old[k]) * sigma_old[k])
      v[, k] <- ifelse(denom > 0, (alpha_old[k] * dnorm(y, mean_k, sigma_old[k])) / denom,0.5)
    }
    
    # --- CM-step ---
    pi_new <- colMeans(z)
    beta_new <- matrix(0, nrow =p, ncol = K)
    sigma_new <- numeric(K)
    eta_new <- numeric(K)
    alpha_new <- numeric(K)
    
    for (k in 1:K) {
      denom_z <- sum(z[, k])
      if (denom_z < .Machine$double.eps) denom_z <- .Machine$double.eps
      
      alpha_new[k] <- sum(z[, k] * v[, k]) / denom_z
      
      W <- diag(z[, k] * (v[, k] + (1 - v[, k]) / eta_old[k]))
      
      XtWX <- t(X) %*% W %*% X
      XtWy <- t(X) %*% W %*% y
      beta_new[,k] <- qr.solve(XtWX, XtWy)
      
      num_sigma <- sum(z[, k] * (v[, k] + (1 - v[, k]) / eta_old[k]) *
                         (y - X %*% beta_new[,k])^2)
      denom_sigma <- sum(z[, k])
      if (denom_sigma < .Machine$double.eps) denom_sigma <- .Machine$double.eps
      sigma_new[k] <- sqrt(num_sigma / denom_sigma)
    }
    for (k in 1:K) {
      denom_eta <- sum(z[, k] * (1 - v[, k]))
      if (denom_eta < .Machine$double.eps) denom_eta <- .Machine$double.eps
      eta_val <- sum(z[, k] * (1 - v[, k]) * (y - X %*%beta_new[,k])^2) /
        (denom_eta * sigma_new[k]^2)
      if (is.na(eta_val) || eta_val < 1) eta_val <- 1.5
      eta_new[k] <- eta_val
    }
    
    # Order by intercept
    ord <- order(beta_new[1, ])
    beta_new <- beta_new[,ord , drop = FALSE]
    sigma_new <- sigma_new[ord]
    alpha_new <- alpha_new[ord]
    eta_new <- eta_new[ord]
    pi_new <- pi_new[ord]
    
    # Convergence check
    if (!any(is.na(c(pi_new, beta_new, sigma_new, alpha_new, eta_new)))) {
      if (max(abs(pi_new - pi_old)) < tol &&
          max(abs(beta_new - beta_old)) < tol &&
          max(abs(sigma_new - sigma_old)) < tol &&
          max(abs(alpha_new - alpha_old)) < tol &&
          max(abs(eta_new - eta_old)) < tol) {
        break
      }
    }
    
    # Update
    pi_old <- pi_new
    beta_old <- beta_new
    sigma_old <- sigma_new
    alpha_old <- alpha_new
    eta_old <- eta_new
  }
  
  # Log-likelihood
  loglik <- sum(log(rowSums(sapply(1:K, function(k) {
    pi_new[k] * (
      alpha_new[k] * dnorm(y, X%*%beta_new[,k], sigma_new[k]) +
        (1 - alpha_new[k]) * dnorm(y, X%*%beta_new[,k], sqrt(eta_new[k]) * sigma_new[k])
    )
  }))))
  
  for (k in 1:K) {
    mean_k <- X%*%beta_new[,k]
    z[, k] <- pi_new[k] * (
      alpha_new[k] * dnorm(y, mean_k, sigma_new[k]) +
        (1 - alpha_new[k]) * dnorm(y, mean_k, sqrt(eta_new[k]) * sigma_new[k])
    )
  }
  z <- z / rowSums(z)
  df=k*(ncol(beta_new)+3)+k-1
  AIC <- -2 * loglik + 2*df
  BIC <- -2 * loglik + log(n) * df

  out=list(pi = pi_new, beta = beta_new, sigma = sigma_new,
           alpha = alpha_new, AIC=AIC, BIC=BIC, eta = eta_new, loglik = loglik, z=z)
  return(out)
}

###GMoLE####
GMoLE.fit=function(x,y,k=NULL,init.pi=NULL,init.beta=NULL,init.sigma=NULL){
  n=length(y)
  if(!is.null(init.beta)){init.beta=init.beta}else{init.beta=matrix(rnorm(k*p),p,k)}  # intercept, slope
  if(!is.null(init.sigma)){init.sigma=init.sigma}else{init.sigma=rgamma(k,1,1)}
  fit_mle <- hmeEM(x = x, y = y, k = k, lambda=init.pi[,1],beta = init.beta,sigma=init.sigma)
  ##Reorder
  num_params_mle <- k*(nrow(fit_mle$beta)+1)+(k-1)*length(fit_mle$w)
  loglik=fit_mle$loglik
  
  # Information criteria
  AIC <- -2 * loglik + 2 * num_params_mle
  BIC <- -2 * loglik + log(n) * num_params_mle
  return(list(fit=fit_mle,AIC=AIC,BIC=BIC))
}

####CGMoLE####
CGMoLE.fit <- function(x, y, K = 2, max_iter = 500, tol = 1e-6, verbose = FALSE,init.pi=NULL,init.beta=NULL,init.sigma=NULL,init.eta=NULL,init.alpha=NULL) {
  n <- length(y)
  X <- cbind(1, x)   # intercept + slope
  p=ncol(X)
  
  ## --- Initialization ---
  gamma <- matrix(0, nrow = p, ncol = K)
  # regression coefficients (K x p)
  if(!is.null(init.beta)){beta=init.beta}else{beta=matrix(rnorm(K*p),p,K)}  # intercept, slope
  if(!is.null(init.sigma)){sigma=init.sigma}else{sigma=rgamma(k,1,1)}
  
  # contamination (per component)
  if(!is.null(init.alpha)){alpha=init.alpha}else{alpha <- rep(0.85, K)}  # proportion "good"
  if(!is.null(init.eta)){eta=init.eta}else{eta<- rep(2.0,  K)}  # variance inflation for "bad"
  
  # storage
  z <- matrix(0, n, K)  # responsibilities
  v <- matrix(0, n, K)  # prob good within component
  loglik_prev <- -Inf
  
  #function for row-softmax via independent logits then normalize
  gating_probs <- function(gamma) {
    G <- matrix(0, n, K)
    for (k in 1:K) {
      G[, k] <- invlogit(X%*%gamma[,k])
    }
    G / rowSums(G)  # normalize to sum 1
  }
  
  # log-likelihood
  comp_density <- function(k) {
    mu <- X %*% beta[,k]
    alpha[k]*dnorm(y, mu, sigma[k]) +
      (1 - alpha[k]) * dnorm(y, mu, sqrt(eta[k]) * sigma[k])
  }
  total_loglik <- function() {
    pix <- gating_probs(gamma)
    dens <- pix[,1]*comp_density(1) + pix[,2]*comp_density(2)
    sum(log(pmax(dens, .Machine$double.xmin)))
  }
  
  ## --- ECM iterations ---
  for (iter in 1:max_iter) {
    # E-step: responsibilities (z) and good-prob (v)
    pix <- gating_probs(gamma)
    dens <- matrix(0, n, K)
    for (k in 1:K) {
      mu_k <- X %*% beta[,k]
      dens_good <- dnorm(y, mu_k, sigma[k])
      dens_bad  <- dnorm(y, mu_k, sqrt(eta[k]) * sigma[k])
      dens[, k] <- alpha[k]*dens_good + (1-alpha[k])*dens_bad
      
      # v_ik = P(good | y_i, k)
      v[, k] <- (alpha[k]*dens_good) / pmax(dens[, k], .Machine$double.xmin)
    }
    z_unnorm <- pix * dens
    z <- z_unnorm / rowSums(z_unnorm)
    
    # CM-step 1: update alpha, beta, sigma (eta, gamma fixed)
    for (k in 1:K) {
      w_k <- z[, k] * (v[, k] + (1 - v[, k]) / eta[k])  # effective weights
      
      # alpha_k
      alpha[k] <- sum(z[, k] * v[, k]) / pmax(sum(z[, k]), .Machine$double.xmin)
      alpha[k] <- min(max(alpha[k], 1e-3), 0.999)       # keep in (0,1)
      
      # beta_k (WLS)
      W <- diag(pmax(w_k, 1e-10), n, n)
      XtW <- t(X) %*% W
      beta[,k] <- solve(XtW %*% X, XtW %*% y)
      
      # sigma_k
      resid_k <- y - X %*% beta[,k]
      sigma[k] <- sqrt( sum(w_k * resid_k^2) / pmax(sum(z[, k]), 1e-10) )
      sigma[k] <- pmax(sigma[k], 1e-6)
    }
    
    # CM-step 2: update gamma (gating) — only for component 2 (comp 1 fixed to 0)
    # Maximize Q wrt gamma[2,]: weighted binomial regression with "soft labels"
    # Target: z_{i2} ~ Binomial with covariate x and weight 1 (or z-sum)
    # We use iteratively reweighted least squares via optim (BFGS)
    Q_gamma2 <- function(g) {
      # g = (g0, g1) for comp 2; comp 1 logits fixed at 0
      # Our gating_probs() uses independent logits then normalize;
      # keep consistent here:
      pi1 <- invlogit(0 + 0*x)         # all 0.5 before normalization
      pi2 <- invlogit(g[1] + g[2]*x)
      denom <- pi1 + pi2
      pi2n <- pi2/denom                # normalized prob for comp 2
      # Expected complete-data loglik for gating
      sum(z[,2] * log(pmax(pi2n, 1e-12)) + z[,1] * log(pmax(1 - pi2n, 1e-12)))
    }
    opt <- optim(gamma[2, ], fn = function(par) -Q_gamma2(par),
                 method = "BFGS", control = list(reltol = 1e-8, maxit = 200))
    gamma[2, ] <- opt$par
    
    # CM-step 3: updae eta (variance inflation) per component
    for (k in 1:K) {
      resid_k <- y - X %*% beta[,k]
      num <- sum(z[, k] * (1 - v[, k]) * (resid_k^2))
      den <- pmax(sum(z[, k] * (1 - v[, k])) * sigma[k]^2, 1e-10)
      eta[k] <- pmax(num / den, 1.01)    # keep eta > 1
      if (!is.finite(eta[k])) eta[k] <- 2.0
    }
    
    # check convergence via log-likelihood
    ll <- total_loglik()
    if (verbose && iter %% 10 == 0) cat(sprintf("Iter %d, loglik=%.4f\n", iter, ll))
    if (abs(ll - loglik_prev) < tol) break
    loglik_prev <- ll
  }
  
  # final responsibilities & gating
  pix <- gating_probs(gamma)
  dens <- cbind(comp_density(1), comp_density(2))
  z <- (pix * dens) / rowSums(pix * dens)
  
  df <- k*(ncol(beta)+3) + (k-1)*ncol(gamma)
  
  AIC <- -2 * ll + 2 * df
  BIC <- -2 * ll + log(n) * df
  
  out=list(
    beta = beta, sigma = sigma, alpha = alpha, eta = eta,
    gamma = gamma, z = z, pi_x = pix, loglik = total_loglik(),
    iters = iter, converged = (iter < max_iter),
    AIC=AIC,BIC=BIC
  )
  return(out)
}

####S-GMoLE####
S_GMoLE.fit <- function(x, y, K = 2, max_iter = 500,init.pi=NULL,init.beta=NULL,init.sigma=NULL) {
  n <- length(y)
  z <- cbind(1, x)  # design matrix
  p <- ncol(z)
  t <- x
  
  # initialize
   if(!is.null(init.pi)){init.pi=init.pi}else{init.pi=matrix(rep(1/K,K),n,K)}
  if(!is.null(init.beta)){beta=init.beta}else{beta=matrix(rnorm(K*p),p,K)}  # intercept, slope
  if(!is.null(init.sigma)){sigma=init.sigma}else{sigma=rgamma(k,1,1)}
  
  loglik_vec <- c()
  pi_x <- init.pi
  for (iter in 1:max_iter) {
    # --- E-step ---
    r0 <- sapply(1:K, function(k) {
      pi_x[, k] * dnorm(y, mean = z %*% beta[, k], sd = sigma[k])
    })
    r_sum <- rowSums(r0)
    # Handle numerical issues: replace zeros with machine epsilon
    r_sum[r_sum == 0] <- .Machine$double.xmin
    r <- r0 / r_sum
    # --- M-step (Beta, Sigma) ---
    Beta1 <- matrix(0, nrow = p, ncol = K)
    sigma21 <- numeric(k)
    for (k in 1:K) {
      # Add ridge regularization for stability
      W <- diag(pmax(r[, k], 1e-10))  # Ensure weights are positive
      XtWX <- t(z) %*% W %*% z + diag(1e-6, p)  # Ridge regularization
      XtWy <- t(z) %*% W %*% y
      Beta1[, k] <- solve(XtWX, XtWy)
      
      # Ensure variance is positive
      sigma21[k] <- max(sum(r[, k] * (y - z %*% Beta1[, k])^2) / sum(r[, k]), 1e-6)
    }
    # --- Update mixing proportions nonparametrically ---
    pi_x=NULL
    for(k in 1:K){
      zn=r[,k]
      temp_data=data.frame(y=zn,x=x)
      mod=nnet(y~x,data=temp_data,size=2,decay=0.01,linout=T,trace=F)
      pi_x=cbind(pi_x,predict(mod,temp_data))
    }
    pi_x <- pi_x / rowSums(pi_x)
    pi_x[is.na(pi_x)] <- 1/K
    
    # --- loglik with numerical stability ---
    loglik_components <- sapply(1:K, function(k) {
      pi_x[, k] * dnorm(y, mean = z %*% Beta1[, k], sd = sqrt(sigma21[k]))
    })
    loglik_sum <- rowSums(loglik_components)
    loglik_sum[loglik_sum == 0] <- .Machine$double.xmin  # Replace zeros
    loglik <- sum(log(loglik_sum))
    
    loglik_vec <- c(loglik_vec, loglik)
    
    # Check convergence (only if we have at least 2 iterations)
    if (iter > 2) {
      if (abs(loglik_vec[iter] - loglik_vec[iter - 1]) < 1e-6) break
    }
    
    beta <- Beta1
    sigma <- sqrt(sigma21)
  }
  
  # Compute BIC
  n_params <- K*(nrow(Beta1)+1)
  
  # Effective nonparametric parameters
  e_params <- (K-1)*length(mod$wts)
  
  # Total degrees of freedom
  df <- n_params + e_params
  
  # Information criteria
  AIC_val <- -2 * loglik + 2 * df
  BIC_val <- -2 * loglik + log(n) * df
  
  sigma1=sigma0
  
  g <- sapply(1:K, function(k) {
    pi_x[, k] * dnorm(y, mean = z %*% Beta1[, k], sd = sigma1[k])
  })
  gn=g/rowSums(g)
  return(list(beta = Beta1, Sigma = sigma1, mix.prop =pi_x,llk=loglik_vec,
              loglik = loglik, AIC = AIC_val, BIC = BIC_val, iterations = iter,z=gn))
}

####S-CG-MoLE####
S_CG_MoLE.fit <- function(x, y, K = 2, max_iter = 500, tol = 1e-6, verbose = FALSE,init.pi=NULL,init.beta=NULL,init.sigma=NULL,init.eta=NULL,init.alpha=NULL) {
  n=length(y)
  X=cbind(1,x)
  p=ncol(X)
  
  # Initialize parameters
  if(!is.null(init.pi)){init.pi=init.pi}else{init.pi=matrix(rep(1/K,K),n,K)}
  if(!is.null(init.beta)){beta=init.beta}else{beta=matrix(rnorm(K*p),p,K)}  # intercept, slope
  if(!is.null(init.sigma)){sigma=init.sigma}else{sigma=rgamma(K,1,1)}
  
  if(!is.null(init.alpha)){alpha=init.alpha}else{alpha <- rep(0.85, K)}  # proportion "good"
  if(!is.null(init.eta)){eta=init.eta}else{eta<- rep(2.0,  K)}  # variance inflation for "bad"
  
  # responsibilities and v
  z <- matrix(1/K, n, K)      # responsibilities (will be updated in E-step)
  v <- matrix(0.8, n, K)      # probability of being a "good" obs within component
  
  # initial non-parametric mixing proportions (equal)
  pi_x <- init.pi
  
  loglik_prev <- -Inf
  iter <- 0
  llvec=NULL
  
  for (iter in 1:max_iter) {
    # ----- E-step -----
    # Compute component densities (good and bad) for each component
    dens_good_mat <- matrix(0, n, K)
    dens_bad_mat  <- matrix(0, n, K)
    comp_dens_mat <- matrix(0, n, K)  # alpha*good + (1-alpha)*bad (per component)
    for (k in 1:K) {
      mu_k <- X%*%beta[, k]
      dens_good_mat[, k] <- dnorm(y, mean = mu_k, sd = sigma[k])
      dens_bad_mat[, k]  <- dnorm(y, mean = mu_k, sd = sqrt(eta[k]) * sigma[k])
      comp_dens_mat[, k] <- alpha[k] * dens_good_mat[, k] + (1 - alpha[k]) * dens_bad_mat[, k]
    }
    
    # Compute weighted mixture for responsibilities using pi_x
    z_unnorm <- pi_x * comp_dens_mat
    row_sums_unnorm <- rowSums(z_unnorm)
    row_sums_unnorm[row_sums_unnorm == 0] <- 1e-300
    z <- z_unnorm / row_sums_unnorm  # normalized responsibilities
    
    # Update v: probability observation is "good" conditional on component and y
    # v[i,k] = P(good | component k, y_i) = alpha[k]*dens_good / comp_dens
    for (k in 1:K) {
      comp_dens_k <- comp_dens_mat[, k]
      v[, k] <- (alpha[k] * dens_good_mat[, k]) / pmax(comp_dens_k, 1e-300)
      v[, k] <- pmin(pmax(v[, k], 1e-10), 1 - 1e-10)
    }
    
    # ----- CM-step -----
    for (k in 1:K) {
      # Update alpha_k (within-component good proportion)
      alpha[k] <- sum(z[, k] * v[, k]) / pmax(sum(z[, k]), 1e-10)
      alpha[k] <- pmin(pmax(alpha[k], 0.01), 0.99)
      
      # Update beta via weighted least squares using weights w_k
      w_k <- z[, k] * (v[, k] + (1 - v[, k]) / pmax(eta[k], 1.01))
      w_k <- pmax(w_k, 1e-10)
      XtWX <- t(X) %*% (w_k * X) + diag(1e-6, p)
      XtWy <- t(X) %*% (w_k * y)
      beta[,k] <- solve(XtWX, XtWy)
      
      # Update sigma_k
      residuals_k <- y - X%*%beta[, k]
      sigma[k] <- sqrt(sum(w_k * residuals_k^2) / pmax(sum(z[, k]), 1e-10))
      sigma[k] <- pmax(sigma[k], 1e-6)
    }
    # Non-parametric update of mixing proportions: smooth the responsibilities z[,k] over x
    pi_x=NULL
    for(k in 1:K){
      zn=z[,k]
      temp_data=data.frame(y=zn,x=x)
      mod=nnet(y~x,data=temp_data,size=2,decay=0.01,linout=T,trace=F)
      pi_x=cbind(pi_x,predict(mod,temp_data))
    }
    pi_x <- pi_x / rowSums(pi_x)
    pi_x[is.na(pi_x)] <- 1/K
    
    # Update eta (contamination factor)
    for (k in 1:K) {
      residuals_k <- y - X%*%beta[, k]
      res2 <- residuals_k^2
      num <- sum(z[, k] * (1 - v[, k]) * res2)
      den <- pmax(sum(z[, k] * (1 - v[, k])) * sigma[k]^2, 1e-10)
      eta[k] <- pmax(num / den, 1.01)  # keep > 1
      eta[k] <- pmin(pmax(eta[k], 1.01), 1e6)
    }
    
    # Compute log-likelihood: sum_i log( sum_k pi_x[i,k] * comp_dens_mat[i,k] )
    for (k in 1:K) {
      mu_k <- X%*%beta[, k]
      dens_good_mat[, k] <- dnorm(y, mean = mu_k, sd = sigma[k])
      dens_bad_mat[, k]  <- dnorm(y, mean = mu_k, sd = sqrt(eta[k]) * sigma[k])
      comp_dens_mat[, k] <- alpha[k] * dens_good_mat[, k] + (1 - alpha[k]) * dens_bad_mat[, k]
    }
    mix_dens <- rowSums(pi_x * comp_dens_mat)
    mix_dens <- pmax(mix_dens, 1e-300)
    loglik <- sum(log(mix_dens))
    
    if (iter > 1 && abs(loglik - loglik_prev) < tol) break
    loglik_prev <- loglik
    llvec=c(llvec,loglik)
  }
  
  n_params <- K*(ncol(beta)+3)
  
  # Effective nonparametric degrees of freedom
  e_params <- (K-1)*length(mod$wts)
  
  # Total degrees of freedom
  df <- n_params + e_params
  
  AIC_val <- -2 * loglik + 2 * df
  BIC_val <- -2 * loglik + log(n) * df
  
  out=list(beta = beta,
           sigma = sigma,
           alpha = alpha,
           eta = eta,
           z = z,llk=llvec,
           v = v,
           pi_x = pi_x,
           loglik = loglik,
           iterations = iter,AIC=AIC_val,BIC=BIC_val)
  return(out)
}




