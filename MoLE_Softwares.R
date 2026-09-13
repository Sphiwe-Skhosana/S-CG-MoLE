library(mixtools)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ContaminatedMixt)
if(!require("nnet")){install.packages("nnet")}else{library(nnet)}


initialise.random=function(x,y,K){
  # Initialization (for K>1)
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
invlogit<-function(z) 1/(1+exp(-z))

###G-MoLE (Note: if u=NULL then we get G-MLR)
mlogit_update <- function(U, R, w, steps = 5) {
  q<-ncol(U)
  k<-ncol(R)
  for (s in 1:steps) {
    eta <- cbind(0, U %*% w)
    P   <- exp(eta - apply(eta, 1, max))
    P   <- P / rowSums(P)
    
    grad <- as.vector(crossprod(U, R[, -1, drop = FALSE] - P[, -1, drop = FALSE]))
    H    <- matrix(0, q * (k - 1), q * (k - 1))
    for (a in 1:(k - 1)) {
      for (b in 1:(k - 1)) {
        d <- P[, a + 1] * ((a == b) - P[, b + 1])
        H[((a - 1) * q + 1):(a * q), ((b - 1) * q + 1):(b * q)] <- crossprod(U, d * U)
      }
    }
    w <- w + matrix(solve(H + diag(1e-6, q * (k - 1)), grad), q, k - 1)
  }
  w
}

GMoLE.fit<-function(x, y, k = 2, u = x, max_iter = 500, tol = 1e-6,
                      init.pi = NULL, init.beta = NULL, init.sigma = NULL) {
  n<-length(y)
  X<-cbind(1, as.matrix(x))
  U<-if (is.null(u)) matrix(1, n, 1) else cbind(1, as.matrix(u))
  p<-ncol(X)
  q<-ncol(U)
  
  ##Initialization
  beta<-if(is.null(init.beta))  matrix(rnorm(k * p), p, k) else init.beta
  sigma<-if(is.null(init.sigma)) rgamma(k, 1, 1)            else init.sigma
  pi_x<-if(is.null(init.pi))    matrix(1 / k, n, k)        else init.pi
  w<-matrix(0, q, k - 1)
  
  loglik_vec<-numeric(0)
  
  for (iter in 1:max_iter) {
    ##E-step
    r0<-sapply(1:k, function(j) pi_x[, j] * dnorm(y, X %*% beta[, j], sigma[j]))
    r<-r0 / pmax(rowSums(r0), .Machine$double.xmin)
    
    ##M-step: experts
    for (j in 1:k) {
      wj<-pmax(r[,j], 1e-10)
      beta[,j]<-solve(crossprod(X, wj * X) + diag(1e-6, p), crossprod(X, wj * y))
      sigma[j]<-sqrt(max(sum(wj * (y - X %*% beta[, j])^2) / sum(wj), 1e-6))
    }
    
    ##M-step: gating
    w<-mlogit_update(U, r, w)
    eta<-cbind(0, U %*% w)
    pi_x<-exp(eta - apply(eta, 1, max))
    pi_x<-pi_x / rowSums(pi_x)
    
    ##Observed-data log-likelihood
    dens<-sapply(1:k, function(j) pi_x[, j] * dnorm(y, X %*% beta[, j], sigma[j]))
    loglik<-sum(log(pmax(rowSums(dens), .Machine$double.xmin)))
    loglik_vec<-c(loglik_vec, loglik)
    
    if (iter > 1 && abs(loglik_vec[iter] - loglik_vec[iter - 1]) < tol) break
  }
  
  ##Information criteria
  ## beta (p) + sigma (1) per expert, plus (k-1) gating vectors of length q
  df<-k*(p + 1) + (k - 1) * q
  AIC<--2 * loglik + 2 * df
  BIC<--2 * loglik + log(n) * df
  
  ##Relabel components (increasing expert intercept)
  ## The gating coefficients are recentred on the new reference component.
  ord<-order(beta[1, ])
  beta<-beta[, ord, drop = FALSE]
  sigma<-sigma[ord]
  pi_x<-pi_x[, ord, drop = FALSE]
  Wf<-cbind(0, w)[, ord, drop = FALSE]
  w<- Wf[,-1, drop = FALSE] - Wf[, 1]
  
  g  <- sapply(1:k, function(j) pi_x[, j] * dnorm(y, X %*% beta[, j], sigma[j]))
  gn<-g/ rowSums(g)
  if(is.null(u)) pi_x=colMeans(pi_x)
  list(beta = beta, sigma = sigma, w = w, mix.prop = pi_x, z = gn,
       loglik = loglik, llk = loglik_vec, df = df, AIC = AIC, BIC = BIC,
       iterations = iter,
       ## hmeEM-style names, for existing downstream code
       fit = list(beta = beta, sigma = sigma, lambda = pi_x,
                  posterior = gn, w = w, loglik = loglik))
}

##CG-MoLE (Note that if u=NULL then we get CG-MLR)
CGMoLE.fit <- function(x, y, K = 2, u = x, max_iter = 500, tol = 1e-6,
                       verbose = FALSE, init.pi = NULL, init.beta = NULL,
                       init.sigma = NULL, init.eta = NULL, init.alpha = NULL) {
  
  n<-length(y)
  X<-cbind(1, as.matrix(x))
  U<-if (is.null(u)) matrix(1, n, 1) else cbind(1, as.matrix(u))
  p<-ncol(X)
  q<-ncol(U)
  
  ##Initialization
  beta<-if (is.null(init.beta))  matrix(rnorm(K * p), p, K) else init.beta
  sigma<-if (is.null(init.sigma)) rgamma(K, 1, 1)            else init.sigma
  alpha<-if (is.null(init.alpha)) rep(0.85, K)               else init.alpha
  eta<-if(is.null(init.eta))  rep(2.0, K)                else init.eta
  pi_x<-if(is.null(init.pi))  matrix(1 / K, n, K)        else init.pi
  gamma<-matrix(0, q, K - 1)               # gating, component 1 as reference
  
  z<-matrix(1 / K, n, K)                   # responsibilities
  v<-matrix(0.8, n, K)                     # P(good | component k)
  llvec<-numeric(0)
  
  ##ECM iterations
  for (iter in 1:max_iter) {
    
    ##E-step
    dens_good<-dens_bad<-comp_dens <- matrix(0, n, K)
    for (k in 1:K) {
      mu <- X %*% beta[, k]
      dens_good[,k]<-dnorm(y,mu,sigma[k])
      dens_bad[,k]<-dnorm(y,mu,sqrt(eta[k])*sigma[k])
      comp_dens[,k]<-alpha[k]*dens_good[,k]+(1-alpha[k])*dens_bad[,k]
    }
    zu<-pi_x*comp_dens
    z<-zu/pmax(rowSums(zu),.Machine$double.xmin)
    
    for (k in 1:K) {
      v[,k]<-alpha[k]*dens_good[, k]/pmax(comp_dens[,k],.Machine$double.xmin)
      v[,k]<-pmax(v[,k], 1e-10)
    }
    
    ##CM-step 1: alpha, beta, sigma, eta
    for (k in 1:K) {
      nk<-max(sum(z[, k]),1e-10)
      
      alpha[k]<-min(max(sum(z[,k] * v[, k]) / nk, 1e-3), 0.999)
      
      w <- pmax(z[, k] * (v[, k] + (1 - v[, k]) / max(eta[k], 1.01)), 1e-10)
      beta[,k] <- solve(crossprod(X, w * X) + diag(1e-6, p), crossprod(X, w * y))
      
      res<-y - X %*% beta[,k]
      sigma[k]<-max(sqrt(sum(w * res^2) / nk), 1e-6)
      
      num<-sum(z[, k] * (1 - v[, k]) * res^2)
      den<-max(sum(z[, k] * (1 - v[, k])) * sigma[k]^2, 1e-10)
      eta[k]<-min(max(num / den, 1.01), 1e6)
    }
    
    ##CM-step 2: gating
    gamma <- mlogit_update(U, z, gamma)
    lin<- cbind(0, U %*% gamma)
    pi_x<- exp(lin - apply(lin, 1, max))
    pi_x<- pi_x / rowSums(pi_x)
    
    ##Observed-data log-likelihood
    for (k in 1:K) {
      mu <- X %*% beta[, k]
      comp_dens[, k] <- alpha[k] * dnorm(y, mu, sigma[k]) +
        (1 - alpha[k]) * dnorm(y, mu, sqrt(eta[k]) * sigma[k])
    }
    ll<-sum(log(pmax(rowSums(pi_x * comp_dens), .Machine$double.xmin)))
    llvec<-c(llvec, ll)
    if (verbose && iter %% 10 == 0) cat(sprintf("Iter %d, loglik=%.4f\n", iter, ll))
    
    if (iter > 1 && abs(llvec[iter] - llvec[iter - 1]) < tol) break
  }
  
  ##Information criteria
  ## beta (p) + sigma + alpha + eta per component, plus (K-1) gating vectors
  df<-K * (p + 3) + (K - 1) * q
  AIC <--2 * ll + 2 * df
  BIC <--2 * ll + log(n) * df
  
  ##Relabel components (increasing expert intercept) ----
  ## Gating coefficients are recentred on the new reference component.
  ord <- order(beta[1, ])
  beta<- beta[, ord, drop = FALSE]
  sigma<- sigma[ord]
  alpha<- alpha[ord]
  eta<- eta[ord]
  pi_x<- pi_x[, ord, drop = FALSE]
  z <- z[, ord, drop = FALSE]
  v<- v[, ord, drop = FALSE]
  Gf<- cbind(0, gamma)[, ord, drop = FALSE]
  gamma <- Gf[,-1, drop = FALSE] - Gf[, 1]
  
  if(is.null(u)) pi_x=colMeans(pi_x)
  
  list(beta = beta, sigma = sigma, alpha = alpha, eta = eta, gamma = gamma,
       pi_x = pi_x, z = z, v = v,
       loglik = ll, llk = llvec, df = df, AIC = AIC, BIC = BIC,
       iterations = iter, converged = (iter < max_iter))
}

## S-G-MoLE
S_GMoLE.fit <- function(x, y, K = 2, u = x, max_iter = 500, tol = 1e-6,
                        init.pi = NULL, init.beta = NULL, init.sigma = NULL,
                        init.eta = NULL, init.alpha = NULL) {
  n<-length(y)
  X<-cbind(1, as.matrix(x))
  p<-ncol(X)
  dat<-data.frame(u)                       # gating covariate(s) for the smoother
  
  pi_x<-if (is.null(init.pi))    matrix(1 / K, n, K)        else init.pi
  beta<-if (is.null(init.beta))  matrix(rnorm(K * p), p, K) else init.beta
  sigma<-if (is.null(init.sigma)) rgamma(K, 1, 1)            else init.sigma
  
  loglik_vec <- numeric(0)
  
  for (iter in 1:max_iter) {
    ##E-step
    r0 <- sapply(1:K, function(k) pi_x[, k] * dnorm(y, X %*% beta[, k], sigma[k]))
    r  <- r0 / pmax(rowSums(r0), .Machine$double.xmin)
    
    ##M-step: beta, sigma
    for (k in 1:K) {
      w <- pmax(r[, k], 1e-10)
      beta[,k]<-solve(crossprod(X, w * X) + diag(1e-6, p), crossprod(X, w * y))
      sigma[k]<- sqrt(max(sum(w * (y - X %*% beta[, k])^2) / sum(w), 1e-6))
    }
    ##Nonparametric update of the mixing proportions (in u)
    for (k in 1:K) {
      d   <- data.frame(z = r[, k], dat)
      mod <- nnet(z ~ ., data = d, size = ceiling((2/3)*q)+1, decay = 0.01,
                  linout = TRUE, trace = FALSE)
      pi_x[, k] <- predict(mod, d)
    }
    pi_x<-pmax(pi_x, 1e-10)                # linout can return negatives
    pi_x<-pi_x/rowSums(pi_x)
    ##Observed-data log-likelihood
    dens<-sapply(1:K, function(k) pi_x[, k] * dnorm(y, X %*% beta[, k], sigma[k]))
    loglik<-sum(log(pmax(rowSums(dens), .Machine$double.xmin)))
    loglik_vec<-c(loglik_vec, loglik)
    
    if (iter > 1 && abs(loglik_vec[iter] - loglik_vec[iter - 1]) < tol) break
  }
  ##Information criteria
  n_params<- K * (p + 1)                    # beta (p per comp.) + sigma
  e_params<- (K - 1) * length(mod$wts)      # effective np. parameters
  df<- n_params + e_params
  AIC_val<- -2 * loglik + 2 * df
  BIC_val<- -2 * loglik + log(n) * df
  
  ##Label switching: Relabel components (increasing intercept)
  ord<-order(beta[1, ])
  beta<-beta[, ord, drop = FALSE]
  sigma<-sigma[ord]
  pi_x<-pi_x[, ord, drop = FALSE]
  
  g <-sapply(1:K, function(k) pi_x[, k] * dnorm(y, X %*% beta[, k], sigma[k]))
  gn<-g / rowSums(g)
  
  list(beta = beta, Sigma = sigma, mix.prop = pi_x, z = gn,
       loglik = loglik, llk = loglik_vec, AIC = AIC_val, BIC = BIC_val,
       df = df, iterations = iter)
}

## S-CG-MoLE
S_CG_MoLE.fit <- function(x, y, K = 2, u = x, max_iter = 500, tol = 1e-6,
                          verbose = FALSE, init.pi = NULL, init.beta = NULL,
                          init.sigma = NULL, init.eta = NULL, init.alpha = NULL) {
  n <-length(y)
  X<-cbind(1, as.matrix(x)); q=as.matrix(u)
  p<-ncol(X)
  q<-ncol(u)
  dat<-data.frame(u)                       # gating covariate(s) for the smoother
  
  ## If initial values are not provided, generate them randomly
  beta<-if (is.null(init.beta))  matrix(rnorm(K * p), p, K) else init.beta
  sigma<-if (is.null(init.sigma)) rgamma(K, 1, 1)            else init.sigma
  alpha<-if (is.null(init.alpha)) rep(0.85, K)               else init.alpha
  eta<-if (is.null(init.eta))   rep(2.0, K)                else init.eta
  pi_x<-if (is.null(init.pi))    matrix(1 / K, n, K)        else init.pi
  
  z<-matrix(1 / K, n, K)                   # responsibilities
  v<-matrix(0.8, n, K)                     # P(good | component k)
  llvec<-numeric(0)
  
  for (iter in 1:max_iter) {
    ##E-step
    dens_good <- dens_bad <- comp_dens <- matrix(0, n, K)
    for (k in 1:K) {
      mu<-X %*% beta[, k]
      dens_good[,k]<- dnorm(y, mu, sigma[k])
      dens_bad[,k]<-dnorm(y, mu, sqrt(eta[k]) * sigma[k])
      comp_dens[,k]<- alpha[k] * dens_good[, k] + (1 - alpha[k]) * dens_bad[, k]
    }
    zu <- pi_x * comp_dens
    z  <- zu / pmax(rowSums(zu), 1e-300)
    
    for (k in 1:K) {
      v[, k]<-alpha[k] * dens_good[, k] / pmax(comp_dens[, k], 1e-300)
      v[, k]<-pmax(v[, k], 1e-10)
    }
    
    ##CM-steps: alpha, beta, sigma, eta
    for (k in 1:K) {
      nk<-max(sum(z[, k]), 1e-10)
      
      alpha[k]<-min(max(sum(z[, k] * v[, k]) / nk, 0.01), 0.99)
      
      w<-pmax(z[, k] * (v[, k] + (1 - v[, k]) / max(eta[k], 1.01)), 1e-10)
      beta[,k]<-solve(crossprod(X, w * X) + diag(1e-6, p), crossprod(X, w * y))
      
      res<-y - X %*% beta[, k]
      sigma[k]<-max(sqrt(sum(w * res^2) / nk), 1e-6)
      
      num<-sum(z[,k]*(1 - v[, k]) * res^2)
      den<-max(sum(z[,k]*(1 - v[, k])) * sigma[k]^2, 1e-10)
      eta[k]<-min(max(num / den, 1.01), 1e6)
    }
    
    ##Nonparametric update of the mixing proportions (in u)
    for (k in 1:K) {
      d <-data.frame(zk = z[, k], dat)
      mod<-nnet(zk ~ ., data = d, size =ceiling((2/3)*q)+1 , decay = 0.01,
                  linout = TRUE, trace = FALSE)
      pi_x[, k] <- predict(mod, d)
    }
    pi_x<-pmax(pi_x, 1e-10)
    pi_x<-pi_x / rowSums(pi_x)
    
    ##Observed-data log-likelihood
    for (k in 1:K) {
      mu<-X %*% beta[, k]
      comp_dens[, k] <- alpha[k] * dnorm(y, mu, sigma[k]) +
        (1 - alpha[k]) * dnorm(y, mu, sqrt(eta[k]) * sigma[k])
    }
    loglik<- sum(log(pmax(rowSums(pi_x * comp_dens), 1e-300)))
    llvec<- c(llvec, loglik)
    if (verbose) cat("iter", iter, " loglik =", loglik, "\n")
    
    if(iter > 1 && abs(llvec[iter] - llvec[iter - 1]) < tol) break
  }
  
  ##Information criteria
  n_params<-K*(p + 3)                    # beta (p) + sigma + alpha + eta
  e_params<-(K-1)*length(mod$wts)
  df<- n_params+e_params
  AIC_val<--2*loglik + 2 * df
  BIC_val<--2*loglik + log(n) * df
  
  ##Label switching: Relabel components (increasing intercept)
  ord<-order(beta[1, ])
  beta<-beta[, ord, drop = FALSE]
  sigma<-sigma[ord]
  alpha<-alpha[ord]
  eta<-eta[ord]
  pi_x<-pi_x[,ord, drop = FALSE]
  z<-z[,ord,drop = FALSE]
  v<-v[,ord,drop = FALSE]
  list(beta=beta, sigma = sigma, alpha = alpha, eta = eta,
       pi_x = pi_x, z = z, v = v,
       loglik = loglik, llk = llvec, AIC = AIC_val, BIC = BIC_val,
       df = df, iterations = iter)
}
