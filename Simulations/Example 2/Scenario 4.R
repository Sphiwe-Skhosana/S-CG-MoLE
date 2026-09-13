# Function to simulate data from mixture of regressions with Gaussian distributed errors
simulate_data <- function(n,x,u, beta1, beta2,sigma1,sigma2) {
  # Calculate mixing proportions using non-monotonic function
  pi_u <- 0.1 + 0.8 * sin(pi * u)
  
  # Generate component membership
  z <- rbinom(n, 1, pi_u)
  
  m=n*0.05
  idx=sample(1:n,m)
  # Generate the response
  y <- numeric(n-m)
  y[idx]<- runif(m,-15,15)
  for (i in setdiff(1:n,idx)) {
    if (z[i] == 1) {
      y[i] <- beta1[1] + beta1[2] * x[i] + rnorm(1,0,sigma1)
    } else {
      y[i] <- beta2[1] + beta2[2] * x[i] + rnorm(1,0,sigma2)
    }
  }
  return(list(x = x,u=u, y = y, true_pi = pi_x, true_z = z))
}

##Simulations
num_samples=100
# Arrays to store results
beta1_estimates <- matrix(0, nrow = num_samples, ncol = 2*6)
beta2_estimates <- matrix(0, nrow = num_samples, ncol = 2*6)
sigma_estimates <- matrix(0, nrow = num_samples, ncol = 2*6)


# True parameters
beta1 <- c(0, 1)
beta2 <- c(4, 1)
sigma1 <- 1
sigma2 <- 1

k=2
for(n in c(1e3)){
  # Array to store all mixing proportion estimates
  pi_est1=pi_est2=pi_est3=pi_est4=pi_est5=pi_est6=array(0, dim = c(n, 2, num_samples))
  # Generate predictor variable
  x <- sort(rnorm(n, 0, 1))
  u <- runif(n)
  pi_u <- 0.1 + 0.8 * sin(pi * u)
  count=0
  # Run simulation study
  while(count<num_samples) {
    model1=model2=model3=model4=model5=model6=NULL
    data <- simulate_data(n, x,u, beta1, beta2, sigma1, sigma2)
    init.model=list(beta0=cbind(beta1,beta2),sigma0=c(sigma1,sigma2),pi0=cbind(data$true_pi,1-data$true_pi))
    try({model1 <- GMoLE.fit(x=data$x,u=NULL, y=data$y,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = T)
    try({model2 <- CGMoLE.fit(x=data$x,u=NULL, y=data$y,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = T)
    try({model3 <- GMoLE.fit(x=data$x,u=u, y=data$y,k=k,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = T)
    try({model4 <- CGMoLE.fit(x=data$x,u=u, y=data$y,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = T)
    try({model5 <- S_GMoLE.fit(x=data$x,u=u,y=data$y, max_iter = 1e3, K=k,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = T)
    try({model6 <- S_CG_MoLE.fit(x=data$x,u=u, y=data$y,max_iter = 1e3,init.beta = init.model$beta0,init.sigma=init.model$sigma0,init.pi=init.model$pi0)},silent = F)
    ##Store results
    if(!is.null(model1)&!is.null(model2)&!is.null(model3)&!is.null(model4)&!is.null(model5)){
      count=count+1;i=count 
      beta1_estimates[i,] <- c(model1$fit$beta[,1],model2$beta[,1],model3$fit$beta[,1],model4$beta[,1],model5$beta[,1],model6$beta[,1])
      beta2_estimates[i,] <- c(model1$fit$beta[,2],model2$beta[,2],model3$fit$beta[,2],model4$beta[,2],model5$beta[,2],model6$beta[,2])
      sigma_estimates[i,] <- c(model1$fit$sigma,model2$sigma,model3$fit$sigma,model4$sigma,model5$Sigma,model6$sigma)
      pi_est1[,,i] <- matrix(model1$fit$lambda,n,k,byrow=T)-cbind(pi_u,1-pi_u); pi_est2[,,i] <- matrix(model2$pi,n,k,byrow=T)-cbind(pi_u,1-pi_u);pi_est3[, , i] <- model3$fit$lambda-cbind(pi_u,1-pi_u)
      pi_est4[,,i] <- model4$pi_x-cbind(pi_u,1-pi_u);pi_est5[, , i] <- model5$mix.prop-cbind(pi_u,1-pi_u);pi_est6[, , i] <- model6$pi_x-cbind(pi_u,1-pi_u)
    }
  }
  write.csv(round(beta1_estimates,4),paste0("beta1_estimates_",n,".csv"))
  write.csv(round(beta2_estimates,4),paste0("beta2_estimates_",n,".csv"))
  write.csv(round(sigma_estimates,4),paste0("sigma_estimates_",n,".csv"))
  out_b_beta1=out_b_beta2=out_b_sigma=NULL
  for(j in 1:2){
    out_b_beta1=rbind(out_b_beta1,colMeans(beta1_estimates)-matrix(cbind(beta1,beta2)[,j],nrow=1,ncol=ncol(beta1_estimates),byrow=T))
    out_b_beta2=rbind(out_b_beta2,colMeans(beta2_estimates)-matrix(cbind(beta2,beta1)[,j],nrow=1,ncol=ncol(beta1_estimates),byrow=T))
    out_b_sigma=rbind(out_b_sigma,colMeans(sigma_estimates)-matrix(c(sigma1,sigma2)[cbind(c(1,2),c(2,1))[,j]],nrow=1,ncol=ncol(sigma_estimates),byrow=T))
  }
  idx1=which.min(apply(abs(out_b_beta1),1,max));idx2=which.min(apply(abs(out_b_beta2),1,max));idx3=which.min(apply(abs(out_b_sigma),1,max))
  ##BIAS
  bias_beta1=out_b_beta1[idx1,]
  bias_beta2=out_b_beta2[idx2,]
  bias_sigma=out_b_sigma[idx3,]
  ##MSE
  out_m_beta1=out_m_beta2=out_m_sigma=NULL
  for(j in 1:2){
    out_m_beta1=rbind(out_m_beta1,colMeans((beta1_estimates-matrix(cbind(beta1,beta2)[,j],nrow=num_samples,ncol=ncol(beta1_estimates),byrow=T))^2))
    out_m_beta2=rbind(out_m_beta2,colMeans((beta2_estimates-matrix(cbind(beta2,beta1)[,j],nrow=num_samples,ncol=ncol(beta1_estimates),byrow=T))^2))
    out_m_sigma=rbind(out_m_sigma,colMeans((sigma_estimates-matrix(c(sigma1,sigma2)[cbind(c(1,2),c(2,1))[,j]],nrow=num_samples,ncol=ncol(sigma_estimates),byrow=T))^2))
  }
  idx1=which.min(apply(out_m_beta1,1,max));idx2=which.min(apply(out_m_beta2,1,max));idx3=which.min(apply(out_m_sigma,1,max))
  mse_beta1=out_m_beta1[idx1,]
  mse_beta2=out_m_beta2[idx2,]
  mse_sigma=out_m_sigma[idx3,]
  mse_pi1=apply(pi_est1,3,function(x) sum(rowSums(x^2))/n)
  mse_pi2=apply(pi_est2,3,function(x) sum(rowSums(x^2))/n)
  mse_pi3=apply(pi_est3,3,function(x) sum(rowSums(x^2))/n)
  mse_pi4=apply(pi_est4,3,function(x) sum(rowSums(x^2))/n)
  mse_pi5=apply(pi_est5,3,function(x) sum(rowSums(x^2))/n)
  mse_pi6=apply(pi_est6,3,function(x) sum(rowSums(x^2))/n)
  mse_pi=cbind(c(mean(mse_pi1),sd(mse_pi1)),c(mean(mse_pi2),sd(mse_pi2)),c(mean(mse_pi3),sd(mse_pi3)),c(mean(mse_pi4),sd(mse_pi4)),c(mean(mse_pi5),sd(mse_pi5)),c(mean(mse_pi6),sd(mse_pi6)))
  colnames(mse_pi)=c("GMLR","CGMLR","GMoLE","CGMoLE","SGMoLE","SCGMoLE");
  rownames(mse_pi)=c("mean","sd")
  bias_beta=rbind(bias_beta1,bias_beta2);colnames(bias_beta)=c(rep("GMLR",2),rep("CGMLR",2),rep("GMoLE",2),rep("CGMoLE",2),rep("SGMoLE",2),rep("SCGMoLE",2))
  names(bias_sigma)=c(rep("GMLR",k),rep("CGMLR",k),rep("GMoLE",k),rep("CGMoLE",k),rep("SGMoLE",k),rep("SCGMoLE",k))
  mse_beta=rbind(mse_beta1,mse_beta2);colnames(mse_beta)=c(rep("GMLR",2),rep("CGMLR",2),rep("GMoLE",2),rep("CGMoLE",2),rep("SGMoLE",2),rep("SCGMoLE",2))
  mse_sigma=rbind(mse_sigma)
  colnames(mse_sigma)=c(rep("GMLR",k),rep("CGMLR",k),rep("GMoLE",k),rep("CGMoLE",k),rep("SGMoLE",k),rep("SCGMoLE",k))
  rownames(bias_beta)=c("comp1_beta","comp2_beta")
  #rownames(bias_sigma)="sigma"
  rownames(mse_beta)=c("comp1_beta","comp2_beta")
  #rownames(mse_sigma)="sigma"
  write.csv(round(bias_beta,6),paste0("bias_beta_",n,".csv"))
  write.csv(round(bias_sigma,6),paste0("bias_sigma_",n,".csv"))
  write.csv(round(mse_beta,6),paste0("mse_beta_",n,".csv"))
  write.csv(round(mse_sigma,6),paste0("mse_sigma_",n,".csv"))
  write.csv(round(mse_pi,6),paste0("mse_pi_",n,".csv"))
}