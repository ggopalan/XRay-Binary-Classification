#Giri Gopalan (Joint work with Luke Bornn and Saku Vrtilek)

#PURPOSE: R Script to perform elliptical slice sampling for the classification of X-Ray binary black hole systems into 1 of 3 Classes: 1)black hole 2)non-pulsar 3)pulsar
library(Rcpp)
library(RcppArmadillo)
library(MASS)

#########################
#Overview of Model Parameters:
# - alpha_1, alpha_2, alpha_3 which correspond to the baseline propensity that each data point corresponds to a BH, non-pulsar, or pulsar.
# - beta_1, beta_2, beta_3 which correspond to the marginal effect that the latent variables have on the propensity to be of each type.
# - Y_1, Y_2, Y_3 are latent real valued variables which govern the propensity of each data point to be of each type.
#Overview of Model:
# - alphas have independent unit normal priors.
# - betas have independent unit normal priors.
# - Y_j are drawn according to a Gaussian process prior with squared exponential kernel.
# - The class of  data point i is drawn according to a categorical distribution where p_j is proportional to exp(alpha_j+beta_j*Y_ij), and j indexes over (3) classes.
# - The classes of the data points are independent conditioning on Y_j, alphas, and betas.
#NOTE: Y refers to the latent variables in the code, but this is Z in the paper.
#######################
#load data
load("rcpp_bbh_inverse.RData")

#log-likelihood function necessary for ESS
log_lik <- function(f){
	Ystar <- matrix(NA,nrow=length(Y),ncol=3)
	Ystar[,1] <- f[1:length(Y)]
	Ystar[,2] <- f[(length(Y)+1):(2*length(Y))]
	Ystar[,3] <- f[(2*length(Y)+1):(3*length(Y))]
	alpha <- f[(3*length(Y)+1):(3*length(Y)+3)]
	beta <- f[(3*length(Y)+4):(3*length(Y)+6)]
	llh_probs <- exp(t(alpha + beta*t(Ystar)))/rowSums(exp(t(alpha + beta*t(Ystar))))
	sum_llh_probs <- 0
	for(k in 1:length(Y))
	{
		sum_llh_probs <- sum_llh_probs+log(llh_probs[k,Y[k]])
	}
	return (sum_llh_probs)
}
#Implementation of the ESS algorithm
for(i in 2:N){
  print(i)
  f <- rep(0,3*length(Y)+6)
  f[1:length(Y)] <- Y_samples[,1,i-1]
  f[(length(Y)+1):(2*length(Y))] <- Y_samples[,2,i-1]
  f[(2*length(Y)+1):(3*length(Y))] <- Y_samples[,3,i-1]
  f[(3*length(Y)+1):(3*length(Y)+3)] <- alpha_samples[i-1,]
  f[(3*length(Y)+4):(3*length(Y)+6)] <- beta_samples[i-1,]
  llh_thresh <- log_lik(f) + log(unif_samples[i])
  f_star <- f*cos(theta[i])+norm_samples[i,]*sin(theta[i])
  while(log_lik(f_star) < llh_thresh)
  {
    if (theta[i] < 0) 
    {
      theta_min[i] <- theta[i]
    }
    else
    {
      theta_max[i] <- theta[i]
    } 
    theta[i] <- runif(n=1,min=theta_min[i],max=theta_max[i])  
    f_star <- f*cos(theta[i])+norm_samples[i,]*sin(theta[i])     
  }
  alpha_samples[i,] <- f_star[(3*length(Y)+1):(3*length(Y)+3)]
  beta_samples[i,] <- f_star[(3*length(Y)+4):(3*length(Y)+6)]
  Y_samples[,1,i] <- f_star[1:length(Y)]
  Y_samples[,2,i] <- f_star[(length(Y)+1):(2*length(Y))]
  Y_samples[,3,i] <- f_star[(2*length(Y)+1):(3*length(Y))]
}
#Save the samples
save.image("bbh_ess_10000reps.RData")
