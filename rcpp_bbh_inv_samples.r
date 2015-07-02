#Giri Gopalan (Joint Work With Saku Vrtilek and Luke Bornn)
#Astrostatistics
#R script to calculate inverses and samples necessary for MCMC (elliptical slice sampling) sampler
#Necessary packages

library(mvtnorm)
library(MASS)
library(inline)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)


############################################
#INLINE RCPP FUNCTIONS
############################################
incl <- 'using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;
inline MatrixXd AtA(const MapMatd& A) {
int n(A.cols());
return MatrixXd(n,n).setZero().selfadjointView<Lower>()
.rankUpdate(A.adjoint());
}'
#Inverting a Matrix With RcppEigen
invCpp <- 'using Eigen::Map;
using Eigen::MatrixXd;
// Map the double matrix AA from R
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
// evaluate and return the inverse of A
const MatrixXd At(A.inverse());
return wrap(At);'
rcpp_inv <- cxxfunction(signature(AA='matrix'),invCpp,"RcppEigen")
#Computing the Moore-Penrose Pseudoinverse using RcppArmadillo
rcpp_arm_pinv_src<- '
arma::mat Z = Rcpp::as<arma::mat>(ZS);
arma::mat Y =pinv(Z);
return wrap(Y);'
rcpp_arm_pinv <- cxxfunction(signature(ZS="numeric"),body=rcpp_arm_pinv_src,plugin="RcppArmadillo")
#Computing the Cholesky decomposition of a matrix with RcppEigen
cholCpp <- 'const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const MatrixXd LT((A.llt().matrixU()));
return wrap(LT);'
rcpp_chol <- cxxfunction(signature(AA='matrix'),cholCpp,"RcppEigen",incl)
#RMV_Norm
rcpp_rmvnorm <- function(n,chol_S)
{
	m <- dim(chol_S)[1]
	mat <- matrix(rnorm(n*m),nrow=n,ncol=m)
	return (mat%*%chol_S)
}
#########################
# Load the data
load('dat.RData')
#Compute inverse and pre-allocate samples
N <- 10000
phi <- .1 #parameter which controls the decay of the exponential kernel
X <- training[,5:7]
Y <- as.numeric(training[,1])
dist_mat <- as.matrix(dist(X))
#calculate and store the covariance matrix for the squared exponential kernel
cov_mat <- exp(-(dist_mat^2)/phi) #GP covariance matrix for the squared exponential kernel
cov_mat_inv <- ginv(cov_mat)#precompute the inverse and store for later use
cov_mat_chol <- rcpp_chol(cov_mat) #precompute Cholesky root and store for later use
alpha_samples <- matrix(NA,nrow=N, ncol=3) # N by 3 matrix to store alpha samples
beta_samples <- matrix(NA,nrow=N, ncol=3) # N by 3 matrix to store beta samples
Y_samples <- array(NA, dim=c(length(Y),3,N)) #N by m by 3 array to store all latent samples
alpha_samples[1,] <- c(0,0,0)  
beta_samples[1,] <- c(10,10,10)
Y_samples[,1,1] <- rcpp_rmvnorm(1,cov_mat_chol)
Y_samples[,2,1] <- rcpp_rmvnorm(1,cov_mat_chol)
Y_samples[,3,1] <- rcpp_rmvnorm(1,cov_mat_chol)
#calculate the block covariance matrix for the MVN vector (Y1,Y2,Y3,alpha,beta)
block_cov_matrix <- matrix(0,nrow=(3*length(Y)+6),ncol=(3*length(Y)+6))
block_cov_matrix[1:length(Y),1:length(Y)] <- cov_mat
block_cov_matrix[(length(Y)+1):(2*length(Y)),(length(Y)+1):(2*length(Y))] <- cov_mat
block_cov_matrix[(2*length(Y)+1):(3*length(Y)),(2*length(Y)+1):(3*length(Y))] <- cov_mat
block_cov_matrix[(3*length(Y)+1):(3*length(Y)+6),(3*length(Y)+1):(3*length(Y)+6)] <- diag(6)
block_chol <-rcpp_chol(block_cov_matrix)
norm_samples <- rcpp_rmvnorm(N, block_chol)
#compute uniform samples for ESS
unif_samples <- runif(n=N)
theta <- runif(n=N,min=0,max=2*pi)
theta_min <- theta-2*pi
theta_max <- theta+2*pi
save.image("rcpp_bbh_inverse.RData")
