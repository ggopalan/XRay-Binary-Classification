#Giri Gopalan (Joint work with Luke Bornn and Saku Vrtilek)


#PURPOSE: The purpose of this R script is to draw from the posterior predictive distribution of the compact star type for a system we are trying to predict the type of, given posterior samples (output of bbh_ess_sampler.r).
#From these posterior predictive draws we allow the scientist to make a decision about what the system type is. (E.g. take the mode of the posterior predictive draws for the compact object type.)


#USAGE: Rscript predictive_distribution_sampler.r input_file_name output_file_name

require(MASS)
require(fields)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(inline)

#Load the test data, which is assumed to be in 'Validation'
load("bbh_ess_10000reps.RData")
args <- commandArgs(trailingOnly = TRUE)
testdata <- read.table(args[1])
names(testdata) <- c("System","Date","CC3","CC1","CC2")

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
cholCpp <- 'const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const MatrixXd LT((A.llt().matrixL()));
return wrap(LT);'
rcpp_chol <- cxxfunction(signature(AA='matrix'),cholCpp,"RcppEigen",incl)

X_test <- testdata[,3:5]
cov_mat_sq <- as.matrix(dist(X_test))
cov_mat_sq <- exp(-(cov_mat_sq^2)/phi)
dist_mat_predict <- as.matrix(dist(rbind(X, X_test)))
cov_mat_pred <- exp(-(dist_mat_predict^2)/phi)
cov_mat_pred <- cov_mat_pred[(length(Y)+1):dim(cov_mat_pred)[1],1:length(Y)]
print(dim(cov_mat_pred))
latent_sample <- matrix(rep(0, 3*dim(cov_mat_pred)[1]),ncol=3)
grid_prediction_mean <- matrix(rep(0,3*dim(cov_mat_pred)[1]),ncol=3)
test <- cov_mat_sq-cov_mat_pred%*%cov_mat_inv%*%t(cov_mat_pred)
new_test <- rcpp_chol(test)
#RcppArmadillo code for calculating the predicted probabilities of compact object types on input points, which are the 3 CCI values at input data. (We use Armadillo to speed up costly matrix multiplication.)
predict_probs_src <-'
Rcpp::NumericVector vecY_samples(Y_samples_);
Rcpp::NumericVector vecp_samples(p_samples_);
Rcpp::IntegerVector dim(dim_);
arma::cube Y_samples(vecY_samples.begin(), dim[3], dim[1], dim[2]);
arma::cube prob_samples(vecp_samples.begin(), dim[2]/100,dim[0],dim[1]);
arma::mat a_samples = Rcpp::as<arma::mat>(a_samples_);
arma::mat b_samples = Rcpp::as<arma::mat>(b_samples_);
arma::mat cov_mat_inv = Rcpp::as<arma::mat>(cov_mat_inv_);
arma::mat cov_mat_pred = Rcpp::as<arma::mat>(cov_mat_pred_);
arma::mat cov_mat_sq = Rcpp::as<arma::mat>(cov_mat_sq_);
arma::mat grid_prediction = Rcpp::as<arma::mat>(grid_prediction_);
arma::mat latent_sample = Rcpp::as<arma::mat>(latent_sample_);
arma::mat grid_prediction_mean = Rcpp::as<arma::mat>(grid_prediction_mean_);
arma::mat chol_grid_prediction_cov = Rcpp::as<arma::mat>(grid_prediction_cov_);
arma::mat Z = Rcpp::as<arma::mat>(Z_);
double sum = 0;
for(int i = 1; i < dim[2]/100; i++)
{
    Rcout << i;
    grid_prediction_mean = cov_mat_pred*cov_mat_inv*Y_samples.slice(100*i);
    for(int k = 0; k < dim[1]; k++)
    {
        Z = arma::randn(dim[0]);
        latent_sample.col(k) = grid_prediction_mean.col(k)+chol_grid_prediction_cov*Z;
    }
	for(int j = 0; j < dim[0]; j++)
	{
        sum = 0;
        for(int k = 0; k < dim[1]; k++)
		{
            prob_samples(i,j,k) = exp(a_samples(100*i,k)+b_samples(100*i,k)*latent_sample(j,k));
			sum += prob_samples(i,j,k);
		}
        for(int k = 0; k < dim[1]; k++)
		{
			prob_samples(i,j,k) = prob_samples(i,j,k)/sum;
		}
	}
}
return wrap(prob_samples);
'
predict_probs <- cxxfunction(signature(Y_samples_="numeric",p_samples_="numeric",dim_="integer",a_samples_="numeric",b_samples_="numeric",cov_mat_inv_="numeric",cov_mat_pred_="numeric",cov_mat_sq_="numeric",grid_prediction_="matrix",latent_sample_="numeric",grid_prediction_mean_="numeric",grid_prediction_cov_="numeric",Z_ = "numeric"),body=predict_probs_src,plugin="RcppArmadillo")
#Cal to auxillary Cpp function
prob_samples <- predict_probs(as.vector(Y_samples),rep(0,dim(cov_mat_pred)[1]*N*3/100),c(dim(cov_mat_pred)[1],3,N,dim(cov_mat_inv)[1]),alpha_samples,beta_samples,cov_mat_inv,cov_mat_pred,cov_mat_sq,matrix(rep(0,3*dim(cov_mat_pred)[1]),ncol=3),latent_sample,grid_prediction_mean,new_test,cov_mat_sq)
#Use the posterior-predictive probabilities to draw from the multinomial distribution
unif_draws <- matrix(runif(dim(cov_mat_pred)[1]*100),ncol=100)
post_pred_samples <- matrix(rep(0,dim(cov_mat_pred)[1]*100),ncol=100)
for(i in 1:100)
{
    for(j in 1:dim(cov_mat_pred)[1])
    {
        if(unif_draws[j,i] < prob_samples[i,j,1])
        {
            post_pred_samples[j,i] = 1
        }
        else if((prob_samples[i,j,1] <= unif_draws[j,i]) & (unif_draws[j,i] < prob_samples[i,j,1]+prob_samples[i,j,2]))
       {
            post_pred_samples[j,i] = 2
        }
        else
        {
            post_pred_samples[j,i] = 3
        }
    }
}
#Remove the burn-in
post_pred_samples <- post_pred_samples[,2:100]
#NOTE: to predict a system type, we take the posterior predictive mode
#Save the predictions
save(post_pred_samples,file=args[2])
