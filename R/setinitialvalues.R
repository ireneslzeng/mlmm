#'RSTAN code for mlmc() function
#'@description Generate initial value for parameters
#'@import MASS Matrix stats stats4 ggplot2 Rcpp rstan
#'
#'@param npred number of predictors for the completed data regression model
#'@param np number of second level unit i.e. proteins/genes
#'@param npred_miss number of predictors for missingness
#'@param npred_sub number of predictors for subjects
#'@param nmiss number of observations with missing responses values 
#'@param nsid number of subjects
#'@return u_mean initial value of second level unit mean i.e. mean of protein abundance for protein 1
#'@return u_std initial value of second level unit standard deviation i.e. std of protein abundance for protein 1
#'@return beta2_mean initial mean value for sampling unit ,i.e. subject
#'@return beta_Std inital std for sampling unit
#'@return beta2_theta_shape inital shape value for gamma distributed beta_std
#'@return beta2_theta_rate inital rate value for gamma distributed beta_std
#'
#'@examples
#'testexmp=setinitvalues(npred=2,np=3,npred_miss=3,npred_sub=2,nmiss=10,nsid=30)
#'@export

setinitvalues=function(npred=npred, np=np, npred_miss=npred_miss, npred_sub=npred_sub, nmiss=nmiss, nsid=nsid)
{return(list(alpha_response=0.008, pVAR=solve(stats::rWishart(1, df=npred+1, Sigma=as.matrix(Matrix::Diagonal(npred)))[,,1]),
ita=stats::rgamma(1,1,1/10), U_latent=MASS::mvrnorm(np,mu=rep(0,npred), Sigma=as.matrix(Matrix::Diagonal(npred))),
g=stats::rnorm(npred,0,1), alpha_mu=stats::rnorm(npred_miss,0,1), alpha_latent=stats::rnorm(npred_miss,0,1),
alpha_theta=stats::rgamma(npred_miss,1,1), beta2_latent=MASS::mvrnorm(nsid,mu=rep(0,npred_sub),
Sigma=as.matrix(Matrix::Diagonal(npred_sub))), beta2_mu=MASS::mvrnorm(nsid,mu=rep(0,npred_sub),
Sigma=as.matrix(Matrix::Diagonal(npred_sub))), beta2_theta=stats::rgamma(npred_sub,1,1),
y_m_latent=stats::rnorm(nmiss,0,1)))
}

