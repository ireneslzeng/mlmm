#'RSTAN code for mlmc() function
#'@description Generate initial values for parameters
#'@import MASS Matrix stats stats4 ggplot2 Rcpp rstan
#'
#'@param npred number of predictors for the regression model
#'@param np number of protein/metabolite units comprised of the
#'response values (i.e. which
#'represents peptides' ion-intensities used to construct
#'protein/metabolite's abundance)   
#'@param npred_miss number of predictors for missingness
#'@param npred_sub number of predictors for the second level units
#'such as subjects
#'@param nmiss number of observations with missing responses values 
#'@param nsid number of second level units i.e. subjects
#'@param censor_lim_upp upper-limit of censored value of responses 
#'@param ita_a shape parameter for gamma distributed prior-ita(std of
#'response value)
#'@param ita_b rate parameter for gamma distributed prior-ita
#'@param g_mu mean of normal distributed location parameter g for 
#'re-parameterising U (regression coefficient of model in completed data)
#'@param g_sig std of normal distributed location parameter g
#'@param alpha_mu_u mean of normal distributed location parameter alpha_mu for 
#'re-parameterising alpha (regression coefficient of logistic
#'regression model for missing prob) 
#'@param alpha_mu_s std of normal distributed location parameter alpha_mu 
#'@param alpha_theta_a shape parameter of gamma distributed dispersion parameter
#'alpha_theta for re-parameterising alpha
#'@param alpha_theta_b rate parameter of gamma distributed alpha_theta
#'@param beta2_theta_a shape parameter of gamma distributed beta2_theta
# for re-parameterising beta2 (regression coefficient for second level 
#'units,i.e. subject)
#'@param beta2_theta_b rate parameter of gamma distributed dispersion
#'parameter beta2_theta
#'@return pVAR precision matrix for predictors in completed data model
#'@return U_latent standardized multinormal distributed latern
#'variable to re-parameterise regression coefficient U 
#'@return g location parameter to re-parameterise U
#'@return alpha_mu mean value for alpha(regression coefficient of
#'model for missing probability)
#'@return alpha_latent standardized normal distributed latent variable
#'to re-parameterise alpha
#'@return beta2_latent standardized multivariate normal distributed
#'latent variable to re-parameterising beta2
#'@return beta2_mu mean of the multivariate normal distributed beta2
#'@return y_m_latent standardized normal distributed latent variable
#'to re-parameterise response variable  
#'
#'@examples
#'testexmp=setinitvalues(npred=2,np=3,npred_miss=3,npred_sub=2,nmiss=10,
#'nsid=30)
#'@export

setinitvalues=function(npred=npred, np=np, npred_miss=npred_miss,
    npred_sub=npred_sub, nmiss=nmiss,nsid=nsid,
    censor_lim_upp=0.008, ita_a=1, ita_b=1/10, g_mu=0, g_sig=1,
    alpha_mu_u=0,alpha_mu_s=1, alpha_theta_a=1,
    alpha_theta_b=1, beta2_theta_a=1,beta2_theta_b=1)

{return(list(alpha_response=censor_lim_upp,    
pVAR=solve(stats::rWishart(1, df=npred+1, 
Sigma=as.matrix(Matrix::Diagonal(npred)))[,,1]),    
ita=stats::rgamma(1,ita_a,ita_b), 
U_latent=MASS::mvrnorm(np,mu=rep(0,npred),    
Sigma=as.matrix(Matrix::Diagonal(npred))),
g=stats::rnorm(npred,g_mu,g_sig), 
alpha_mu=stats::rnorm(npred_miss,alpha_mu_u,alpha_mu_s),    
alpha_latent=stats::rnorm(npred_miss,0,1),    
alpha_theta=stats::rgamma(npred_miss,alpha_theta_a,alpha_theta_b),    
beta2_latent=MASS::mvrnorm(nsid,mu=rep(0,npred_sub),    
Sigma=as.matrix(Matrix::Diagonal(npred_sub))),     
beta2_mu=MASS::mvrnorm(nsid,mu=rep(0,npred_sub),
Sigma=as.matrix(Matrix::Diagonal(npred_sub))),    
beta2_theta=stats::rgamma(npred_sub,1,1),    
y_m_latent=stats::rnorm(nmiss,0,1)))
}
