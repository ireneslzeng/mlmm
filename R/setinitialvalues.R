#'The function to set initial values for parameters: setinitvalues().
#'@description Generate initial values for parameters
#'@import Mass Matrix stats stats4 ggplot2 Rcpp rstan
#'@importFrom stats rWishart rgamma rnorm
#'@importFrom Matrix Diagonal 
#'@importFrom MASS mvrnorm    
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
#'@param censor_lim_upp upper-limit of censored value of responses.
#'The default value is 0.001 according to an experiment device. User can
#'change it according to the data. 
#'@param ita_a shape parameter for gamma distributed prior-ita (std of
#'response value).  The default is set to 1. 
#'@param ita_b rate parameter for gamma distributed prior-ita. The default
#'is set to 1/10. The default values of shape and rate parameters provide
#'a reasonable wide range of initial value for ita. Users can change it
#'accordingly. 
#'@param g_mu mean of normal distributed location parameter g for 
#'re-parameterising U (regression coefficient of model in the
#'completed data). The default value is set to 0 for the mean of standard
#'normal distribution.   
#'@param g_sig std of normal distributed location parameter g. The default
#'is set to 1 for the std of normal distribution.
#'@param alpha_mu_u mean of normal distributed location parameter alpha_mu
#'for re-parameterising alpha (regression coefficient of logistic
#'regression model for missing prob). The default is set to 0.  
#'@param alpha_mu_s std of normal distributed location parameter alpha_mu 
#'@param alpha_theta_a shape parameter of gamma-distributed dispersion
#'parameter alpha_theta for re-parameterising alpha. Default value is set
#'to 1 as a natural starting value. 
#'@param alpha_theta_b rate parameter of gamma-distributed alpha_theta.
#'Default value uses 1/10 as for ita_b. Both default values of shape(_a) and
#'rate(_b) of alpha_theta can be changed to give a wider range (_b=1/10)
#'or a narrower range (_b=0.5).  
#'@param beta2_theta_a shape parameter of gamma-distributed beta2_theta
#'for re-parameterizing beta2 (regression coefficient for second level 
#'units,i.e. subject). Default value uses 1. 
#'@param beta2_theta_b rate parameter of gamma distributed dispersion
#'parameter beta2_theta. Default value used 1/10, same as for ita_b.
#'@return pVAR precision matrix for predictors in completed data model
#'@return U_latent standardized multinormal distributed latern
#'variable to re-parameterise regression coefficient U. 
#'@return g location parameter to re-parameterise U.
#'@return alpha_mu mean value for alpha(regression coefficient of
#'model for missing probability).
#'@return alpha_latent standardized normal distributed latent variable
#'to re-parameterize alpha.
#'@return beta2_latent standardized multivariate normal distributed
#'latent variable to re-parameterising beta2.
#'@return beta2_mu mean of the multivariate normal distributed beta2
#'@return y_m_latent standardized normal distributed latent variable
#'to re-parameterise response variable.  
#'    
#'@examples
#'testexmp=setinitvalues(npred=2,np=3,npred_miss=3,npred_sub=2,nmiss=10,
#'nsid=30)
#'@export

setinitvalues=function(npred, np, npred_miss,
npred_sub, nmiss, nsid,
censor_lim_upp=0.008, ita_a=1, ita_b=1/10, g_mu=0, g_sig=1,
alpha_mu_u=0,alpha_mu_s=1, alpha_theta_a=1,
alpha_theta_b=1/10, beta2_theta_a=1,beta2_theta_b=1/10)
    
{return(list(alpha_response=censor_lim_upp,    
pVAR=solve(rWishart(1, df=npred+1, 
Sigma=as.matrix(Diagonal(npred)))[,,1]),    
ita=rgamma(1,ita_a,ita_b), 
U_latent=mvrnorm(np,mu=rep(0,npred),    
Sigma=as.matrix(Diagonal(npred))),
g=rnorm(npred,g_mu,g_sig), 
alpha_mu=rnorm(npred_miss,alpha_mu_u,alpha_mu_s),    
alpha_latent=rnorm(npred_miss,0,1),    
alpha_theta=rgamma(npred_miss,alpha_theta_a,alpha_theta_b),    
beta2_latent=mvrnorm(nsid,mu=rep(0,npred_sub),    
Sigma=as.matrix(Diagonal(npred_sub))),     
beta2_mu=mvrnorm(nsid,mu=rep(0,npred_sub),
Sigma=as.matrix(Diagonal(npred_sub))),    
beta2_theta=rgamma(npred_sub,1,1),    
y_m_latent=rnorm(nmiss,0,1)))
}
