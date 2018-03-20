//RSTAN code for mlmc() function
//include "license.stan" // GPL2+
data
 {
  	int<lower=0> nobs;
	int ncensor;
	int nmiss;
	int npred;
	int npred_miss;
	int npred_sub;
	int nsid;
  	int np;
  	int respond_dep;
 	int<lower=0> sid[nobs];
	int<lower=0> sid_m[nmiss];
 	int<lower=0> pid[nobs];
	int<lower=0> pid_m[nmiss];

	matrix[nobs,npred] pred;
	matrix[nmiss,npred] pred_m;
  	matrix[nobs,npred_miss] pred_miss;
  	matrix[nmiss,npred_miss] pred_miss_m;
	matrix[nobs,npred_sub] pred_sub;
  	matrix[nmiss,npred_sub] pred_sub_m;
  	real y[nobs];
  	int miss_m[nmiss];
  	int miss_obs[nobs];
	corr_matrix[npred] R;
	vector[npred] Sigma_sd;
	cov_matrix[npred] prec;
	vector[npred]  mn;

  //Define data for censored responses
  	real censor_lim;
  	int<lower=0> sid_c[ncensor];
  	int<lower=0> pid_c[ncensor];
  	matrix[ncensor,npred] pred_c;
  	matrix[ncensor,npred_sub] pred_sub_c;
 }

transformed data
   {	cov_matrix[npred] T;
    	cov_matrix[npred] invprec;
    	T=diag_matrix(Sigma_sd)*R*diag_matrix(Sigma_sd);
    	invprec=inverse(prec);
   }

parameters
  { 	matrix[np,npred] U_latent;
    	row_vector[npred] g;
    	cov_matrix[npred] pVAR;
    	real<lower=0> ita;

    	vector[npred_miss] alpha_latent;
    	vector[npred_miss] alpha_mu;
    	vector[npred_miss] alpha_theta;
    	real alpha_response;

    	matrix[nsid,npred_sub] beta2_latent;
    	row_vector[npred_sub] beta2_theta;
    	matrix[nsid,npred_sub] beta2_mu;

    	real y_m_latent[nmiss];
  }

transformed parameters
  { 	matrix[np,npred] U;
    	matrix[nsid,npred_sub] beta2;
    	vector[npred_miss] alpha;
    	real mu_m[nmiss];
    	real mu[nobs];
    	real mu_c[ncensor];
    	real y_m[nmiss];
    	real<lower=0,upper=1> pmiss[nobs];
    	real<lower=0,upper=1> pmiss_m[nmiss];

    for (sub in 1:nsid)
    beta2[sub]=beta2_mu[sub]+dot_product(beta2_theta,beta2_latent[sub]);

    for (prot in 1:np)
    U[prot]=g+U_latent[prot]*pVAR;

    alpha=alpha_mu+dot_product(alpha_theta,alpha_latent);  		

    for (pep in 1:nobs)
      {mu[pep]=dot_product(pred_sub[pep],beta2[sid[pep]])+dot_product(pred[pep],U[pid[pep]]);
       if (respond_dep==1) pmiss[pep]=inv_logit(dot_product(alpha,pred_miss[pep])+alpha_response*mu[pep]);
       else pmiss[pep]=inv_logit(dot_product(alpha,pred_miss[pep]));
       if (pmiss[pep]==0)  pmiss[pep]=0.001;
       if (pmiss[pep]==1)  pmiss[pep]=0.999;
	    }

    for (pep2 in 1:nmiss)
      {mu_m[pep2]=dot_product(beta2[sid_m[pep2]],pred_sub_m[pep2])+dot_product(pred_m[pep2],U[pid_m[pep2]]);
       y_m[pep2]=mu_m[pep2]+y_m_latent[pep2]*ita;
       if (respond_dep==1) pmiss_m[pep2]=inv_logit(dot_product(alpha,pred_miss_m[pep2])+alpha_response*mu_m[pep2]);
	     else pmiss_m[pep2]=inv_logit(dot_product(alpha,pred_miss[pep2]));
       if (pmiss_m[pep2]==0)  pmiss_m[pep2]=0.001;
       if (pmiss_m[pep2]==1)  pmiss_m[pep2]=0.999;
	    }
    for (pep3 in 1:ncensor)
      {mu_c[pep3]=dot_product(beta2[sid_c[pep3]],pred_sub_c[pep3])+dot_product(U[pid_c[pep3]],pred_c[pep3]);
 	     if (mu_c[pep3]> censor_lim) mu_c[pep3]=censor_lim;
       }
   }

model
  {
  	for (sub in 1:nsid)
  	{beta2_latent[sub]~normal(0,1);
	 beta2_mu[sub]~normal(0,1);}
	 beta2_theta~gamma(1,1);
   	 g~multi_normal(mn,T);
	 pVAR~inv_wishart(npred,invprec);

	 for (prot in 1:np)
	 U_latent[prot]~multi_normal(mn,R);

	 ita~gamma(1,1);
   	 y~normal(mu,ita);

	 alpha_latent~normal(0,1);
	 alpha_mu~normal(0,1);
   	 alpha_theta~gamma(1,1);
	 for (pep in 1:nobs)
   		miss_obs[pep]~bernoulli(pmiss[pep]);

	 for (pep2 in 1:nmiss)
	  {miss_m[pep2]~bernoulli(pmiss_m[pep2]);
    	   y_m_latent[pep2]~normal(0,1);}

  if (ncensor>0)
    	{for (pep in 1:ncensor)
   	target +=log(Phi((censor_lim-mu_c[pep])/ita)+0.001);}
  	}
