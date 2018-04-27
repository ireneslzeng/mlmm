#'mlmm() function. 
#'
#'@description mlmm() handles Bayesian multilevel model with response
#'variable 'that has missing values that depends on the response value
#'itself.  'Apart from the response value, the missingness is also
#'known to 'associate with the other variables.  The method is created
#'for analysing 'mass-spectrometry data when it has abundance-dependant
#'missing and censored 'values, and there are prior information
#'available for the associations 'between the probability of missing
#'and the known variables.  'The function mlmm is written for response
#'variable has no censored values 'while mlmc function include imputing
#'censored values.
#'
#'@import MASS Matrix Rcpp stats4 ggplot2
#'
#'@param formula_completed The main regression model formula. 
#'It has the same formula format as lmr() and it is used to 
#'define the first level response and its explanatory variables.
#'
#'@param formula_missing The logistic regression model formula. 
#'It has the same formula as formula_completed.
#'
#'@param formula_subject The second level formula in the multilevel
#'model 'which is used to define second level unit such as subject and
#'explanatory 'variables.
#'
#'@param pdata The dataset contains response and predictors in a long
#'format.  Response is a vector with an indictor variable to define
#'the corresponding unit. The data needs to have the following
#'rudimental variables: the indicator variable for first level
#'response, the indicator variable for second level unit such as
#'subject or a sampling unit, an indicator for missingness and
#'indictor of censoring.  Missingness and censor are two different
#'classification, there should not have any overlap between
#'missingness and censored.  Data structure can be referenced from the
#'example and reference papers.
#'
#'@param respond_dep_missing An indicator of whether response value is 
#'missing-dependant.
#'
#'@param pidname Variable name to define the multilevel response unit, 
#'i.e. protein name or gene name
#'
#'@param sidname Vriable name to define the subject unit,
#'i.e. patient id or sampling id
#'
#'@param prec_prior prior precision matrix of the explanatory
#'variables at the first level unit in the multilevel model, for
#'example, the variables to predict ion intensity. Its dimension will
#'be q x q, where q is the number of explanatory variables at the
#'right-hand side of formula_completed. The default is a matrix with
#'diagonal value of 0.01 and off-diagonal value of 0.005.
#'
#'@param alpha_prior prior for coefficients of predictors to missing
#'probability in the logistic regression. Its length will be equal to
#'the number of variables at the right-hand side of the
#'formula_missing.  Default is a vector of zeros.
#'
#'@param iterno Number of iterations for the posterior samplings
#'
#'@param chains rstan() parameter to define number of chains of
#'posterior samplings.
#'
#'@param thin rstan() parameter to define the frequency of iterations
#'saved.

#'@param seed random seed for rstan() function
#'
#'@param algorithm rstan() parameter which has three options c(NUTS,HMC,
#'Fixed_param).
#'
#'@param warmup Number of iterations for burn-out in stan.
#'
#'@param adapt_delta_value Adaptive delta value is an adaptation
#'parameters 'for sampling algorithms,default is 0.85, value between
#'0-1.

#'@param savefile A logical variable to indicate if the sampling files
#'are to be saved.
#'
#'@param usefit A logical variable to indicate if the model use the existing 
#'fit.
#'
#'@return Return of the function is the result fitted by stan(). It
#'will have 'the summarized parameters from all chains and summary
#'results for each 'chain.
#'Plot() function will return the visualization of the mean and parameters.
#'@examples
#'#\dontrun{
#'library(MASS)
#'set.seed(150)
#'var2=abs(rnorm(1000,0,1));treatment=c(rep(0,500),rep(1,500))
#'geneid=rep(seq_len(20),50);sid=c(rep(seq_len(25),20),rep(seq_len(25)+25,20))
#'cov1=rWishart(1,df=100,Sigma=diag(rep(1,100)))
#'u=rnorm(100,0,1)
#'mu=mvrnorm(n=1,mu=u,cov1[,,1])
#'sdd=rgamma(1,shape=1,scale=1/10)
#'var1=(1/0.85)*var2+2*treatment
#'for (i in seq_len(1000)) {var1[i]=var1[i]+rnorm(1,mu[geneid[i]],sdd)}
#'miss_logit=var2*(-0.9)+var1*(0.01)
#'probmiss=exp(miss_logit)/(exp(miss_logit)+1)
#'miss=rbinom(1000,1,probmiss);table(miss)
#'pdata=data.frame(var1,var2,treatment,miss,geneid,sid)
#'for ( i in seq_len(1000)) if (pdata$miss[i]==1) pdata$var1[i]=NA;
#'pidname="geneid";sidname="sid";
#'#copy and paste the following formulas to the mmlm() function respectively
#'formula_completed=var1~var2+treatment
#'formula_missing=miss~var2
#'formula_censor=censor~1
#'formula_subject=~treatment
#'model3=mlmm(formula_completed=var1~var2+treatment,formula_missing=miss~var2,
#'formula_subject=~treatment,pdata=pdata,respond_dep_missing=TRUE,
#'pidname="geneid",sidname="sid",iterno=10,chains=2,usefit=TRUE,
#'savefile=FALSE)
#'#}
#'@export

mlmm=function(formula_completed, formula_missing, formula_subject, 
pdata, respond_dep_missing=TRUE, pidname, sidname, prec_prior=NULL, 
alpha_prior=NULL, iterno=100, chains=3, thin=1, seed=125, 
algorithm="NUTS", warmup=floor(iterno/2), adapt_delta_value=0.90, 
savefile=FALSE, usefit=FALSE)

{current.na.action=options('na.action');options(na.action='na.pass')

    t=stats::terms(formula_completed)
    mf=stats::model.frame(t,pdata,na.action='na.pass')
    mm=stats::model.matrix(mf,pdata);t2=stats::terms(formula_missing)
    mf2=stats::model.frame(t2,pdata,na.action='na.pass')
    mm2=stats::model.matrix(mf2,pdata);missing=stats::model.response(mf2)

    if (!is.null(formula_subject)) {tt3=stats::terms(formula_subject)
    mf3=stats::model.frame(tt3,pdata,na.action='na.pass')
    mm3=stats::model.matrix(mf3,pdata)}
    y_all=stats::model.response(mf)

options('na.action' = current.na.action$na.action)

#######################Prepare data;
    ns=length(y_all)
    nmiss=table(missing)[2]
    nobs=ns-nmiss
    datamiss=subset(pdata,(missing==1))
    dataobs=subset(pdata,(missing!=1))
    npred=dim(mm)[2];npred_miss=dim(mm2)[2];npred_sub=dim(mm3)[2]
    sid=dataobs[,sidname];sid_m=datamiss[,sidname]
    nsid=length(table(pdata[,sidname]))
    pid=dataobs[,pidname]
    pid_m=datamiss[,pidname]
    np=length(table(pdata[,pidname]))
    pred=mm[(missing!=1),]
    pred_miss=mm2[(missing!=1),]
    pred_sub=mm3[(missing!=1),]
    pred_m=mm[(missing==1),]
    pred_miss_m=mm2[(missing==1),]
    pred_sub_m=mm3[(missing==1),]
    miss_m=datamiss[,colnames(mf2)[1]]
    miss_obs=dataobs[,colnames(mf2)[1]]
    if (respond_dep_missing) respond_dep=1 else respond_dep=0

#data for prior;
    R=as.matrix(Matrix::Diagonal(npred))
#Assign default prior for precision matrix
#of explanatory variables at the first leve;
    if (is.null(prec_prior)){ prec_prior=matrix(nrow=npred,ncol=npred)
    for (i in seq_len(npred))
        for (j in seq_len(npred)) {
        if (i==j) prec_prior[i,j]=0.1 else prec_prior[i,j]=0.005}
        }
    mn=rep(0,npred)
    Sigma_sd=rep(10,npred)
#Assign default prior value for assoication of missing prob and variables

    if (is.null(alpha_prior)) alpha_prior=rep(1,npred_miss)
    
    prstan_data=list(
    nobs=nobs,nmiss=nmiss,nsid=nsid,np=np,npred=npred,
    npred_miss=npred_miss,npred_sub=npred_sub,
    respond_dep=respond_dep,y=y_all[(missing!=1)],
    sid=sid,pid=pid,pred=pred,pred_miss=pred_miss,pred_sub=pred_sub,
    miss_obs=miss_obs,
    sid_m=sid_m,pid_m=pid_m,pred_m=pred_m,
    pred_miss_m=pred_miss_m,pred_sub_m=pred_sub_m,miss_m=miss_m,
    R=R,Sigma_sd=Sigma_sd,mn=mn,
    prec_prior=prec_prior,alpha_prior=alpha_prior)
    if (respond_dep==1) 
    parsstr=c("U","beta2","alpha","alpha_response") else {
    parsstr=c("U","beta2","alpha")
    }
    mlmm_path=.libPaths("mlmm")
    initvalue1=function () {
    setinitvalues(
    npred=npred,np=np,
    npred_miss=npred_miss,
    npred_sub=npred_sub,nmiss=nmiss,
    nsid=nsid
    )
    }
    
    if (usefit==TRUE){
        stanfit <- stanmodels$mmlm_code
        if (savefile==TRUE) fitmlmm=rstan::sampling(stanfit,data=prstan_data,
        iter=iterno, pars=parsstr, seed=seed, thin=thin,
        algorithm=algorithm, warmup=warmup, chains=chains,
        control=list(adapt_delta=adapt_delta_value),
        sample_file=file.path(getwd(),"samples")) else {
        fitmlmm=rstan::sampling(stanfit,data=prstan_data,iter=iterno,
        pars=parsstr, seed=seed, thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),
        sample_file=NULL)}
        }  else {
        if (savefile==TRUE) 
        fitmlmm=rstan::stan(file=file.path(mlmm_path,
        "mlmm/exec/mmlm_code.stan"),
        data=prstan_data,iter=iterno,
        init=initvalue1,pars=parsstr,seed=seed,thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),
        sample_file=file.path(getwd(),"samples"),save_dso=FALSE) else {
        fitmlmm=rstan::stan(file=file.path(mlmm_path,
        "mlmm/exec/mmlm_code.stan"),
        data=prstan_data,iter=iterno,
        init=initvalue1,pars=parsstr,seed=seed,thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),sample_file=NULL,
        save_dso=FALSE)
        }
        }
    print(fitmlmm)
    if (savefile==TRUE) 
    utils::write.csv(as.array(fitmlmm),file=file.path(getwd(),
    "outsummary.csv"),
    row.names=TRUE)
    return(fitmlmm)
}
