#'mlmc() function.
#'@useDynLib rstanarm, .registration = TRUE
#'@description mlmc() handles Bayesian multilevel model with response
#'variable 'that has left-censored values, and missing values that
#'depends on the 'response value itself. Apart from the response value,
#'the missingness is 'also known to associate with the other variables.
#'The method is created for 'analysing mass-spectrometry data when it
#'has abundance-dependant missing and 'censored values, and there are
#'prior information available for the 'associations between the
#'probability of missing and the known variables.  'The imputed values
#'for the censored response are outputed as part of the 'parameters.
#'
#'@import methods MASS Matrix Rcpp stats4 ggplot2 
#'
#'@param formula_completed The main regression model formula; It has
#'the same 'formula format as lmr() and it is used to define the first
#'level response 'and its explanatory variables.
#'
#'@param formula_missing  The logistic regression model formula; 
#'It has the same formula as formula_completed.
#'
#'@param formula_censor The formula used in the program to define the 
#'observations with censored values.
#'
#'@param formula_subject The second level formula in the multilevel
#'model 'which is used to define responses such as subject and 'its
#'explanatory variables.
#'
#'@param pdata The dataset contains response and predictors in a long
#'format.  Response is a vector with an indictor variable to define
#'the corresponding 'unit. The data needs to have the following
#'rudimental variables: 'the indicator variable for first level
#'response, 'second level indicator variable for subject 'such as
#'subject id or a sampling unit, 'an indicator for missingness and
#'indictor of censoring.  'Missingness and censor are two different
#'classifications, 'when these two variables are tabulated, there must
#'not have any observation 'defined as censored and missing.  'Data
#'structure can be referred from the example and vignette.
#'
#'@param respond_dep_missing A logical variable to indicate whether 
#'response value is missing-dependant.
#'
#'@param response_censorlim The detectable limit for the response value, 
#'i.e. 1 mg per Liter for intensity value.
#'
#'@param pidname Variable name to define the multilevel response unit, 
#'i.e. protein name or gene name.
#'
#'@param sidname Variable name to define the subject unit, 
#'i.e. patient id or sampling id.
#'
#'@param prec_prior prior precision matrix of the explanatory
#'variables 'at the first level unit in the multilevel model, for
#'example, the variables 'to predict ion intensity. Its dimension will
#'be q x q, where q is the number 'of explanatory variables at the
#'right-hand side of formula_completed.The 'default is a matrix with
#'diagonal values of 0.01 and off-diagonal values of '0.005.
#'
#'@param alpha_prior prior for coefficients of predictors to missing
#'probability in the logistic regression. Its length will be equal to
#'the 'number of variables at the right-hand side of the
#'formula_missing.  Default is a vector of zeros.
#'
#'@param iterno Number of iterations for the posterior samplings.
#'
#'@param chains rstan parameter to define number of chains of
#'posterior 'samplings.
#'
#'@param thin rstan parameter to define the frequency of iterations
#saved.
#'
#'@param seed random seed for rstan function.
#'
#'@param algorithm rstan parameter which has three options NUTS, HMC, 
#'Fixed param.
#'
#'@param warmup Number of iterations for burn-out in stan.
#'
#'@param adapt_delta_value Adaptive delta value is an adaptation
#'parameters for sampling algorithms, default is 0.85, value between
#'0-1.
#'
#'@param savefile A logical variable to indicate if the sampling files
#'are to 'be saved.
#'
#'@param usefit A logical variable to indicate if the model use the
#'existin gfit.
#'
#'@return Return of the function is the result fitted by stan. It will
#'have 'the summarized parameters from all chains and summary results
#'for each chain.
#'
#'
#'@examples
#'\dontshow{
#'testfun=function()
#'{set.seed(100)
#'var2=abs(rnorm(50,0,1))
#'var1=(1/0.85)*var2;
#'geneid=rep(c(1:10),5);sid=c(rep(c(1:25),2),rep(c(26:50),2));
#'censor=rep(0,50)
#'pidname="geneid";sidname="sid";
#'miss_logit=var2*(-0.9);
#'miss=rbinom(50,1,exp(miss_logit)/(exp(miss_logit)+1));
#'censor=rep(0,50)
#'for (i in seq_len(50)) {if (var1[i]<0.05) censor[i]=1}
#'for ( i in seq_len(50)) {if ((miss[i]==1) & (censor[i]==1)) miss[i]=0};
#'for ( i in seq_len(50)) {if (miss[i]==1) var1[i]=NA;
#'if (censor[i]==1) var1[i]=0.05}
#'pdata=data.frame(var1,var2,miss,censor,geneid,sid)
#'model1=mlmc(
#'formula_completed=var1~var2,formula_missing=miss~var2,
#'formula_censor=censor~1,formula_subject=~var2,
#'pdata=pdata,response_censorlim=0.05,
#'respond_dep_missing=TRUE,pidname="geneid",sidname="sid",
#'iterno=10,chains=2,savefile=FALSE,usefit=TRUE)
#'}
#'system.time(testfun())
#'}
#'\dontrun{
#'set.seed(150)
#'library(MASS)
#'var2=abs(rnorm(800,0,1));treatment=c(rep(0,400),rep(1,400));
#'var1=(1/0.85)*var2+2*treatment;
#'geneid=rep(seq_len(50),16);
#'sid=c(rep(seq_len(25),16),rep(seq_len(50)+25,16))
#'cov1=rWishart(1,df=50,Sigma=diag(rep(1,50)))
#'u=rnorm(50,0,1);mu=mvrnorm(n=1,mu=u,cov1[,,1])
#'sdd=rgamma(1,shape=1,scale=1/10);
#'for (i in seq_len(800)) {var1[i]=var1[i]+rnorm(1,mu[geneid[i]],sdd)}
#'miss_logit=var2*(-0.9)+var1*(0.001);
#'miss=rbinom(800,1,exp(miss_logit)/(exp(miss_logit)+1));
#'censor=rep(0,800)
#'for (i in seq_len(800)) {if (var1[i]<0.002) censor[i]=1}
#'pdata=data.frame(var1,var2,treatment,miss,censor,geneid,sid);
#'for ( i in seq_len(800)) 
#'{if ((pdata$miss[i]==1) & (pdata$censor[i]==1)) pdata$miss[i]=0};
#'for ( i in seq_len(800)) {
#'if (pdata$miss[i]==1) pdata$var1[i]=NA;
#'if (pdata$censor[i]==1) pdata$var1[i]=0.002};
#'pidname="geneid";sidname="sid";
#'#copy and paste the following formulas to the mlmm() function respectively
#'formula_completed=var1~var2+treatment;
#'formula_missing=miss~var2;
#'formula_censor=censor~1;
#'formula_subject=~treatment;
#'response_censorlim=0.002;
#'
#'model1=mlmc(formula_completed=var1~var2+treatment,
#'formula_missing=miss~var2,
#'formula_censor=censor~1,formula_subject=~treatment,pdata=pdata,
#'response_censorlim=0.002,
#'respond_dep_missing=TRUE,pidname="geneid",sidname="sid",
#'iterno=50,chains=2,savefile=FALSE,usefit=FALSE)
#'}
#'@export

mlmc=function(formula_completed,formula_missing,
    formula_censor=NULL,formula_subject,
    pdata,respond_dep_missing=TRUE,
    response_censorlim=NULL,
    pidname,sidname,
    prec_prior=NULL,alpha_prior=NULL,
    iterno=100,chains=3,thin=1,seed=125,
    algorithm="NUTS",warmup=floor(iterno/2),
    adapt_delta_value=0.90,
    savefile=FALSE,usefit=FALSE)

{current.na.action=options('na.action');options(na.action='na.pass')
    
    t=stats::terms(formula_completed)
    mf=stats::model.frame(t,pdata,na.action='na.pass')
    mm=stats::model.matrix(mf,pdata)
    
    t2=stats::terms(formula_missing)
    mf2=stats::model.frame(t2,pdata,na.action='na.pass')
    mm2=stats::model.matrix(mf2,pdata)
    missing=stats::model.response(mf2)
    
    if (!is.null(formula_subject))
    {tt3=stats::terms(formula_subject);
    mf3=stats::model.frame(tt3,pdata,na.action='na.pass');
    mm3=stats::model.matrix(mf3,pdata)} else {
    stop("No subject formula defined")}
    if (!is.null(formula_censor)){
    tt4=stats::terms(formula_censor);
    mf4=stats::model.frame(tt4,pdata,na.action='na.pass');
    mm4=stats::model.matrix(mf4,pdata);
    censor=stats::model.response(mf4)
    } else ncensor=0;
    
    if (!is.null(response_censorlim))
    censor_lim=response_censorlim else {
    stop("No response censor limit defined")
    }
    y_all=stats::model.response(mf)
    options('na.action' = current.na.action$na.action)
    #unit test detect errors for censor and missing
    checkt=table(missing,censor)
    if (checkt[4]>0) 
    stop("responses have overlapped definition in censor and missing");
    #######################Prepare data
    ns=length(y_all)
    if (!is.null(formula_censor)) ncensor=table(censor)[2]
    nmiss=table(missing)[2]
    nobs=ns-ncensor-nmiss
    datamiss=subset(pdata,(missing==1))
    dataobs=subset(pdata,(missing!=1 & censor!=1))
    datacen=subset(pdata,(censor==1))
    if (is.null(dim(mm))) npred=1 else npred=dim(mm)[2]
    if (is.null(dim(mm2))) npred_miss=1 else npred_miss=dim(mm2)[2]
    if (is.null(dim(mm3))) npred_sub=1 else npred_sub=dim(mm3)[2]
    sid=dataobs[,sidname]
    sid_m=datamiss[,sidname]
    nsid=length(table(pdata[,sidname]))
    pid=dataobs[,pidname]
    pid_m=datamiss[,pidname]
    np=length(table(pdata[,pidname]))
    pred=as.matrix(mm[(missing!=1 & censor!=1),]) 
    pred_miss=as.matrix(mm2[(missing!=1 & censor!=1),])
    pred_sub=as.matrix(mm3[(missing!=1 & censor!=1),])
    pred_m=as.matrix(mm[(missing==1),])
    pred_miss_m=as.matrix(mm2[(missing==1),])
    pred_sub_m=as.matrix(mm3[(missing==1),])
    if ((!is.null(formula_censor)))
    {npred_c=dim(mm4)[2]
    sid_c=datacen[,sidname]
    pid_c=datacen[,pidname]
    np_c=length(table(pid_c))
    pred_c=as.matrix(mm[(censor==1),])
    pred_sub_c=as.matrix(mm3[(censor==1),])
    }
    miss_m=datamiss[,colnames(mf2)[1]]
    miss_obs=dataobs[,colnames(mf2)[1]]
    if (respond_dep_missing) respond_dep=1 else respond_dep=0
    #data for prior
    R=as.matrix(Matrix::Diagonal(npred))
    #Assign default prior for precision matrix 
    #of explnatory variables at the first leve
    if (is.null(prec_prior)){
    prec_prior=matrix(nrow=npred,ncol=npred)
    for ( i in seq_len(npred))
    for (j in seq_len(npred))
    {if (i==j) prec_prior[i,j]=0.1 else prec_prior[i,j]=0.005}
    }
    mn=rep(0,npred)
    Sigma_sd=as.vector(rep(10,npred))
    #Assign default prior value for assoication between missing prob 
    #& variables
    if (is.null(alpha_prior)) alpha_prior=rep(0,npred_miss)
    
    if ((!is.null(formula_censor)))
    {prstan_data=list(
    nobs=nobs, nmiss=nmiss, ncensor=ncensor, nsid=nsid, np=np,
    npred=npred, npred_miss=npred_miss, npred_sub=npred_sub,
    censor_lim=censor_lim, respond_dep=respond_dep,
    y=y_all[(missing!=1 & censor!=1)], sid=sid, pid=pid, 
    pred=pred, pred_miss=pred_miss, pred_sub=pred_sub, miss_obs=miss_obs,
    sid_m=sid_m, pid_m=pid_m, pred_m=pred_m, 
    pred_miss_m=pred_miss_m, pred_sub_m=pred_sub_m, 
    miss_m=miss_m, sid_c=sid_c, pid_c=pid_c, pred_c=pred_c, 
    pred_sub_c=pred_sub_c, npred_c=npred_c, 
    R=R,Sigma_sd=Sigma_sd, mn=mn,
    prec_prior=prec_prior,alpha_prior=alpha_prior)
    } else {
    prstan_data=list(
    nobs=nobs, nmiss=nmiss, ncensor=ncensor, nsid=nsid, np=np,
    npred=npred, npred_miss=npred_miss, npred_sub=npred_sub, 
    censor_lim=censor_lim, respond_dep=respond_dep, 
    y=y_all[(missing!=1 & censor!=1),], sid=sid, pid=pid,
    pred=pred, pred_miss=pred_miss, pred_sub=pred_sub, miss_obs=miss_obs, 
    sid_m=sid_m, pid_m=pid_m, pred_m=pred_m, 
    pred_miss_m=pred_miss_m, pred_sub_m=pred_sub_m, miss_m=miss_m,
    sid_c=0, pid_c=0,
    pred_c=matrix(nrow=1,ncol=1), pred_sub_c=matrix(nrow=1,ncol=1), npred_c=0,
    R=R,Sigma_sd=Sigma_sd,mn=mn,
    prec_prior=prec_prior,alpha_prior=alpha_prior)
    }

    if (respond_dep==1) parsstr=c("U","beta2","alpha","alpha_response") else {
    parsstr=c("U","beta2","alpha")
    }
    
    mlmm_path=.libPaths("mlmm")
    initvalue1=function () {
    setinitvalues(npred=npred,np=np,npred_miss=npred_miss,npred_sub=npred_sub,
    nmiss=nmiss,nsid=nsid)}
    if (usefit==TRUE) { 
        stanfit=stanmodels$mlmc_code
        if (savefile==TRUE)  
        fitmlmc=rstan::sampling(stanfit,data=prstan_data,iter=iterno,
        init=initvalue1,
        pars=parsstr,seed=seed,thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),
        sample_file=file.path(getwd(),"samples")) else {
        fitmlmc=rstan::sampling(stanfit,data=prstan_data,iter=iterno,
        init=initvalue1,
        pars=parsstr, seed=seed, thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),sample_file=NULL)
        }
        }  else {  
        if (savefile==TRUE) 
        fitmlmc=rstan::stan(file=file.path(mlmm_path,
        "mlmm/exec/mlmc_code.stan"),
        data=prstan_data, iter=iterno,init=initvalue1,pars=parsstr,
        seed=seed,thin=thin,algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),
        sample_file=file.path(getwd(),"samples"),save_dso=FALSE) else 
        {
        stanfit=stanmodels$mlmc_code
        fitmlmc=rstan::stan(file=file.path(mlmm_path,
        "mlmm/exec/mlmc_code.stan"),
        data=prstan_data,
        iter=iterno, init=initvalue1,pars=parsstr, seed=seed, thin=thin,
        algorithm=algorithm,warmup=warmup,chains=chains,
        control=list(adapt_delta=adapt_delta_value),sample_file=NULL,
        save_dso=FALSE)}
        }

    print(fitmlmc); 
    if (savefile==TRUE) 
    utils::write.csv(as.array(fitmlmc),file=file.path(getwd(),
    "outsummary.csv"),
    row.names=TRUE)
    return(fitmlmc)
}

