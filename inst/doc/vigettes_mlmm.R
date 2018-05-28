## ----setup, include=TRUE---------------------------------------------------
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo = TRUE)

## ----installation1,eval=FALSE----------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  BiocInstaller::biocLite("mlmm")
#  library(mlmm)

## ----installation2,eval=FALSE----------------------------------------------
#  devtools::install_github("https://github.com/ireneslzeng/mlmm")
#  library(mlmm)

## ----read data,eval=TRUE---------------------------------------------------
data(pdata,package="mlmm")
pdata[1:2,]

## ----table for miss and censor,eval=TRUE-----------------------------------
table(pdata$miss,pdata$censor)

## ----re-assign miss and censor,eval=TRUE-----------------------------------
n=dim(pdata)[1]
for (i in seq_len(n))
if (pdata$miss[i]==1 && pdata$censor[i]==1) pdata$miss[i]=0

## ----set formular examples,eval=TRUE---------------------------------------
formula_completed=var1~var2+treatment;
formula_missing=miss~var2;
formula_censor=censor~1;
formula_subject=~treatment;
response_censorlim=0.002;

## ----example 1a,eval=TRUE--------------------------------------------------
library(mlmm)
model2=mlmc(formula_completed=var1~var2+treatment,
formula_missing=miss~var2, 
formula_censor=censor~1, formula_subject=~sid, 
pdata=pdata, response_censorlim=0.002, respond_dep_missing=TRUE, 
pidname="geneid", sidname="sid", 
iterno=100, chains=2, savefile=TRUE, usefit=TRUE)

## ----example 1b-1c,eval=FALSE----------------------------------------------
#  model1=mlmc(formula_completed=var1~var2+treatment,
#  formula_missing=miss~var2,
#  formula_censor=censor~1, formula_subject=~sid,
#  pdata=pdata, response_censorlim=0.002,
#  respond_dep_missing=TRUE, pidname="geneid", sidname="sid",
#  alpha_prior=c(0,0.001), iterno=100, chains=2,
#  savefile=TRUE, usefit=FALSE)
#  
#  prec_example=matrix(c(0.01,0.001,0.001,0.01),nrow=2,ncol=2)
#  
#  model3=mlmc(formula_completed=var1~var2,
#  formula_missing=miss~var2,
#  formula_censor=censor~1, formula_subject=~sid+treatment,
#  pdata=pdata, response_censorlim=0.002, respond_dep_missing=TRUE,
#  pidname="geneid", sidname="sid", prec_prior=prec_example,
#  iterno=100, chains=2, savefile=TRUE, usefit=FALSE)
#  
#  model3b=mlmc(formula_completed=var1~var2,
#  formula_missing=miss~var2,
#  formula_censor=censor~1, formula_subject=~sid+treatment,
#  pdata=pdata, response_censorlim=0.002, respond_dep_missing=TRUE,
#  pidname="geneid", sidname="sid", prec_prior=prec_example,
#  iterno=1000, chains=2, savefile=TRUE, usefit=TRUE)
#  

## ----example 2a,eval=TRUE--------------------------------------------------
model5=mlmm(formula_completed=var1~var2+treatment, 
formula_missing=miss~var2, 
formula_subject=~sid, pdata=pdata,
respond_dep_missing=FALSE, pidname="geneid", sidname="sid",
iterno=100, chains=2, savefile=TRUE, usefit=TRUE)

## ----example 2b prior,eval=TRUE--------------------------------------------
prec_example=matrix(c(0.01,0.001,0.001,0.001,0.01,0.001,0.001,0.001,0.01),
nrow=3,ncol=3)

## ----example 2b, 2c and 2d,eval=FALSE--------------------------------------
#  model4=mlmm(formula_completed=var1~var2+treatment,
#  formula_missing=miss~var2,
#  formula_subject=~sid, pdata=pdata,
#  respond_dep_missing=FALSE, pidname="geneid", sidname="sid",
#  alpha_prior=c(0,0.001), prec_prior=prec_example,
#  iterno=100, chains=2, savefile=TRUE, usefit=FALSE)
#  
#  #Example 2c. Use default priors
#  model5=mlmm(formula_completed=var1~var2+treatment,
#  formula_missing=miss~var2,
#  formula_subject=~sid+treatment, pdata=pdata,
#  respond_dep_missing=FALSE, pidname="geneid",
#  sidname="sid", iterno=100, chains=2, savefile=TRUE, usefit=FALSE)
#  
#  #Example 2d. set usefit=TRUE to call the compiled model from the previous run
#  model5b=mlmm(formula_completed=var1~var2+treatment,
#  formula_missing=miss~var2,
#  formula_subject=~sid+treatment, pdata=pdata,
#  respond_dep_missing=FALSE, pidname="geneid",
#  sidname="sid", iterno=1000, chains=2, savefile=TRUE, usefit=TRUE)

## ----plot one,echo=TRUE,eval=TRUE------------------------------------------
summaryreader=read.csv(file=file.path(getwd(),"outsummary.csv"),header=TRUE,
sep=",",skip=0)

iterno=dim(summaryreader)[1];burnin=iterno/2

U.1.1=rowMeans(matrix(c(summaryreader$chain.1.U.1.1,
summaryreader$chain.2.U.1.1),nrow=iterno,ncol=2))[burnin:iterno]

meanU=mean(U.1.1)
qU=quantile(U.1.1,p=seq(0,1,by=0.025))
scale=seq(0,1,by=0.025)

plot(scale,qU,pch=19,ylab="quantiles of estimate",xlab="quantiles")

segments(0,qU[names(qU)=="50%"],1,qU[names(qU)=="50%"],lwd=2,col="red")
segments(0,qU[names(qU)=="2.5%"],1,qU[names(qU)=="2.5%"],lty=2,lwd=2,col="red")
segments(0,qU[names(qU)=="97.5%"],1,qU[names(qU)=="97.5%"],lty=2,lwd=2,
col="red")

legend(0.5,qU[names(qU)=="50%"],"median",cex=0.8,bty="n")
legend(0.03,qU[names(qU)=="2.5%"],"2.5%",cex=0.8,bty="n")
legend(0.90,qU[names(qU)=="97.5%"],"97.5%",cex=0.8,bty="n")
qU

## ----plot two,echo=TRUE,eval=TRUE------------------------------------------
sample1reader=read.csv(file=file.path(getwd(),"samples_1.csv"),header=TRUE,
sep=",",skip=25)

sample2reader=read.csv(file=file.path(getwd(),"samples_2.csv"),header=TRUE,
sep=",",skip=25)

#plot variable U.1.1 - the intercept of first unit
trajectory_length=dim(sample1reader)[1]

plot(seq(1,trajectory_length,by=1),sample1reader$U.1.1,xlab="trajectory
number",ylab="U.1.1",type="n",
ylim=c(min(sample1reader$U.1.1,sample2reader$U.1.1,na.rm=TRUE),
max(sample1reader$U.1.1,sample2reader$U.1.1,na.rm=TRUE)))

trajectory=seq(1,trajectory_length,by=1)

lines(trajectory,sample1reader$U.1.1)
lines(trajectory,sample2reader$U.1.1,col="red")

## ----sessionInfo,echo=TRUE,eval=TRUE---------------------------------------
sessionInfo()

