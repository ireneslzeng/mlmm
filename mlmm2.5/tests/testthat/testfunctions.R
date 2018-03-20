data(pdata)
test_that("Test that data and formula match1",
{
 #copy and paste the following formulas to the mmlm() function respectively
 expect_match(tryCatch({model_test=mlmc(formula_completed=var1~var2+treatment,formula_missing=miss~var2,
 formula_censor=censor~1,formula_subject=,pdata=pdata,response_censorlim=0.002,
 respond_dep_missing=FALSE,pidname="geneid",sidname="sid",
 iterno=50,chains=2,savefile=FALSE,usefit=FALSE)},error=function(e) {print("errors in data")}),"errors in data")
 
 expect_match(tryCatch({model_test=mlmm(formula_completed=var1~var2+treatment,formula_missing=miss~var2,
 formula_subject=,pdata=pdata,respond_dep_missing=TRUE,
 pidname="geneid",sidname="sid",pathname=pathdir,iterno=10,chains=2,usefit=FALSE,savefile=FALSE)
 },error=function(e) {print("errors in formula")}),"errors in formula")

})

test_that("Test that data and formula match2",
{
 #copy and paste the following formulas to the mmlm() function respectively
 expect_match(tryCatch({model_test=mlmc(formula_completed=var1~var2+treatment,formula_missing=,
 formula_censor=censor~1,formula_subject=~treatment,pdata=pdata,response_censorlim=0.002,
 respond_dep_missing=FALSE,pidname="geneid",sidname="sid",
 iterno=50,chains=2,savefile=FALSE,usefit=FALSE)},error=function(e) {print("errors in data")}),"errors in data")
 
 expect_match(tryCatch({model_test=mlmm(formula_completed=var1~var2+treatment,formula_missing=,
 formula_subject=~sid,pdata=pdata,respond_dep_missing=TRUE,
 pidname="geneid",sidname="sid",pathname=pathdir,iterno=10,chains=2,usefit=FALSE,savefile=FALSE)
 },error=function(e) {print("errors in formula")}),"errors in formula")

})


test_that("Test that to use set initial value",
{
expect_equal(setinitvalues(npred=2,np=3,npred_miss=3,npred_sub=2,nmiss=10,nsid=30)$`alpha_response`,0.008)
})
