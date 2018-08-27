####------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### Prints output from modelSampler()
###  . Runs under R
### 
### ------------------------------------------------------------
###  Written by:
###
###  Tanujit Dey                        tanujit.dey@gmail.com
###  Department of Quantitative Health Sciences		
###  Cleveland Clinic
###  Cleveland, OH
###  -------------------------------------------------------------
###  THIS PROGRAM SHOULD NOT BE COPIED, USED, MODIFIED, OR 
###  DISSEMINATED IN ANY WAY WITHOUT SPECIFIC WRITTEN PERMISSION 
###  FROM THE AUTHOR.
###################################################################


print.modelSampler = function (x, ...) 
{

if (!inherits(x, "modelSampler"))
  stop("x must be of class 'modelSampler'")

o.r=rev(order(abs(apply(x$beta.all,1,mean,na.rm=TRUE))))
n.cov=dim(x$modelTracker)[1]
n.sample=dim(x$modelTracker)[2]

cat("-------------------------------------------------------------------","\n")
cat("No. predictors               :",n.cov,"\n")
cat("No. sampled values           :",n.sample,"\n")
cat("Estimated complexity         :",round(mean(x$complexity),3),"+/-",
                                     round(sd(x$complexity),3),"\n")
cat("Prob. visiting new model     :",round(x$coverage[n.sample],3),"\n")
cat("\n")
cat("Model selection results:","\n")
print(x$FPE[o.r,],justify="right",print.gap=3)
cat("\n")
cat("Top models stratified by size:","\n")
print(x$FPEstrat[o.r,],justify="right",print.gap=3)
print(x$FPEstrat.pen,justify="right",print.gap=3)
cat("-------------------------------------------------------------------","\n")

}

