###-------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### Prints output from boot.modelSampler()
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


print.boot.modelSampler = function (x, ...) 
{
  if (!inherits(x, "boot.modelSampler"))
  stop("x must be of class 'boot.modelSampler'")
  n.cov=length(x$oob.pe.hard)
  beta.names=colnames(x$beta.count)
  o.r=rev(order(abs(x$beta.ensemble)))

  khat=min((1:n.cov)[x$oob.pe.hard==min(x$oob.pe.hard,na.rm=TRUE)],na.rm=TRUE) 
  if (khat!=Inf)
 {
    	
	cat("-------------------------------------------------------------------","\n")
	cat("Optimal model obtained via ensemble out-of-bagging:")
	cat("\n")
	print(beta.names[o.r[1:khat]])
	cat("-------------------------------------------------------------------","\n")
}
  else
 {
	cat("-------------------------------------------------------------------","\n")
	cat("There is no model to select via ensemble out-of-bagging")
	cat("\n")
	cat("-------------------------------------------------------------------","\n")
}

}

