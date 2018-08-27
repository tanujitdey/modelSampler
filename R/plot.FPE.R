
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### Plot for visualizing instability of FPE model selction
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

plot.FPE = function(x, ...)
{

  if (!inherits(x, "boot.modelSampler")) stop("x is not of class boot.modelSampler")
  def.par=par(no.readonly = TRUE) # save default, for resetting...
  par(mfrow = c(1, 2), bg="cornsilk")
  
  
  aic.full = as.numeric(x$aicbic.full[,1]==TRUE)
  bic.full = as.numeric(x$aicbic.full[,2]==TRUE)

 n.cov =length(aic.full)
 B = dim(x$track.aic)[1]

aic.adj = aic.full - apply(x$track.aic,2,sum)/B
bic.adj = bic.full - apply(x$track.bic,2,sum)/B

aicbic.dat =cbind(aic.adj,bic.adj)
rownames(aicbic.dat)=c(1:n.cov)

barplot(aicbic.dat[,1],ylim=c(-1,1),border = "dark blue",
	  plot.title = title(paste("AIC models")),
    	  xlab = "variable index",
    	  ylab = "Misclassification Rate",las=2)


barplot(aicbic.dat[,2],ylim=c(-1,1),border = "dark blue",
	  plot.title = title(paste("BIC models")),
    	  xlab = "variable index",
    	  ylab = "",las=2)
	 

par(def.par)#- reset to default

}
