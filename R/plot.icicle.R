###-------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ###------------------------------------------------------------- 
### Plot for visualizing model space
###  . Runs under R
### 
### ------------------------------------------------------------
###  Written by:
###
###  Tanujit Dey                        tanujit.dey@gmail.com
###  Department of Quantitative Health Sciences		
###  Cleveland Clinic
###  Cleveland, OH
###  ###-------------------------------------------------------------
###  THIS PROGRAM SHOULD NOT BE COPIED, USED, MODIFIED, OR 
###  DISSEMINATED IN ANY WAY WITHOUT SPECIFIC WRITTEN PERMISSION 
###  FROM THE AUTHOR.
################################################################

plot.icicle = function(x, main=deparse(substitute(x)), ...)
{
  
  if (!inherits(x, "boot.modelSampler")) stop("x is not of class boot.modelSampler")
  
  
  ### --------------------------------------------------------------
  ### Required computations for graphics
  ### --------------------------------------------------------------
  
  f.boot=x$beta.count
  modelSize=(1:dim(f.boot)[1])[apply(f.boot,1,sum)>0]
  f.boot=as.matrix(f.boot[modelSize,])
  n.cov=dim(f.boot)[2]
  
  if(n.cov==1) return(NULL)
  
  n.models=length(modelSize)  
  beta.names=colnames(f.boot)
  
  aa=apply(f.boot,1,sum)
  ab=apply(cbind(1:n.models),1,function(j){j*f.boot[j,]/aa[j]})
  
  
  if(n.cov >=5)
  {
    x.at=pretty(1:n.cov,min(50,n.cov))
    x.at=x.at[x.at>=1 & x.at<=n.cov]
    y.at=pretty(1:n.models,min(50,n.models))
    y.at=x.at[x.at>=1 & x.at<=n.cov]
    
    def.par=par(no.readonly = TRUE) # save default, for resetting...
    par(bg="cornsilk")
    
    ### contour plot
    filled.contour(1:n.cov, 1:n.models, ab, color = terrain.colors,
                   plot.title = title(paste(main,"Icicle Plot")),
                   xlab = "",
                   ylab = "Model size",
                   plot.axes = { axis(1, at=x.at, labels=beta.names[x.at], tick=TRUE, las=2,cex.axis=0.7)
                     axis(2, at=y.at, labels=modelSize[y.at],cex.axis=0.7)},
                   key.axes = axis(4, cex.axis=0.7, seq(0, 1, by = .1)))
    
  }
  
  else
  {
    x.at=1:n.cov
    y.at=1:n.models
    
    def.par=par(no.readonly = TRUE) # save default, for resetting...
    par(bg="cornsilk")
    
    
    ### contour plot
    filled.contour(1:n.cov, 1:n.models, ab, color = terrain.colors,
                   plot.title = title(paste(main,"Icicle Plot")),
                   xlab = "",
                   ylab = "Model size",
                   plot.axes = { axis(1, at=x.at, labels=beta.names, tick=TRUE, las=2,cex.axis=0.7)
                     axis(2, at=y.at, labels=modelSize,cex.axis=0.7)},
                   key.axes = axis(4, cex.axis=0.7, seq(0, 1, by = .1)))
    
  }
  
  
  par(def.par)#- reset to default
  
  
}
