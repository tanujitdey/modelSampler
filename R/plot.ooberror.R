###--------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### Plot for OOB error rate
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



plot.ooberror = function(x, main=deparse(substitute(x)), ...)
{

  if (!inherits(x, "boot.modelSampler")) stop("x is not of class boot.modelSampler")
  n.cov = length(x$oob.pe.hard)
  modelSize = (1:n.cov)[is.na(x$oob.pe.hard)==FALSE]

  if(n.cov >=2)
{
  def.par=par(no.readonly = TRUE) # save default, for resetting...
  par(bg="cornsilk")
  
  
  x.err = modelSize
  y.err = na.omit(x$oob.pe.hard)
  err.se = x$oob.se[x$oob.se >0]

  err.l = y.err - err.se
  err.u = y.err + err.se
  
  ylim.l = min(err.l, na.rm =TRUE)
  ylim.u = max(err.u, na.rm =TRUE)
 

  if (n.cov < 50) pch.symbol=16 else pch.symbol=20
     
  plot(x.err,y.err,ylim = c(ylim.l, ylim.u), type = "n", xlab="Dimension",       
       xaxt = "n", ylab="Out-of-Bag Prediction Error")
  axis(1,x.err,las=2,cex.axis=0.7)	

  title(main = main)

  
  	if(n.cov >=3)
	{     		

		mat <- cbind(x.err,err.l,err.u)
            mat <- mat[sort.list(mat[,1]),]
		
            x.grid <- mat[,1]
            lower <- mat[,2]
            upper <- mat[,3]
    			
            polygon(c(x.grid,rev(x.grid)),c(lower,rev(upper)),
            col="gray", border=FALSE)
	
		points(x.err,y.err, col = "red", lwd=2, cex.axis=0.7, pch = pch.symbol)
  		y.loess <- loess(y.err ~ x.err, data.frame(x.err, y.err))
  		y.predict <- predict(y.loess, data.frame(x.err))
  		lines(x.err,y.predict)
	}

}

else
{
return(NULL)
}
par(def.par)#- reset to default
}


