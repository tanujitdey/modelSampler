###--------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### Plots from the output of FPE modelSampler()
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


plot.modelSampler= function(x,...){
  
  if (!inherits(x, "modelSampler")) stop("x is not of class modelSampler")
  if(is.null(x$FPEstrat.pen) == "TRUE") stop("This plot wrapper is not useful for this analysis")

  o.r=rev(order(abs(apply(x$beta.all,1,mean,na.rm=T))))
  n.cov=dim(x$modelTracker)[1]
  n.sample=dim(x$modelTracker)[2]
  beta.names=rownames(x$modelTracker)
  imagePlot=vector("list",n.cov)
  modelFreq=as.double(x$FPEstrat.pen[1,])
  modelSize=as.double(colnames(x$FPEstrat.pen))
  for (k in 1:n.cov) {
    imagePlot[[k]]=(1:n.sample)[x$modelTracker[k,]==1]
  }

  
def.par=par(no.readonly = TRUE) # save default, for resetting...
par(bg="cornsilk")

layout(cbind(c(2,1),c(4,1),c(5,3)), c(4,4,2), c(2,6), TRUE)
par(mar=c(3,5,1,1))

plot(c(1,n.cov),
     c(1,n.sample),
     type="n",
     xaxt="n",
     xlab="",
     ylab="Iteration")
axis(1,at=(1:n.cov),labels=beta.names[o.r],tick=TRUE,las=2)  

for (k in 1:n.cov) {
  points(jitter(rep(k,length(imagePlot[[o.r[k]]]))),
         imagePlot[[o.r[k]]],
         pch=21,cex=0.30,
         bg="green3")
}
par(mar=c(2,3,2,1))
c.hist=hist(x$complexity,
            breaks=seq(min(x$complexity),max(x$complexity),length=10),
            plot=FALSE)
barplot(c.hist$counts,
        axes=FALSE,
        names.arg=round(c.hist$breaks[-1],2),
        main="Complexity",
        space=0)
par(mar=c(4,0,1,1))
plot(x$coverage,
     1:n.sample,
     type="l",
     yaxt="n",
     ylab="",
     xlab="Coverage")


par(mar=c(2,1,2,1))
matplot(modelSize,t(x$FPEstrat.pen[-1,]),
        type="l",
        xaxt="n",
        lty=1,
        col=c(1:3),
        main="Penalization")
axis(1,at=modelSize,labels=modelSize,tick=TRUE)


par(mar=c(2,3,2,1))
barplot(modelFreq,
        axes=FALSE,
        main="Model Freq.",
        names.arg=round(modelSize,2),
        space=0,
        horiz=FALSE)

par(def.par)#- reset to default

}

