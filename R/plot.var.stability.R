###-------------------------------------------------------------
#################################################################
### Last update: Aug 2018
### ------------------------------------------------------------- 
### variable stability plot over model space
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


plot.var.stability=function(x, #boot object 
filter.flag=FALSE, #use best subset of variables?
high.dim=25, #special treatment if n.cov>high.dim
...)

{
if (!inherits(x, "boot.modelSampler")) stop("x is not of class boot.modelSampler")

def.par=par(no.readonly = TRUE) # save default, for resetting...
par(bg="gray")


### Get beta coefficients
### Check that n.cov>1
### Sort beta using the BMA

beta.v=as.matrix(x$beta.stability)
n.cov=dim(beta.v)[2]
if(n.cov==1) return(NULL)
modelSize=(1:n.cov)[apply(beta.v,1,sum)!=0]
n.models=length(modelSize)
beta.v=beta.v[,rev(order(abs(x$beta.ensemble)))]

### Get prediction error
### Check NA status
### Check filter flag status
### Generate beta label

pe=x$oob.pe.hard
if (all(is.na(pe))) return(NULL)
min.pe=min(na.omit((1:n.cov)[pe == min(pe,na.rm=TRUE)]))
if (min.pe==1) filter.flag=FALSE
beta.label=1:n.cov
xt=beta.label[beta.v[min.pe,]==0]
beta.label[xt] = ''

### Running median smoother for high dimensional case
if (n.cov > high.dim) {
khat <- 1 + 2 * min((n.models-1)%/% 2, ceiling(0.25*n.models))
beta.v[modelSize,]=apply(rbind(beta.v[modelSize,]),2,function(x){runmed(x,k=khat)})
}


if (filter.flag==FALSE)

{
beta.v = rbind(beta.v[modelSize,])
o.r=rev(order(pe,na.last=NA))

ab.line.min = o.r[n.models]
ab.line.rest=o.r[-n.models]

if (n.cov <= high.dim) pch.sym = 16 else pch.sym = "*"
if (n.cov <= high.dim) 
{ 
plot.type = "b"
plot.lwd = 2
}
else 
{
plot.type = "l"
plot.lwd = 1
}


i=1:n.models
k=1
j=1.2
dy=0.2 

split.screen(rbind(c(0, .8,0,1), c(.8,.95,0,1))) 
screen(1)
cus.col=rgb((1:n.models)/n.models, g=0,b=0)

matplot(modelSize, beta.v, 
xlab = "Model Size", type = plot.type, xaxt = "n", col = rainbow(n.models),
pch = pch.sym ,lwd=plot.lwd, lty=1, ylab = "Coefficients", cex.axis =0.7,
main = "Variable Stability Plot",
ylim=c(min(beta.v),max(beta.v))) 
axis(1,modelSize,las=2,cex.axis=0.7) 
axis(4, at = beta.v[nrow(beta.v), ], label = beta.label, las=2, cex.axis = 0.7)

if (n.cov<=high.dim)
{

abline(v = ab.line.rest,lty=1,lwd=6,col ="light gray")
abline(v = ab.line.rest,lty=1,lwd=2,col =cus.col[-n.models])
abline(v = ab.line.min,lty=1,lwd=2,col =cus.col[n.models])
} 

else
{
abline(v = ab.line.rest,lty=3,lwd=1,col =cus.col[-n.models])
abline(v = ab.line.min,lty=1,lwd=2,col =cus.col[n.models])
}




screen(2)

plot.window(xlim = c(dy, j), ylim = c(min(i-.5), max(i+.4)), xaxs = "i", 
yaxs = "i")
rect((k-1)*j+ dy, i-0.5, k*j, i+0.4, col = cus.col)
mtext("Color key for vertical lines",side=1, cex=.7) 
mtext("Prediction error value: from largest to smallest",side=4,cex=.9)

close.screen(all=TRUE) 
}


else
{
beta.v = rbind(beta.v[modelSize,-xt])
beta.label=beta.label[-xt]

o.r=rev(order(pe,na.last=NA))
n.models = length(o.r)

ab.line.min = o.r[n.models]
ab.line.rest=o.r[-n.models]

if (n.cov <= high.dim) pch.sym = 16 else pch.sym = "*"
if (n.cov <= high.dim) plot.type = "b" else plot.type = "l"

i=1:n.models
k=1
j=1.2
dy=0.2 

split.screen( rbind(c(0, .8,0,1), c(.8,.95,0,1)) ) 
screen(1)
cus.col=rgb((1:n.models)/n.models, g=0,b=0)

matplot(modelSize, beta.v, 
xlab = "Model Size", type = plot.type, col = rainbow(min.pe), xaxt = "n",
pch = pch.sym ,lwd=1, lty=1,ylab = "Coefficients", cex.axis =0.7,
main = "Variable Stability Plot",
ylim=c(min(beta.v),max(beta.v)))
axis(1,modelSize,las=2,cex.axis=0.7) 
axis(4, at = beta.v[nrow(beta.v), ], label = beta.label, las=2, cex.axis = 0.7)


if (n.cov<=high.dim)
{

abline(v = ab.line.rest,lty=1,lwd=6,col ="light gray")
abline(v = ab.line.rest,lty=1,lwd=2,col =cus.col[-n.models])
abline(v = ab.line.min,lty=1,lwd=2,col =cus.col[n.models])
} 

else
{
abline(v = ab.line.rest,lty=3,lwd=1,col =cus.col[-n.models])
abline(v = ab.line.min,lty=1,lwd=2,col =cus.col[n.models])
}

screen(2)

plot.window(xlim = c(dy, j), ylim = c(min(i-.5), max(i+.4)), xaxs = "i", 
yaxs = "i")
rect((k-1)*j+ dy, i-0.5, k*j, i+0.4, col = cus.col)
mtext("Color key for vertical lines",side=1, cex=.7) 
mtext("Prediction error value: from largest to smallest",side=4,cex=.9)

close.screen( all=TRUE) 

}

par(def.par)#- reset to default

}


