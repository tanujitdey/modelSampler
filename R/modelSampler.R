#################################################################
### Last update: Aug 2018
### ##------------------------------------------------------------- 
### Spike and slab model sampler
###  . Runs under R
### 
#################################################################
### Last update: Aug 2018
###------------------------------------------------------------- 
### Spike and slab model sampler
###  . Runs under R
###--------------------------------------------------------------
### DESCRIPTION:
###   Stochastic algorithm for determining best subset of models
###   in linear regression
### ------------------------------------------------------------
###  Written by:
###
###  Tanujit Dey                        tanujit.dey@gmail.com
###  Department of Quantitative Health Sciences		
###  Cleveland Clinic
###  Cleveland, OH
###-------------------------------------------------------------
###  THIS PROGRAM SHOULD NOT BE COPIED, USED, MODIFIED, OR 
###  DISSEMINATED IN ANY WAY WITHOUT SPECIFIC WRITTEN PERMISSION 
###  FROM THE AUTHOR.
###################################################################

options(object.size=Inf,expressions=100000,memory=Inf,width=150)

modelSampler=function(
                formula,          #formula (required)
                data,             #data (required)
                V.small=0.05,     #small hypervariance atom 
                V.big=NULL,       #big hypervariance atom (NULL--> V.big=n)
                n.iter1=2500,     #no. of burn-in iterations
                n.iter2=2500,     #no. iterations sampled after burn-in
                fast=FALSE,       #break beta update into 'beta.blocks' chunks
                beta.blocks=100,  #size of beta updates (applies for fast=T)
                complexity=NULL,  #if null, estimate complexity
                verbose=TRUE,     #print iterations + other friendly facts?
                seed=NULL,        #set random generator seed
                ...)
{

###-------------------------------------------------------------
#  Some useful functions
### --------------------------------------------------------------

varianceSample=function(b,w,v0,V){
   w1=(1-w)*exp(-b^2/(2*v0))/sqrt(v0)
   w2=w*exp(-b^2/(2*V))/sqrt(V)
   if (w1==0) w1=1e-4  #numerical issue
   if (w2==0) w2=1e-4  #numerical issue
   sample(c(v0,V),1,prob=c(w1,w2))
}
mse=function(x){mean(x^2,na.rm=TRUE)}
resample <- function(x, size, ...)
       if(length(x) <= 1) { if(!missing(size) && size == 0) x[FALSE] else x
       } else sample(x, size, ...)

### --------------------------------------------------------------
###  Read data & preprocess
### --------------------------------------------------------------

if (verbose) cat("Running model sampler...","\n")	

#  Set seed if appropriate
#  Read data + set dimensions
#  Remove intercept (if it exists)
#  Watch out for singular X matrices
#  Set some tolerance values for large p and large n
#  Rescale X and Y
#  Compute mse
#  Compute XX, x.y


### preliminary checks for formula and data
if (!inherits(formula, "formula")) stop("'formula' is not a formula object.")
if (is.null(data)) stop("'data' is missing.")
if (!is.data.frame(data)) stop("'data' must be a data frame.")


if (is.null(seed)) seed=-1.0*abs(round(rnorm(1,sd=1e5)))
mf=match.call(expand.dots = FALSE)
m=match(c("formula", "data"),names(mf),0)
mf.org=mf=mf[c(1, m)]
mf$drop.unused.levels=TRUE
mf[[1]]=as.name("model.frame")
mf=eval(mf,parent.frame())
mt=attr(mf,"terms")
Y=model.response(mf,"numeric")
X.org=as.matrix(model.matrix(mt,mf))

beta.names=unlist(dimnames(X.org)[2])
if(any(beta.names=="(Intercept)")){
 if (length(beta.names) == 1) stop("Model contains only an intercept.  Add predictors.")
 int.pt=(beta.names=="(Intercept)")
 X.org=as.matrix(X.org[,!int.pt])
 beta.names=beta.names[!int.pt]
}


n.data=nrow(X.org)
n.cov=length(beta.names)
Y.org=c(scale(Y,scale=FALSE))
max.n=2500
max.p=1500
max.i=25000*100
max.xx=250
n.adj=n.data-n.cov

      
if (length(Y.org) <= 1)
      stop("Less than one observation in the data. Analysis is not meaningful.")


if ((n.data*n.cov) >= max.i)
     svs.qr.rank=0 
else{
     svs.qr=qr(X.org,tol=1e-10)
     svs.qr.rank=svs.qr$rank
}
if (svs.qr.rank!=n.cov || n.cov>=max.p){
      cat("*** X or dimension issues: computing constrained mse ***","\n")
      library(randomForest)
      rf.fit=randomForest(X.org,Y.org,ntree=1000,importance=TRUE)
      cov.order=order(rf.fit$importance[,1])[1:pmin(n.data-2,n.cov,100)]
      res.fit=lsfit(as.matrix(X.org[,cov.order]),Y.org)$res
      mse.org=sum(res.fit^2)/(n.data-pmin(n.data-2,n.cov,100))
}
else{
      svs.res=qr.resid(svs.qr,Y.org)
      mse.org=sum(svs.res^2)/n.adj
}
sf=sqrt(n.data/mse.org)
Y=Y.org*sf
X=scale(X.org)*sqrt(n.data/(n.data-1))  ##n/(n.data-1) corrects variance
X[is.na(X)]=0  ##get rid of NA's caused by scaling degenerate X

if (n.cov>max.n){# reduce computations for large dimensions
    nrep=ceiling(n.data/max.xx)
    XX=0
    for (j in 1:nrep){
       b.temp=(max.xx*(j-1)+1):(min(n.data,max.xx*j))
       XX=XX + as.matrix(t(X[b.temp,])%*%X[b.temp,])
    }
}
else{
    XX=as.matrix(t(X)%*%X)
}
sum.xy=t(X)%*%Y



### --------------------------------------------------------------
###  Initialize Gibbs parameters
###  Set up output matrices and vectors

burnin.flag=TRUE
n.sample=0
beta.all=modelTracker=matrix(0,n.cov,n.iter2)
hpm=rep(0,n.cov)
mss=aic=bic=coverage=rep(0,n.iter2)
coverage[1]=1
beta=rnorm(n.cov)
if (is.null(V.big)) V.big=n.data
hypervariance=resample(c(V.small,V.big),size=n.cov,replace=TRUE)
if (!is.null(complexity))
    complexity.flag=FALSE
else{
   complexity.flag=TRUE
   complexity=runif(1)
}
complexity.vec=rep(0,n.iter2)
### --------------------------------------------------------------
###  Gibbs Routine
### --------------------------------------------------------------
  
	
### -------------------------------------------------------------
# burn in followed by sampled values

for (i in 1:(n.iter1+n.iter2)){

        ### ---------------------------------	
        # sample beta
        # fast update for large dimensions available

        if (fast==TRUE & n.cov>beta.blocks){
             nrep=ceiling(n.cov/beta.blocks)
	     b.seq=c(1:n.cov)
	     for (j in 1:nrep){
	       b.size=min(beta.blocks,length(b.seq))
	       b.sample=resample(1:length(b.seq),size=b.size)
	       b.temp=b.seq[b.sample]
	       sum.xy.temp=sum.xy[b.temp]-XX[b.temp,-b.temp]%*%beta[-b.temp]
	       bvar.c=diag(1/(hypervariance[b.temp]),b.size) +
	                      XX[b.temp,b.temp]/n.data
	       qr.var=qr(bvar.c)
	       chol.var=t(chol(bvar.c))
	       beta[b.temp]=c(qr.coef(qr.var,
	            chol.var%*%matrix(rnorm(b.size),nrow=b.size)+sum.xy.temp/n.data))
	       if (j<nrep) b.seq=b.seq[-b.sample]
	     }
	}
	else{
  	     bvar.c=diag(1/hypervariance,n.cov)+XX/n.data
	     qr.var=qr(bvar.c)
             chol.var=t(chol(bvar.c))
	     beta=c(qr.coef(qr.var,
	              chol.var%*%matrix(rnorm(n.cov),nrow=n.cov)+sum.xy/n.data))
	}
	
        ### ---------------------------------	
        # sample the hypervariance

        hypervariance=apply(cbind(beta),1,varianceSample,w=complexity,v0=V.small,V=V.big)
        
        ### ---------------------------------	
        # sample the complexity parameter

        nonzero.pt=(hypervariance==V.big)
        if (complexity.flag==TRUE)
          complexity=rbeta(1,(1+sum(1*(nonzero.pt))),(1+sum(1*(!nonzero.pt))))
        if (1-complexity<1e-3) complexity=1-1e-3 #numerical issue
        if (complexity<1e-3)   complexity=1e-3   #numerical issue
	  		
	  		
        ### ---------------------------------			
        # verbose details
        # save output values (model tracker, mss, aic, bic)

        if (i<=n.iter1 & verbose & i%%250==0)
           cat(paste("burn-in iteration ",i,",",sep=""),
                     "complexity=",signif(complexity,3),"\n")
        if (i>n.iter1){
           if (burnin.flag==TRUE & verbose){
                cat("burn-in completed, now outputing values...","\n")
                burnin.flag=FALSE
           }
        }
        if (i>n.iter1){
         n.sample=n.sample+1
         beta.all[,n.sample]=beta/sf
         complexity.vec[n.sample]=complexity
         modelTracker[nonzero.pt,n.sample]=1
         hpm[nonzero.pt]=hpm[nonzero.pt]+1
         if (sum(nonzero.pt)>0){
            b.m=solve(diag(1/V.big,sum(nonzero.pt))
                +as.matrix(XX[nonzero.pt,nonzero.pt])/n.data)%*%sum.xy[nonzero.pt]/n.data
           mss.alpha=mse(Y-as.matrix(X[,nonzero.pt])%*%b.m)/sf^2
         }
         else
           mss.alpha=mse(Y)/sf^2
         mss[n.sample]=mss.alpha
         aic[n.sample]=mss.alpha+2*mse.org*sum(nonzero.pt)/n.data
         bic[n.sample]=mss.alpha+mse.org*log(n.data)*sum(nonzero.pt)/n.data
         if (n.sample>1){ 
           if (!any(mss[1:(n.sample-1)]==mss.alpha)) coverage[n.sample]=1
         }
	 if (verbose & n.sample%%250==0)
              cat(paste("sampled value ",n.sample,",",sep=""),
                  "complexity=",signif(complexity,3),"\n")
        }
      }
if (verbose) cat("done","\n")


### --------------------------------------------------------------
# . Gibbs finished.
# . Final Details
# . Find hpm.  Find top AIC and BIC model
# . Find top R-squared model, by model size
# . Record size of models visited and frequencies
### --------------------------------------------------------------

hpm=hpm/n.sample
aic.top=order(aic)[1]
bic.top=order(bic)[1]
aic.model=modelTracker[,aic.top]==1
bic.model=modelTracker[,bic.top]==1
colnames(modelTracker)=1:n.sample
rownames(modelTracker)=beta.names
coverage=cumsum(coverage)/(1:n.sample)
FPE=as.data.frame(cbind(hpm,
                         as.character(aic.model),
                         as.character(bic.model)))
colnames(FPE)=c("hpm","aic","bic")
rownames(FPE)=beta.names
m.size=apply(modelTracker,2,sum)
FPEstrat=FPEstrat.mss=FPEstrat.aic=FPEstrat.bic=FPEstrat.pen=modelSize=modelFreq=NULL
if (sum(m.size) != 0) {
  for (k in 1:n.cov){
    if (sum(m.size==k)>0){
      modelSize=c(modelSize,k)
      modelFreq=c(modelFreq,sum(m.size==k))
      mss.pt=resample((1:n.sample)[m.size==k & mss==min(mss[m.size==k])],1)
      FPEstrat=cbind(FPEstrat,as.character(modelTracker[,mss.pt]==1))
      FPEstrat.mss=c(FPEstrat.mss,mss[mss.pt])
      FPEstrat.aic=c(FPEstrat.aic,aic[mss.pt])
      FPEstrat.bic=c(FPEstrat.bic,bic[mss.pt])
    }
  }
  FPEstrat=as.data.frame(FPEstrat)
  FPEstrat.pen=as.data.frame(rbind(modelFreq,FPEstrat.mss,FPEstrat.aic,FPEstrat.bic))
  colnames(FPEstrat)=colnames(FPEstrat.pen)=modelSize
  rownames(FPEstrat)=rownames(beta.all)=beta.names
  rownames(FPEstrat.pen)=c("freq","mss","aic","bic")
}

### --------------------------------------------------------------
### Return the goodies
### --------------------------------------------------------------

out <- list(
            formula=mf.org,                #original formula
            modelTracker=modelTracker,     #models visited
            beta.all=beta.all,             #sampled beta values
            FPE=FPE,                       #FPE models (hpm, aic, bic)
            FPEstrat=FPEstrat,             #FPE models stratified by model size
            FPEstrat.pen=FPEstrat.pen,     #Penalty for stratified FPE models
            hpm=hpm,                       #posterior inclusion probability
            mss=mss,                       #mean rss for each model
            aic=aic,                       #aic for each model
            bic=bic,                       #bic for each model
            coverage=coverage,             #prob model visited by iteration
            complexity=complexity.vec,     #complexity values
		sf=sf					 # sigma-square hat
       )

class(out) <- "modelSampler"
return(out)

}


