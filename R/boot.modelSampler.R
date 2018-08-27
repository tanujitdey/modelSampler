#################################################################
### Last update: Aug 2018
### ###------------------------------------------------------------- 
### Bootstrap wrapper
### Calls modelSampler
###
### Bagged ensemble yields best model of size k, for any k.
### Hard shrunk ensemble predictor used to determine optimal k.
### Optimal model determined using prediction error.
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
###################################################################

boot.modelSampler= function(
                      formula,          #formula (required)
                      data,             #data (required)
                      n.iter1 = 2500,   #no. of burn-in iterations
                      n.iter2 = 2500,   #no. iterations sampled after burn-in
                      B = 20,           #no. bootstraps
                      verbose = TRUE,   #verbose output?
                      ...)
{

### --------------------------------------------------------------
#  Some useful functions
### --------------------------------------------------------------

mse=function(x){mean(x^2,na.rm=TRUE)}


### --------------------------------------------------------------
###  Read data & preprocess
### --------------------------------------------------------------

if (verbose) cat("Running out-of-bag model sampler...","\n")	

### preliminary checks for formula and data
if (!inherits(formula, "formula")) stop("'formula' is not a formula object.")
if (is.null(data)) stop("'data' is missing.")
if (!is.data.frame(data)) stop("'data' must be a data frame.")



#  Read data + set dimensions
#  Remove intercept (if it exists)
#  Watch out for singular X matrices

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

      
if (length(Y) <= 1)
      stop("Less than one observation in the data.  Analysis is not meaningful.")

### --------------------------------------------------------------
### Set up output matrices and vectors
### Set random seeds
### --------------------------------------------------------------

beta.ensemble=rep(0,n.cov)
beta.count=beta.stability=matrix(0,n.cov,n.cov)
pred.ensemble=matrix(NA,n.data,B)
pred.hard=oob.hard=matrix(0,n.data,n.cov)
oob.pe.ensemble=save.seed=save.sf=rep(NA,B)
track.aic=track.bic=matrix(0,B,n.cov)


##--------------------------
## FPE models from full data
##--------------------------

full.ms = modelSampler(formula = formula,
                     data = data,
                     n.iter1 = n.iter1,
                     n.iter2 = n.iter2,
                     verbose = FALSE)

aicbic.full = full.ms$FPE[,-1]
 

### --------------------------------------------------------------
###  BOOTSTRAP LOOP
### --------------------------------------------------------------

for (b in 1:B){

if (verbose) cat("Bootstrap draw: ",b,"\n")

### Draw bootstrap indices 
### Create boot and oob data
### Save the bootstrap seed for later re-creation
### Call model sampler

save.seed[b]=round(runif(1,-1e6,1e6))
set.seed(save.seed[b])
boot.sample=sample(1:n.data,n.data,replace=TRUE)
oob.sample=setdiff(1:n.data,boot.sample)
boot.ms=modelSampler(formula = formula,
                     data = data[boot.sample,],
                     n.iter1 = n.iter1,
                     n.iter2 = n.iter2,
                     verbose = FALSE)
save.sf[b]=boot.ms$sf


###
### calculate FPE matrix
###

track.aic[b,] = as.numeric(boot.ms$FPE[,2]==TRUE)
track.bic[b,] = as.numeric(boot.ms$FPE[,3]==TRUE)


### Calculate full ensemble
### Standardize bagged and OOB X-predictors
### Note OOB standardization (X and Y) must be based on bagged data only

X.mean=apply(as.matrix(X.org[boot.sample,]),2,mean)
X.sd=apply(as.matrix(X.org[boot.sample,]),2,sd)*sqrt((n.data-1)/n.data)
Y.mean=mean(Y[boot.sample])
X.boot=scale(as.matrix(X.org[boot.sample,]),center=X.mean,scale=X.sd)
X.boot[,X.sd==0]=0
X.oob=scale(as.matrix(X.org[oob.sample,]),center=X.mean,scale=X.sd)
X.oob[,X.sd==0]=0
beta.bma=apply(boot.ms$beta.all,1,mean,na.rm=TRUE)
beta.ensemble=beta.ensemble+beta.bma
pred.ensemble[oob.sample,b]=Y.mean+X.oob%*%beta.bma
oob.pe.ensemble[b]=mse(Y-apply(pred.ensemble,1,mean,na.rm=TRUE))
if (verbose) {
  cat("Ensemble PE:                ",
       round(mean(oob.pe.ensemble,na.rm=TRUE),3),"\n")
}

### Condition on model size:
### -- Compute conditional beta
### -- Use conditional beta to find best model of size k
### -- Track number of times a variable is in best model of size k
### -- Compute PE for hard shrunk ensemble predictor
### -- Ridge regression used.  !! Note ridge parameter=1 !!

m.size=apply(boot.ms$modelTracker,2,sum)
if (sum(m.size) != 0) {
  for (k in 1:n.cov) {
    if (sum(m.size==k)>0) {
      if (n.cov == 1) {
        beta.cond=mean(boot.ms$beta.all[,m.size==k],na.rm=TRUE)
      }
      else { 
        beta.cond=apply(as.matrix(boot.ms$beta.all[,m.size==k]),1,mean,na.rm=TRUE)
      }
      hard.pt=rev(order(abs(beta.cond)))[1:k]
      beta.count[k,hard.pt]=beta.count[k,hard.pt]+1
      X.hard=as.matrix(X.boot[,hard.pt])
      beta.hard=solve(t(X.hard)%*%X.hard+diag(1,k))%*%t(X.hard)%*%(Y[boot.sample]-Y.mean)
      pred.hard[oob.sample,k]=pred.hard[oob.sample,k]+
           Y.mean+as.matrix(X.oob[,hard.pt])%*%beta.hard
      oob.hard[oob.sample,k]=oob.hard[oob.sample,k]+1
      if (verbose) {
        cat("(Model size ",k,", PE):     (",
            sum(m.size==k),", ",
            round(mse(Y-(pred.hard/oob.hard)[,k]),3),")",
            sep="","\n")
      }
    }
  }
 }


}

### Repeat bootstrap loop
### Use same seeds to draw bootstrap data
### Compute beta stability values
beta.ensemble=beta.ensemble/B
modelSize=(1:n.cov)[apply(beta.count,1,sum)>0]
o.pt=rev(order(abs(beta.ensemble)))
if (verbose) cat("Re-seeding and computing stability values...\n")
for (b in 1:B) {
  if (verbose) cat("Iteration:",b,"\n")
  set.seed(save.seed[b])
  boot.sample=sample(1:n.data,n.data,replace=TRUE)
  X.mean=apply(as.matrix(X.org[boot.sample,]),2,mean)
  X.sd=apply(as.matrix(X.org[boot.sample,]),2,sd)*sqrt((n.data-1)/n.data)
  Y.mean=mean(Y[boot.sample])
  X.boot=scale(as.matrix(X.org[boot.sample,]),center=X.mean,scale=X.sd)
  X.boot[,X.sd==0]=0
  for (k in modelSize) {
    hard.pt=o.pt[1:k]
    X.hard=as.matrix(X.boot[,hard.pt])
    beta.hard=solve(t(X.hard)%*%X.hard+diag(1,k))%*%t(X.hard)%*%(Y[boot.sample]-Y.mean)
    beta.stability[k,hard.pt]=beta.stability[k,hard.pt]+beta.hard*save.sf[b]
  }
}
beta.stability=beta.stability/B


### --------------------------------------------------------------
###  Final details
### --------------------------------------------------------------

### Set row and column names
### Final details for PE
### Use Jackknife estimate for standard error
rownames(beta.count)=rownames(beta.stability)=1:n.cov
colnames(beta.count)=colnames(beta.stability)=beta.names
rownames(track.aic)=rownames(track.bic)=1:B
colnames(track.aic)=colnames(track.bic)=beta.names
oob.pe.hard=apply(Y-pred.hard/oob.hard,2,mse)
oob.se=sqrt(apply(((Y-pred.hard/oob.hard)^2-oob.pe.hard)^2,2,sum,na.rm=T))/n.data


### --------------------------------------------------------------
### Return the goodies
### --------------------------------------------------------------

out <- list(
            beta.count=beta.count,            #used for describing model space
            beta.ensemble=beta.ensemble,      #bagged ensemble estimator
            oob.pe.hard=oob.pe.hard,          #OOB PE for hard shrunk ensemble
            oob.pe.ensemble=oob.pe.ensemble,  #OOB PE for full ensemble
            track.aic=track.aic,	      #aic instability matrix 
            track.bic=track.bic,	      #bic instability matrix 
	    aicbic.full=aicbic.full,	      #aicbic models from full data
	    oob.se=oob.se,		      #standard error for OOB PE for hard shrunk ensemble
            beta.stability=beta.stability     #beta coefficient for stability plot
           ) 

class(out) <- "boot.modelSampler"
return(out)

}

