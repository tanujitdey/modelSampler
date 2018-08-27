\name{modelSampler}
\alias{modelSampler}
\title{A unique tool for variable selection and model exploration in linear regression}
\description{
    Core function for the package. A call to \code{modelSampler} initiates a Gibbs 
sampler for drawing values from the posterior of a rescaled spike and slab model. 
Results from the Gibbs sampler are used to derive optimal AIC, BIC and highest 
posterior models from a restricted model search. The core function can also be called 
from its bootstrap wrapper, \code{boot.modelSampler} which can be used to assess the stability of 
AIC and BIC model selection as well as providing a more stable set of variables.
  
}
\usage{
    modelSampler(formula,         
       data,            
       V.small=0.05,     
       V.big=NULL,       
       n.iter1=2500,     
       n.iter2=2500,     
       fast=FALSE,           
       beta.blocks=100,  
       complexity=NULL,  
       verbose=TRUE,   
       seed=NULL,        
             ...)
   }
\arguments{
   
\item{formula}{A symbolic description of the model that is
        to be fit.}
\item{data}{Data frame containing the predictors (variables) in the
        model.}
\item{V.small}{Small hypervariance set to implement selective shrinkage.}
\item{V.big}{Big hypervariance. If null, V.big = n, sample size}
\item{n.iter1}{Number of burn-in iterations.}
\item{n.iter2}{Number of iterations sampled after burn-in.}
\item{fast}{Break up beta update into 'beta.blocks' chunks. Typically set to 'FALSE'.}
\item{beta.blocks}{Size of beta updates (only used when fast=TRUE).}
\item{complexity}{Model complexity parameter, which is estimated by Gibbs sampler.}
\item{verbose}{Print iterations and other user friendly outputs.}
\item{seed}{Set random generator seed.}
\item{...}{Further arguments passed to or from other methods.}
}

\details{
The specially designed Bayesian rescaled spike and slab model is 
designed to induce a type of regularization called selective shrinkage (for details see, 
reference). Selective shrinkage is due to the type of two-point prior used for the hypervariance in the 
prior as well as the choice of \code{V.big}, which by default is set to the sample size. 

}
\value{
     An object of class \code{modelSampler}, which is a list with the
     following components: 
\item{formula}{The original formula used in calling \code{modelSampler}.}
\item{modelTracker}{Total models visited after burn-in sampling.}
\item{beta.all}{Sampled beta values after burn-in sampling.}
\item{FPE}{Lists of variables selected by AIC and BIC. Also returns 
 posterior inclusion probability of each variable.}
\item{FPEstrat}{Returns top models stratified by size. 
 Selection criterion is minimum residual sum of squares (RSS).}
\item{FPEstart.pen}{Returns FPE values of the models stratified by model size. Also returns 
frequencies of models visited by \code{modelSampler}.}
\item{hpm}{Returns the posterior inclusion probability of each variable.}
\item{mss}{Returns minimum RSS values of each model visited by \code{modelSampler}.}
\item{aic}{Returns AIC values of each model visited by \code{modelSampler}.}
\item{bic}{Returns BIC values of each model visited by \code{modelSampler}.}
\item{coverage}{Returns a vector of probability of visiting a new model at each iteration
 visited by \code{modelSampler}.}
\item{complexity}{Returns a vector of estimated complexity parameters at each iteration by \code{modelSampler}.}

}


\author{
    Tanujit Dey \email{tanujit.dey@gmail.com} 

}

\references{
   
Ishwaran, H. and Rao, J. S. (2003). Detecting differentially expressed genes in
microarrays using Bayesian model selection. \emph{J. Amer. Stat. Assoc.}, \bold{98},
438 -- 455.

Ishwaran, H. and Rao, J. S. (2005).  Spike and slab gene selection for multigroup
microarray data. \emph{J. Amer. Stat. Assoc.}, \bold{100}, 764 -- 780.


Ishwaran, H. and Rao, J. S. (2005).  Spike and slab variable selection: frequentist
and Bayesian strategies. \emph{Ann. Statist.}, \bold{33}, 730 -- 773.

Dey, T. (2013).  modelSampler: An R Tool for Variable Selection and Model Exploration in Linear Regression. Journal of Data Science, \bold{11(2)}, 371--387.   
}
\seealso{
  \code{boot.modelSampler},
 \code{print.boot.modelSampler},
 \code{print.modelSampler},
 \code{plot.modelSampler},
  \code{plot.icicle},
  \code{plot.FPE},
\code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{

 # Example 1:

  data(Pollute, package = "modelSampler") 
  ms.out <- modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, verbose=TRUE)

  # Print several outputs from modelSampler. 
  
  print(ms.out)
  
  # Returns a collection of graphics which includes a complexity plot; 
  # a penalization plot which depicts model size specific estimated 
  # minimum residual sum of squares, AIC, BIC values; a dimensionality plot 
  # of several model sizes visited by modelSampler; 
  # an image plot to visualize variable importance
  # and a coverage plot depicting the probability of visiting new model by Gibbs sampler. 
  # For details of each plot, see plot.modelSampler.

  plot.modelSampler(ms.out)

  # Based on preliminary analysis, an out-of-bag technique is used
  # estimate prediction error (PE). Based on estimated PE, 
  # the best model of size "k" is being selected.
  
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  
  # Prints selected subset of variables, based on estimated prediction error.
  
  print(ms.boot)
  
  # This plot will give an idea about instability of FPE model selection criteria.
    
  plot.FPE(ms.boot)
  
  # This plot will depict the model space.

  plot.icicle(ms.boot, main="The Pollute data")
 
  # Graphical visualization for selecting "the" best model based on estimated 
  # prediction error of hard shrunk predictors.

  plot.ooberror(ms.boot, main="The Pollute data")
  }
\keyword{file}
