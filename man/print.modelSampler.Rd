\name{print.modelSampler}
\alias{print.modelSampler}
\title{Print summary output from modelSampler}
\description{
    \code{modelSampler} returns an \code{object} to get several outputs. 
}
\usage{
    print.modelSampler(x,...)
}
\arguments{
    \item{x}{An object of class \code{modelSampler},
    created by the function \code{modelSampler}. See below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
 Prints model size specific FPE values, namely; minimum residual sum of squares, minimum AIC and minimum BIC. 
Prints variables selected by FPE criteria AIC and BIC. Also prints posterior inclusion probabilities of 
each variable and top models (based on \eqn{R^2}{R^2} value) stratified by size. 
}
\author{
    Tanujit Dey \email{tanujit.dey@gmail.com} 

}

\references{
    Dey, T. (2013).  modelSampler: An R Tool for Variable Selection and Model Exploration in Linear Regression. Journal of Data Science, \bold{11(2)}, 371-387.
}
\seealso{
  \code{boot.modelSampler},
  \code{modelSampler},
  \code{plot.modelSampler},
 \code{print.boot.modelSampler},
  \code{plot.icicle},
  \code{plot.FPE},
  \code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.out <- modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, verbose=TRUE)
  print.modelSampler(ms.out)
}
\keyword{file}
