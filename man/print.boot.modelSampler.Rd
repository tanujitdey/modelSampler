\name{print.boot.modelSampler}
\alias{print.boot.modelSampler}
\title{Print summary output from boot.modelSampler}
\description{
 At the end of bootstrap run, it returns an \code{object} which prints the "best" subsets of variables. 
}
\usage{
    print.boot.modelSampler(x,...)
}
\arguments{
    \item{x}{An object of class \code{boot.modelSampler},
    created by the function \code{boot.modelSampler}. See below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Outputs from \code{boot.modelSampler}. 
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
 \code{print.modelSampler},
  \code{plot.icicle},
  \code{plot.FPE},
  \code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  print.boot.modelSampler(ms.boot)
}
\keyword{file}
