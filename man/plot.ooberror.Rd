\name{plot.ooberror}
\alias{plot.ooberror}
\title{Plot of OOB Error Rate}
\description{
    Plot out-of-bag (OOB) error rate from the
    conditional ensemble estimation.
}
\usage{
    plot.ooberror(x, main = deparse(substitute(x)), ...)
}
\arguments{
    \item{x}{An object of class \code{boot.modelSampler},
    created by the function \code{boot.modelSampler}. See below.}
    \item{main}{Title of the plot.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  An out-of-bag estimated prediction error (PE) plot will give a
visualization of the "best" dimension of model. Red points are model size
specific estimated PE. Dotted black lines are \eqn{PE \pm s.e.(PE)}{PE +/- s.e(PE)}. 
Smoothed black line represents loess estimated line of PE.

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
  \code{print.boot.modelSampler},
  \code{print.modelSampler},
  \code{plot.modelSampler},
  \code{plot.icicle},
  \code{plot.var.stability},
  \code{plot.FPE}.
  
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  plot.ooberror(ms.boot, main="Pollute Data")
}
\keyword{file}
