\name{plot.FPE}
\alias{plot.FPE}
\title{AIC-BIC instability plot}
\description{
    Plot to visualize the instability of FPE models.
 }
\usage{
    plot.FPE(x, ...)
}
\arguments{
    \item{x}{An object of class \code{boot.modelSampler},
    as that created by the function \code{boot.modelSampler}. See below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  When \code{modelSampler} is being called in the \code{boot.modelSampler}, 
  each FPE model is kept track of in order to study the instability of these procedures. After
  finishing total bootstrap iteration, all FPE models are compared to the FPE models
  produced from full data set. This total comparison is visualized through a simple bar plot.
  Positive values less than \code{1} indicate that a given variable was selected by the FPE criteria in both the full data set
and bootstrapped data set, but was not selected at all times during bootstrap
iterations. Negative values greater than \code{-1} imply that a particular variable is
selected several times by the FPE criteria for bootstrapped data but was not
selected in the full data set.
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
  \code{plot.ooberror},
  \code{plot.icicle},
  \code{plot.var.stability},
  \code{plot.modelSampler}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  plot.FPE(ms.boot)
}
\keyword{file}
