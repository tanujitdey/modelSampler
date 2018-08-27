\name{plot.icicle}
\alias{plot.icicle}
\title{Icicle plot}
\description{
    Graphical analysis of the \code{object} created by \code{boot.modelSampler}.
}
\usage{
    plot.icicle(x,main = deparse(substitute(x)), ...)
}
\arguments{
    \item{x}{An object of class \code{boot.modelSampler},
    as that created by the function \code{boot.modelSampler}. See below.}
	\item{main}{Title of the plot.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  An Icicle plot will help user to assess the stability of variables and provide a
stable set of selected variables.
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
  \code{print.modelSampler},
 \code{print.boot.modelSampler},
  \code{plot.modelSampler},
  \code{plot.FPE},
  \code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  plot.icicle(ms.boot, main="Pollute Data")

}
\keyword{file}
