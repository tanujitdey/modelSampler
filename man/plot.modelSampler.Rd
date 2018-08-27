\name{plot.modelSampler}
\alias{plot.modelSampler}
\title{Graphical analysis of modelSampler}
\description{
    Graphical analysis of the \code{object} created by \code{modelSampler}.
}
\usage{
    plot.modelSampler(x,...)
}
\arguments{
    \item{x}{An object of class \code{modelSampler},
    created by the function \code{modelSampler}. See below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Five plots are produced.  Going from top to bottom, left to right: 

(1)  \code{modelSampler} estimates a complexity parameter. A Complexity Plot 
indicating the range of estimated complexity parameters,
thus allowing the user to interpret the dimensionality of model
space. A higher complexity value characterizes a larger model.

(2)  A Penalization Plot depicting FPE values as a function of dimension.

(3) A Dimensionality Plot depicting model size visited by \code{modelSampler}.

(4)  An Image Plot for visualization of models sampled as a function of the
number of Monte Carlo iterations. The importance of a variable can be
indirectly assessed by this plot.

(5) A Coverage Plot depicting the probability of visiting a new model at
each iteration.
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
  \code{plot.icicle},
  \code{plot.FPE},
  \code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.out <- modelSampler(MortRate~., Pollute, n.iter1=2500, 
  n.iter2=2500, verbose=TRUE)
  plot.modelSampler(ms.out)
}
\keyword{file}
