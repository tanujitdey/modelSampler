\name{plot.var.stability}
\alias{plot.var.stability}
\title{Variable stability plot}
\description{
    Variable stability plot over the model space.
}
\usage{
    plot.var.stability(
		x, 
		filter.flag = FALSE,
		high.dim = 25,
		...)
}
\arguments{
    \item{x}{An object of class \code{boot.modelSampler},
    created by the function \code{boot.modelSampler}. See below.}
    \item{filter.flag}{It is an optional argument. Default value is \code{FALSE}. If the value is \code{TRUE}, 
	then it plots only the best subset of variables over the entire model space. It's useful
	in high dimensional data. When number of variables are very high, it's hard to visualize the path
      important variables over the model space. This argument will help user to see the path of the best
	subset of variables.}
   \item{high.dim}{An optional argument. \code{25} is default setting. Usful in high dimensional setting. 
   This argument implements median smoothing to the estimators to enchance better graphical depiction.}
 \item{...}{Further arguments passed to or from other methods.}
}
\details{
 Coming soon...
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
  \code{plot.ooberror},
  \code{plot.FPE}.
  
}
\examples{

#---------------------------------------------------------------------------
# Example 1:
# In the following example Ozone data set is used to draw the plot. Total number
# of variables in the data set is 12. Therefore optional arguments are not so
# useful to use.
#----------------------------------------------------------------------------

data(Ozone, package = "modelSampler") 
ms.boot <- boot.modelSampler(ozone~., 
			 Ozone, 
			 n.iter1=2500, 
  			 n.iter2=2500, 
			 B=20, 
			 verbose = TRUE)

plot.var.stability(ms.boot, 
			filter.flag=FALSE, 
			high.dim = 25)

#--------------------------------------------------------------------------------------
# Example 2:
# In this example Ozone interaction data set is used. Total number of variables in this
# data set is 89. So, it's now useful to use optional arguments. This example will plot 
# only the best subset of variables over the modified model space.
#---------------------------------------------------------------------------------------

data(OzoneI, package = "modelSampler") 
ms.boot <- boot.modelSampler(ozone~., 
			 OzoneI, 
			 n.iter1=2500, 
  			 n.iter2=2500, 
			 B=20, 
			 verbose = TRUE)

plot.var.stability(ms.boot, 
			filter.flag=TRUE, 
			high.dim = 25)

}
\keyword{file}
