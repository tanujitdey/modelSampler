\name{boot.modelSampler}
\alias{boot.modelSampler}
\title{An ensemble technique to select the best subset of variables}
\description{
    Using Out-of-bag technique, the best subset of variables is found. 
}
\usage{
    boot.modelSampler(formula, 
		data, 
		n.iter1=2500,
		n.iter2=2500,
		B=20,
		verbose=TRUE,
			...)
}
\arguments{
   	\item{formula}{A symbolic description of the model that is
        to be fit. This argument is required for \code{modelSampler}.}
	\item{data}{Data frame containing the predictors (variables) in the
       model. This argument is required for \code{modelSampler}.}
	\item{n.iter1}{Number of burn-in samples. This argument is required for \code{modelSampler}.}
  \item{n.iter2}{Number of samples after burn-in. This argument is required for \code{modelSampler}.}
   \item{B}{Total number of bootstrap iterations.}
\item{verbose}{Prints pretty summary from each bootstrap iteration.}	
    \item{...}{Further arguments passed to or from other methods.}
}

\details{
  \code{boot.modelSampler} is a bootstrap wrapper that calls
the primary function \code{modelSampler}. The user specifies \code{B}, the number of
bootstrap draws to use, and the wrapper then makes \code{B} calls to the primary
function. Each call uses a bootstrap draw of the original data.
For each bootstrap draw, a hard shrunk posterior mean
is computed for each model size visited by the  \code{modelSampler}.  Hard shrunk
estimators are then combined over the \code{B} draws to form an ensemble for each
given model size.  Out-of-bagging is then used to estimate the
prediction error for each of these ensemble hard shrunk predictors and the
predictor with the smallest prediction error is determined. The dimension
of this predictor is defined to be the optimal model size. The optimal
model is chosen by then selecting the first ordered \code{k} variables. Ordering is
based on an ensemble BMA predictor formed by averaging the Bayesian model
averaged (BMA) estimator over the \code{B} bootstrap draws.

}

\value{
     An object of class \code{boot.modelSampler}, which is a list with the
     following components:
\item{beta.count}{Returns a matrix whose each row corresponds to model size and each column
 corresponds to variable name. The entries corresponds to each cell is the number of times each variable
 being identified as the most significant variable according to the model size.}
\item{beta.ensemble}{Bagged ensemble estimators.}
\item{oob.pe.hard}{Out-of-bag estimated prediction error for hard shrunk predictors for each model.}
\item{oob.pe.ensemble}{Out-of-bag estimated prediction error for the ensemble predictors.}
\item{track.aic}{Returns a matrix of AIC models during each bootstrap iteration.}
\item{track.bic}{Returns a matrix of BIC models during each bootstrap iteration.}
\item{aicbic.full}{AIC-BIC models from full data.}
\item{oob.se}{Standard error of out-of-bag estimated prediction error for hard shrunk
 predictors for each model.}

}


\author{
    Tanujit Dey \email{tanujit.dey@gmail.com} 

}

\references{
    Dey, T. (2013).  modelSampler: An R Tool for Variable Selection and Model Exploration in Linear Regression. Journal of Data Science, \bold{11(2)}, 371-387.
}


   

\seealso{
  \code{modelSampler},
 \code{print.boot.modelSampler},
 \code{print.modelSampler},
 \code{plot.modelSampler},
 \code{plot.FPE},
  \code{plot.icicle},
\code{plot.var.stability},
  \code{plot.ooberror}. 
}
\examples{
  data(Pollute, package = "modelSampler") 
  ms.boot <- boot.modelSampler(MortRate~., Pollute, n.iter1 = 2500, 
  n.iter2=2500, B=20, verbose = TRUE)
  print(ms.boot)
  }
\keyword{file}
