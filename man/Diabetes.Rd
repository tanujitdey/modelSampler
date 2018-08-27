\name{Diabetes}
\docType{data}
\alias{Diabetes}
\title{Diabetes data}
\description{
In this data set, 442 diabetes patients were measured on 10 baseline variables. For details see reference below.

}

\format{
 The data frame containing:  
\tabular{ll}{

age			\tab 	  Age in years\cr
sex			\tab 	  categorical variable\cr
bmi			\tab 	  Body mass index\cr
bp			\tab    Average blood pressure\cr
s1 - s6		\tab    Six blood serum measurements\cr
Y			\tab    A quantitative measure of disease progression recorded one year after baseline measurement.
}
}

\source{
The data set have been taken from \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/}
}

\references{
Efron, B. , Hastie, T. ,  Johnstone, I. and Tibshirani, R. (2003). Least angle regression. \emph{Annals of Statistics}, \bold{32(2)}, 
407 -- 499.
}

\examples{data(Diabetes, package = "modelSampler")}
\keyword{datasets}

