\name{Pollute}
\docType{data}
\alias{Pollute}
\title{Pollute data}
\description{
 This data set consist of mortality rates against socio-economic, 
meteorological and pollution variables for 60 statistical areas in the USA. For details of the 
data set, see references below.

}

\format{ The data set have 60 observations with 15 variables. \dQuote{MortRate} is the outcome variable.
 The data frame containing:  
\tabular{ll}{

Rainfall	\tab  Mean annual precipitation in inches\cr
JanTemp	\tab  Mean January temperature in degrees Fahrenhiet\cr
JulTemp 	\tab  Mean July temperature in degrees Fahrenhiet\cr
X.Over65 	\tab  Percent of 1960 SMSA population which is 65 years of age or over\cr
AvHHSize	\tab  Population per household\cr
YrsEduc 	\tab  Median school years completed for those over 25 in 1960 SMSA\cr
X.GdHouse	\tab  Percent of housing units which are sound with all facilities\cr
PopDens 	\tab  Population per square mile in urbanized area in 1960\cr
X.NWhite	\tab  Percent of 1960 urbanized area population which is non-white\cr
WCollar 	\tab  Percent employment in white-collar occupations in 1960 urbanized area\cr
X.Poor 	\tab  Percent families with income under \$3000 in 1960 urbanized area(cm)\cr
HydCarb 	\tab  Relative population potential of hydrocarbons\cr
NOx	 	\tab  Relative population potential of oxides of nitrogen\cr
SO2		\tab  Relative population potential of sulpher dioxides\cr
RelHumid	\tab  Percent relative humidity, annual average at 1 PM\cr
MortRate	\tab  Total age adjusted mortality rate, expressed as death per 100000 population.
}
}

\source{
Alan Miller at \url{http://users.bigpond.net.au/amiller/}
}

\references{
Gunst, R.F. and Mason, R.L. (1980). \emph{Regression Analysis and its Application, Mercel Decker: New York},  
370 -- 371.


McDonald, G.C. and Schwing, R.C. (1973). Instabilities of regression coefficient estimates relating air pollution to mortality.\emph{Technometrics}, \bold{15},
 463 -- 482.
}

\examples{data(Pollute, package = "modelSampler")}
\keyword{datasets}

