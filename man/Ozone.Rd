\name{Ozone}
\docType{data}
\alias{Ozone}
\title{Ozone data}
\description{
A data frame with 203 observations on 13 variables, each observation
is one day. Original data is modified to delete all the missing
values.
}

\format{
\dQuote{ozone} being the outcome variable, other 12 variables are
independent variables: 
\tabular{ll}{

month		\tab Month: 1 = January, ..., 12 = December\cr
daymonth	\tab Day of month\cr
dayweek	\tab Day of week: 1 = Monday, ..., 7 = Sunday\cr
ozone		\tab Daily maximum one-hour-average ozone reading\cr
press		\tab 500 millibar pressure height (m) measured at Vandenberg AFB\cr
wind		\tab Wind speed (mph) at Los Angeles International Airport (LAX)\cr
humid		\tab Humidity (\%) at LAX\cr
tempS		\tab Temperature (degrees F) measured at Sandburg, CA\cr
tempE		\tab Temperature (degrees F) measured at El Monte, CA\cr
ibht		\tab Inversion base height (feet) at LAX\cr
pressg	\tab Pressure gradient (mm Hg) from LAX to Daggett, CA\cr
ibtemp	\tab Inversion base temperature (degrees F) at LAX\cr
vis		\tab Visibility (miles) measured at LAX
}
}
\source{
Leo Breiman, Department of Statistics, UC Berkeley.  Data used in
Leo Breiman and Jerome H. Friedman (1985), Estimating optimal
transformations for multiple regression and correlation, \emph{J. Amer. Stat. Assoc.}, \bold{80}, 
580 -- 598.
}

\examples{data(Ozone, package = "modelSampler")}
\keyword{datasets}

