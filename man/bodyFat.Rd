\name{bodyFat}
\docType{data}
\alias{bodyFat}
\title{Body Fat data}
\description{
Lists estimates of the percentage of body fat determined by
underwater weighing and various body circumference measurements for
251 men. (see references below).


A variety of popular health books suggest that the readers assess
their health, at least in part, by estimating their percentage of
body fat. In Bailey (1994), for instance, the reader can estimate
body fat from tables using their age and various skin-fold
measurements obtained by using a caliper. Other texts give
predictive equations for body fat using body circumference
measurements (e.g. abdominal circumference) and/or skin-fold
measurements. See, for instance, Behnke and Wilmore (1974), pp.
66-67; Wilmore (1976), p. 247; or Katch and McArdle (1977), pp.
120-132).

}

\format{bodyFat & Density determined from underwater weighing
 Percent body fat from Siri's (1956) equation. \dQuote{bodyFat} is the outcome variable.
 The data frame containing:  
\tabular{ll}{

age		\tab  Age (years)\cr
weight	\tab  Weight (lbs)\cr
height 	\tab  Height (inches)\cr
neck 		\tab  Neck circumference (cm)\cr
chest		\tab  Chest circumference (cm)\cr
abdomen 	\tab  Abdomen 2 circumference (cm)\cr
hip 		\tab  Hip circumference (cm)\cr
thigh 	\tab  Thigh circumference (cm)\cr
knee 		\tab  Knee circumference (cm)\cr
ankle 	\tab  Ankle circumference (cm)\cr
biceps 	\tab  Biceps (extended) circumference (cm)\cr
forearm 	\tab  Forearm circumference (cm)\cr
wrist 	\tab  Wrist circumference (cm)
}

(Measurement standards are apparently those listed in Benhke and
Wilmore (1974), pp. 45-48 where, for instance, the abdomen 2
circumference is measured "laterally, at the level of the iliac
crests, and anteriorly, at the umbilicus".)

}

\source{
The data set have been taken from Statlib at  \url{http://lib.stat.cmu.edu/datasets/}
}


\references{

Bailey, C.(1994). Smart Exercise: Burning Fat, Getting Fit.
\emph{Houghton-Mifflin Co., Boston}, 179 -- 186.

Behnke, A.R. and Wilmore, J.H. (1974). \emph{Evaluation and
Regulation of Body Build and Composition, Prentice-Hall, Englewood
Cliffs, N.J.}

Siri, W.E. (1956), Gross composition of the body, \emph{Advances in
Biological and Medical Physics, vol. IV, edited by J.H. Lawrence
and C.A. Tobias, Academic Press, Inc., New York.}

Katch, Frank and McArdle, William (1977). Nutrition, Weight Control,
and Exercise, \emph{Houghton Mifflin Co., Boston.}

Wilmore, Jack (1976). Athletic Training and Physical Fitness:
Physiological Principles of the Conditioning Process, \emph{Allyn
and Bacon, Inc., Boston.}
}

\examples{data(bodyFat, package = "modelSampler")}
\keyword{datasets}

