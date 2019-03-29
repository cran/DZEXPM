\name{DZEXPM-package}
\alias{DZEXPM-package}
\alias{DZEXPM}
\docType{package}
\title{
Estimation and Prediction of Skewed Spatial Processes 
}
\description{
DZEXPM contains the function to estimate and predict skewed spatial processes. 
}
\details{
Asymmetric spatial processes arise naturally in finance, economics, hydrology
and ecology. For such processes, Double Zero Expectile (DZEXP) normal process is proposed in Majumdar and Paul (2015). By using a Bayesian methodology, Wang, Yang and Majumdar (2018) show that by adding measurement error to the DZEXP model (DZEXPM), a reasonably flexible model is obtained, which is also computationally tractable in the literature. 

As an example, a skewed data set on maximum annual temperature obtained from weather stations in Louisiana and Texas of year 2003 is attached to test the DZEXPM package.
}
\author{
Jiangyan Wang, Miao Yang and Anandamayee Majumdar
}
\references{
Majumdar, A and Paul, D. (2015). Zero expectile processes and Bayesian spatial regression. Journal of Computational and Graphical Statistics, 25(3).

Wang, J., Yang, M. and Majumdar, A. (2018). Comparative Study and Sensitivity Analysis of Skewed Spatial Processes. Computational Statistics, 33, 75-98.
}
%%\examples{
%%dzexpm(y_ful, x_ful, n_ful, n, u1, u2, theta, p, iter, loopiter)
%%}
\keyword{spatial statistics, MCMC, Bayesian method, spatial regression, zero expectile}
