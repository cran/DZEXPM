\name{dzexpm}
\alias{dzexpm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The estimation and prediction function for the skewed spatial processes
}
\description{
This function provides a way to combine the estimation and prediction procedure for the skewed spatial processes.
}
\usage{
dzexpm(y_ful, x_ful, n_ful, n, u1, u2, theta, p, iter, loopiter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y_ful}{the full response variable of a dataset}
\item{x_ful}{the full covariates of a dataset}
\item{n_ful}{the full number of the spatial locations in a dataset}
\item{n}{the number of spatial locations being employed to estimate and evaluate the spatial process}
\item{u1}{horizontal coordinate of the spatial locations}
\item{u2}{vertical coordinate of the spatial locations}
\item{p}{is a useful parameter to capture the "skewness" of the process, which belongs in (0,1)}
\item{theta}{the range parameter of the covariance kernel in the spatial process}
\item{iter}{the MCMC iterations to update the parameters}
\item{loopiter}{the burn in size}
}
\details{
This function is generated according to the spatial model, Double Zero Expectile Normal Process with measurement error, which has good prediction performance. See the references for details.
}
\value{
A date.frame is returned which includes quantities "DIC", "coverage", "bias.med", "mpe.med" and "SD.med", where "DIC" stands for the Deviance Information Criterion (DIC), see Gelman et al. (2004), and "coverage", "bias.med", "mpe.med" and "SD.med" represent the median of prediction coverage, bias, mean
prediction error and standard error over MCMC scheme iterations, respectively, for the predictive performance of the skewed spatial process.
}
\references{
Wang, J., Yang, M. and Majumdar, A. (2017). Comparative Study and Sensitivity Analysis of Skewed Spatial Processes.


Gelman, A., Carlin, J., Stern, H. and Rubin, D. (2004). Bayesian Data Analysis: Second Edition. Texts in Statistical Science. CRC Press. ISBN 1-58488-388-X. LCCN 2003051474. MR 2027492.
}
\author{
Jiangyan Wang, Miao Yang and Anandamayee Majumdar
}
%%\note{
%%  ~~further notes~~
%%}

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
  real<- as.matrix(do.call(cbind, MaxTemp03))
  n_ful<- nrow(real)
  y_ful<- real[,2]
  x_ful<- real[,7]
  x_ful<- log(x_ful)
  x_ful<- scale(x_ful)
  u1<- real[,3]
  u2<- real[,4]
  theta<- 2/150000
  n<- 56
  p<- 0.1
 dzexpm(y_ful, x_ful, n_ful, n, u1, u2, theta, p, 50, 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
