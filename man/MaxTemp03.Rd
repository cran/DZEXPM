\name{MaxTemp03}
\alias{MaxTemp03}
\docType{data}
\title{
Maximum annual temperature in the year 2003 of 76 spatially located stations in the states of Louisiana and Texas
}
\description{
This data set contains the maximum annual temperature in the year 2003
of 76 spatially located stations in the states of Louisiana and Texas, as the response
variable, and elevation as the covariate.
}
\usage{data("MaxTemp03")}
\format{
  A data frame with 76 observations on the following 7 variables.
  \describe{
    \item{\code{stn}}{a numeric vector for the station number}
    \item{\code{temMax}}{a numeric vector for the maximum temperature}
    \item{\code{utmn}}{a numeric vector for the North UTM (UNIVERSAL TRANSVERSE MERCARTOR GRID SYSTEM)}
    \item{\code{utme}}{a numeric vector for the South UTM (UNIVERSAL TRANSVERSE MERCARTOR GRID SYSTEM)}
    \item{\code{lat}}{a numeric vector for latitude}
    \item{\code{long}}{a numeric vector for longitude}
    \item{\code{elevate}}{a numeric vector for elevation}
  }
}
\details{
This is an illustrative data
set obtained for a collection of stations in Texas and Louisiana. The collection was from
the NOAA Satellite and Information Service involving global climate and weather. 
}
\source{
The data link is: http://www.ncdc.noaa.gov /cgi-bin/res40.pl?page=gsod.html
}
\references{
Wang, J., Yang, M. and Majumdar, A. (2017). Comparative Study and Sensitivity Analysis of Skewed Spatial Processes.
}
%%\examples{
%%ata(data)
%%## maybe str(data) ; plot(data) ...
%%}
\examples{
data(MaxTemp03)
}
\keyword{datasets}
