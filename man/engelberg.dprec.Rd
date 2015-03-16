\name{engelberg.dprec}
\alias{engelberg.dprec}
\docType{data}
\title{ Example datasets for gpd }
\description{
  Daily precipitation at Engelberg Switzerland from Jan 1, 1901 till 
  Dec 31, 2005.
  Data are in mm.

}
\usage{
  data(engelberg.dprec)
}
\format{
  A data frame with two columns: First column (prec) is the daily 
  precipitation (in mm per day), second column (date) is the date of 
  the measurement. Format of the date is "yyyy.mm.dd".
}
\details{

}
\source{
  Data from MeteoSwiss Zurich, Switzerland
}
\references{
  
}
\examples{
data(engelberg.dprec)
engelberg.dprec[1:368,]

}
\keyword{datasets}
