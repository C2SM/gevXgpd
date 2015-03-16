\name{extract.peaks}
\alias{extract.peaks}
\title{ Extract peaks from a dataset }
\description{
  This is a common entry point to declustering methods for extracting
  peaks from a dataset for later use in a peaks-over-threshold analysis.
  Declustering methods currently available "minimum separation" and
  "all exceedances" (i.e. no declustering). 
}
\usage{
extract.peaks(dat, threshold, method = "min.sep", ...)
}
\arguments{
  \item{dat}{ A data frame with the dataset from which to extract 
     peaks. The first column
     comprises the measurements taken at regular time intervals. 
     The remaining columns don't influence but they are output with
     the peaks. I.e. the date in one of the remaining columns is useful. }
  \item{threshold}{ A single value for the threshold. Peaks are only
     returned if they exceed (or equal) the threshold. }
  \item{method}{ The method of declustering to use. A character string.
     Methods available are: \code{method="all.exceeds"} which returns all
     measurements (rows of dat) that exceed the 
     threshold. \code{method="all.exceeds"} returns all measurements that
     exceed the threshold AND are at least a minimum time period separated.
     Additional parameters to the methods can be looked up in the help
     pages of the methods (i.e. \code{help("min.sep")}.}
  \item{\dots}{ Additional (specific) parameters to the methods. }
}
\details{
  
}
\value{
  A subset of rows of the input data frame \code{dat}, i.e. threshold
  exceedances satisfying the conditions of the declustering.
}
\references{  }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{all.exceeds}}, \code{\link{min.sep}} }
\examples{
data(engelberg.dprec)
pks <- extract.peaks(engelberg.dprec,threshold=30,
                          method="min.sep",min.sep=5)
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
