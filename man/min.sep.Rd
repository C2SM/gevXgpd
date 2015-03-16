\name{min.sep}
\alias{min.sep}
\title{ Extract peaks from a dataset }
\description{
  A simple declustering technique for subsequent peaks-over-threshold
  analysis. The method extracts all exceedances above (or equal) a threshold
  which are at least a certain time period apart. This function is not 
  usually called by the user directly. Use \code{extract.peaks} instead.
}
\usage{
min.sep(dat, threshold, min.sep = 1)
}
\arguments{
  \item{dat}{ A dataframe with the measurements in the first column. 
     The measurements need to be listed in regular time increments.}
  \item{threshold}{ The threshold. }
  \item{min.sep}{ The minimum separation between consecutive exceedances.
     If \code{min.sep=5} then two consecutive peaks need to be at least
     5 observation increments appart, i.e. 4 no-peaks measurements in
     between. If there are more than one threshold exceedances in a
     period of min.sep duration then the largest is returned.}
}
\details{
}
\value{
  A subset of the rows of \code{dat}.}
\references{  }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{extract.peaks}}, \code{\link{all.exceeds}} }
\examples{
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
