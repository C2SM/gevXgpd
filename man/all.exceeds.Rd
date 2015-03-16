\name{all.exceeds}
\alias{all.exceeds}
\title{ Extract peaks from dataset }
\description{
  A trivial method of declustering for subsequent peaks-over-threshold
  analysis. Extracts those items of a dataset exceeding (or beeing equal) 
  to a dataset. Usually this function is not called by the user directly.
  Use \code{extract.peaks} instead.
}
\usage{
all.exceeds(dat, threshold)
}
\arguments{
  \item{dat}{ A dataframe with the measurements in the first column. }
  \item{threshold}{ The threshold. }
}
\details{
  No declustering is performed here at all. Only exceedances extracted.
}
\value{
  A subset of rows of the input data frame \code{dat}, i.e. threshold
  exceedances
}
\references{  }
\author{ Christoph Frei }
\note{  }
\seealso{ \code{\link{extract.peaks}}, \code{\link{min.sep}} }
\examples{
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
