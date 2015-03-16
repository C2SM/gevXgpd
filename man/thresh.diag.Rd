\name{thresh.diag}
\alias{thresh.diag}
\title{ Diagnostics for GPD threshold choice }
\description{
  Calculates diagnostics for choosing a suitable threshold in a 
  peak over threshold analysis. This function scans through a series
  of thresholds and fits a GPD to the threshold exceedances. Mean
  exceedances, modified scale and shape are returned (for each threshold)
  together with corresponding confidence intervals. Results can be
  plotted in \code{plot.thresh.diag}.
}
\usage{
thresh.diag(pks, threshs, probs = c(0.025, 0.5, 0.975), min.exceeds = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pks}{ A data frame with independent peaks. This can be obtained
      from a data set through \code{extract.peaks} }
  \item{threshs}{ A vector of thresholds for which to calculate 
      the diagnostics. }
  \item{probs}{ The probs defining the confidence interval and estimated
      value of the diagnostics. It is necessary to include 0.5 in the vector.
      The first and last element define the probs for the lower and upper
      quantile of the confidence. }
  \item{min.exceeds}{ Minimum number of exceedances necessary to calculate
      the diagnostics. If there are less exceedances NA is returned in
      the diagnostics. }
}
\details{
  Confidence intervals are maximum likelihood intervals for the shape
  and modified scale. For mean exceedances the confidence is estimated
  from the asymptotic normality of means directly from the sample.
}
\value{
  A list with
  \item{threshs }{The vector of thresholds as in the input.}
  \item{num.exceed }{A vector with the number of threshold exceedances
       for each threshold.}
  \item{q.scale.star }{A matrix with quantiles of the modified scale. 
       Results are for all threshs in the different rows. Different
       quantiles are in different columns.}
  \item{q.shape }{A matrix with quantiles of the shape. Structure as
       above.}
  \item{q.mean.exceed }{A matrix with quantiles of the mean exceedance. 
       Structure as above.}
}
\references{ Coles 2001 }
\author{ Christoph Frei }
\note{ }

\seealso{ \code{\link{plot.thresh.diag}} }
\examples{
data(engelberg.dprec)
threshold <- 5
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
threshs <- seq(5,85,by=2.5)
diag <- thresh.diag(pks,threshs)

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
