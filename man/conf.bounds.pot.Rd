\name{conf.bounds.pot}
\alias{conf.bounds.pot}
\title{ Confidence bounds for GPD analysis }
\description{
  Calculates confidence intervals for return values (quantiles) 
  and for distribution parameters of a GPD fit to data. 
  Three types of confidence intervals can be
  calculated: asymptotic maximum likelihood intervals, intervals
  estimated by parameteric resampling and intervals from likelihood
  profiles. The asymptotic mle is the default. 
  Resampling and likelihood profiles are more accurate but take 
  longer to calculate. Asymptotic mle intervals are always symetric 
  with respect to the fitted value.
}
\usage{
conf.bounds.pot(pot, ret = c(5, 10, 50, 100, 500), 
                probs = c(0.025, 0.975), size = 500, cal = "MLE",
                profs.out = FALSE)
}
\arguments{
  \item{pot}{ An object from \code{fitGPD} containing the data and
         fitted GPD distribution. }
  \item{ret}{ A vector of return periods for which the confidence 
         intervals shall be calculated. }
  \item{probs}{ A vector of probs specifying the quantiles of the 
         confidence intervals. Default is the 95% confidence interval. 
         It is possible (recommended) to calculate several confidence
         bounds at once, e.g. probs = c(0.025,0.25,0.5,0.75,0.975) }
  \item{size}{ The number of samples to simulate when \code{cal="SIM"} i.e. 
         the interval is calculated by parameteric resampling. }
  \item{cal}{ The method to use for calculating the confidence bounds. 
         A character string. Either \code{cal="SIM"} (parametric resampling),
         \code{cal="MLE"} (asymptotic maximum likelihood) or 
         \code{cal="LPROF"} (likelihood profile). }
  \item{profs.out}{ TRUE or FALSE, if likelihood profiles should be provided
         in the output. Is only relevant for \code{cal="LPROF"}.}
}
\details{
  The calculation with parametric resampling can take long, e.g. nearly
  one minute. Also \code{cal="LPROF"} can take some time to calculate. For
  \code{cal="SIM"} the time to calculate depends on \code{size} only, while
  for \code{cal="LPROF"} the time increases linearly with the number of
  return values \code{ret} specified.

  This function is also used to calculate the confidence intervals that
  will be plotted in \code{plot.conf.pot}. In this case it is meaningful
  to specify several values of return periods in \code{ret} for nice 
  plotting afterwards.
}
\value{
  A list with elements
  \item{conf.retv }{A matrix with length(ret) rows and length(probs) columns. 
           Column k represents the quantile of the confidence distribution
           for the k-th probability in \code{probs}.
           Row i represents the quantiles of the confidence distribution
           of return values for the i-th return period in \code{ret}. }
  \item{conf.paras}{A matrix with three rows and length(probs) columns.
           Quantiles
           of the confidence distribution for each parameter are given
           in the different rows. row 1 for lamda (number of 
           exceedances per time unit), row 2 
           for scale parameter, row 3 for shape parameter.}
  \item{paras.sim}{Only provided for \code{cal="SIM"}. An array with 3
           columns containing the parameter estimates for the simulated
           data sets. This could be used to determine confidence
           distributions for functions of the parameters. }
  \item{paras.cov}{Only provided for \code{cal="MLE"}. A 3x3 (or 2x2 in
           case of exponential fits) matrix giving
           the variance-covariance matrix of the parameter estimation.
           This is a direct result of the fit. It can be used to 
           determine asymptotic MLE confidence distributions for functions 
           of the parameters. The ordering of parameters is lamda (number of 
           events per time unit), alpha (scale parameter), k (shape 
           parameter).}
  \item{lamda }{A single number: Number of exceedances (peaks) per time 
           unit. The same as that in argument \code{xval}}
  \item{cal }{A character string giving the method of calculation for
           the confidence interval. This is the same like in input.}
  \item{prof.chi }{Only if \code{cal=="LPROF"} and \code{profs.out==TRUE}. 
           The likelihood profile 
           for parameter \code{chi}, i.e. position. The output is a 
           three column matrix, in turn the chi parameter, the conditional 
           likelihood maximum and the probs of the sampling distribution.}
  \item{prof.alpha }{Only if \code{cal=="LPROF"} and \code{profs.out==TRUE}. 
           The likelihood profile 
           for parameter \code{alpha}, i.e. scale. See also above.}
  \item{prof.k }{Only if \code{cal=="LPROF"} and \code{profs.out==TRUE}. 
           The likelihood profile 
           for parameter \code{k}, i.e. shape. See also above.}

}
\references{ Coles 2001 for the maximum likelihood confidence interval
       and likelihood profile confidence interval. }
\author{ Christoph Frei }
\note{  }
\seealso{ \code{\link{fitGPD}}, \code{\link{plot.conf.pot}} }
\examples{
data(engelberg.dprec)
threshold <- 30
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
npks <- length(pks[,1])
nyrs <- as.numeric(substr(
                  engelberg.dprec[dim(engelberg.dprec)[1],"date"],1,4))-
             as.numeric(substr(engelberg.dprec[1,"date"],1,4))
dat.pot <- fitGPD(pks,lamda=npks/nyrs,chi=threshold)
conf.sim <- conf.bounds.pot(dat.pot,cal="SIM")
conf.sim$conf.retv
conf.mle <- conf.bounds.pot(dat.pot,cal="MLE")
conf.mle$conf.retv
conf.lprof <- conf.bounds.pot(dat.pot,cal="LPROF")
conf.lprof$conf.retv

conf.sim$conf.paras
conf.mle$conf.paras
conf.lprof$conf.paras
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
