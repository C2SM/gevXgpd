\name{fitGPD}
\alias{fitGPD}
\title{ Fit a GPD distribution }
\description{
  Estimate parameters of a GPD or (more specifically) a EXPONENTIAL
  distribution to a data set. This is a common task of Extreme
  Value Analysis when the data samples are independent exceedances of
  a large threshold. (Peak over Threshold Method.) Two methods
  of estimation are currently implemented: Maximum Likelihood Estimation, 
  and L-moments estimation. 
}
\usage{
fitGPD(data.df,chi=min(data.df[,1],na.rm=TRUE),
       lamda=1,ret=c(5,10,50,100,500), 
       dist = "GPD", estim = "mlik", method="Nelder-Mead", maxit=10000)
}
\arguments{
  \item{data.df}{ A dataframe with the peak values in the first 
       column. Additional columns are possible but are not used in the
       function.  }
  \item{chi}{ The threshold used to determine the peaks. Default is the
       minimum peak value. It is recommended to specify the threshold
       explicitly. A value of \code{chi=0} can be used if the data are
       exceedances. }
  \item{lamda}{ The number of peaks (exceedances) per time unit.
       The default \code{lamda=1} makes that return periods 
       in \code{ret} are given in units of the average recurrence
       between peaks. }
  \item{ret}{ A vector of return periods for which fitted return levels 
     should be provided in the resulting list (element \code{fit}. }
  \item{dist}{ With \code{dist="EXPON"}, the fit is restricted to the family 
     of Exponential distributions. \code{dist="GPD"} (default) is for fitting
     GPD in general.}
  \item{estim}{ The method to use for parameter estimation: 
     \code{estim="mlik"}: Maximum Likelihood, 
     \code{estim="lmom"}: L-moments }
  \item{method}{ A character string specifying the method of numerical 
     optimisation to be used in the likelihood based estimation. 
     This parameter is passed on to optim(). }
  \item{maxit}{ The maximum number of iterative steps in optim() for 
     likelihood estimation. }

}
\details{
  The estimation of a GPD distribution is considered here as a three
  parameter problem where the frequency of threshold exceedences 
  \code{lamda} is a parameter too. The estimated \code{lamda} is 
  usually not exactly the \code{lamda} provided as input argument which
  is due to numerical inaccuracies. 

  As starting values for the iterative maximum likelihood 
  estimations the method of moments estimates for an exponential 
  fit are taken, i.e.
  shape=0 and scale= mean of exceedances. 

  If there is no convergence in \code{maxit} iterative steps a warning
  is issued. Moreover the parameter values as well as the values of 
  the fitted return levels (list element \code{fit} in the result list
  are set to NA. Also the covariance matrix 
  (list element \code{cov}) will be left as a matrix filled with NAs. 
}
\value{
  A list with the following elements
  \item{data }{The input data.frame for which the distribution was fitted.}
  \item{lamda }{The frequency of peaks (exceedances) per time unit.}
  \item{alpha }{The scale parameter of the fitted distribution 
     (single value).}
  \item{chi }{The threshold used.}
  \item{k }{The shape parameter of the fitted distribution (single value). 
     Note that we use a sign convention, which
     differs from that in Coles (2001) to keep consistency with the 
     convention in block maximum analysis.}
  \item{fit }{A two column matrix. First column: Selection of return
     periods (from the input). Second column: corresponding extreme 
     values (return levels).}
  \item{cov }{A 3x3 matrix (in the general GPD case) or a 2x2 matrix (in 
     EXPON case) giving the variance-covariance of the parameter 
     estimates. (Order is \code{lamda},\code{alpha},\code{k}.) 
     This is determined within optim() - the inverse of the hessian
     matrix - and is only available for maximum likelihood estimates. It
     is used for calculating confidence intervals.}
  \item{dist }{A character string indicating the distribution family that 
     was fitted to the data. Can be "GPD" (Generalized Extreme Value) or
     "EXPON" (Exponential Distribution). }
  \item{ks.pval }{The p-value of the Kolmogorov-Smirnov Goodness-of-fit
     test. This is the result of applying \code{\link{ks.test}}. This should
     not be interpreted as confirming goodness-of-fit, because the 
     the test presumes that the data has not been used for estimating the
     parameters of the distribution. The result may help to identify serious
     misfits that may be due to data quality problems or misconception of
     the usage of GPD. }
  \item{estim }{A character string indicating the method of estimation 
     used (i.e. the method specified in the inputs to \code{fitGPD} }
  \item{method}{ A character string specifying the method of numerical 
     optimisation. Same as input. Only returned for \code{estim="mlik"}. }
  \item{maxit}{ The maximum number of iterative steps in optim(). 
     Same as input. Only returned for \code{estim="mlik"}. }
}
\references{ Palutikov et al. 1999 for a general description. 
   Hosking 1990, (Eq. 2.3 and p. 114 and Table 2) for L-moments 
   estimation. Coles 2001 for likelihood estimation.  }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{plot.pot}}, \code{\link{conf.bounds.pot}}, 
          \code{\link{plot.conf.pot}} }
\examples{
data(engelberg.dprec)
threshold <- 30
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
npks <- length(pks[,1])
nyrs <- as.numeric(substr(
                  engelberg.dprec[dim(engelberg.dprec)[1],"date"],1,4))-
             as.numeric(substr(engelberg.dprec[1,"date"],1,4)) + 1
dat.pot <- fitGPD(pks,lamda=npks/nyrs,chi=threshold)
dat.pot$alpha
dat.pot$k
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
