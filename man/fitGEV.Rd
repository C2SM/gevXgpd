\name{fitGEV}
\alias{fitGEV}
\title{ Fit a GEV distribution }
\description{
  Estimate parameters of a GEV or (more specifically) a GUMBEL 
  distribution to a set of data samples. This is a common task of Extreme
  Value Analysis for the case where the data samples are maxima taken from
  equally sized blocks of original data. (Block Maximum Method.) Three methods
  of estimation are currently implemented: Maximum Likelihood Estimation, 
  Modified Maximum Likelihood using a Geophysical Prior Distribution for
  the shape parameter, and L-moments estimation. 
}
\usage{
fitGEV(data.df,lamda=1,ret=c(5,10,50,100,500)/lamda, 
       dist = "GEV", estim = "mlik", method="Nelder-Mead", maxit=10000)
}
\arguments{
  \item{data.df}{ A data frame with the data to be fitted in the first 
     column. The data frame
     can have more than one column. Additional columns are not used.}
  \item{lamda}{ The number of observations (blocks) per time unit. 
     Values other than
     lamda=1 can be used to express the return period in other units than 
     in units of the block size.}
  \item{ret}{ A vector of return periods for which fitted return levels 
     should be provided in the resulting list (element \code{fit}.}
  \item{dist}{ With \code{dist="GUMBEL"}, the fit is restricted to the family 
     of GUMBEL Distributions. \code{dist="GEV"} (default) is for fitting
     GEV in general.}
  \item{estim}{ The method to use for parameter estimation: 
     \code{estim="mlik"}: Maximum Likelihood, 
     \code{estim="lmom"}: L-moments,
     \code{estim="pmlik"}: Modified Maximum Likelihood (Geophysical Prior)}
  \item{method}{ A character string specifying the method of numerical 
     optimisation to be used in the maximum likelihood (\code{estim="mlik"}) 
     and the modified maximum likelihood (\code{estim="pmlik"}) estimations. 
     This parameter is passed on to optim(). }
  \item{maxit}{ The maximum number of iterative steps in optim() for 
     maximum likelihood (\code{estim="mlik"}) and modified maximum 
     likelihood (\code{estim="pmlik"}) estimation. }
}
\details{
  As starting values for the iterative maximum likelihood and modified
  likelihood estimations the L-moment estimate is used. 

  If there is no convergence in \code{maxit} iterative steps a warning
  is issued. Moreover the parameter values as well as the values of 
  the fitted return levels (list element \code{fit} in the result list
  are set to NA. Also the covariance matrix 
  (list element \code{cov}) will be left as a matrix filled with NAs. 

  The modified maximum likelihood estimation is exactly like that in 
  Martins and Stedinger (2000). The prior distribution for shape is
  implemented hard wired, i.e. there is no control on that from the
  functions parameters.  
}
\value{
  A list with the following elements
  \item{data }{The input data.frame for which the distribution was fitted.}
  \item{lamda }{The number of observations per time unit. Values other than
     lamda=1 can be used to express the return period in other units than 
     in units of the block size.}
  \item{alpha }{The scale parameter of the fitted distribution 
     (single value).}
  \item{chi }{The location parameter of the fitted distribution 
     (single value).}
  \item{k }{The shape parameter of the fitted distribution (single value). 
     Note that we use sign convention of Zwiers and Kharin (1998), which
     differs from that in Coles (2001)}
  \item{fit }{A two column matrix. First column: Selection of return
     periods (from the input). Second column: corresponding extreme 
     values (return levels).}
  \item{cov }{A 3x3 matrix (in the general GEV case) or a 2x2 matrix (in 
     GUMBEL case) giving the variance-covariance of the parameter 
     estimates. This is determined within optim() and is only available
     for maximum likelihood and modfied maximum likelihood estimates.}
  \item{dist }{A character string indicating the distribution family that 
     was fitted to the data. Can be "GEV" (Generalized Extreme Value) or
     "GUMBEL" (Gumbel Distribution). }
  \item{ks.pval }{The p-value of the Kolmogorov-Smirnov Goodness-of-fit
     test. This is the result of applying \code{\link{ks.test}}. This should
     not be interpreted as confirming goodness-of-fit, because the 
     the test presumes that the data has not been used for estimating the
     parameters of the distribution. The result may help to identify serious
     misfits that may be due to data quality problems or misconception of
     the usage of GEV. }
  \item{estim }{A character string indicating the method of estimation 
     used (i.e. the method specified in the inputs to \code{fitGEV}). }
  \item{method}{ A character string specifying the method of numerical 
     optimisation. Same as input. Only returned for \code{estim="mlik"} and
     \code{estim="pmlik"}. }
  \item{maxit}{ The maximum number of iterative steps in optim(). 
     Same as input. Only returned for \code{estim="mlik"} and
     \code{estim="pmlik"}. }
}
\references{ Zwiers and Kharin (1998) for parameter convention. 
   Hosking 1990, (Eq. 2.3 and p. 114 and Table 2) for L-moments 
   estimation. Coles 2001 for likelihood estimation. Martins and 
   Stedinger (2000, p. 740) for modified maximum likelihood estimation.}
\author{ Christoph Frei }
\note{ 

}
\seealso{ \code{\link{DGEV}}, \code{\link{plot.xval}}, 
          \code{\link{conf.bounds.xval}} }
\examples{
data(engelberg.annmax)
# Fitting l-moments and maximum likelihood
xval.lmom <- fitGEV(engelberg.annmax,estim="lmom")
xval.mlik <- fitGEV(engelberg.annmax,estim="mlik")
xval.lmom[c("alpha","chi","k")]  # results for L-moments estimation
xval.mlik[c("alpha","chi","k")]  # results for maximum likelihood estimation

# fit a gumbel distribution (i.e. force shape to be zero)
xval.gumbel.mlik <- fitGEV(engelberg.annmax,estim="mlik",dist="GUMBEL")

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
