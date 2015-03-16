\name{lik.prof}
\alias{lik.prof}
\title{ Calculate profile of the likelihood function }
\description{
  For a specified likelihood function caculate a profile near the
  best fit by varying one of the parameters. (Maximum of the likelihood
  function with one parameter kept fix. This is used to calculate
  profile likelihood confidence intervals for that parameter. 
  The confidence distribution 
  is calculated for the selected probabilities. The sampling of the
  profile stops as
  soon as the outermost confidence probability is reached.
  This function is usually not called directly by the user. Likelihood
  profile confidence intervals should be determined from
  \code{conf.bounds.xval.lprof}.
}
\usage{
lik.prof(data, LikFun, index, start, 
         probs = c(0.025, 0.975), frac.incr = 1/10, 
         method = "BFGS", maxit = 10000, ...)
}
\arguments{
  \item{data}{ A single vector with the dataset. }
  \item{LikFun}{ The likelihood function. This must be a prespecified
     function with arguments \code{par} - a vector with the parameters -
     and \code{data} - the dataset. There may be additional arguments
     to \code{LikFun} which can be specified through the \code{\dots} 
     argument in \code{lik.prof}.}
  \item{index}{ An integer number giving the index of the parameter in
     argument \code{par} to \code{LikFun} for which to calculate confidence
     quantiles, i.e. the parameter for which to profile the likelihood 
     function. }
  \item{start}{ A vector giving the starting value to use in the 
     optimisation. This is a vector of the same length (and order) like
     argument \code{par} in \code{LikFun}. It is recommended to use the
     maximum likelihood best fit as starting value. }
  \item{probs}{ The probabilities for which confidence quantiles should
     be calculated. }
  \item{frac.incr}{ A real valued increment to use in
     profiling the likelihood function. Is specified in units of the
     sampling distribution standard deviation of the respective parameter.
     The sampling distribution standard deviation is calculated from
     the hessian matrix of the likelihood surface at its maximum. }
  \item{method}{ The method of optimisation to use. This argument is
     simply passed on to \code{optim}. Here we do not use the default
     method ("Nelder-Mead") of \code{optim} because this is unreliable for 
     one dim optimisation and issues warnings. One-dim optimisation is  
     performed when a likelihood profile is calculated for a GUMBEL or
     an EXPON distribution. }
  \item{maxit}{ The maximum number of iterations to use in maximising
     the likelihood surface. }
  \item{\dots}{ Additional arguments to \code{LikFun} if necessary. }
}
\details{
  The profiling is implemented as follows:

  (a) First the maximum likelihood best fit is calculated. Moreover the 
  covariance matrix of the parameters at the likelihood maximum is 
  determined.

  (b) An increment "dpar" of the parameter is calculated as 
  the \code{frac.incr}
  fraction of one standard deviation. The increment determines the 
  resolution with which the likelihood profile is sampled.

  (c) The profile is calculated for parameter values 
  par.best.fit + j * dpar with j=1,2,3,... and with j=-1,-2,-3,...
  The starting value of the optimisation is taken to be the parameter
  set that maximises the likelihood at the particular preset of the
  parameter under question. The profiling is stopped as soon as the
  most extreme quantile (according to the inputs in \code{probs}) is
  reached.

  (d) All the quantiles of the sampling distribution are subsequently 
  calculated by interpolation of the likelihood profile with a spline.
}
\value{
  A list with elements:
  \item{conf.par }{A vector with the confidence quantiles of the 
     parameter for the specified probs. Formally this is a matrix
     so that the columns can be named with the \code{probs}. }
  \item{profile }{A matrix with 3 columns representing the profile.
     Column 1: parameter value, column 2: the maximum of the 
     likelihood function conditional on the parameter value, column 3: 
     the probability of the parameter in the sampling distribution 
     (i.e. using the chi-squared distribution for the deviance statistic.}
}
\references{ See Coles 2001 for max likelihood confidence intervals.  }
\author{ Christoph Frei }
\note{  }

\seealso{ The likelihood function for GEV \code{\link{gev.lik}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
