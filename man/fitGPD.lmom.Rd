\name{fitGPD.lmom}
\alias{fitGPD.lmom}
\alias{fitGPD.mlik}
\alias{fitEXPON.lmom}
\alias{fitEXPON.mlik}
\title{ Fit GPD / Exponential distribution }
\description{
  The functions for fitting GPD or Eexponential distributions to data. These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{fitGPD} instead.
}
\usage{
fitGPD.lmom(data.df, chi = min(data.df[,1],na.rm=TRUE), 
            lamda = 1, ret = c(5, 10, 50, 100, 500), ...)
fitEXPON.lmom(data.df, chi = min(data.df[,1],na.rm=TRUE), 
            lamda = 1, ret = c(5, 10, 50, 100, 500), ...)
fitGPD.mlik(data.df, chi = min(data.df[,1],na.rm=TRUE),
            lamda = 1, ret = c(5,10,50,100,500),
            method="Nelder-Mead", maxit=10000, ...)
fitEXPON.mlik(data.df, chi = min(data.df[,1],na.rm=TRUE),
            lamda = 1, ret = c(5,10,50,100,500),
            method="Nelder-Mead", maxit=10000, ...) 
}
\arguments{
  \item{data.df}{ A data frame with the data to be fitted in the first 
     column. The data frame
     can have more than one column. Additional columns are not used.}
  \item{chi}{ The threshold used to determine the peaks. Default is the
       minimum peak value. It is recommended to specify the threshold
       explicitly. A value of \code{chi=0} can be used if the data are
       exceedances. }
  \item{lamda}{ The number of observations (blocks) per time unit. 
     Values other than
     lamda=1 can be used to express the return period in other units than 
     in units of the block size.}
  \item{ret}{ A vector of return periods for which fitted return levels 
     should be provided in the resulting list (element \code{fit}.}
  \item{method}{ A character string specifying the method of numerical 
     optimisation to be used in the maximum likelihood (\code{estim="mlik"}) 
     and the modified maximum likelihood (\code{estim="pmlik"}) estimations. 
     This parameter is passed on to optim(). }
  \item{maxit}{ The maximum number of iterative steps in optim() for 
     maximum likelihood (\code{estim="mlik"}) and modified maximum 
     likelihood (\code{estim="pmlik"}) estimation. }
  \item{\dots}{ Placeholder for parameters not used by the method. }
}
\details{
  These are the functions for individual methods and should not be used by
  the user. The user should invoke \code{fitGPD} instead.
}
\value{
}
\references{  }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{fitGPD}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
