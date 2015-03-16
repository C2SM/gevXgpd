"fitEXPON.lmom" <-
function(data.df,chi=min(data.df[,1],na.rm=TRUE),
         lamda=1,ret=c(5,10,50,100,500),...) {
# ==================================================================
#  estimates parameters of Exponential Distr. fit to data.
#  the estimation is based on the method of L-moments (Hosking 1990)
#  also estimates values with a given return period
#  
#  ARGUMENTS:
#  data.df  : one column dataframe with data (those data exceeding
#             the threshold but in original units)
#  OPTIONAL ARGUMENTS:
#  chi      : threshold used to select the data. Default is chi=0
#             such that threshold excess values y=x-chi are
#             treated without specification of chi.
#  lamda    : frequency of threshold exceedence per year. 
#             This is needed to determine return periods in years.
#             Default value is 1. In this case return periods are
#             in time units of average event return periods.
#  ret      : return periods for which extreme values will be 
#             estimated; can be a list of values > 1.
#             
#  OUTPUT: (a list of type xval with the following elements)
#  data:    original dataframe
#  lamda:   mean exceedence per year
#  alpha:   scale parameter of EXPON
#  chi  :   fixed threshold value used for event definition
#  k    :   shape parameter of GPD (k=0, forced for EXPON)
#  fit  :   extreme values for predefined return periods ret
#  dist :   what distribution was fitted (GPD or EXPON)
#  estim:   what estimator was used (mlik or lmom)

   # PREPARATIONS
   nam <- names(data.df)[1]
   len  <- dim(data.df)[1]
   sor <- sort(as.array(data.df[,nam]))
   jm1 <- seq(1,len)-1
   jm2 <- seq(1,len)-2

   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # CHECK IF LARGER THAN THRESHOLD
   if (any(sor<chi)) {
       stop("** ERROR ** data must be larger than threshold chi")}
   
   # CALCULATE PSTAR AND B VALUES 
   # (see Hosking 1990, Eq. 2.3 and p. 114.)
   p00 <-  1
   p10 <- -1
   p11 <-  2
   p20 <-  1
   p21 <- -6
   p22 <-  6
   b0 <- mean(sor)
   b1 <- mean(jm1/(len-1)*sor)
   b2 <- mean(jm1*jm2/(len-1)/(len-2)*sor)

   # ESTIMATE L-MOMENTS (actually only l1 and l2 needed for GPD)
   l1 <- p00*b0
   l2 <- p10*b0+p11*b1
   l3 <- p20*b0+p21*b1+p22*b2

   # ESTIMATE PARAMETERS FROM L-MOMENTS  
   # (see Hosking 1990, Table 2, notice that l1 actually is l1-chi)
   #k     <- (l1-chi)/l2 - 2.0
   k     <- 0.0
   alpha <- (1+k)*(l1-chi)

   # Conduct Kolmogorov-Smirnov Goodness-of-Fit Test
   # Note: This is a misuse of the KS-Test because the distribution paremeters
   # were determined from the data itself. The result may nonetheless pinpoint to
   #Êserious misfits.
   ks.pval <- suppress.Warnings(
           ks.test(x=data.df[,nam]-chi,y="FGPD",alpha=alpha,k=k)$p.value)

   # DETERMINE EXTREME VALUES FOR GIVEN RETURN PERIODS
   fitted <- XGPD(foft(ret),alpha,chi,k)  
 
   # CONSTRUCT pot OBJECT FOR OUTPUT
   list(
       data=data.df,
       lamda=lamda,
       alpha=alpha,
       chi=chi,
       k=k,
       fit=cbind(ret,fitted),
       dist="EXPON",
       ks.pval=ks.pval,
       estim="lmom"
   )
}

