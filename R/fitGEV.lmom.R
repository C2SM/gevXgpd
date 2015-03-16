"fitGEV.lmom" <-
function(data.df,lamda=1,ret=c(5,10,50,100,500)/lamda, ...) {
# ==================================================================
#  Estimates parameters of GEV-fit to sample data.
#  The estimation is based on the method of L-moments (Hosking 1990).
#  Also estimates values with a given return period.
#  
#  ARGUMENTS:
#  data.df  : one column dataframe with data
#  OPTIONAL ARGUMENTS:
#  ret      : return periods for which extreme values will be 
#             estimated; can be a list of values > 1.
#             
#  OUTPUT: (a list of type xval with the following elements)
#  data:    original dataframe
#  alpha:   scale parameter of GEV
#  chi  :   location parameter of GEV
#  k    :   shape parameter of GEV
#  fit  :   extreme values for predefined return periods ret
#  dist :   what distribution was fitted (GEV or GUMBEL)
#  estim:   what estimator was used (mlik or lmom)

   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # PREPARATIONS
   nam <- names(data.df)[1]
   len  <- dim(data.df)[1]
   sor <- sort(as.array(data.df[,nam]))
   jm1 <- seq(1,len)-1
   jm2 <- seq(1,len)-2
   
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

   # ESTIMATE L-MOMENTS
   l1 <- p00*b0
   l2 <- p10*b0+p11*b1
   l3 <- p20*b0+p21*b1+p22*b2

   # ESTIMATE PARAMETERS FROM L-MOMENTS  
   # (see Hosking 1990, Table 2)
   zz    <- 2/(3+l3/l2) - log(2)/log(3)
   k     <- 7.8590*zz + 2.9554*zz^2
   gk    <- gamma(1+k)
   alpha <- l2*k/((1-2^(-k))*gk)
   chi   <- l1 + alpha*(gk-1)/k

   # Conduct Kolmogorov-Smirnov Goodness-of-Fit Test
   # Note: This is a misuse of the KS-Test because the distribution paremeters
   # were determined from the data itself. The result may nonetheless pinpoint to
   #Êserious misfits.
   ks.pval <- suppressWarnings(
           ks.test(x=data.df[,nam],y="FGEV",alpha=alpha,chi=chi,k=k)$p.value)

   # DETERMINE EXTREME VALUES FOR GIVEN RETURN PERIODS
   fitted <- XGEV(foft(ret),alpha,chi,k)  
   
   # CONSTRUCT xval OBJECT FOR OUTPUT
   list(
       data=data.df,
       lamda=lamda,
       alpha=alpha,
       chi=chi,
       k=k,
       fit=cbind(ret,fitted),
       dist="GEV",
       ks.pval=ks.pval,
       estim="lmom"
   )
}

