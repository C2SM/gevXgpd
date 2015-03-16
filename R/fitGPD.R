"fitGPD" <-
function(..., dist="GPD",estim="mlik") {
# ==================================================================
#  Estimates parameters of GPD-fit to sample data.
#  The fit can be restricted to the family of EXPONENTIAL Distributions
#  when dist="EXPON".
#  Options for estimation are Maximum Likelihood (estim="mlik")
#  and L-moments (estim="lmom").
#  For definition of other parameters see sub functions.
#  

if ((dist != "GPD") & (dist != "EXPON")) {
    stop("** ERROR ** dist must be either EXPON or GPD")}

if ((estim != "mlik") & (estim != "lmom")) {
    stop("** ERROR ** estim must be either lmom or mlik")}

switch(dist,
    "EXPON" = switch(estim,
	"lmom" = fitEXPON.lmom(...),
	"mlik" = fitEXPON.mlik(...)),
    "GPD" = switch(estim,
	"lmom" = fitGPD.lmom(...),
	"mlik" = fitGPD.mlik(...)))
}

