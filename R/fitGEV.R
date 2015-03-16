"fitGEV" <-
function(..., dist="GEV",estim="mlik") {
# ==================================================================
#  Estimates parameters of GEV-fit to sample data.
#  The fit can be restricted to the family of GUMBEL Distributions
#  when dist="GUMBEL".
#  Options for estimation are Maximum Likelihood (estim="mlik")
#  and L-moments (estim="lmom").
#  This is only the wrapper to distribute between different estimation
#  techniques.

if ((dist != "GEV") & (dist != "GUMBEL")) {
    stop("** ERROR ** dist must be either GUMBEL or GEV")}

if ((estim != "mlik") & (estim != "pmlik") & (estim != "lmom")) {
    stop("** ERROR ** estim must be either lmom or mlik")}

switch(dist,
    "GUMBEL" = switch(estim,
	"lmom" = fitGUMBEL.lmom(...),
	"mlik" = fitGUMBEL.mlik(...)),
    "GEV" = switch(estim,
	"lmom" = fitGEV.lmom(...),
	"mlik" = fitGEV.mlik(...),
	"pmlik" = fitGEV.pmlik(...)))
}

