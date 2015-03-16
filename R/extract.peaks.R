"extract.peaks" <-
function(dat,threshold,method="min.sep",...) {
       FUN <- match.fun(method)
       FUN(dat,threshold,...)
}

