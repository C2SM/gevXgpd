"min.sep" <-
function(dat,threshold,min.sep=1) {
      #preps
      x <- dat[,1]
      nn <- length(x)
      mn <- min(x,na.rm=TRUE)
      x[is.na(x)] <- mn
      fill.value <- mn-1

      # checks
      if (threshold <= mn) {
        stop("Threshold smaller than smallest value of dataset.")
      }

      "peak.indices" <- 
      function(x,threshold,min.sep,fill.value) {
         # get indices of peaks (maxima) which are at least min.sep
         #units apart.

         # build matrix with shifted arrays in columns
         nn <- length(x)
         width <- 2*(min.sep-1)+1
         xmat <- matrix(fill.value,nrow=nn+width-1,ncol=width)
         for (k in 1:width) {
            xmat[(width-k)+(1:nn),k] <- x
         }
         # drop leading and ending rows
         ii <- min.sep:(nn+min.sep-1)
         xmat <- xmat[ii,]
         # determine max in matrix rows
         row.max <- apply(xmat,FUN=max,MARGIN=c(1),na.rm=TRUE)

         # select peaks exceeding threshold
         ind <- ((row.max == xmat[,min.sep]) & (row.max >= threshold))
         if (sum(ind) == 0) {
           ind <- c()
         } else {
           ind <- which(ind)
         }

         #clean up
         rm(xmat,row.max)
 
         return(ind)
      }

      #first go (indices of peaks that are at least min.sep apart
      ind1 <- peak.indices(x,threshold=threshold,
                           min.sep=min.sep,fill.value=fill.value)
      if (length(ind1)==0) {
        return(dat[ind1,])
      }

      # mask all values within min.sep-1 around the peak
      if (min.sep > 1) {
         ind.mask <- c()
         for (dd in 1:(min.sep-1)) {
            ind.mask <- c(ind.mask,ind1-dd,ind1+dd)
         }
         ind.mask <- ind.mask[!((ind.mask > length(x)) | (ind.mask < 1))]
         x[ind.mask] <- fill.value
      }

      # second go (indices of peaks with possibly hidden secondary peaks)
      ind <- peak.indices(x,threshold=threshold,
                          min.sep=min.sep,fill.value=fill.value)

      # clean up
      rm(ind1)

      # return peaks
      return(dat[ind,])
}

