##########################################
# Spectra modification functions used in #
#          workflow function             #
##########################################

#' Remove fragments below x% of base peak intensity
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (int_tresh / 100 )
}

#' Normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

#' Remove precursor ion
removePrecursor <- function(window = 1) {
  
  .between <- function(x, a, b) {
    (x >= a) & (x <= b)
  }  
  
  function(x, precursorMz, ...) {
    x[!.between(x[,1], precursorMz - 0.5*window, precursorMz+0.5*window),,drop=FALSE]
  }
    
}
