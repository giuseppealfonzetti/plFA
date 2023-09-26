#' @export
setGeneric('getThetaPath', function(x, ..., verbose = TRUE) standardGeneric('getThetaPath'), signature = 'x')

#' @export
setMethod('getThetaPath', 'StoFit', function(x) x@pathTheta)

#' @export
setMethod('getThetaPath', 'PlFaFit', function(x) getThetaPath(x@stoFit))
