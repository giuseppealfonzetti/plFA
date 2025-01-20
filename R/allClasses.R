#' list with constraints
setClass('Constraints', slots = c(loadings = 'matrix', corrflag = 'logical', stdlv="logical", loglatsd="vector"))

#' List collecting problem dimensions
setClass('Dimensions', slots = c(n     = 'numeric',
                                 p     = 'numeric',
                                 q     = 'numeric',
                                 cat   = 'vector',
                                 pairs = 'numeric',
                                 nthr  = 'numeric',
                                 nload = 'numeric',
                                 ncorr = 'numeric',
                                 nvar  = 'numeric',
                                 npar  = 'numeric'))

#' Quantities specific for stochastic estimation
setClass('StoFit', slots = c( trajSubset  = 'vector',
                              projIters   = 'vector',
                              pathTheta   = 'matrix',
                              pathAvTheta = 'matrix',
                              pathGrad    = 'matrix',
                              pathValNll  = 'matrix',
                              control     = 'list',
                              lastIter    = 'numeric',
                              convergence = 'numeric',
                              cppTime     = 'ANY'))


# Object returned by fit_plFA
setClass('PlFaFit', slots = c(
  freq     = 'matrix',
  valfreq  = 'matrix',
  freqTime = 'numeric',
  cnstr    = 'Constraints',
  dims     = 'Dimensions',
  method   = 'character',
  init     = 'vector',
  numFit   = 'ANY',
  stoFit   = 'StoFit',
  theta    = 'vector',
  RTime    = 'numeric',
  cores    = 'numeric'))

