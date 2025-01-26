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
setClass('StoFit', slots = c( path_iters   = 'vector',
                              path_theta   = 'list',
                              path_avtheta = 'list',
                              path_nll     = 'vector',
                              nll          = "numeric",
                              control      = 'list',
                              last_iter     = 'numeric',
                              convergence  = 'logical',
                              cppTime      = 'ANY'))


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

