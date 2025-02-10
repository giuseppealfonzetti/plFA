#' list with constraints
setClass('Constraints', slots = c(loadings = 'matrix', corrflag = 'logical', stdlv="logical", loglatsd="vector", llc="ANY"))

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

#' List of convergence diagnostics
setClass('Convergence', slots = c(convergence_full  = 'numeric',
                                  convergence_burn  = 'logical',
                                  proj_after_burn   = 'logical',
                                  finite_gr         = 'logical',
                                  neg_pdiff         = 'logical'))

#' Quantities specific for stochastic estimation
setClass('StoFit', slots = c( path         = 'list',
                              nll          = "numeric",
                              control      = 'list',
                              last_iter    = 'numeric',
                              burnt        = 'numeric',
                              convergence  = 'Convergence'))


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
  cores    = 'numeric',
  version  = 'character'))

