#' list with constraints
setClass('Constraints', slots = c(loadings = 'matrix', corrflag = 'numeric'))

#' List collecting problem dimensions
setClass('Dimensions', slots = c(n = 'numeric', p = 'numeric', q = 'numeric', cat = 'vector', pairs = 'numeric'))

#' Quantities specific for stochastic estimation
setClass('StoFit', slots = c(
  trajSubset = 'vector',
  projIters = 'vector',
  pathTheta = 'matrix',
  pathAvTheta = 'matrix',
  pathGrad = 'matrix',
  control = 'list',
  cppTime = 'ANY'
))
# Object returned by fit_plFA
setClass('PlFaFit', slots = c(
  cnstr = 'Constraints',
  dims = 'Dimensions',
  method = 'character',
  init = 'vector',
  numFit = 'ANY',
  stoFit = 'StoFit',
  theta = 'vector',
  RTime = 'numeric'))

