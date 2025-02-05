#' plFAlavaan Class
#'
#' This is a class that extends the lavaan class.
#'
#' @importFrom lavaan lavaan
#' @export
setClass(
  Class = "plFAlavaan",
  contains = "lavaan"
)

print_before <- function(object) {
  class(object) <- "lavaan"
  cat("plFA", object@external$plFA@version, "\n")
  cat("  \u2A09 ")
  callNextMethod()
}

print_stuff <- function(what, res) {
  x <-
    stringr::str_pad(
      res, 54 - nchar(paste0("  ", what, " ")), side = "left"
    )
  x <- sprintf(paste0("  ", what," %s"), x)
  cat(x)
}

setMethod("show", "plFAlavaan", function(object) {
  print_before(object)
  print_stuff("Total time", paste0(round(object@timing$total, 2), " s"))
  cat("\n")
  print_stuff(" - Data reduction", paste0(round(object@external$plFA@freqTime, 2), " s"))
  cat("\n")
  print_stuff(" - Optimization", paste0(round(object@external$plFA@RTime, 2), " s"))
  cat("\n")
  print_stuff(" - Computing variance", paste0(round(object@timing$vcov, 2), " s"))
})

setMethod("summary", "plFAlavaan", print_before)

setMethod("coef", "plFAlavaan", function(object, ...) {
  class(object) <- "lavaan"
  callNextMethod()
})

setMethod("logLik", "plFAlavaan", function(object, ...) {
  class(object) <- "lavaan"
  object@Options$estimator <- "ML"
  callNextMethod()
})

setMethod("AIC", "plFAlavaan", function(object, ...) {
  object@loglik$AIC
})

setMethod("BIC", "plFAlavaan", function(object, ...) {
  object@loglik$BIC
})
