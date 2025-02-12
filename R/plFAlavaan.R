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

# setMethod("show", "plFAlavaan", function(object) {
#   print_before(object)
#   print_stuff("Total time", paste0(round(object@timing$total, 2), " s"))
#   cat("\n")
#   print_stuff(" - Data reduction", paste0(round(object@external$plFA@freqTime, 2), " s"))
#   cat("\n")
#   print_stuff(" - Optimization", paste0(round(object@external$plFA@RTime, 2), " s"))
#   cat("\n")
#   print_stuff(" - Computing variance", paste0(round(object@timing$vcov, 2), " s"))
# })

setMethod("show", "plFAlavaan", print_before)

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

setMethod("plotTraj", "plFAlavaan", function(OBJ) {
  plotTraj(OBJ@external$plFA)
})

setMethod(
  f = "update",
  signature = signature(object = "plFAlavaan"),
  definition = function(object, ..., evaluate = TRUE) {
    # 1. Extract the original call from the object
    cl <- object@call

    # 2. Update the call with new arguments provided in '...'
    extra_args <- list(...)
    if (length(extra_args) > 0) {
      for(arg in names(extra_args)) {
        cl[[arg]] <- extra_args[[arg]]
      }
    }
    if (!("se" %in% names(extra_args))) {
      # If vars already computed, use it! Saves time
      cat("Skipping computation of vars!")
      extra_args$vars <- object@external$vars
    }

    # 3. Replace the function name in the call from lavaan::lavaan to plFA::cfa
    cl[[1]] <- quote(plFA::cfa)

    # 4. Evaluate the updated call (if requested) or return the call for inspection
    if (evaluate) {
      return(eval(cl, parent.frame()))
    } else {
      return(cl)
    }
  }
)
