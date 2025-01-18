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

setMethod("show", "plFAlavaan", function(object) {

  class(object) <- "lavaan"
  cat("plFA x ")
  callNextMethod()

})

setMethod("summary", "plFAlavaan", function(object, ...) {

  class(object) <- "lavaan"
  cat("plFA x ")
  callNextMethod()

})

setMethod("coef", "plFAlavaan", function(object, ...) {
  class(object) <- "lavaan"
  callNextMethod()
})
