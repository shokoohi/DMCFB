.methLevels <- function(object) {
    return(assays(object)$methLevels)
}

.replace.methLevels <- function(object, value) {
    assays(object)$methLevels <- value
    return(object)
}

#' @rdname methLevels-method
#' @aliases methLevels-method methLevels
setMethod("methLevels", signature(object = "BSDMC"), .methLevels)

#' @rdname methLevels-method
#' @aliases methLevels-method methLevels<-
setReplaceMethod("methLevels", signature(object = "BSDMC", value = "matrix"),
    .replace.methLevels)
