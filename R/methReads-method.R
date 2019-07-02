.methReads <- function(object) {
    return(assays(object)$methReads)
}

.replace.methReads <- function(object, value) {
    assays(object)$methReads <- value
    return(object)
}

#' @rdname methReads-method
#' @aliases methReads-method methReads
setMethod("methReads", signature(object = "BSDMC"), .methReads)

#' @rdname methReads-method
#' @aliases methReads-method methReads<-
setReplaceMethod("methReads", signature(object = "BSDMC", value = "matrix"),
    .replace.methReads)

