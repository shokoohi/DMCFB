.totalReads <- function(object) {
    return(assays(object)$totalReads)
}

.replace.totalReads <- function(object, value) {
    assays(object)$totalReads <- value
    return(object)
}

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads
setMethod("totalReads", signature(object = "BSDMC"), .totalReads)

#' @rdname totalReads-method
#' @aliases totalReads-method totalReads<-
setReplaceMethod("totalReads", signature(object = "BSDMC", value = "matrix"),
    .replace.totalReads)
