.cBSDMC <- function(methReads, totalReads, methLevels, rowRanges,
                    colData = DataFrame(row.names = colnames(methReads)),
                    metadata = list(), ...) {
    new("BSDMC", SummarizedExperiment(
    assays = SimpleList(methReads = methReads, totalReads = totalReads,
    methLevels = methLevels),
    rowRanges = rowRanges, colData = colData, metadata = list()))
    }

#' @rdname cBSDMC-method
#' @aliases cBSDMC-method cBSDMC
setMethod("cBSDMC", signature(methReads = "matrix", totalReads = "matrix",
    methLevels = "matrix", rowRanges = "GRanges"), .cBSDMC)
