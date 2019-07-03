.combine.BSDMC <- function(obj1, obj2) {
    if (any(is.element(colnames(obj1), colnames(obj2)))) {
    stop("The BSDMC objects to combine should not have samples in common!")
    }
    colData.new <- rbind(colData(obj1), colData(obj2))
    rowRanges.new <- sort(unique(c(rowRanges(obj1), rowRanges(obj2))))
    ind.match.obj1 <- findOverlaps(rowRanges.new, rowRanges(obj1),
    select = "first"
    )
    ind.match.obj2 <- findOverlaps(rowRanges.new, rowRanges(obj2),
    select = "first")
    nr <- length(rowRanges.new)
    nc <- nrow(colData.new)
    methReads.new <- matrix(integer(length = nr * nc),
    ncol = nc, nrow = nr,
    dimnames = list(names(rowRanges.new), rownames(colData.new)))
    methReads.new[, ] <- cbind(
    methReads(obj1)[ind.match.obj1, ],
    methReads(obj2)[ind.match.obj2, ]
    )
    methReads.new[is.na(methReads.new)] <- 0L
    totalReads.new <- matrix(integer(length = nr * nc),
    ncol = nc, nrow = nr,
    dimnames = list(names(rowRanges.new), rownames(colData.new))
    )
    totalReads.new[, ] <- cbind(
    totalReads(obj1)[ind.match.obj1, ],
    totalReads(obj2)[ind.match.obj2, ]
    )
    totalReads.new[is.na(totalReads.new)] <- 0L
    methLevels.new <- matrix(
    double(length = nr * nc),
    ncol = nc, nrow = nr, dimnames = list(names(rowRanges.new
    ), rownames(colData.new))
    )
    methLevels.new[, ] <- cbind(
    methLevels(obj1)[ind.match.obj1, ],
    methLevels(obj2)[ind.match.obj2, ]
    )
    methLevels.new[is.na(totalReads.new)] <- 0.0

    z <- cBSDMC(
    colData = colData.new, rowRanges = rowRanges.new,
    methReads = methReads.new, totalReads = totalReads.new,
    methLevels = methLevels.new
    )
    return(z)
    }

#' @rdname combine-method
#' @aliases combine-method combine
setMethod("combine", signature(obj1 = "BSDMC", obj2 = "BSDMC"), .combine.BSDMC)
