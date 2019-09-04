#' @title cBSDMC method
#' @description Creates a \code{\link{BSDMC-class}} object
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name cBSDMC-method
#' @import GenomicRanges
#' @import S4Vectors
#' @inheritParams params
#' @details The rows of a \code{BSDMC} object represent ranges (in genomic
#' coordinates) of interest. The ranges of interest are described by a
#' \code{GRanges} or a \code{GRangesList} object,
#' accessible using the \code{rowRanges} function.
#' The \code{GRanges} and \code{GRangesList} classes
#' contains sequence (e.g., chromosome) name, genomic coordinates, and strand
#' information. Each range can be annotated with additional data; this data
#' might be used to describe the range or to
#' summarize results (e.g., statistics of differential abundance) relevant to
#' the range. Rows may or may not have row names; they often will not.
#' @return A \code{\link{BSDMC-class}}
#' @examples
#' set.seed(1980)
#' nr <- 150
#' nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc / metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ2 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, methStates = meths, methVars = methv, colData = cd1
#' )
#' OBJ2
#' @exportMethod cBSDMC
setGeneric("cBSDMC", function(
    methReads, totalReads, methLevels, rowRanges,
    colData = DataFrame(row.names = colnames(methReads) ),
    metadata = list(), ...) standardGeneric("cBSDMC"),
    signature = c("methReads", "totalReads", "methLevels", "rowRanges"))

#' @title methReads method
#' @description Returns \code{methReads} stored in \code{\link{BSDMC-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @import S4Vectors
#' @name methReads-method
#' @inheritParams params
#' @return A matrix
#' @examples
#' nr <- 150
#' nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, colData = cd1
#' )
#' methReads(OBJ1)
#' @exportMethod methReads
setGeneric("methReads", function(object) standardGeneric("methReads"))

#' @title methReads method
#' @description Assigns \code{methReads} to \code{\link{BSDMC-class}}
#' @name methReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMC-class}} object
#' @examples
#' methReads(OBJ1) <- methc
#' @exportMethod methReads<-
setGeneric("methReads<-", function(object,value) standardGeneric("methReads<-"))

#' @title totalReads method
#' @description Returns \code{totalReads} stored in \code{\link{BSDMC-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' nr <- 150
#' nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, colData = cd1
#' )
#' totalReads(OBJ1)
#' @exportMethod totalReads
setGeneric("totalReads", function(object) standardGeneric("totalReads"))

#' @title totalReads method
#' @description Assigns \code{totalReads} to \code{\link{BSDMC-class}}
#' @name totalReads-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMC-class}} object
#' @examples
#' totalReads(OBJ1) <- metht
#' @exportMethod totalReads<-
setGeneric("totalReads<-", function(object, value)
    standardGeneric("totalReads<-"))

#' @title methLevels method
#' @description Returns \code{methLevels} stored in \code{\link{BSDMC-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name methLevels-method
#' @import S4Vectors
#' @inheritParams params
#' @return A matrix
#' @examples
#' nr <- 150
#' nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, colData = cd1
#' )
#' methLevels(OBJ1)
#' @exportMethod methLevels
setGeneric("methLevels", function(object) standardGeneric("methLevels"))

#' @title methLevels method
#' @description Assigns \code{methLevels} to \code{\link{BSDMC-class}}
#' @name methLevels-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMC-class}} object
#' @examples
#' methLevels(OBJ1) <- methl
#' @exportMethod methLevels<-
setGeneric("methLevels<-", function(object, value)
    standardGeneric("methLevels<-"))

#' @title combine method
#' @description combine two \code{\link{BSDMC-class}} or
#' two \code{\link{BSDMC-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name combine-method
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMC-class}} or \code{\link{BSDMC-class}}
#' @examples
#' set.seed(1980)
#' nr <- 150
#' nc <- 8
#' metht <- matrix(as.integer(runif(nr * nc * 2, 0, nr)), nr)
#' methc <- matrix(
#'   rbinom(n = nr * nc, c(metht), prob = runif(nr * nc * 2)),
#'   nr, nc * 2
#' )
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group = rep("G1", each = nc), row.names = LETTERS[1:nc])
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc[, 1:nc], totalReads = metht[, 1:nc],
#'   methLevels = methl[, 1:nc], colData = cd1
#' )
#' cd2 <- DataFrame(
#'   Group = rep("G2", each = nc),
#'   row.names = LETTERS[nc + 1:nc]
#' )
#' OBJ2 <- cBSDMC(
#'   rowRanges = r1, methReads = methc[, nc + 1:nc], totalReads =
#'     metht[, nc + 1:nc], methLevels = methl[, nc + 1:nc], colData = cd2
#' )
#' OBJ3 <- combine(OBJ1, OBJ2)
#' OBJ3
#' @exportMethod combine
setGeneric("combine", function(obj1, obj2) standardGeneric("combine"))

#' @title readBismark method
#' @description reads BS-Seq data
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name readBismark-method
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @inheritParams params
#' @return A \code{\link{BSDMC-class}} object
#' @examples
#' fn <- list.files(system.file("extdata", package = "DMCHMM"))
#' fn.f <- list.files(system.file("extdata", package = "DMCHMM"),
#'   full.names = TRUE
#' )
#' OBJ <- readBismark(fn.f, fn)
#' cdOBJ <- DataFrame(Cell = factor(c("BC", "TC", "Mono"),
#'   labels = c("BC", "TC", "Mono")
#' ), row.names = c("BCU1568", "BCU173", "BCU551"))
#' colData(OBJ) <- cdOBJ
#' OBJ
#' @exportMethod readBismark
setGeneric("readBismark", function(files, colData)
    standardGeneric("readBismark"))

#' @title findDMCFB method
#' @description DMC identification via Bayesian functional regression models
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name findDMCFB-method
#' @inheritParams params
#' @return \code{\link{BSDMC-class}} object
#' @import BiocParallel
#' @import GenomicRanges
#' @import S4Vectors
#' @import MASS
#' @importFrom stats as.formula binomial coefficients coef contrasts glm logLik
#' @importFrom benchmarkme get_ram
#' @importFrom speedglm speedglm
#' @importFrom data.table is.data.table
#' @importFrom utils combn
#' @importFrom splines ns
#' @importFrom tibble is_tibble as_tibble
#' @importFrom arm bayesglm
#' @importFrom matrixStats rowVars
#' @importFrom fastDummies dummy_cols
#' @importFrom matrixStats colQuantiles
#' @importFrom utils memory.size
#' @examples
#' set.seed(1980)
#' nr <- 1000
#' nc <- 4
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, colData = cd1
#' )
#' OBJ2 <- findDMCFB(OBJ1,
#'   bwa = 10, bwb = 10, nBurn = 50, nMC = 50, nThin = 1,
#'   alpha = 0.05, nCores = 2, pSize = 500, sfiles = FALSE
#' )
#' OBJ2
#' @exportMethod findDMCFB
setGeneric("findDMCFB", function(object, bwa, bwb, nBurn, nMC, nThin, alpha,
    sdv, nCores, pSize, sfiles)
    standardGeneric("findDMCFB"))

#' @title plotDMCFB method
#' @description Plotting the results of DMC identifation stored in
#' a \code{\link{BSDMC-class}} object
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name plotDMCFB-method
#' @inheritParams params
#' @import BiocParallel
#' @import GenomicRanges
#' @import S4Vectors
#' @import grDevices
#' @import rtracklayer
#' @import graphics
#' @importFrom graphics abline axis layout legend lines par points segments
#' @return Plot
#' @examples
#' set.seed(1980)
#' nr <- 1000
#' nc <- 4
#' metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
#' methc <- matrix(rbinom(n = nr * nc, c(metht), prob = runif(nr * nc)), nr, nc)
#' methl <- methc / metht
#' r1 <- GRanges(rep("chr1", nr), IRanges(1:nr, width = 1), strand = "*")
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(
#'   Group = rep(c("G1", "G2"), each = nc / 2),
#'   row.names = LETTERS[1:nc]
#' )
#' OBJ1 <- cBSDMC(
#'   rowRanges = r1, methReads = methc, totalReads = metht,
#'   methLevels = methl, colData = cd1
#' )
#' OBJ2 <- findDMCFB(OBJ1,
#'   bwa = 10, bwb = 10, nBurn = 50, nMC = 50, nThin = 1,
#'   alpha = 0.05, nCores = 2, pSize = 500, sfiles = FALSE
#' )
#' plotDMCFB(OBJ2)
#' @exportMethod plotDMCFB
setGeneric("plotDMCFB", function(object, region, nSplit,parList)
    standardGeneric("plotDMCFB"))
