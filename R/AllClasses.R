#' @title params
#' @description parameters name and their descriptions
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name params
#' @param methReads The matrix \code{methReads} contains the number of
#' methylated reads spanning a CpG-site. The rows represent the CpG sites in
#' \code{rowRanges} and the columns represent the samples in \code{colData}.
#' @param totalReads The matrix \code{totalReads} contains the number of reads
#' spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
#' and the columns represent the samples in \code{colData}.
#' @param methLevels The matrix \code{methLevels} contains the predicted
#' methylation level spanning a CpG-site using Bayesian functional regression
#' models. The rows represent the CpG sites in \code{rowRanges} and the columns
#' represent the samples in \code{colData}.
#' @param rowRanges A \code{\link{GRanges}} or \code{\link{GRangesList}}
#' object describing the ranges of interest. Names, if present, become the row
#' names of the \code{\link{SummarizedExperiment}} object. The length of the
#' \code{\link{GRanges}} or \code{\link{GRangesList}} must equal the number of
#' rows of the matrices in \code{assays}. If \code{rowRanges} is missing, a
#' \code{\link{SummarizedExperiment}} instance is returned.
#' @param colData Object of class \code{"DataFrame"} containing information on
#' variable values of the samples
#' @param metadata A \code{list} of storing MCMC samples or DMCs
#' @param object A \code{\link{BSDMC-class}}
#' object
#' @param value An integer matrix
#' @param name A character list
#' @param obj1 A \code{\link{BSDMC-class}}
#' @param obj2 A \code{\link{BSDMC-class}}
#' @param files A character list
#' @param file A character
#' @param nCores An integer value specifying the number of machine cores for
#' parallel computing
#' @param pSize An integer value specifying the number of cytosines in a regrion
#' to be used in a Bayesian functiona regression model for DMC detection
#' @param bwa An integer value specifying the band-width size of B-spline basis
#' matrix for a natural cubic spline for the group-specific effects of the
#' Bayesian functional regression model
#' @param bwb An integer value specifying the band-width size of B-spline basis
#' matrix for a natural cubic spline for the individual-specific effects of the
#' Bayesian functional regression model
#' @param nBurn An integer value specifying the number of burn-in samples
#' @param nThin An integer value specifying the thining number in MCMC
#' @param nMC An integer value specifying the number of MCMC samples after
#' burn-in
#' @param sdv An double value specifying the standard deviation of priors
#' @param alpha A numeric value specifying the level of \eqn{\alpha} in credible
#' interval \eqn{(1-\alpha)\%}
#' @param col A character vector indicating which colors to alternate.
#' @param sfiles A logical value indicating whether files to be saved or not.
#' @param region An integer vector of length two specifying which subset of the
#' object to be plotted
#' @param nSplit A integer value specifying the number of subsets must be done
#' for plotting the results of DMC identification
#' @param parList A list specifying plots parameters, see \code{\link{par}}
#' @param ... other possible parameters
#' @docType NULL
NULL

#' @title BSDMC object
#' @description The \code{BSDMC} object is an S4 class that represents
#' differentially methylated CpG sites (DMCs) in BS-Seq Data.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name BSDMC-class
#' @import methods
#' @import SummarizedExperiment
#' @slot methReads An integer matrix
#' @slot totalReads An integer matrix
#' @slot methLevels A numeric matrix
#' @param methReads The matrix \code{methReads} contains the number of
#' methylated reads spanning a CpG-site. The rows represent the CpG sites in
#' \code{rowRanges} and the columns represent the samples in \code{colData}.
#' @param totalReads The matrix \code{totalReads} contains the number of reads
#' spanning a CpG-site. The rows represent the CpG sites in \code{rowRanges}
#' and the columns represent the samples in \code{colData}.
#' @param methLevels The matrix \code{methLevels} contains the predicted
#' methylation level spanning a CpG-site using Bayesian functional regression
#' models. The rows represent the CpG sites in \code{rowRanges} and the columns
#' represent the samples in \code{colData}.
#' @docType class
#' @keywords object
#' @return A \code{\link{BSDMC-class}} object
#' @seealso \code{\link{RangedSummarizedExperiment-class}}
#' \code{\link{GRanges-class}}
#' @examples
#' nr <- 500; nc <- 16
#' metht <- matrix(as.integer(runif(nr * nc, 0, nr)), nr)
#' methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
#' meths <- matrix(as.integer(runif(nr * nc, 0, 10)), nr)
#' methl <- methc/metht
#' methv <- matrix((runif(nr * nc, 0.1, 0.5)), nr)
#' r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
#' names(r1) <- 1:nr
#' cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
#' OBJ2 <- cBSDMC(rowRanges=r1,methReads=methc,totalReads=metht,
#' methLevels=methl,methStates=meths,methVars=methv,colData=cd1)
#' OBJ2
#' @exportClass BSDMC
BSDMC <- setClass("BSDMC", representation(
    methReads = "matrix", totalReads = "matrix", methLevels = "matrix"
    ), prototype (methReads = matrix(), totalReads = matrix(), methLevels =
    matrix()), contains = "RangedSummarizedExperiment")

setValidity("BSDMC", function(object){
    if(length(assays(object)) != 3)
        return("The assays slot in BSDMC object must be of length three.")
    if(!(all( is.element(names(assays(object)), c(
        "methReads", "totalReads", "methLevels")))))
        return("The assays slot in BSDMC object must contain totalReads,
        methReads and methStates.")
    if(!all( vapply(assays(object), class, character(1)) == "matrix" ))
        return("The methReads, totalReads and methLevels slots
        of a BSDMC object must be matrices.")
    if(!all(vapply(assays(object), typeof, character(1)) ==
            c("integer", "integer", "double")))
        return("The methReads, totalReads and methLevels slots
        of a BSDMC object must be integer, integer and double,
        repectively.")}
        )
