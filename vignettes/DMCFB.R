## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=TRUE, message=FALSE-------------------------------------------
library(DMCFB)
fn <- list.files(system.file("extdata",package = "DMCFB"))
fn.f <- list.files(system.file("extdata",package="DMCFB"), full.names=TRUE)
OBJ <- readBismark(fn.f, fn)
cdOBJ <- DataFrame(Cell = factor(c("BC", "TC","Mono"),
levels = c("BC", "TC", "Mono")), row.names = c("BCU1568","BCU173","BCU551"))
colData(OBJ) <- cdOBJ
OBJ

## ---- eval=TRUE, message=FALSE-------------------------------------------
library(DMCFB)
set.seed(1980)
nr <- 1000; nc <- 8
metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
methl <- methc/metht
r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
names(r1) <- 1:nr
cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
OBJ2 <- cBSDMC(rowRanges=r1,methReads=methc,totalReads=metht,
methLevels=methl,colData=cd1)
OBJ2

## ---- eval=FALSE---------------------------------------------------------
#  library(DMCFB)
#  start.time <- Sys.time()
#  path0 <- "..//BCData/" # provide the path to the files
#  namelist.new <- list.files(path0,pattern="blk",full.names=F)
#  namelist.new.f <- list.files(path0,pattern="blk",full.names=T)
#  type <- NULL
#  for(i in seq_along(namelist.new)){
#      type[i] <- unlist(strsplit(namelist.new[i], split=c('_'), fixed=TRUE))[2]
#  }
#  type
#  table(type)
#  indTC <- which(type=="TC")
#  indBC <- which(type=="BC")
#  indMono <- which(type=="Mono")
#  namelist.new <- namelist.new[c(indBC,indMono,indTC)]
#  namelist.new.f <- namelist.new.f[c(indBC,indMono,indTC)]
#  BLKDat <- readBismark(namelist.new.f, namelist.new)
#  colData1 <- DataFrame(Group = factor(
#      c(rep("BC",length(indBC)), rep("Mono",length(indMono)),
#        rep("TC", length(indTC))), levels = c("BC", "Mono", "TC")),
#      row.names = colnames(BLKData))
#  colData(BLKDat) <- colData1
#  BLK.BC.Mono.TC <- sort(BLKDat)
#  DMC.obj = findDMCFB(object = BLKDat, bwa = 30, bwb = 30, nBurn = 300, nMC = 300,
#                      nThin = 1, alpha = 5e-5, pSize = 500, sfiles = FALSE)

## ---- fig.width=7--------------------------------------------------------
library(DMCFB)
set.seed(1980)
nr <- 1000; nc <- 8
metht <- matrix(as.integer(runif(nr * nc, 0, 100)), nr)
methc <- matrix(rbinom(n=nr*nc,c(metht),prob = runif(nr*nc)),nr,nc)
methl <- methc/metht
r1 <- GRanges(rep('chr1', nr), IRanges(1:nr, width=1), strand='*')
names(r1) <- 1:nr
cd1 <- DataFrame(Group=rep(c('G1','G2'),each=nc/2),row.names=LETTERS[1:nc])
OBJ1 <- cBSDMC(rowRanges=r1,methReads=methc,totalReads=metht,
methLevels=methl,colData=cd1)
OBJ2 = findDMCFB(object = OBJ1, bwa = 30, bwb = 30, nBurn = 10, nMC = 10,
                    nThin = 1, alpha = 0.05, pSize = 500, sfiles = FALSE)
plotDMCFB(OBJ2, region = c(1,400), nSplit = 2)

