.plotDMCFB <- function(object, region, nSplit, parList){
    if(missing(object) | length(rowData(object))<6){
        stop("Provide a DMCFB object with DMCs in 'rowData'.")
    }
    if(!missing(region)){
        if(length(region)!=2){
            stop("Provide the start and the end of region to be plotted.")
        }else{
            object = object[(region[1]):(region[2]),]
        }
    }

    if(missing(nSplit)) nSplit = 1
    if(missing(parList)) parList = list()

    pars = list(mfrow = min(3, (dim(rowData(object))[2]-1)/5), type = "p",
                main = NULL, xlab = "Position", ylab = "Methylation",
                col = c("blue", "red"), ticks = NULL,  nticks=7, pch=c(4,20),
                mgp = c(2,1,0), oma = c(0,0,0,0) + 0, mar = c(3,3,1,3) + 0,
                cex = 0.5, cex.axis = 0.9, cex.lab = 0.9, cex.main = 0.9,
                lty = 1, lwd = 1, legend = TRUE, ablines = NULL,
                legend.adj=c(1,0.5))


    pars[intersect(names(pars),names(parList))] <- parList[
        intersect(names(pars),names(parList))]

    if(!is.numeric(pars$mfrow)){
        stop("Provide an integer greater than 0 for 'mfrow'.")
    }
    if(!any(pars$type%in%c("l","p"))){
        stop("Specify 'l' for lines or 'p' for points in 'type'.")
    }
    if(length(pars$col)==1){
        pars$col[2] = pars$col[1]
    }
    if(length(pars$pch)==1){
        pars$pch[2] = pars$pch[1]
    }

    lcont = (dim(rowData(object))[2]-1)/5
    ContranstNames = names(rowData(object))[2:(lcont+1)]

    sidelevels = NULL
    for(i in seq_along(ContranstNames)){
        sidelevels = cbind(sidelevels, c(unlist(strsplit(unlist(strsplit(
            ContranstNames[i], split=c('vs'), fixed=TRUE))[1], split=c(
                'Cs'), fixed=TRUE))[2], unlist(strsplit(ContranstNames[
                    i], split=c('vs'), fixed=TRUE))[2]))
    }

    if(is.null(pars$main)){
        pars$main = apply(sidelevels, 2, function(x){paste(x[1],"vs",x[2])})
    }else if(length(pars$main)!=lcont){
        stop(paste0("Provide ",lcont, " names for the titles of pair-wise
                    comparision plots."))
    }

    Pos.list = split(seq_along(object), ceiling(seq_along(object)/length(
        object)*nSplit))
    aa = length(Pos.list)
    if(length(Pos.list[[aa]])<7){
        Pos.list[[aa-1]] = c(Pos.list[[aa-1]] , Pos.list[[aa]])
        Pos.list = Pos.list[-aa]
    }

    for (il in seq_along(Pos.list)) {

        object.sub = object[unlist(Pos.list[il]),]

        xmax = ceiling(max(start(object.sub)+1-min(start(object.sub))) * 1.03)
        xmin = floor(max(start(object.sub)+1-min(start(object.sub))) * -0.03)

        def.arg <- list()

        if(lcont>1){
            for (i in seq_len(lcont)) {
                if(i%%pars$mfrow==0){
                    def.arg[[i]] = list(xaxt='n', yaxt='n', bty='n', xaxs='i',
                                        yaxs='i',las=1,pch=20,xlim=c(xmin,xmax),
                                        ylim=c(-0.1,1.1), xlab=pars$xlab,
                                        ylab=pars$ylab,mgp = pars$mgp, main =
                                            pars$main[i],cex.lab=pars$cex.lab,
                                        cex.main = pars$cex.main)
                }else{
                    def.arg[[i]]<-list(
                        xaxt='n', yaxt='n', bty='n', xaxs='i',yaxs='i', las=1,
                        pch=20,xlim=c(xmin,xmax), ylim=c(-0.1,1.1),xlab="",
                        ylab=pars$ylab, mgp=pars$mgp,main= pars$main[i],
                        cex.lab=pars$cex.lab,cex.main = pars$cex.main)
                }
            }
            NameBP <- start(object.sub) - min(start(object.sub)) + 1
            if(is.null(pars$nticks)){
                ticks = as.integer(seq(1, max(start(object.sub))-min(start(
                    object.sub))+1,length.out = 7))
                pars$nticks = 7
            }else{
                ticks = as.integer(seq(1, max(start(object.sub))-min(start(
                    object.sub))+1,length.out = pars$nticks))
            }
            nlabs=c()
            for(kk in seq_along(ticks)){
                nlabs[kk] = sum(NameBP<= ticks[kk])
            }
            nlabs = (start(object.sub))[nlabs]
        }else{
            def.arg[[1]]<-list(
                xaxt='n', yaxt='n', bty='n', xaxs='i',yaxs='i',las=1,pch=20,
                xlim=c(xmin,xmax), ylim=c(-0.1,1.1),xlab=pars$xlab,
                ylab=pars$ylab, mgp = pars$mgp, main=pars$main[1],
                cex.lab=pars$cex.lab,cex.main = pars$cex.main)
            NameBP <- start(object.sub) - min(start(object.sub)) + 1
            if(is.null(pars$nticks)){
                ticks = as.integer(seq(1, max(start(object.sub))-min(start(
                    object.sub))+1,length.out = 7))
                pars$nticks = 7
            }else{
                ticks = as.integer(seq(1, max(start(object.sub))-min(start(
                    object.sub))+1,length.out = pars$nticks))
            }
            nlabs=c()
            for(kk in seq_along(ticks)){
                nlabs[kk] = sum(NameBP<= ticks[kk])
            }
            nlabs = (start(object.sub))[nlabs]
        }

        old.par <- par(no.readonly = TRUE)
        if(length(ContranstNames)>1){
            nf <- layout( matrix(seq_len(pars$mfrow), pars$mfrow,1, byrow=TRUE))
        }else{
            nf <- layout( matrix(c(1), 1, 1, byrow = TRUE ))
        }
        par( oma = pars$oma,  mar = pars$mar)
        if(lcont>1){
            for (i in seq_len(lcont)) {

                d2 <- data.frame(dir=rep(0,length(object.sub)))
                ind.o = rowData(object.sub)[,1+lcont+i]=="Hypo"
                ind.e = rowData(object.sub)[,1+lcont+i]=="Equal"
                ind.r = rowData(object.sub)[,1+lcont+i]=="Hyper"
                d2$NameBP <- start(object.sub) - min(start(object.sub)) + 1
                d2$col = NA
                d2$col[ind.o]="black"
                d2$col[ind.e]="black"
                d2$col[ind.r]="black"
                d2$dmcH = NA
                d2$dmcH[ind.r]= 1.08
                d2$dmcH[ind.o]= NA
                d2$dmcH[ind.e]= NA
                d2$dmcH[rowData(object.sub)[,1+i]==0] = NA
                d2$dmcL = NA
                d2$dmcL[ind.r]= NA
                d2$dmcL[ind.e]= NA
                d2$dmcL[ind.o]= -0.08
                d2$dmcL[rowData(object.sub)[,1+i]==0] = NA

                do.call("plot",c(NA,def.arg[[i]]))
                if(lcont>1){
                    if(i%%pars$mfrow==0){
                        axis(1, at=ticks, labels=nlabs)
                    }else{
                        axis(1, at=ticks, labels=NA)
                    }
                }else{
                    axis(1, at=ticks, labels=nlabs)
                }
                axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(
                    0,0.2,0.4,0.6,0.8,1),cex.axis=pars$cex.axis)
                axis(4, at=c(-0.05), labels=paste0(sidelevels[
                    1,i]," < ",sidelevels[2,i]), las=1, cex.axis=pars$cex.axis,
                    line = -2)
                axis(4, at=c(1.05), labels=paste0(sidelevels[
                    1,i]," > ", sidelevels[2,i]), las=1, cex.axis=pars$cex.axis,
                    line = -2)

                with(d2, segments(x0=NameBP, y0=rep(1.02,length(
                    d2$dmcH)),x1=NameBP, y1 = d2$dmcH, lwd=pars$lwd, col = col))
                with(d2, segments(x0=NameBP, y0=rep(-0.02,length(
                    d2$dmcL)), x1=NameBP, y1 = d2$dmcL, lwd=pars$lwd, col=col))
                if(pars$type=="l"){
                    lines(d2$NameBP, rowMeans(methLevels(object.sub)[
                        ,colData(object.sub)[[1]]== sidelevels[
                            1,i]]), col=pars$col[1],cex=pars$cex)
                    lines(d2$NameBP, rowMeans(methLevels(object.sub)[
                        ,colData(object.sub)[[1]]== sidelevels[
                            2,i]]), col=pars$col[2],cex=pars$cex)
                }else{
                    points(d2$NameBP, rowMeans(methLevels(object.sub)[
                        ,colData(object.sub)[[1]]== sidelevels[
                            1,i]]), col=pars$col[1],cex=pars$cex,
                        pch=pars$pch[1])
                    points(d2$NameBP, rowMeans(methLevels(object.sub)[
                        ,colData(object.sub)[[1]]== sidelevels[
                            2,i]]), col=pars$col[2],cex=pars$cex, pch=pars$pch[
                                2])
                }

                if(!is.null(pars$ablines)){
                    abline(v=pars$ablines, col="gray50", lwd=pars$lwd)
                }
                if(pars$legend){
                    legend(pars$legend.adj[1],pars$legend.adj[2], c(sidelevels[
                        ,i]), lwd = pars$lwd, col = pars$col, cex = 0.95,
                        text.col = pars$col, pch = pars$pch, bty= 'n')
                }
            }
        }else{
            par( oma = pars$oma,  mar = pars$mar)
            d2 <- data.frame(dir=rep(0,length(object.sub)))
            ind.o = rowData(object.sub)[,1+lcont+1]=="Hypo"
            ind.e = rowData(object.sub)[,1+lcont+1]=="Equal"
            ind.r = rowData(object.sub)[,1+lcont+1]=="Hyper"
            d2$NameBP <- start(object.sub) - min(start(object.sub)) + 1
            d2$col = NA
            d2$col[ind.o]="black"
            d2$col[ind.e]="black"
            d2$col[ind.r]="black"
            d2$dmcH = NA
            d2$dmcH[ind.r]= 1.08
            d2$dmcH[ind.o]= NA
            d2$dmcH[ind.e]= NA
            d2$dmcH[rowData(object.sub)[,1+1]==0] = NA
            d2$dmcL = NA
            d2$dmcL[ind.r]= NA
            d2$dmcL[ind.e]= NA
            d2$dmcL[ind.o]= -0.08
            d2$dmcL[rowData(object.sub)[,1+1]==0] = NA

            do.call("plot",c(NA,def.arg[[1]]))
            axis(1, at=ticks, labels=nlabs)
            axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c(
                0,0.2,0.4,0.6,0.8,1),cex.axis=pars$cex.axis)
            axis(4, at=c(-0.05), labels=paste0(sidelevels[1,i]," < ",sidelevels[
                2,i]), las=1, cex.axis=pars$cex.axis , line = -2)
            axis(4, at=c(1.05),  labels=paste0(sidelevels[1,i]," > ",sidelevels[
                2,i]), las=1, cex.axis=pars$cex.axis , line = -2)

            with(d2, segments(x0=NameBP, y0=rep(1.02,length(
                d2$dmcH)), x1 = NameBP, y1 = d2$dmcH, lwd=pars$lwd, col = col))
            with(d2, segments(x0=NameBP, y0=rep(-0.02,length(
                d2$dmcL)), x1 = NameBP, y1 = d2$dmcL, lwd=pars$lwd, col = col))
            if(pars$type=="l"){
                lines(d2$NameBP, rowMeans(methLevels(object.sub)[
                    ,colData(object.sub)[[1]]== sidelevels[1,i]]
                ), col=pars$col[1],cex=pars$cex)
                lines(d2$NameBP, rowMeans(methLevels(object.sub)[
                    ,colData(object.sub)[[1]]== sidelevels[2,i]]), col=pars$col[
                        2],cex=pars$cex)
            }else{
                points(d2$NameBP, rowMeans(methLevels(object.sub)[
                    ,colData(object.sub)[[1]]== sidelevels[1,i]]), col=pars$col[
                        1],cex=pars$cex, pch=pars$pch[1])
                points(d2$NameBP, rowMeans(methLevels(object.sub)[,colData(
                    object.sub)[[1]]== sidelevels[2,i]]), col=pars$col[
                        2],cex=pars$cex, pch=pars$pch[2])
            }
            if(!is.null(pars$ablines)){
                abline(v=pars$ablines, col="gray50", lwd=pars$lwd)
            }
            if(pars$legend){
                legend(pars$legend.adj[1],pars$legend.adj[2], c(sidelevels[
                    ,i]), lwd = 1, col = pars$col,cex = 0.95, text.col=pars$col,
                    pch = pars$pch,bty= 'n')
            }
        }
        par(old.par)

    }
    }

#' @rdname plotDMCFB-method
#' @aliases plotDMCFB-method plotDMCFB
setMethod("plotDMCFB", signature(
    object="BSDMC", region="ANY", nSplit="ANY", parList="ANY"), .plotDMCFB)
