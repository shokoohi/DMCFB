
.MCMCFB <- function(object, bwa, bwb, nBurn, nMC, nThin, sdv, nCores, pSize,
    sfiles) {
    if (missing(object)) {
        stop("No BSData object is specified.\n")
    }

    if (missing(pSize)) pSize <- 500
    if (missing(bwa)) bwa <- 30
    if (missing(bwb)) bwb <- 30
    if (missing(nBurn)) nBurn <- 500
    if (missing(nMC)) nMC <- 500
    if (missing(nThin)) nThin <- 1
    if (missing(sfiles)) sfiles <- TRUE

    if (missing(sdv)) {
        sdv <- sqrt(mean(rowVars(methReads(
        object
    ) / totalReads(object)), na.rm = TRUE))
    message("The priors's SD = ", round(sdv,4), ", estimated from data ... ")
    }

    xp <- start(object)
    nPos <- length(xp)
    if(.Platform$OS.type=="unix"){
    maxram <- (get_ram())[[1]] / 1e+9
    }else{
        maxram <- memory.size(max = FALSE) / 1e+9
    }

    if (missing(nCores)) {
    if (maxram / multicoreWorkers() > 5) {
        nCores <- multicoreWorkers()
    } else {
        nCores <- floor(maxram / (5 * nBurn * nThin * pSize * bwa *
        bwb / 10000000))
        nCores <- max(1, min(nCores, multicoreWorkers()))
    }
    }

    message("Number of assigned cores: ", nCores, " ... ")

    if ((maxram / nCores) > 10) {
    aaa <- nCores * 2 * pSize
    Pos.list1 <- list()
    if (nPos <= (floor(5 * nCores / 4) * pSize)) {
        Pos.list1 <- list(xp = xp)
    } else {
        Pos.list1 <- split(xp, ceiling(seq_along(xp) / (nCores * 2 * pSize)))
        aa <- length(Pos.list1)
        if (length(Pos.list1[[aa]]) < (floor(nCores / 4) * pSize)) {
        Pos.list1[[aa - 1]] <- c(Pos.list1[[aa - 1]], Pos.list1[[aa]])
        Pos.list1 <- Pos.list1[-aa]
        }
    }
    } else {
        aaa <- nCores * pSize
        Pos.list1 <- list()
        if (nPos <= (floor(5 * nCores / 4) * pSize)) {
        Pos.list1 <- list(xp = xp)
        } else {
        Pos.list1 <- split(xp, ceiling(seq_along(xp) / (nCores * pSize)))
        aa <- length(Pos.list1)
        if (length(Pos.list1[[aa]]) < (floor(nCores / 4) * pSize)) {
        Pos.list1[[aa - 1]] <- c(Pos.list1[[aa - 1]], Pos.list1[[aa]])
        Pos.list1 <- Pos.list1[-aa]
    }
    }
    }

    formula <- as.formula(paste("logit(MethRead/ReadDepth)", paste(
    c(names(colData(object))),
    collapse = " + "
    ), sep = " ~ "))

    message("------------------------------------------------------------")
    message("Fitted model: ")
    if (length(c(names(colData(object)))) > 1) {
    formula2 <- as.formula(paste0("logit(MethRead/ReadDepth) ~ ", paste0(
    "F(", c(names(colData(object)))[1], ") + "
    ), paste(
    c(names(colData(object)))[-1],
    collapse = " + "
    )))
    message(formula2[2]," ", formula2[1]," ", formula2[3])
    } else {
    formula2 <- as.formula(paste0("logit(MethRead/ReadDepth) ~ ", paste0(
    "F(", c(names(colData(object)))[1], ") "
    )))
    message(formula2[2]," ", formula2[1]," ", formula2[3])
    }

    message("------------------------------------------------------------")
    message("Creating ", length(Pos.list1), " batches of genomic positions ...")


    object.df <- as.data.frame(colData(object))
    object.df <- object.df[all.vars(formula)[
    all.vars(formula) %in% names(colData(object))
    ]]
    ind1 <- vapply(object.df, is.character, logical(1))
    object.df[ind1] <- lapply(object.df[ind1], as.factor)
    ind1 <- vapply(object.df, is.factor, logical(1))
    object.df.factor <- as.data.frame(object.df[, ind1])
    names(object.df.factor) <- names(object.df)[ind1]
    ind1 <- vapply(object.df.factor, nlevels, numeric(1)) > 1
    fnames <- names(object.df.factor)[ind1]
    object.df.factor <- as.data.frame(object.df.factor[, ind1])
    names(object.df.factor) <- fnames

    if (ncol(object.df) < 1) stop("There is no covarite in the data.")

    tmpfiles <- list()
    for (indK in seq_along(Pos.list1)) {
    object.par <- object[start(object) %in% Pos.list1[[indK]], ]
    strand(object.par) <- "*"
    object.par <- sort(object.par)
    message(
    "Running batch ", sprintf(paste0("%", nchar(trunc(
        length(Pos.list1)
        )), "d"), indK), "/", length(
        Pos.list1
        ), "; ", paste0((seqnames(object.par))@values), "; ",
        sprintf(paste0("%", nchar(trunc(aaa)), "d"), dim(
        object.par
        )[1]), " positions; Region [",
        sprintf(
        paste0("%", nchar(trunc(
        tail(Pos.list1[[length(Pos.list1)]])[6]
        )), "d"),
        start(object.par)[1]
        ), ", ", sprintf(paste0("%", nchar(trunc(
        tail(Pos.list1[[length(Pos.list1)]])[6]
        )), "d"), tail(
        start(object.par)
        )[6]), "]; Date ", paste0(
        Sys.time()
        ))

    xp <- start(object.par)
    nPos <- nrow(object.par)

    Pos.list <- list()
    if (nPos <= pSize) {
        Pos.list <- list(xp = xp)
    } else {
        Pos.list <- split(xp, ceiling(seq_along(xp) / pSize))
        aa <- length(Pos.list)
        if (length(Pos.list[[aa]]) < pSize * 0.6) {
        Pos.list[[aa - 1]] <- c(Pos.list[[aa - 1]], Pos.list[[aa]])
        Pos.list <- Pos.list[-aa]
        }
    }

    mylist <- NULL

    for (i in seq_along(Pos.list)) {
        mylist[[i]] <- list(
        object = object.par[
        (start(object.par)) %in% Pos.list[[i]],
        ],
        bwa = bwa, bwb = bwb, nBurn = nBurn, nMC = nMC, nThin = nThin, sdv = sdv
        )
    }

    # options(warn=-1)
    .fitpar <- function(i) {
        mmlist <- i
        object.par2 <- mmlist[[1]]
        bwa <- mmlist[[2]]
        bwb <- mmlist[[3]]
        nBurn <- mmlist[[4]]
        nMC <- mmlist[[5]]
        nThin <- mmlist[[6]]
        sdv <- mmlist[[7]]
        nseq <- length(colnames(object.par2))
        formula <- as.formula(paste("logit(MethRead/ReadDepth)", paste(
        c(names(colData(object.par2))),
        collapse = " + "
        ), sep = " ~ "))
        object.par2.df <- as.data.frame(colData(object.par2))
        object.par2.df <- object.par2.df[all.vars(
        formula
        )[all.vars(formula) %in% names(colData(object.par2))]]
        ind1 <- vapply(object.par2.df, is.character, logical(1))
        object.par2.df[ind1] <- lapply(object.par2.df[ind1], as.factor)
        ind1 <- vapply(object.par2.df, is.factor, logical(1))
        object.par2.df.factor <- as.data.frame(object.par2.df[, ind1])
        names(object.par2.df.factor) <- names(object.par2.df)[ind1]
        ind1 <- vapply(object.par2.df.factor, nlevels, numeric(1)) > 1
        fnames <- names(object.par2.df.factor)[ind1]
        object.par2.df.factor <- as.data.frame(object.par2.df.factor[, ind1])
        names(object.par2.df.factor) <- fnames

        xp <- start(object.par2)
        Y0 <- methReads(object.par2)
        M0 <- totalReads(object.par2)

        nits <- nBurn + nMC * nThin

        G <- length(unique(object.par2.df.factor[, 1]))

        L <- length(xp)

        nz <- 2 + sum(vapply(
        object.par2.df.factor, nlevels,
        numeric(1)
        )[-1] - 1) +
        dim(object.par2.df)[2] - dim(object.par2.df.factor)[2]


        Xns1 <- ns(xp, df = bwa)
        Xns2 <- ns(xp, df = bwb)

        ni <- as.vector(table(object.par2.df.factor[, 1]))

        Subject <- rep(seq_len(nseq), each = L)

        Beta.old <- matrix(0, nrow = nseq, ncol = bwb)

        GammaDeltaEta.old <- rep(0, G * bwa + nz)

        Beta.store <- array(0, c(nMC, nseq, bwb))
        GammaDeltaEta.store <- matrix(0, nrow = nMC, ncol = G * bwa + nz)

        Xb <- rep(1, nseq) %x% Xns2
        Xg <- rep(1, nseq) %x% Xns1

        ConT <- contrasts(object.par2.df[, 1])

        Xd.mat <- NULL
        for (i in seq_len(dim(ConT)[2])) {
        Xd.mat1 <- NULL
        for (j in seq_len(dim(ConT)[1])) {
            Xd.mat1 <- rbind(Xd.mat1, rep(ConT[j, i], ni[j]) %x% Xns1)
        }
        Xd.mat <- cbind(Xd.mat, Xd.mat1)
        }

        Xe <- cbind(rep(1, nseq * L), log(as.vector(data.matrix(M0)) + 1))

        if (dim(object.par2.df.factor)[2] > 1) {
        results <- as.data.frame(object.par2.df.factor[, -1])
        names(results) <- names(object.par2.df.factor)[-1]
        results <- fastDummies::dummy_cols(results, remove_first_dummy = TRUE)
        results <- subset(results, select = names(results)[
            !(names(results) %in% names(object.par2.df.factor))
        ])
        for (jj in names(results)) {
            Xe <- cbind(Xe, rep(results[, jj], each = L))
        }
        for (jj in colnames(object.par2.df)[
            !(colnames(object.par2.df) %in% fnames)
        ]) {
            Xe <- cbind(Xe, rep(object.par2.df[, jj], each = L))
        }
        }

        Xgde <- cbind(Xg, Xd.mat, Xe)
        data.all <- cbind(as.vector(data.matrix(Y0)), as.vector(
        data.matrix(M0)
        ))

        Beta.mean <- as.vector(Xns2 %*% t(Beta.old))

        GammaDeltaEta.mean <- Xgde %*% GammaDeltaEta.old

        #options(show.error.messages = FALSE)
        suppressWarnings(fit0 <- try(
        glm(cbind(data.all[, 1], data.all[, 2] - data.all[
            , 1
        ]) ~ offset(GammaDeltaEta.mean + Beta.mean) - 1,
        family =
            "binomial"
        ),
        silent = TRUE
        ))
        #options(show.error.messages = TRUE)
        if ((is(fit0)[1] == "try-error")) {
        #options(warn = -1)
        fit0 <- bayesglm(cbind(data.all[, 1], data.all[, 2] - data.all[, 1]) ~
        offset(GammaDeltaEta.mean + Beta.mean) - 1,
        family = binomial(logit), maxit = 1
        )
        #options(warn = 0)
        } else if (sum(is.na(coefficients(fit0))) > 0) {
        #options(warn = -1)
        fit0 <- bayesglm(cbind(data.all[, 1], data.all[, 2] - data.all[, 1]) ~
        offset(GammaDeltaEta.mean + Beta.mean) - 1,
        family = binomial(logit), maxit = 1
        )
        #options(warn = 0)
        }

        # log.like = as.numeric(logLik(fit0))

        BPrec <- diag(rep(1 / sdv, bwb))
        GDEPrec <- diag(rep(1 / sdv, bwa * G + nz))

        ico <- 0
        for (myiter in seq_len(nits)) {
        gdata <- data.frame(Y = data.all[, 1], m = data.all[
            , 2
        ], Beta.mean, Xgde)

        #options(show.error.messages = FALSE)
        suppressWarnings(fitg <- try(
            speedglm(cbind(Y, m - Y) ~ offset(
            Beta.mean
            ) - 1 + Xgde, data = gdata, family = binomial(
            logit
            )),
            silent = TRUE
        ))
        #options(show.error.messages = TRUE)

        if ((is(fitg)[1] == "try-error")) {
            #options(warn = -1)
            fitg <- bayesglm(cbind(Y, m - Y) ~ offset(
            Beta.mean
            ) - 1 + Xgde, data = gdata, family = binomial(
            logit
            ), maxit = 1)
            #options(warn = 0)
            } else if (sum(is.na(coefficients(fitg))) > 0) {
            #options(warn = -1)
            fitg <- bayesglm(cbind(Y, m - Y) ~ offset(
            Beta.mean
            ) - 1 + Xgde,
            data = gdata,
            family = binomial(logit), maxit = 1
            )
            #options(warn = 0)
            }

            mle.Sig <- summary(fitg)$cov.unscaled

            #options(show.error.messages = FALSE)
            suppressWarnings(mle.Prec <- try(
            ginv(mle.Sig),
            silent = TRUE
            ))
            suppressWarnings(Post.Cov <- try(
            ginv(mle.Prec + GDEPrec),
            silent = TRUE
            ))
            #options(show.error.messages = TRUE)

            myw <- 0
            while (myw < 100 & ((is(mle.Prec)[1] == "try-error") | (
            is(Post.Cov)[1] == "try-error"))) {
            myw <- myw + 1
            #options(warn = -1)
            fitg <- bayesglm(cbind(Y, m - Y) ~ offset(
            Beta.mean
            ) - 1 + Xgde,
            data = gdata,
            family = binomial(logit), maxit = 1
            )
            #options(warn = 0)

            #options(show.error.messages = FALSE)
            suppressWarnings(mle.Prec <- try(
            ginv(mle.Sig),
            silent = TRUE
            ))

            suppressWarnings(Post.Cov <- try(
            ginv(mle.Prec + GDEPrec),
            silent = TRUE
            ))
            #options(show.error.messages = TRUE)
            }

            Post.mean <- Post.Cov %*% (mle.Prec %*% coef(fitg))
            GammaDeltaEta.old <- mvrnorm(1, mu = Post.mean, Sigma = Post.Cov)
            GammaDeltaEta.mean <- Xgde %*% GammaDeltaEta.old

            for (i in seq_len(nseq)) {
            Y <- Y0[, i]
            m <- M0[, i]
            ivec <- Subject == i

            idata <- data.frame(
            Y = Y, m = m, GammaDeltaEta.mean = GammaDeltaEta.mean[
            ivec
            ], Xns2
            )

            #options(show.error.messages = FALSE)
            suppressWarnings(fiti <- try(
            speedglm(cbind(Y, m - Y) ~ offset(
            GammaDeltaEta.mean
            ) - 1 + Xns2,
            data = idata,
            family = binomial(logit)
            ),
            silent = TRUE
            ))
            #options(show.error.messages = TRUE)

            if ((is(fiti)[1] == "try-error")) {
            #options(warn = -1)
            fiti <- bayesglm(cbind(Y, m - Y) ~ offset(
            GammaDeltaEta.mean
            ) - 1 + Xns2,
            data = idata,
            family = binomial(logit), maxit = 1
            )
            options(warn = 0)
            } else if (sum(is.na(coefficients(fiti))) > 0) {
            options(warn = -1)
            fiti <- bayesglm(cbind(Y, m - Y) ~ offset(
            GammaDeltaEta.mean
            ) - 1 + Xns2,
            data = idata,
            family = binomial(logit), maxit = 1
            )
            #options(warn = 0)
            }

            mle.Sig <- summary(fiti)$cov.unscaled

            #options(show.error.messages = FALSE)
            suppressWarnings(mle.Prec <- try(
            ginv(mle.Sig),
            silent = TRUE
            ))

            suppressWarnings(Post.Cov <- try(
            ginv(mle.Prec + BPrec),
            silent = TRUE
            ))
            #options(show.error.messages = TRUE)

            myw <- 0
            #options(warn = -1)
            while (((is(mle.Prec)[1] == "try-error") | (
            is(Post.Cov)[1] == "try-error")) & myw < 100) {
            myw <- myw + 1

            fiti <- bayesglm(cbind(Y, m - Y) ~ offset(
            GammaDeltaEta.mean
            ) - 1 + Xns2,
            data = idata,
            family = binomial(logit), maxit = 1
            )

            mle.Sig <- summary(fiti)$cov.unscaled

            #options(show.error.messages = FALSE)
            suppressWarnings(mle.Prec <- try(
            ginv(mle.Sig),
            silent = TRUE
            ))

            suppressWarnings(Post.Cov <- try(
            ginv(mle.Prec + BPrec),
            silent = TRUE
            ))
            #options(show.error.messages = TRUE)
            }
            #options(warn = 0)

            Post.mean <- Post.Cov %*% (mle.Prec %*% coef(fiti))
            Beta.old[i, ] <- mvrnorm(1, mu = Post.mean, Sigma = Post.Cov)
        }
        Beta.mean <- as.vector(Xns2 %*% t(Beta.old))

        #options(show.error.messages = FALSE)
        suppressWarnings(fit0 <- try(
            glm(cbind(data.all[, 1], data.all[, 2] - data.all[, 1]) ~
            offset(GammaDeltaEta.mean + Beta.mean) - 1,
            family = "binomial"
            ),
            silent = TRUE
        ))
        #options(show.error.messages = TRUE)

        if ((is(fit0)[1] == "try-error")) {
            #options(warn = -1)
            fit0 <- bayesglm(cbind(data.all[, 1], data.all[, 2] - data.all[
            , 1
            ]) ~ offset(GammaDeltaEta.mean + Beta.mean) - 1,
            family = binomial(logit), maxit = 1
            )
            #options(warn = 0)
        } else if (sum(is.na(fit0$coefficients)) > 0) {
            #options(warn = -1)
            fit0 <- bayesglm(cbind(data.all[, 1], data.all[
            , 2
            ] - data.all[
            , 1
            ]) ~ offset(GammaDeltaEta.mean + Beta.mean) - 1,
            family = binomial(logit), maxit = 1
            )
            #options(warn = 0)
        }

        # log.like = as.numeric(logLik(fit0))

        if (myiter > nBurn & (myiter %% nThin == 0)) {
            ico <- ico + 1
            Beta.store[ico, , ] <- Beta.old
            GammaDeltaEta.store[ico, ] <- GammaDeltaEta.old
        }
        }
        rm(gdata, Xb, Xg, Xd.mat)

        Cmb <- combn(G, 2)
        Cmb <- cbind(Cmb[, seq_len(G - 1)], matrix(c(1, G + 1), 2, 1), Cmb[
        , -c(seq_len(G - 1))
        ] + 1)

        GammaDeltaEta.hold <- matrix(0, nrow = nMC, ncol = bwa)
        GammaDeltaEta.hold <- cbind(GammaDeltaEta.hold, GammaDeltaEta.store)

        GMAT <- list()
        for (k in seq_len(dim(Cmb)[2])) {
        gfit.mat <- matrix(0, nrow = nMC, ncol = L)
        for (i in seq_len(nMC)) {
            gfit.mat[i, ] <- Xns1 %*% (GammaDeltaEta.hold[i, ((
            Cmb[2, k] - 1) * bwa + 1):(Cmb[2, k] * bwa)] - GammaDeltaEta.hold[
            i, ((Cmb[1, k] - 1) * bwa + 1):(Cmb[1, k] * bwa)
            ])
        }
        GMAT[[k]] <- gfit.mat
        }

        bhat <- matrix(0, nrow = nseq, ncol = L)
        bstar <- matrix(0, nrow = nseq, ncol = L)
        for (j in seq_len(nseq)) {
        bfit.mat <- matrix(0, nrow = nMC, ncol = L)
        bfit.mat1 <- matrix(0, nrow = nMC, ncol = L)
        for (i in seq_len(nMC)) {
            bfit.mat[i, ] <- Xns2 %*% Beta.store[i, j, ]
            bfit.mat1[i, ] <- bfit.mat[i, ] +
            Xgde[((j - 1) * L + 1):(j * L), ] %*% GammaDeltaEta.store[i, ]
        }
        bhat[j, ] <- apply(bfit.mat, 2, median)
        bstar[j, ] <- apply(bfit.mat1, 2, median)
        }
        rm(Xgde)
        return(list(bhat = bhat, GMAT = GMAT, bstar = bstar))
    }

    # optbp <- MulticoreParam(workers = nCores, progressbar = FALSE)
    optbp <- SnowParam(workers = nCores)
    # registerDoParallel(nCores)
    # optbp <- DoparParam()
    register(optbp)
    .myfun <- function(mylist, .fitpar) {
        bplapply(mylist, .fitpar, BPPARAM = optbp)
    }
    totres <- .myfun(mylist, .fitpar)

    tmpf <- tempfile(
        pattern = paste0("object.part", indK, "-"),
        fileext = ".rds"
    )
    saveRDS(totres, tmpf)
    tmpfiles[[indK]] <- tmpf
    rm(totres)
    }

    message("------------------------------------------------------------")
    message("Combining ", length(tmpfiles), " objects ...")
    object.new <- NULL
    object.new <- readRDS(tmpfiles[[1]])
    if (length(tmpfiles) > 1) {
    for (indK in seq_along(tmpfiles)[-1]) {
        object.hold <- readRDS(tmpfiles[[indK]])
        object.new <- c(object.new, object.hold)
        rm(object.hold)
    }
    }
    message("Objects are combined.")

    metadata(object)$MClist <- object.new

    SNAM <- paste0(
    "BSMCMC-pSize", pSize, "bwa", bwa, "bwb", bwb, "nBurn", nBurn, "nMC",
    nMC, "nThin", nThin, ".rds"
    )
    if (sfiles) {
    saveRDS(object, paste0(SNAM))
    message(
        "BSDMC object (includes MCMC samples) stored as '", paste0(
        SNAM
        ), "'."
    )
    }
    return(object)
}

.findFB <- function(object, alpha) {
    if (missing(object)) {
    stop("No BSDMC object is specified.\n")
    }

    if (missing(alpha)) alpha <- 0.05

    formula <- as.formula(paste("logit(Methylation)", paste(
    c(names(colData(object))),
    collapse = " + "
    ), sep = " ~ "))

    object.df <- as.data.frame(colData(object))
    object.df <- object.df[all.vars(formula)[all.vars(
    formula
    ) %in% names(colData(object))]]
    ind1 <- vapply(object.df, is.character, logical(1))
    object.df[ind1] <- lapply(object.df[ind1], as.factor)
    ind1 <- vapply(object.df, is.factor, logical(1))
    object.df.factor <- as.data.frame(object.df[, ind1])
    names(object.df.factor) <- names(object.df)[ind1]
    ind1 <- vapply(object.df.factor, nlevels, numeric(1)) > 1
    fnames <- names(object.df.factor)[ind1]
    object.df.factor <- as.data.frame(object.df.factor[, ind1])
    names(object.df.factor) <- fnames

    if (ncol(object.df) < 1) stop("There is no covarite in the data.\n")

    G <- length(unique(object.df.factor[, 1]))
    nPos <- length(start(object))

    Cmb <- combn(G, 2)
    Cmb <- cbind(Cmb[, seq_len(G - 1)], matrix(c(1, G + 1), 2, 1), Cmb[
    , -c(seq_len(G - 1))] + 1)

    MClist <- metadata(object)$MClist

    for (mm in seq_along(MClist)) {
    gfit <- list()
    for (k in seq_len(dim(Cmb)[2])) {
    gfit[[k]] <- t(colQuantiles(MClist[[mm]]$GMAT[[k]],
        probs = c(alpha / 2, 0.5, 1 - alpha / 2)))
    }
    MClist[[mm]]$gfit <- gfit
    }

    if (length(MClist[[1]]$gfit) > 2) {
    gfit.pairs <- lapply(MClist, function(elt) {
        elt$gfit
    })
    gfit.pairs <- lapply(gfit.pairs, function(x) {
        do.call("rbind", x[-1])
    })
    gfit.pairs <- do.call("cbind", gfit.pairs)
    DMCs.pairs <- matrix(rep(c(0, 1), (dim(gfit.pairs)[1] / 3) * (
        dim(gfit.pairs)[2])), nrow = dim(gfit.pairs)[1] / 3, ncol = dim(
        gfit.pairs
    )[2])
    DMCs.dir <- matrix(rep(
        c("Equal", "Hyper", "Hypo"),
        (dim(gfit.pairs)[1] / 3) * (
        dim(gfit.pairs)[2])
    ), nrow = dim(gfit.pairs)[1] / 3, ncol = dim(
        gfit.pairs
    )[2])
    DMCs <- matrix(rep(c(0, 1), dim(gfit.pairs)[2]), nrow = dim(gfit.pairs)[
        2
    ], ncol = 1)
    for (j in seq_len(dim(gfit.pairs)[1] / 3)) {
        DMCs.pairs[j, ] <- as.numeric(gfit.pairs[(j - 1) * 3 + 1, ] *
        gfit.pairs[(j - 1) * 3 + 3,] > 0)
        DMCs.dir[j, gfit.pairs[(j - 1) * 3 + 2, ] > 0] <- "Hyper"
        DMCs.dir[j, gfit.pairs[(j - 1) * 3 + 2, ] < 0] <- "Hypo"
        DMCs.dir[j, DMCs.pairs[j, ] == 0] <- "Equal"
    }
    DMCs[, 1] <- (colMeans(DMCs.pairs) > 0) * 1
    DMCs.dir <- t(DMCs.dir)
    gfit.pairs <- t(gfit.pairs)
    gfit.indiv <- gfit.pairs[, seq(2, dim(gfit.pairs)[2], by = 3)]
    gfit.indiv <- gfit.indiv[, as.numeric(object.df.factor[, 1])]
    DMCs.pairs <- t(DMCs.pairs)
    } else {
    gfit.pairs <- lapply(MClist, function(elt) {
        elt$gfit
    })
    gfit.pairs <- lapply(gfit.pairs, function(x) {
        do.call("rbind", x[-1])
    })
    gfit.pairs <- do.call("cbind", gfit.pairs)
    DMCs <- matrix(rep(c(0, 1), (dim(gfit.pairs)[2])), nrow = dim(gfit.pairs)[
        2
    ], ncol = 1)
    DMCs.dir <- matrix(rep(c("Equal", "Hyper", "Hypo"), (dim(
        gfit.pairs
    )[2])), nrow = dim(gfit.pairs)[2], ncol = 1)
    j <- 1
    DMCs[, j] <- as.numeric(gfit.pairs[(j - 1) * 3 + 1, ] *
        gfit.pairs[(j - 1) * 3 + 3, ] > 0)
    DMCs.dir[gfit.pairs[(j - 1) * 3 + 2, ] > 0, j] <- "Hyper"
    DMCs.dir[gfit.pairs[(j - 1) * 3 + 2, ] < 0, j] <- "Hypo"
    DMCs.dir[DMCs == 0, j] <- "Equal"
    DMCs.pairs <- DMCs
    gfit.pairs <- t(gfit.pairs)
    gfit.indiv <- gfit.pairs[, seq(2, dim(gfit.pairs)[2], by = 3)]
    }
    np <- paste(combn(levels(object.df.factor[, 1]), 2)[2, ], combn(
    levels(object.df.factor[, 1]), 2
    )[1, ], sep = "vs")
    namesPairCon <- paste0(
    rep(c(paste("Cont", np, sep = "")), each = 3),
    colnames(gfit.pairs)
    )
    namesPairDMC <- c(paste("DMCs", np, sep = ""))
    namesPairDIR <- c(paste("", np, sep = ""))
    colnames(DMCs.pairs) <- namesPairDMC
    colnames(DMCs.dir) <- namesPairDIR
    colnames(gfit.pairs) <- namesPairCon
    colnames(DMCs) <- "DMCs"

    bhat <- t(do.call("cbind", lapply(MClist, function(elt) {
    elt$bhat
    })))
    logitbstar <- t(do.call("cbind", lapply(MClist, function(elt) {
    elt$bstar
    })))
    expit <- function(x) {
    return(1 / (1 + exp(-x)))
    }
    names(bhat) <- colnames(object)
    bstar <- expit(logitbstar)
    colnames(bstar) <- colnames(object)

    object.new <- cBSDMC(
    rowRanges = rowRanges(object),
    methReads = methReads(object),
    totalReads = totalReads(object),
    methLevels = bstar,
    colData = colData(object)
    )
    rowData(object.new) <- DataFrame(DMCs, DMCs.pairs, DMCs.dir, gfit.pairs)

    return(object.new)
}

.findDMCFB <- function(object, bwa, bwb, nBurn, nMC, nThin, alpha, sdv,
    nCores, pSize, sfiles) {
    if (missing(object)) {
    stop("No BSData object is specified.\n")
    }
    if (missing(bwa)) bwa <- 30
    if (missing(bwb)) bwb <- 30
    if (missing(nBurn)) nBurn <- 500
    if (missing(nMC)) nMC <- 500
    if (missing(nThin)) nThin <- 1
    if (missing(pSize)) pSize <- 500
    if (missing(sfiles)) sfiles <- TRUE

    if (!(c("MClist") %in% names(metadata(object)))) {
    message("------------------------------------------------------------")
    message("Running Bayesian functional regression model ...")
    object <- .MCMCFB(
        object, bwa, bwb, nBurn, nMC, nThin, sdv, nCores, pSize, sfiles
    )
    }
    message("------------------------------------------------------------")
    message("Identifying DMCs ...")
    object <- .findFB(object, alpha)
    message("DMCs are identified.")
    SNAM <- paste0(
    "BSDMC-pSize", pSize, "bwa", bwa, "bwb", bwb, "nBurn", nBurn, "nMC",
    nMC, "nThin",
    nThin, "alpha", alpha, ".rds"
    )
    if (sfiles) {
    saveRDS(object, paste0(SNAM))
    message("BSDMC object (includes DMCs) stored as ", SNAM)
    }
    message("------------------------------------------------------------")
    res <- table(factor(rowData(object)[, 1], levels = c(0, 1))) /
    dim(object)[1] * 100
    names(res) <- c("Equal(%)", "DMC(%)")
    message("Percentage of non-DMCs and DMCs:")
    print(res)
    message("------------------------------------------------------------")
    message("Percentage of hyper-, hypo-, and equal-methylated positions:")
    if ((dim(rowData(object))[2] - 1) / 5 == 1) {
    res <- data.matrix(100 * table(factor(unlist(rowData(object)[(2 + (dim(
        rowData(object)
    )[
        2
    ] - 1) / 5):(1 + 2 * (dim(rowData(object))[
        2
    ] - 1) / 5)]), levels = c("Equal", "Hyper", "Hypo"))) / dim(object)[1])
    rownames(res) <- c("Equal(%)", "Hyper(%)", "Hypo(%)")
    colnames(res) <- names(rowData(object))[2 + (dim(rowData(object))[2]-1)/5]
    } else {
        res = rowData(object)[(2 + (dim(rowData(object))[2] - 1) / 5):(1 +
            2 * (dim(rowData(object))[2] - 1) / 5)]

        res[sapply(res, is.character)] <- lapply(res[sapply(res, is.character)],
            function(x){factor(x, levels = c("Equal", "Hyper", "Hypo") )})

        res <- 100*vapply(res, table, FUN.VALUE = numeric(3L)) / dim(object)[1]
        rownames(res) <- c("Equal(%)", "Hyper(%)", "Hypo(%)")
    }
    print(t(res))
    message("------------------------------------------------------------")

    return(object)
}


#' @rdname findDMCFB-method
#' @aliases findDMCFB-method findDMCFB
setMethod("findDMCFB", signature = c(
    object = "BSDMC", bwa = "ANY", bwb = "ANY", nBurn = "ANY", nMC = "ANY",
    nThin = "ANY", alpha = "ANY", sdv = "ANY", nCores = "ANY", pSize = "ANY",
    sfiles = "ANY"
), .findDMCFB)
