##############################################################################
## Functions for different methods
## stdErr is estimated with resampling if method = sc or am.xc,
##        bootstrap otherwise
##############################################################################

doREFit.am.XCHWY <- function(DF, DF0, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    ## Reset PE and SE
    alpha <- beta <- gamma <- rep(0, p)
    ## Start am.xc
    outA <- dfsane(alpha, alphaEq, X = X, Y = Y, T = T, cluster = cluster, mt = mt,
                   weights = NULL, alertConvergence = FALSE, quiet = TRUE, 
                   control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    alpha <- outA$par
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    lambda <- npMLE(Ystar[event == 0], Tstar, Ystar)
    zHat <- as.numeric(mt * npMLE(log(max(Y)), Tstar, Ystar) / lambda)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ## outB <- dfsane(beta, betaEq, X = X, Y = Y, T = T, cluster = cluster,
    ##                delta = status[event == 0], mt = mt,
    ##                alpha = outA$par, zHat = zHat, weights = NULL,
    ##                alertConvergence = FALSE, quiet = TRUE, 
    ##                control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    outB <- BBsolve(beta, betaEq, X = X, Y = Y, T = T, cluster = cluster,
                   delta = status[event == 0], mt = mt,
                   alpha = outA$par, zHat = zHat, weights = NULL, quiet = TRUE)
    list(alpha = outA$par, aconv = outA$convergence, beta = outB$par, bconv = outB$convergence,
         muZ = mean(zHat), zHat = zHat)
}

doREFit.am.GL <- function(DF, DF0, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))  
    ## Reset PE and SE
    alpha <- beta <- gamma <- rep(0, p)
    aSE <- bSE <- da <- va <- db <- vb <- NA
    ## Start cox.GL
    outB <- aftsrr(Surv(Y, status) ~ X, subset = event == 0, B = 0,
                   rankWeights = "logrank", method = "nonsm")
    outA <- dfsane(alpha, ghoshU2, beta = outB$beta, T = ifelse(T == Y, 1e5, T),
                   Y = Y[event == 0], X = as.matrix(X[event == 0, ]),
                   cl = mt + 1, ## unlist(lapply(split(id, id), length)),
                   alertConvergence = FALSE, quiet = TRUE, 
                   control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    outB$par <- -1 * outB$beta
    outA$par <- -1 * outA$par
    list(alpha = outA$par, aconv = outA$convergence,
         beta = outB$par, bconv = outB$convergence, muZ = NA)
}

doREFit.cox.HW <- function(DF, DF0, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    ## Reset PE and SE
    alpha <- beta <- gamma <- rep(0, p)
    aSE <- bSE <- da <- va <- db <- vb <- NA
    gamma <- c(0, gamma)
    X <- cbind(1, X[event == 0,])
    ## outA <- dfsane(gamma, HWeq, X = X, Y = Y, T = T, cluster = cluster, mt = mt,
    ##                alertConvergence = FALSE,
    ##                quiet = TRUE, control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    outA <- BBsolve(gamma, HWeq, X = X, Y = Y, T = T, cluster = cluster, mt = mt, quiet = TRUE)
    alpha <- outA$par <- outA$par[-1]
    muZ <- outA$par[1]
    lambda <- npMLE(Y[event == 0], T, Y)
    ## zHat <- as.numeric(mt * npMLE(max(Y), T, Y) / (lambda * exp(X[, -1] %*% alpha)))
    zHat <- as.numeric(mt / (lambda * exp(as.matrix(X[, -1]) %*% alpha)))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    ## muZ <- mean(zHat)
    outB <- dfsane(beta, HWeq2, X = as.matrix(X[,-1]), Y = Y[event == 0],
                   delta = status[event == 0], zHat = zHat/muZ,
                   alertConvergence = FALSE, quiet = TRUE,
                   control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    list(alpha = outA$par, aconv = outA$convergence,
         beta = outB$par, bconv = outB$convergence, muZ = muZ)
}

doREFit.cox.LWYY <- function(DF, DF0, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))   
    out <- dfsane(par = double(p), fn = LWYYeq, X = as.matrix(X[event == 0, ]),
                  Y = Y[event == 0], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                  ## cl = unlist(lapply(split(id, id), length)), 
                  alertConvergence = FALSE, quiet = TRUE,
                  control = list(NM = FALSE, M = 100, noimp = 50, trace = FALSE))
    list(alpha = out$par, beta = rep(0, p), muZ = NA)
}

doREFit.sc.XCYH <- function(DF, DF0, engine, stdErr) {
    ## New data structure:
    ## id <- DF$id
    ## event <- DF$event
    ## status <- DF$status
    ## X <- as.matrix(DF[,-c(1:4)])    
    ## n <- length(unique(id))
    ## p <- ncol(X)
    ## T <- DF$Time
    ## mt <- aggregate(event ~ id, data = DF, sum)$event
    ## Y <- rep(DF$Time[event == 0], mt + 1)
    ## cluster <- unlist(sapply(mt + 1, function(x) 1:x))   
    ## Need to update the following, 10/31/2017
    id <- DF$id
    clsz <- unlist(lapply(split(id, id), length))
    ind <- cumsum(clsz)
    df <- DF
    m <- rep(clsz, clsz) - 1
    y <- rep(df$Time[ind], clsz)
    rmv <- sort(unique(rep(ind, clsz) * (m > 0)))[-1]
    df$Time <- df$Time * (m > 0)
    df <- df[-rmv, ]
    m <- m[-rmv]
    y <- y[-rmv]
    out <- with(df, sarmRV(id, Time, y, df[,-(1:4)], m, seq(0, max(Time), length.out = 100)))
    list(alpha = out$ahat, beta = out$bhat)
}

##############################################################################
# Variance estimation 
##############################################################################
doREFit.Engine.Bootstrap <- function(DF, DF0, engine, stdErr) {
    res <- doREFit(DF, DF0, engine, NULL)
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    clsz <- mt + 1
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    B <- stdErr@B
    betaMatrix <- matrix(0, B, p * 2)
    convergence <- rep(0, B)
    uID <- unique(DF$ID)
    for (i in 1:B) {
        sampled.id <- sample(unique(id), n, TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:n, clsz[sampled.id])
        betaMatrix[i,] <- with(doREFit(DF2, DF0, engine, NULL), c(alpha, beta))
        convergence[i] <- 1 * (betaMatrix[i,] %*% betaMatrix[i,] >
                               1e3 * with(res, c(alpha, beta) %*% c(alpha, beta)))
    }
    converged <- which(convergence == 0)
    if (sum(convergence != 0) > 0) {
        print("Warning: Some bootstrap samples failed to converge")
        tmp <- apply(betaMatrix, 1, function(x) x %*% x)
        converged <- (1:B)[- which(tmp %in% boxplot(tmp, plot = FALSE)$out)]        
    }
    if (all(convergence != 0) || sum(convergence == 0) == 1) {
        print("Warning: some bootstrap samples failed to converge")
        converged <- 1:B
    }
    betaVar <- var(betaMatrix[converged, ], na.rm = TRUE)
    betaSE <- sqrt(diag(as.matrix(betaVar)))
    c(res, list(alphaSE = betaSE[1:p], betaSE = betaSE[1:p + p],
                alphaVar = betaVar[1:p, 1:p], betaVar = betaVar[1:p + p, 1:p + p],
                SEmat = betaMatrix, B = length(converged)))
}

sdOut <- function(dat) {
    ol <- boxplot(dat, plot = FALSE)$out
    if (length(ol) > 0)
        dat <- dat[-which(dat %in% ol)]
    sd(dat, na.rm = TRUE)
}

doREFit.am.XCHWY.resampling <- function(DF, DF0, engine, stdErr) {
    res <- doREFit(DF, DF0, engine, NULL)
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))        
    B <- stdErr@B
    aSE <- bSE <- da <- va <- db <- vb <- NA
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    ua <- matrix(apply(Z, 2, function(x) alphaEq(res$alpha + n ^ (-0.5) * x, X, Y, T, cluster, mt)),
                 nrow = p)
    da <- t(apply(ua, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
    ua2 <- apply(E, 2, function(x) alphaEq(res$alpha, X, Y, T, cluster, mt, weights = x))
    va <- var(t(matrix(ua2, nrow = p)))
    if (qr(da)$rank == p)
        aVar <- solve(da) %*% va %*% t(solve(da))
    if (qr(da)$rank != p)
        aVar <- ginv(da) %*% va %*% t(ginv(da))
    aSE <- sqrt(diag(aVar))
    ub <- matrix(apply(Z, 2, function(x)
        betaEq(X, Y, T, cluster, mt, status[event == 0],
               res$zHat, res$alpha, res$beta + n ^ (-0.5) * x)), nrow = p)
    db <- t(apply(ub, 1, function(x) lm(n ^ (0.5) * x ~ t(Z))$coef[-1]))
    ub2 <- apply(E, 2, function(x) betaEq(X, Y, T, cluster, mt, status[event == 0],
                                          NULL, res$alpha, res$beta, weights = x))
    vb <- var(t(matrix(ub2, nrow = p)))
    if (qr(db)$rank == p)
        bVar <- solve(db) %*% vb %*% t(solve(db))
    if (qr(db)$rank != p)
        bVar <- ginv(db) %*% vb %*% t(ginv(db))
    bSE <- sqrt(diag(bVar))
    c(res, list(alphaSE = aSE, betaSE = bSE, da = da, va = va, db = db, vb = vb, B = stdErr@B))
}

doREFit.sc.XCYH.resampling <- function(DF, DF0, engine, stdErr) {
    id <- DF$id
    clsz <- unlist(lapply(split(id, id), length))
    ind <- cumsum(clsz)
    df <- DF
    m <- rep(clsz, clsz) - 1
    y <- rep(df$Time[ind], clsz)
    rmv <- sort(unique(rep(ind, clsz) * (m > 0)))[-1]
    df$Time <- df$Time * (m > 0)
    df <- df[-rmv, ]
    m <- as.numeric(m[-rmv])
    y <- as.numeric(y[-rmv])
    p <- ncol(df[,-(1:4)])
    out <- with(df, sarmRV(id, Time, y, df[,-(1:4)], m, seq(0, max(Time), length.out = 100)))
    outSE <- with(df, sarmRV.sand(id, Time, y, df[,-(1:4)], m, 
                                  a = out$ahat, b = out$ghat, Bootstrap = stdErr@B))
    list(alpha = out$ahat, beta = out$bhat, alphaSE = outSE$alphaSE, betaSE = outSE$betaSE,
         gamma = out$ghat, gammaSE = sqrt(diag(outSE$ase[(p + 2):(2 * p + 1), (p + 2):(2 * p + 1)])),
         varMat = outSE$ase)
}

##############################################################################
# Nonparametric (~1)
##############################################################################

doNonpara.am.XCHWY <- function(DF, alpha, beta, formula, engine, stdErr) {
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y) + X %*% alpha
    Ta <- log(T) + X %*% alpha
    lambda <- npMLE(Ya[which(event == 0)], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / lambda)
    ly <- ly / max(ly)
    win.ly <- max(ly)
    muZ <- mean(zHat)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[which(cluster == 1)]
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb), zHat / muZ, status[event == 0]))
    win.hy <- max(hy)
    list(t0 = t0, lam = ly * muZ, lamU = rep(NA, ng), lamL = rep(NA, ng), 
         haz = hy, hazU = rep(NA, ng), hazL = rep(NA, ng))
}

doNonpara.cox.NA <- function(DF, alpha, beta, formula, engine, stdErr) {
    ## t0 <- seq(0, max(DF$Time), length.out = 5 * nrow(DF))
    T <- DF$Time
    id <- DF$id
    event <- DF$event
    status <- DF$status
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)  
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    event <- DF$event
    status <- DF$status
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    cluster <- unlist(lapply(split(id, id), function(z) 1:length(z)))
    X <- as.matrix(DF[,-(1:4)])
    DF$T0 <- with(DF, unlist(lapply(split(Time, id), function(x) c(0, x)[1:length(x)])))
    if (!all(X == 0)) {
        tmp <- coxph(as.formula(paste("Surv(T0, Time, event)~",
                                      paste(colnames(DF)[-c(1:4, ncol(DF))], collapse = "+"))),
                     data = DF)
        base <- data.frame(matrix(0, nrow = 1, ncol = ncol(X)))
        names(base) <- names(coef(tmp))
        tmp0 <- survfit(tmp, newdata = base)
        ly <- with(tmp0, approx(time, cumhaz, t0)$y)
        lyU <- -log(with(tmp0, approx(time, upper, t0)$y))
        lyL <- -log(with(tmp0, approx(time, lower, t0)$y))
        tmp <- coxph(as.formula(paste("Surv(Time, status)~",
                                      paste(colnames(DF)[-c(1:4, ncol(DF))], collapse = "+"))),
                     data = DF[event == 0,])
        tmp0 <- survfit(tmp, newdata = base)
        hy <- with(tmp0, approx(time, cumhaz, t0)$y)
        hyU <- -log(with(tmp0, approx(time, upper, t0)$y))
        hyL <- -log(with(tmp0, approx(time, lower, t0)$y))
    }
    if (all(X == 0)) {
        tmp <- coxph(Surv(T0, Time, event) ~ 1, data = DF)
        ly <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        lyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        lyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
        tmp <- coxph(Surv(Time, status) ~ 1, data = DF[event == 0,])
        hy <- with(basehaz(tmp), approx(time, hazard, t0)$y)
        hyU <- -log(with(survfit(tmp), approx(time, upper, t0)$y))
        hyL <- -log(with(survfit(tmp), approx(time, lower, t0)$y))
    }
    list(t0 = t0, lam = ly, lamU = lyU, lamL = lyL,
         haz = hy, hazU = hyU, hazL = hyL)
}

## doNonpara.SE.cox.NA <- function(DF, alpha, beta, engine, stdErr) {
##     ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
##     id <- DF$id
##     B <- stdErr@B
##     cluster <- unlist(lapply(split(id, id), function(z) 1:length(z)))
##     clsz <- unlist(lapply(split(id, id), length))
##     Y <- rep(DF$Time[cumsum(clsz)], clsz)
##     T <- DF$Time[-cumsum(clsz)]
##     nT <- length(T)
##     nY <- length(clsz)
##     tmp <- basehaz(coxph(Surv(T, rep(1, nT)) ~ 1))
##     t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
##     ng <- length(t0)
##     ly <- with(tmp, approx(time, hazard, t0)$y)
##     tmp <- basehaz(coxph(Surv(Time, event) ~ 1, data = DF[cumsum(clsz),]))    
##     hy <- with(tmp, approx(time, hazard, t0)$y)
##     bootH <- bootL <- matrix(NA, length(t0), B) 
##     for (i in 1:B) {
##         bootL[,i] <- with(basehaz(coxph(Surv(T[sample(1:nT, nT, TRUE)], rep(1, nT)) ~ 1)),
##                           approx(time, hazard, t0, yleft = min(hazard), yright = max(hazard))$y)
##         bootH[,i] <- with(basehaz(coxph(Surv(Time, event) ~ 1,
##                                         data = DF[cumsum(clsz),][sample(1:nY, nY, TRUE),])),
##                           approx(time, hazard, t0, yleft = min(hazard), yright = max(hazard))$y)
##     }
##     lyU <- apply(bootL, 1, function(z) quantile(z, .975))
##     lyL <- apply(bootL, 1, function(z) quantile(z, .025))
##     hyU <- apply(bootH, 1, function(z) quantile(z, .975))
##     hyL <- apply(bootH, 1, function(z) quantile(z, .025))
##     list(t0 = t0, ly = ly, lyU = lyU, lyL = lyL, hy = hy, hyU = hyU, hyL = hyL)
## }

doNonpara.cox.HW <- function(DF, alpha, beta, formula, engine, stdErr) {
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    T <- DF$Time
    X <- as.matrix(DF[,-c(1:4)])
    event <- DF$event
    status <- DF$status
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    if (all(X == 0)) alpha <- beta <- 0
    delta <- DF$event
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y)
    Ta <- log(T)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha)))
    ly <- ly / max(ly)
    win.ly <- max(ly)
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    muZ <- mean(zHat)
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb),
                                         exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                                         status[event == 0]))
    win.hy <- max(hy)
    list(t0 = t0, lam = ly * muZ, lamU = rep(NA, ng), lamL = rep(NA, ng), 
         haz = hy, hazU = rep(NA, ng), hazL = rep(NA, ng))
}

doNonpara.SE.am.XCHWY <- function(DF, alpha, beta, formula, engine, stdErr) {
    B <- stdErr@B
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    t0 <- sort(unique(T, Y))
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    ng <- length(t0)
    Ya <- log(Y) + X %*% alpha
    Ta <- log(T) + X %*% alpha
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <-  as.numeric(mt * max(ly) / lambda)
    ly <- ly / max(ly)
    E <- matrix(rexp(length(t0) * B), nrow = length(t0))
    lytmp <- apply(E, 2, function(x) npMLE(t0, exp(Ta), exp(Ya), x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb), zHat / muZ, status[event == 0]))
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z)
        sapply(t0, function(y)
            baseHaz(y, exp(Yb), zHat / muZ, status[event == 0], z)))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(t0 = t0, lam = ly * muZ, lamU = lyU * muZ, lamL = lyL * muZ, haz = hy, hazU = hyU, hazL = hyL)
}

doNonpara.SE.cox.HW <- function(DF, alpha, beta, formula, engine, stdErr) {
    B <- stdErr@B
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    id <- DF$id
    event <- DF$event
    status <- DF$status
    X <- as.matrix(DF[,-c(1:4)])    
    n <- length(unique(id))
    p <- ncol(X)
    T <- DF$Time
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$Time[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    if (all(X == 0)) alpha <- beta <- 0
    delta <- DF$event
    ## t0 <- seq(0, max(Y), length.out = 5 * nrow(DF))
    t0 <- sort(unique(T, Y))
    ng <- length(t0)
    Ya <- log(Y)
    Ta <- log(T)
    lambda <- npMLE(Ya[event == 0], Ta, Ya)
    ly <- npMLE(t0, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / (lambda * exp(as.matrix(X[event == 0,]) %*% alpha)))
    ly <- ly / max(ly)
    E <- matrix(rexp(length(t0) * B), nrow = length(t0))
    lytmp <- apply(E, 2, function(x) npMLE(t0, exp(Ta), exp(Ya), x))
    lytmp <- apply(lytmp, 2, function(z) z / max(z))
    lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
    lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) ## + X %*% beta
    Yb <- Yb[event == 0]
    muZ <- mean(zHat)
    hy <- sapply(t0, function(z) baseHaz(z, exp(Yb),
                                         exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                                         status[event == 0]))
    E <- matrix(rexp(n * B), nrow = n)
    hytmp <- apply(E, 2, function(z)
        sapply(t0, function(y)
            baseHaz(y, exp(Yb), exp(as.matrix(X[event == 0,]) %*% beta) * zHat / muZ,
                    status[event == 0], z)))
    hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
    hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
    list(t0 = t0, lam = ly * muZ, lamU = lyU * muZ, lamL = lyL * muZ,
         haz = hy, hazU = hyU, hazL = hyL)
}

##############################################################################
# Class Definition
##############################################################################

setClass("Engine", representation(tol = "numeric"),
         prototype(tol = 1e-7), contains="VIRTUAL")

setClass("cox.LWYY", contains = "Engine")
setClass("cox.HW", contains = "Engine")
setClass("am.XCHWY", contains = "Engine")
setClass("am.GL", contains = "Engine")
setClass("sc.XCYH", contains = "Engine")

setClass("stdErr")
setClass("bootstrap", representation(B = "numeric"),
         prototype(B = 100), contains="stdErr")
setClass("resampling", representation(B = "numeric"),
         prototype(B = 100), contains="stdErr")


##############################################################################
# Method Dispatch
##############################################################################
setGeneric("doREFit", function(DF, DF0, engine, stdErr) {standardGeneric("doREFit")})

setMethod("doREFit", signature(engine = "cox.LWYY", stdErr = "NULL"), doREFit.cox.LWYY)
setMethod("doREFit", signature(engine = "cox.HW", stdErr = "NULL"), doREFit.cox.HW)
setMethod("doREFit", signature(engine = "am.XCHWY", stdErr = "NULL"), doREFit.am.XCHWY)
setMethod("doREFit", signature(engine = "am.GL", stdErr = "NULL"), doREFit.am.GL)
setMethod("doREFit", signature(engine = "sc.XCYH", stdErr = "NULL"), doREFit.sc.XCYH)

setMethod("doREFit", signature(engine="Engine", stdErr="bootstrap"),
          doREFit.Engine.Bootstrap)

setMethod("doREFit", signature(engine="am.XCHWY", stdErr="resampling"),
          doREFit.am.XCHWY.resampling)
setMethod("doREFit", signature(engine="sc.XCYH", stdErr="resampling"),
          doREFit.sc.XCYH.resampling)

setGeneric("doNonpara", function(DF, alpha, beta, formula, engine, stdErr) {standardGeneric("doNonpara")})
setMethod("doNonpara", signature(engine = "cox.LWYY", stdErr = "NULL"), doNonpara.cox.NA)
setMethod("doNonpara", signature(engine = "cox.HW", stdErr = "NULL"), doNonpara.cox.HW)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "NULL"), doNonpara.am.XCHWY)

setMethod("doNonpara", signature(engine = "cox.LWYY", stdErr = "bootstrap"), doNonpara.cox.NA)
setMethod("doNonpara", signature(engine = "cox.HW", stdErr = "bootstrap"), doNonpara.SE.cox.HW)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "resampling"), doNonpara.SE.am.XCHWY)
setMethod("doNonpara", signature(engine = "am.XCHWY", stdErr = "bootstrap"), doNonpara.SE.am.XCHWY)

## GL method?
setMethod("doNonpara", signature(engine = "am.GL", stdErr = "bootstrap"), doNonpara.SE.am.XCHWY)
setMethod("doNonpara", signature(engine = "am.GL", stdErr = "NULL"), doNonpara.am.XCHWY)

## general model
setMethod("doNonpara", signature(engine = "sc.XCYH", stdErr = "bootstrap"), doNonpara.SE.am.XCHWY)
setMethod("doNonpara", signature(engine = "sc.XCYH", stdErr = "NULL"), doNonpara.SE.am.XCHWY)

##############################################################################
## User's Main Function
##############################################################################

reReg <- function(formula, data, B = 200, 
                  method = c("cox.LWYY", "cox.HW", "am.GL", "am.XCHWY", "sc.XCYH"),
                  se = c("NULL", "bootstrap", "resampling"), plot.ci = FALSE,
                  contrasts = NULL, control = list()) {
    method <- match.arg(method)
    se <- match.arg(se)
    Call <- match.call()
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    if (!is.reSurv(obj)) stop("Response must be a reSurv object")
    formula[[2]] <- NULL
    if (formula == ~ 1) {
        DF <- cbind(obj$reDF[, -5], zero=0)
    } else {
        ## remove intercept
        if (!missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, data))
        if (missing(data)) DF <- cbind(obj$reDF[,-5], model.matrix(formula, parent.frame()))
        DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- DF[order(DF$id, DF$Time), ]
    ## reset ID
    idTmp <- with(DF, as.numeric(unlist(lapply(split(id, id), length))))
    DF$id <- rep(1:length(idTmp), idTmp)
    DF0 <- as.matrix(ddply(DF, "id", head, n=1)[, -c(1:4)])
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class=method), engine.control))
    if (se == "NULL")
        stdErr <- NULL
    else {
        stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class=se), stdErr.control))
        stdErr@B <- B
    }
    if (formula == ~1) {
        fit <- NULL
        fit$alpha <- fit$beta <- rep(NA, ncol(DF) - 4)
        fit$muZ <- NA
        p <- ncol(DF) - 4
        if (plot.ci) {
            stdErr.np.control <- control[names(control) %in% names(attr(getClass("bootstrap"), "slots"))]
            stdErr.np <- do.call("new", c(list(Class = "bootstrap"), stdErr.np.control))
            stdErr.np@B <- B
            fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta, formula = formula, 
                                    engine = engine, stdErr = stdErr.np))
        } else {
            fit <- c(fit, doNonpara(DF = DF, alpha = 0, beta = 0, formula = formula,
                                    engine = engine, stdErr = stdErr))
        }
    } else {
        fit <- doREFit(DF = DF, DF0 = DF0, engine = engine, stdErr = stdErr)
        if (method != "sc.XCYH") {
            if (plot.ci) {
                stdErr.np.control <- control[names(control) %in% names(attr(getClass("bootstrap"), "slots"))]
                stdErr.np <- do.call("new", c(list(Class = "bootstrap"), stdErr.np.control))
                stdErr.np@B <- B
            fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta, formula = formula,
                                    engine = engine, stdErr = stdErr.np))
            } else {
                fit <- c(fit, doNonpara(DF = DF, alpha = fit$alpha, beta = fit$beta, formula = formula,
                                        engine = engine, stdErr = NULL))
            }
        }
            
    }
    class(fit) <- "reReg"
    fit$DF <- DF
    fit$call <- Call
    fit$varNames <- names(DF)[-(1:4)]
    fit$method <- method
    fit$se <- se
    fit
}

##############################################################################
##############################################################################

npMLE <- function(t, tij, yi, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(yi))
    ttmp <- tij[tij != yi]
    ord <- order(ttmp)
    sl <- unique(ttmp[ord])
    l <- ifelse(min(t) < max(sl), which(sl > min(t))[1], length(sl))
    ## res <- vector("double", 1)
    ## tmp <- sl[l:length(sl)]
    tmp <- sl[l:length(sl)]
    tmp <- rev(tmp)
    tij <- rev(tij)
    yi <- rev(yi)
    ## yi <- ifelse(is.infinite(yi), max(yi[!is.infinite(yi)]), yi)
    ## tij <- ifelse(is.infinite(tij), max(tij[!is.infinite(tij)]), tij)
    res <- vector("double", length(tmp)) + 1
    res <- .C("plLambda", as.double(tmp), as.double(tij), as.double(yi), as.double(weights), 
              as.integer(length(tmp)), as.integer(length(yi)),
              out = as.double(res), PACKAGE = "reReg")$out
    out <- rev(res)[sapply(t, function(x) which(rev(tmp) >= x)[1])]
    out <- ifelse(is.na(out), 0, out)
    out <- exp(-out)
}


baseHaz <- function(t, Y, zHat, delta, weights  = NULL) {
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    ind <- which(delta == 1 & Y <= t)
    temp2 <- tmp <- weights[order(Y)]
    ## temp2 <- c(tmp[1], diff(cumsum(tmp)))
    ## temp2[order(Y)] <- temp2
    temp2[order(Y)] <- tmp
    if (length(ind) > 0) {
        out <- sapply(ind, function(x) temp2[x] / sum(zHat * weights * (Y >= Y[x])))
    }
    if (length(ind) == 0)
        out <- 0
    sum(out)
}

alphaEq <- function(alpha, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    ## Ystar <- Y * exp(X %*% alpha)
    ## Tstar <- T * exp(X %*% alpha)
    Ystar <- log(Y) + X %*% alpha
    Tstar <- log(T) + X %*% alpha
    Lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    ## Lambda <- npMLE(Ystar[which(cluster == 1)], log(T), log(Y),
    ## weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("alphaEq", as.double(X[which(cluster == 1), ]), as.double(Lambda),
              as.integer(mt), as.integer(n), as.integer(p),
              ## as.integer(length(weights) %/% n), as.double(weights), 
              out = as.double(res), PACKAGE = "reReg")$out
    res / rep(n * unlist(lapply(split(weights, rep(1:(length(weights) %/% n), each = n)), sum)), each = p)
}

betaEq <- function(X, Y, T, cluster, mt, delta, zHat = NULL, alpha, beta, weights = NULL) {
    p <- ncol(X)
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(zHat)) {
        Ystar <- log(Y) + X %*% alpha
        Tstar <- log(T) + X %*% alpha
        lambda <- npMLE(Ystar[which(cluster == 1)], Tstar, Ystar,
                        weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
        zHat <- as.numeric(weights * mt / lambda)
        zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    }
    Y <- log(Y) + X %*% beta
    Y <- Y[which(cluster == 1)]
    X <- X[which(cluster == 1), ]
    ## delta <- delta[which(cluster == 1)]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("betaEst", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(weights), as.integer(n), as.integer(p),
              as.integer(length(weights) %/% n), 
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}

ghoshU2 <- function(alpha, beta, T, Y, X, cl) {
    ## dim(X) = n by p, dim(Y) = n by 1, dim(T) > n by 1
    d <- max(X %*% (alpha - beta), 0)
    TT <- log(T) - rep(X %*% alpha, cl)
    TY <- log(Y) - X %*% beta - d
    p <- ncol(X)
    .C("ghosh", as.double(TT), as.double(TY), as.double(X), as.integer(cl),
       as.integer(c(0, cumsum(cl)[-length(cl)])),
       as.integer(nrow(X)), as.integer(p), 
       out = as.double(double(p)), PACKAGE = "reReg")$out
}

LWYYeq <- function(beta, X, Y, T, cl) {
    p <- ncol(X)
    res <- vector("double", p)
    wgt <- exp(X %*% beta)
    .C("lwyy", as.double(T), as.double(Y), as.double(X), as.double(wgt), as.integer(cl),
       as.integer(c(0, cumsum(cl)[-length(cl)])), as.integer(nrow(X)), as.integer(p),        
       out = as.double(double(p)), PACKAGE = "reReg")$out       
}


##########################################################################################
## Paper 2: More general models
##########################################################################################

sarm <- function(X, Y, T, id, cluster, method, B = 200) {
    n <- sum(cluster == 1)
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    p <- ncol(X)
    alpha <- beta <- gamma <- rep(0, p)
    muZ <- NULL
    if (method == "M1") {
        gamma <- c(0, gamma)
        out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                       cluster = cluster, mt = mt, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        muZ <- out$par[1]
        alpha <- out$par[-1]
    }
    if (method == "M3") {
        alpha <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                         quiet = TRUE, control = list(M = c(1, 10)))$par
        gamma <- c(0, gamma)
        if (alpha %*% alpha > 100) {
            beta <- c(0, alpha)
        }
        else {
            out <- BBsolve(gamma, coefEq, alpha = alpha, X = X, Y = Y, T = T,
                           cluster = cluster, mt = mt, weights = NULL,
                           quiet = TRUE, control = list(M = c(1, 10)))
            muZ <- out$par[1]
            beta <- out$par[-1] + alpha
        }
    }
    if (method == "M2") {
        ## gamma <- c(0, gamma)
        ## out <- BBsolve(gamma, coefEq, alpha = NULL, X = X, Y = Y, T = T,
        ##                cluster = cluster, mt = mt, weights = NULL,
        ##                quiet = TRUE, control = list(M = c(1, 10)))
        ## muZ <- out$par[1]
        ## alpha <- out$par[-1]
        ## out <- dfsane(alpha, M1eq, X = X, Y = Y, T = T,
        ##                 cluster = cluster, weights = NULL,
        ##                 control = list(NM = TRUE, M = 1, noimp = 50, trace = FALSE))
        out <- BBsolve(alpha, M1eq, X = X, Y = Y, T = T, cluster = cluster, weights = NULL,
                       quiet = TRUE, control = list(M = c(1, 10)))
        alpha <- out$par
    }
    list(alpha = alpha, beta = beta, muZ = muZ)
}

HWeq <-function(gamma, X, Y, T, cluster, mt) {
    n <- sum(cluster == 1)
    Lambda <- npMLE(Y[cluster == 1], T, Y)
    res <- vector("double", length(gamma))
    p <- ncol(X)
    .C("sarm1", as.double(X), as.double(Lambda), as.double(rep(1, n)),
       as.double(gamma), as.integer(mt), as.integer(n), as.integer(p), as.integer(1),
       out = as.double(rep(0, p)), PACKAGE = "reReg")$out                            
}


HWeq2 <-function(beta, X, Y, delta, zHat) {
    n <- nrow(X)
    p <- ncol(X)
    res <- vector("double", p)
    res <- .C("HWb", as.double(Y), as.double(X), as.double(delta), as.double(zHat),
              as.double(X %*% beta), as.integer(n), as.integer(p), as.integer(1),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n
}
    
coefEq <- function(alpha, gamma, X, Y, T, cluster, mt, weights = NULL) {
    n <- sum(cluster == 1)
    if (is.null(weights))
        weights <- rep(1, n)
    if (is.null(alpha))
        alpha <- -1 * gamma[-1]
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha    
    Lambda <- npMLE(Ytmp[cluster == 1], Ttmp, Ytmp,
                    weights = rep(weights, diff(c(which(cluster ==1), length(cluster)+1))))
    Lambda <- Lambda / max(Lambda)
    X <- X[cluster == 1,]
    p <- ncol(X)
    res <- vector("double", (p + 1) * length(weights) %/% n)
    res <- .C("sarm1", as.double(cbind(1, X)), as.double(Lambda), as.double(weights),
              as.double(gamma), as.integer(mt), as.integer(n), as.integer(p+1),
              as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res / n    
}


## Martingal approach
M1eq <- function(alpha, X, Y, T, cluster, weights = NULL) {
    n <- sum(cluster == 1)
    ## mi <- unlist(lapply(split(cluster, id), length)) - 1
    if (is.null(weights))
        weights <- rep(1, n)
    p <- ncol(X)
    Ytmp <- log(Y) + X %*% alpha
    Ttmp <- log(T) + X %*% alpha
    ## Ytmp <- Y * exp(X %*% alpha)
    ## Ttmp <- T * exp(X %*% alpha)
    ind <- which(T != Y)
    Ttmp <- Ttmp[ind]
    Ytmp <- Ytmp[ind]
    X <- X[ind,]
    res <- vector("double", p * length(weights) %/% n)
    res <- .C("sarm2", as.double(X), as.double(Ttmp), as.double(Ytmp), as.double(weights), 
              as.integer(length(Ttmp)), as.integer(p), as.integer(length(weights) %/% n),
              out = as.double(res), PACKAGE = "reReg")$out
    res
}

## Method 1: Z\lambda(t)e^xa  CY's 2004 JASA
## Method 2: Z\lambda(te^xa)
## Method 3: Z\lambda(te^xa)e^xb
## Method 4: Z\lambda(te^xa)e^xa

#########################################################
## General modes in R, need to move this to C sometimes
#########################################################

sarmRV <- function(id, Tij, Yi, X, M, lamEva = NULL, aIni = NULL, bIni = NULL) {
    n <- length(unique(id))
    X <- as.matrix(X)
    p <- ncol(X)
    mm <- matrix((M > 0), nrow(X), ncol(X))
    U1RV <- function(a) {
        tx <- as.vector(log(Tij) + X %*% a)
        yx <- as.vector(log(Yi) + X %*% a)
        tx2 <-  outer(tx, tx,">=")
        txy <-  outer(tx, yx, "<=")
        A <- (tx2 * txy) %*% (X * mm)
        B <- (tx2 * txy) %*% mm
        B[B == 0] <- 1e10
        apply((X - A / B) * mm, 2, sum)
    }
    if (is.null(aIni)) ahat <- dfsane(par = double(p), fn = U1RV, quiet = TRUE, alertConvergence = FALSE,
                                      control = list(NM = TRUE, trace = FALSE))$par
    if (!is.null(aIni)) ahat <- dfsane(par = aIni, fn = U1RV, quiet = TRUE, alertConvergence = FALSE,
                                       control = list(NM = TRUE, trace = FALSE))$par
    ## if (is.null(aIni)) ahat <- BBsolve(par = double(p), fn = U1RV, quiet = TRUE)$par
    ## if (!is.null(aIni)) ahat <- BBsolve(par = aIni, fn = U1RV, quiet = TRUE)$par
    tx <- as.vector(log(Tij) + X %*% ahat)
    yx <- as.vector(log(Yi) + X %*% ahat)
    tx2 <-  outer(tx,tx,">=")
    txy <-  outer(tx, yx, "<=")
    vv <- matrix((M > 0), nrow(X), n)
    yx0 <- as.numeric(unlist(lapply(split(yx, id), unique)))
    txy0 <- outer(tx, yx0, ">=")
    Rn <- (tx2 * txy) %*% (M > 0)
    Rn[Rn == 0] <- 1e10
    ## Lam <- apply(1 - (txy0 * vv) / matrix(Rn, nrow(X), n), 2, prod)
    Lam <- exp(-colSums((txy0 * vv) / matrix(Rn, nrow(X), n)))
    Lam[Lam == 0] <- 1e10
    ind <- cumsum(unlist(lapply(split(id, id), length)))
    U2RV <- function(b) {
        tmp <- ifelse(M[ind] / Lam > 1e5, (M[ind] + .01) / (Lam + .01), M[ind] / Lam)
        as.numeric(t(cbind(1, X[ind,])) %*% (tmp - exp(cbind(1, X[ind,]) %*% b))) / n
        ## as.numeric(t(cbind(1, X[ind,])) %*% (M[ind] / Lam - exp(cbind(1, X[ind,]) %*% b))) / n
    }
    ## if (is.null(bIni)) ghat <- BBsolve(par = double(p + 1), fn = U2RV, quiet = TRUE)$par
    ## if (!is.null(bIni)) ghat <- BBsolve(par = bIni, fn = U2RV, quiet = TRUE)$par
    if (is.null(bIni))
        ghat <- dfsane(par = double(p + 1), fn = U2RV, alertConvergence = FALSE,
                       control = list(NM = TRUE, trace = FALSE))$par
    if (!is.null(bIni))
        ghat <- dfsane(par = bIni, fn = U2RV, alertConvergence = FALSE,
                       control = list(NM = TRUE, trace = FALSE))$par
    list(ahat = ahat, bhat = ghat[-1] + ahat, ghat = ghat, LamTau = ghat[1],
         lamY = approx(x = yx[ind], y = Lam, xout = lamEva, "constant",
                       yright = max(Lam), yleft = min(Lam))$y)
}


sarmRV.sand <- function(id, Tij, Yi, X, M, a = NULL, b = NULL, Bootstrap = 200) {
  n <- length(unique(id))
  X <- as.matrix(X)
  p <- ncol(X)
  mm <- matrix((M > 0), nrow(X), ncol(X))
  clut <- as.numeric(unlist(lapply(split(id, id), length)))
  tmpE <- matrix(rexp(n * Bootstrap), ncol = n)
  tmpN <- matrix(rnorm((2 * p + 1) * Bootstrap), ncol = 2 * p + 1)
  Sn <- function(a, b, e) {
    ## Part S1
    e1 <- rep(e, clut)
    tx <- as.vector(Tij * exp(X %*% a))
    yx <- as.vector(Yi * exp(X %*% a))
    ee1 <- matrix(e1, nrow(X), ncol(X))
    tx2 <-  outer(tx, tx,">=")
    txy <-  outer(tx, yx, "<=")
    A <- (tx2 * txy) %*% (X * mm * ee1)
    B <- (tx2 * txy) %*% (mm * ee1)
    B[B == 0] <- 1e15
    s1 <- apply((X - A / B) * (mm * ee1), 2, sum) / n
    ## lambda
    vv <- matrix((M > 0), nrow(X), n)
    yx0 <- as.numeric(unlist(lapply(split(yx, id), unique)))
    txy0 <- outer(tx, yx0, ">=")
    Rn <- (tx2 * txy) %*% (e1 * (M > 0))
    Rn[Rn == 0] <- 1e15
    ## Lam <- apply(1 - (txy0 * vv * e1) / matrix(Rn, nrow(X), n), 2, prod)
    Lam <- exp(-colSums((txy0 * vv * e1) / matrix(Rn, nrow(X), n)))
    Lam[Lam == 0] <- 1e15
    ind <- cumsum(unlist(lapply(split(id, id), length)))
    ee2 <- matrix(e, nrow(X[ind,]), ncol(X) + 1)
    s2 <- as.numeric(t(cbind(1, X[ind,]) * ee2) %*% (M[ind] / Lam - exp(cbind(1, X[ind,]) %*% b))) / n
    ## list(s1 = s1, s2 = s2)
    return(c(s1, s2))
  }
  V <- var(t(apply(tmpE, 1, function(x) Sn(a, b, x)))) ## / sqrt(n)
  tmp <- t(apply(tmpN, 1, function(x)
                 sqrt(n) * Sn(a + x[1:p] / sqrt(n), b + x[-(1:p)] / sqrt(n), rep(1, n))))
  J0 <- t(coef(lm(tmp[,1:p] ~ tmpN[,1:p] - 1)))
  Jtmp <- t(coef(lm(tmp[,-c(1:p)] ~ tmpN - 1)))
  ## J0 <- coef(lm(tmp[,1:p] ~ tmpN[,1:p] - 1))
  ## Jtmp <- coef(lm(tmp[,-c(1:p)] ~ tmpN - 1))
  J1 <- Jtmp[,1:p]
  J2 <- Jtmp[,-c(1:p)]
  J <- rbind(cbind(J0, matrix(0, ncol = p + 1, nrow = nrow(J0))), cbind(Jtmp))
  if (qr(J)$rank == nrow(J)) J <- solve(J) else J <- ginv(J)
  if (qr(J0)$rank == nrow(J0)) J0 <- solve(J0) else J0 <- ginv(J0)
  ase <- J %*% V %*% t(J)
  list(ase = ase, J = J, V = V,
       alphaSE = sqrt(diag(J0 %*% V[1:p, 1:p] %*% t(J0))),## sqrt(diag(ase)[1:p]),
       betaSE = sqrt(diag(ase[1:p, 1:p] + ase[(p + 2):(2 * p + 1), (p + 2):(2 * p + 1)] +
         2 * ase[1:p, (p + 2):(2 * p + 1)])))
}
