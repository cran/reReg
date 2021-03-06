## globalVariables("DF") ## global variables for reReg

##############################################################################
## Functions for different models
## stdErr is estimated with resampling if method = gsc or am.xc,
##        bootstrap otherwise
##############################################################################

regFit.am.GL <- function(DF, engine, stdErr) {
    DF0 <- DF[DF$event == 0,]
    p <- ncol(DF0) - 6
    Y <- log(DF0$time2)
    X <- as.matrix(DF0[,-(1:6)])
    status <- DF0$terminal
    n <- nrow(DF0)
    log.est <- function(b) {
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(status),
           as.integer(rep(1, n)), as.integer(n), as.integer(p), as.integer(n),
           as.double(rep(1, n)), as.double(rep(1, n)), 
           result = double(p), PACKAGE = "reReg")$result
    }
    fit.b <- eqSolve(engine@par2, log.est, engine@solver)
    bhat <- fit.b$par
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    ghoshU2 <- function(a) {
        d <- max(X %*% (a - bhat))
        tij <- log(DF$time2) - as.matrix(DF[,-(1:6)]) %*% a
        tij <- tij[DF$event == 1]
        ## yi <- Y - X %*% a - d ## Correct version
        yi <- Y - X %*% bhat - d ## Paper version
        if (sum(tij < rep(yi, m)) == 0) return(1e5)
        else
            .C("glU2", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
               as.double(yi), as.double(tij), as.double(X), as.double(rep(1, n)), result = double(p),
               PACKAGE = "reReg")$result
    }
    fit.a <- eqSolve(engine@par1, ghoshU2, engine@solver)
    fit.b$par <- -fit.b$par
    fit.a$par <- -fit.a$par
    out <- list(par1 = fit.a$par, par1.conv = fit.a$convergence,
                par3 = fit.b$par, par3.conv = fit.b$convergence, muZ = NA)
    out$typeRec <- engine@typeRec
    out$typeTem <- engine@typeTem
    return(out)
}

regFit.am.GL.sand <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    DF0 <- DF[DF$event == 0,]
    p <- ncol(DF0) - 6
    Y <- log(DF0$time2)
    X <- as.matrix(DF0[,-(1:6)])
    status <- DF0$terminal
    n <- nrow(DF0)
    m <- aggregate(event ~ id, data = DF, sum)[,2]
    index <- c(1, cumsum(m)[-n] + 1)
    B <- stdErr@B
    E <- matrix(rexp(n * B), nrow = n)
    Z <- matrix(rnorm(p * B), nrow = p)
    Sn <- function(a, b, w, r = "both") {
        s2 <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(status),
                 as.integer(rep(1, n)), as.integer(n), as.integer(p), as.integer(n),
                 as.double(w), as.double(rep(1, n)), 
                 result = double(p), PACKAGE = "reReg")$result
        d <- max(X %*% (a - b))
        tij <- log(DF$time2) - as.matrix(DF[,-(1:6)]) %*% a
        tij <- tij[DF$event == 1]
        ## yi <- Y - X %*% a - d ## Correct version
        yi <- Y - X %*% b - d ## Paper version
        if (sum(tij < rep(yi, m)) == 0) s1 <- NA
        else s1 <- .C("glU2", as.integer(n), as.integer(p), as.integer(index - 1), as.integer(m),
                      as.double(yi), as.double(tij), as.double(X), as.double(w), result = double(p),
                      PACKAGE = "reReg")$result
        if (r == "s1") return(s1)
        if (r == "s2") return(s2)
        return(c(s1, s2))
    }
    V <- var(t(apply(E, 2, function(x) Sn(-res$par1, -res$par3, x))))
    V1 <- V[1:p, 1:p]
    V2 <- V[1:p + p, 1:p + p]
    lmfit1 <- t(apply(Z, 2, function(x) Sn(-res$par1 + x / sqrt(n), -res$par3, rep(1, n), "s1")))
    lmfit2 <- t(apply(Z, 2, function(x) Sn(-res$par1, -res$par3 + x / sqrt(n), rep(1, n), "s2")))
    if (p == 1) {
        J1 <- coef(lm(sqrt(n) * c(lmfit1) ~ c(Z)))[-1]
        J2 <- coef(lm(sqrt(n) * c(lmfit2) ~ c(Z)))[-1]
    } else {        
        J1 <- coef(lm(sqrt(n) * lmfit1 ~ t(Z)))[-1,]
        J2 <- coef(lm(sqrt(n) * lmfit2 ~ t(Z)))[-1,]
    }
    if (qr(J1)$rank == p) aVar <- solve(J1) %*% V1 %*% t(solve(J1))
    else aVar <- ginv(J1) %*% V1 %*% t(ginv(J1))
    if (qr(J2)$rank == p) bVar <- solve(J2) %*% V2 %*% t(solve(J2))
    else bVar <- ginv(J2) %*% V2 %*% t(ginv(J2))    
    aSE <- sqrt(diag(aVar))
    bSE <- sqrt(diag(bVar))
    out <- c(res, list(par1.se = aSE, par3.se = bSE, par1.vcov = aVar, par3.vcov = bVar))
    out$typeRec <- engine@typeRec
    out$typeTem <- engine@typeTem
    return(out)
}

#' @importFrom survival cluster
#' @importFrom stats vcov
regFit.cox.LWYY <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:6)])
    p <- ncol(X)
    T <- DF$time2
    T0 <- unlist(lapply(split(T, id), function(x) c(0, x[-length(x)])))
    fit.coxph <- coxph(Surv(T0, T, event) ~ X + cluster(id))
    out <- list(par1 = coef(fit.coxph), par1.se = sqrt(diag(vcov(fit.coxph))))
    out$typeRec <- engine@typeRec
    out$typeTem <- engine@typeTem
    return(out)
}

#' This is also the ARF in Luo
#' @importFrom survival basehaz
#' @noRd
regFit.cox.GL <- function(DF, engine, stdErr) {
    id <- DF$id
    event <- DF$event
    X <- as.matrix(DF[,-c(1:6)])
    p <- ncol(X)
    T <- DF$time2
    mt <- aggregate(event ~ id, data = DF, sum)$event
    Y <- rep(DF$time2[event == 0], mt + 1)
    X0 <- X[event == 0,,drop = FALSE]
    fit.coxph <- coxph(Surv(T[event == 0], DF$terminal[event == 0]) ~ X0)
    cumHaz <- basehaz(fit.coxph)
    ## cumHaz$hazard <- cumHaz$hazard / max(cumHaz$hazard)
    wgt <- sapply(exp(X0 %*% coef(fit.coxph)), function(x)
        approxfun(cumHaz$time, exp(-cumHaz$hazard * x), yleft = 1,
                  yright = min(exp(-cumHaz$hazard * x)),
                  method = "constant")(T))
    wgt <- 1 / wgt ## ifelse(wgt == 0, 1 / sort(c(wgt))[2], 1 / wgt)
    wgt <- ifelse(wgt > 1e5, 1e5, wgt)
    out <- dfsane(par = engine@par1, fn = coxGLeq, wgt = wgt, 
                  X = as.matrix(X[!event, ]),
                  Y = Y[!event], T = ifelse(T == Y, 1e5, T), cl = mt + 1,
                  alertConvergence = FALSE, quiet = TRUE, control = list(trace = FALSE))
    out <- list(par1 = out$par,
                par3 = coef(fit.coxph),
                par3.se = sqrt(diag(vcov(fit.coxph))))
    out$typeRec <- engine@typeRec
    out$typeTem <- engine@typeTem
    return(out)}

## #' @importFrom rlang is_empty
regFit.general <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    out <- s1(engine@typeRec, DF, engine@eqType, engine@solver, engine@par1, engine@par2)
    if (engine@typeTem != ".") 
        out <- c(out, s2(engine@typeTem, DF, engine@eqType, engine@solver,
                         engine@par3, engine@par4, out$zi))
    out$typeRec <- engine@typeRec
    out$typeTem <- engine@typeTem
    return(out)
}

s1 <- function(type, DF, eqType, solver, par1, par2, Lam0 = NULL, w1 = NULL, w2 = NULL) {
    if (type == "gsc") return(reSC(DF, eqType, solver, par1, par2, Lam0, w1, w2))
    if (type == "cox") return(reCox(DF, eqType, solver, par1, Lam0, w1))
    if (type == "am") return(reAM(DF, eqType, solver, par1, Lam0, w1))
    if (type == "ar") return(reAR(DF, eqType, solver, par1, Lam0, w1))
    return(NULL)
}

s2 <- function(type, DF, eqType, solver, par3, par4, zi, wgt = NULL) {
    if (type == "gsc") return(temSC(DF, eqType, solver, par3, par4, zi, wgt))
    if (type == "cox") return(temCox(DF, eqType, solver, par3, zi, wgt))
    if (type == "am") return(temAM(DF, eqType, solver, par3, zi, wgt))
    if (type == "ar") return(temAR(DF, eqType, solver, par3, zi, wgt))
    return(NULL)
}

regFit.general.sand <- function(DF, engine, stdErr) {
    if (is.na(match(engine@solver, c("dfsane", "BBsolve", "optim", "BBoptim")))) {
        print("Warning: Unidentified solver; BB::dfsane is used.")
        engine@solver <- "dfsane"
    }
    res <- regFit(DF, engine, NULL)
    n <- length(unique(DF$id))
    B <- stdErr@B
    p <- ncol(DF) - 6
    tmpV <- replicate(B,
                      c(s1(engine@typeRec, DF, engine@eqType, NULL, res$par1, res$par2,
                           res$Lam0, rexp(n), rexp(n))$value,
                        s2(engine@typeTem, DF, engine@eqType, NULL, res$par3, res$par4,
                           rexp(n) * res$zi, rexp(n))))    
    V <- var(t(tmpV))
    Z <- matrix(rnorm(ncol(V) * B), B)
    len1 <- length(res$par1)
    len2 <- length(res$par2)
    len3 <- length(res$par3)
    len4 <- length(res$par4)
    na <- len1 + len2
    nb <- len3 + len4
    L <- apply(Z, 1, function(zz) {
    c(s1(engine@typeRec, DF, engine@eqType, NULL,
         res$par1 + zz[1:len1] / sqrt(n), res$par2 + zz[1:len2 + len1] / sqrt(n), res$Lam0)$value,
      s2(engine@typeTem, DF, engine@eqType, NULL,
         res$par3 + zz[1:len3 + len1 + len2] / sqrt(n),
         res$par4 + zz[1:len4 + len1 + len2 + len3] / sqrt(n), res$zi))})
    L <- t(L)
    J <- solve(t(Z) %*% Z) %*% t(Z) %*% (sqrt(n) * L)
    recVar <- solve(J[1:na, 1:na]) %*% V[1:na, 1:na] %*% t(solve(J[1:na, 1:na]))
    par1.vcov <- recVar[1:len1, 1:len1]
    par1.se <- sqrt(diag(par1.vcov))
    res <- c(res, list(par1.vcov = par1.vcov, par1.se = par1.se))
    if (len2 > 0) {
        par2.vcov <- recVar[1:len2 + len1, 1:len2 + len1]
        par2.se <- sqrt(diag(par2.vcov))
        res <- c(res, list(par2.vcov = par2.vcov, par2.se = par2.se))
    }
    if (nb > 0) {
        ind2 <- tail(1:nrow(J), nb)
        J2 <- solve(t(Z[,ind2]) %*% Z[,ind2]) %*% t(Z[,ind2]) %*% (sqrt(n) * L[,ind2])
        temVar <- solve(J2) %*% V[ind2, ind2] %*% t(solve(J2))
        par3.vcov <- temVar[1:len3, 1:len3]
        par3.se <- sqrt(diag(par3.vcov))
        res <- c(res, list(par3.vcov = par3.vcov, par3.se = par3.se))
        if (!is.null(res$par4)) {
            par4.vcov <- temVar[1:len4 + len3, 1:len4 + len3]
            par4.se <- sqrt(diag(par4.vcov))
            res <- c(res, list(par4.vcov = par4.vcov, par4.se = par4.se))
        }
        res$vcovTem <- temVar
    }
    res$vcovRec <- recVar
    return(res)
}

##############################################################################
# Variance estimation 
##############################################################################
regFit.Engine.boot <- function(DF, engine, stdErr) {
    res <- regFit(DF, engine, NULL)
    id <- DF$id
    event <- DF$event
    status <- DF$terminal
    X <- as.matrix(DF[,-c(1:6)])    
    n <- length(unique(id))
    T <- DF$time2
    mt <- aggregate(event ~ id, data = DF, sum)$event
    clsz <- mt + 1
    Y <- rep(DF$time2[event == 0], mt + 1)
    cluster <- unlist(sapply(mt + 1, function(x) 1:x))
    B <- stdErr@B
    uID <- unique(DF$id) # unique(DF$ID)
    bound <- c(res$par1, res$par2, res$par3, res$par4)
    if (stdErr@parallel) {
        cl <- makeCluster(stdErr@parCl)
        clusterExport(cl = cl,
                      varlist = c("DF", "engine"),
                      envir = environment())
        out <- parSapply(cl, 1:B, function(x) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            tmp <- regFit(DF2, engine, NULL)
            return(c(tmp$par1, tmp$par2, tmp$par3, tmp$par4))
        })
        stopCluster(cl)
        bCoef <- t(out)
        convergence <- apply(bCoef, 1, function(x) 1 * (x %*% x > 1e3 * bound %*% bound))
    } else {
        bCoef <- matrix(0, B, length(bound))
        convergence <- rep(0, B)
            for (i in 1:B) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            tmp <- regFit(DF2, engine, NULL)
            bCoef[i,] <- c(tmp$par1, tmp$par2, tmp$par3, tmp$par4)
            convergence[i] <- 1 * (bCoef[i,] %*% bCoef[i,] > 1e3 * bound %*% bound)
        }
    }
    converged <- which(convergence == 0)
    if (sum(convergence != 0) > 0) {
        print("Warning: Some bootstrap samples failed to converge")
        tmp <- apply(bCoef, 1, function(x) x %*% x)
        converged <- (1:B)[- which(tmp %in% boxplot(tmp, plot = FALSE)$out)]        
    }
    if (all(convergence != 0) || sum(convergence == 0) == 1) {
        print("Warning: some bootstrap samples failed to converge")
        converged <- 1:B
    }
    bVar <- var(bCoef[converged, ], na.rm = TRUE)
    bSE <- sqrt(diag(as.matrix(bVar)))
    len1 <- length(res$par1)
    len2 <- length(res$par2)
    len3 <- length(res$par3)
    len4 <- length(res$par4)
    res <- c(res, list(par1.vcov = bVar[1:len1, 1:len1], par1.se = bSE[1:len1], B = length(converged)))
    if (len2 > 0)
        res <- c(res, list(par2.vcov = bVar[1:len2 + len1, 1:len2 + len1],
                           par2.se = bSE[1:len2 + len1]))
    if (len3 > 0)
        res <- c(res, list(par3.vcov = bVar[1:len3 + len1 + len2, 1:len3 + len1 + len2],
                           par3.se = bSE[1:len3 + len1 + len2]))
    if (len4 > 0)
        res <- c(res, list(par4.vcov = bVar[1:len4 + len1 + len2 + len3, 1:len4 + len1 + len2 + len3],
                           par4.se = bSE[1:len4 + len1 + len2 + len3]))
    if (engine@typeRec == "gsc") {
        res$par2.vcov <- bVar[1:len1, 1:len1] + bVar[2:len2 + len1, 2:len2 + len1] + 2 * bVar[1:len1, 2:len2 + len1]
        res$par2.se <- sqrt(diag(res$par2.vcov))
    }
    return(res)
}

##############################################################################
# Class Definition
##############################################################################

setClass("Engine",
         representation(tol = "numeric",
                        par1 = "numeric", par2 = "numeric",
                        par3 = "numeric", par4 = "numeric",
                        baseSE = "logical", 
                        solver = "character", eqType = "character", 
                        typeRec = "character", typeTem = "character"),
         prototype(eqType = "logrank", tol = 1e-7,
                   par1 = 0, par2 = 0, par3 = 0, par4 = 0, 
                   baseSE = FALSE, solver = "dfsane"),
         contains = "VIRTUAL")
setClass("general", contains = "Engine")
setClass("cox.LWYY", contains = "Engine")
setClass("cox.HW", contains = "Engine")
setClass("am.XCHWY", contains = "Engine")
setClass("am.GL", contains = "Engine")
setClass("gsc.XCYH", representation(muZ = "numeric"),
         prototype(muZ = 0), contains = "Engine")
setClass("cox.GL",
         representation(wgt = "matrix"), prototype(wgt = matrix(0)), contains = "Engine")

setClass("stdErr",
         representation(B = "numeric", parallel = "logical", parCl = "numeric"),
         prototype(B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2L),
         contains = "VIRTUAL")

setClass("boot", contains = "stdErr")
setClass("sand", contains = "stdErr")


##############################################################################
# Method Dispatch
##############################################################################
setGeneric("regFit", function(DF, engine, stdErr) {standardGeneric("regFit")})

setMethod("regFit", signature(engine = "general", stdErr = "NULL"), regFit.general)
setMethod("regFit", signature(engine = "general", stdErr = "sand"), regFit.general.sand)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "NULL"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "boot"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.LWYY", stdErr = "sand"), regFit.cox.LWYY)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "NULL"), regFit.cox.GL)
setMethod("regFit", signature(engine = "cox.GL", stdErr = "sand"), regFit.cox.GL)
setMethod("regFit", signature(engine = "am.GL", stdErr = "NULL"), regFit.am.GL)
setMethod("regFit", signature(engine = "Engine", stdErr = "boot"),
          regFit.Engine.boot)
setMethod("regFit", signature(engine = "am.GL", stdErr = "sand"),
          regFit.am.GL.sand)


#' Fits Semiparametric Regression Models for Recurrent Event Data
#'
#' Fits a general (joint) semiparametric regression model for the recurrent event data,
#' where the rate function of the underlying recurrent event process and
#' the hazard function of the terminal event can be specified as a Cox-type model,
#' an accelerated mean model, an accelerated rate model, or a generalized scale-change model.
#' See details for model specifications.
#'
#'
#' Suppose the recurrent event process and the failure events are
#' observed in the time interval \eqn{t\in[0,\tau]},
#' for some constant \eqn{\tau}.
#' We formulate the recurrent event rate function, \eqn{\lambda(t)},
#' and the terminal event hazard function, \eqn{h(t)}, 
#' in the form of
#' \deqn{\lambda(t) = Z \lambda_0(te^{X^\top\alpha}) e^{X^\top\beta}, h(t) = Z h_0(te^{X^\top\eta})e^{X^\top\theta},}
#' where \eqn{\lambda_0(t)} is the baseline rate function,
#' \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{X} is a \eqn{n} by \eqn{p} covariate matrix and \eqn{\alpha},
#' \eqn{Z} is an unobserved shared frailty variable, and
#' \eqn{(\alpha, \eta)} and \eqn{(\beta, \theta)} correspond to the shape and size parameters,
#' respectively.
#' The model includes several popular semiparametric models as special cases,
#' which can be specified via the \code{model} argument with the rate function
#' and hazard function separated by "\code{|}".
#' For examples,
#' Wang, Qin and Chiang (2001) (\eqn{\alpha = \eta = \theta = 0})
#' can be called with \code{model = "cox"};
#' Huang and Wang (2004) (\eqn{\alpha = \eta = 0})
#' can be called with \code{model = "cox|cox"};
#' Xu et al. (2017) (\eqn{\alpha = \beta} and \eqn{\eta = \theta})
#' can be called with \code{model = "am|am"};
#' Xu et al. (2019) (\eqn{\eta = \theta = 0}) can be called with \code{model = "gsc"}.
#' Users can mix the models depending on the application. For example,
#' \code{model = "cox|ar"} postulate a Cox proportional model for the
#' recurrent event rate function and an accelerated rate model for
#' the terminal event hazard function (\eqn{\alpha = \theta = 0}).
#' If only one model is specified without an "\code{|}",
#' it is used for both the rate function and the hazard function.
#' For example, specifying \code{model = "cox"} is equivalent to \code{model = "cox|cox"}.
#' Some models that assumes \code{Z = 1} and requires independent
#' censoring are also implemented in \code{reReg};
#' these includes \code{model = "cox.LWYY"} for Lin et al. (2000),
#' \code{model = "cox.GL"} for Ghosh and Lin (2002),
#' and \code{model = "am.GL"} for Ghosh and Lin (2003).
#'
#' The available methods for variance estimation are:
#' \describe{
#'   \item{boot}{performs nonparametric bootstrap.}
#'   \item{sand}{performs the efficient resampling-based variance estimation.}
#' }
#'
#' The \code{control} list consists of the following parameters:
#' \describe{
#'   \item{tol}{absolute error tolerance.}
#'   \item{alpha, beta, eta, theta}{initial guesses used for root search.}
#'   \item{solver}{the equation solver used for root search. The available options are \code{BB::BBsolve}, \code{BB::dfsane}, \code{BB:BBoptim}, and \code{optim}.}
#'   \item{eqType}{a character string indicating whether the log-rank type estimating equation or the Gehan-type estimating equation (when available) will be used. }
#'   \item{boot.parallel}{an logical value indicating whether parallel computation will be applied when \code{se = "boot"} is called.}
#'   \item{boot.parCl}{an integer value specifying the number of CPU cores to be used when \code{parallel = TRUE}. The default value is half the CPU cores on the current host.}
#' }
#' 
#' @param formula a formula object, with the response on the left of a "~" operator, and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{Recur}.
#' @param data  an optional data frame in which to interpret the variables occurring in the \code{"formula"}.
#' @param subset n optional logical vector specifying a subset of observations to be used
#' in the fitting process.
#' @param B a numeric value specifies the number of bootstraps for variance estimation.
#' When \code{B = 0}, variance estimation will not be performed.
#' @param model a character string specifying the underlying model. See \bold{Details}.
#' @param se a character string specifying the method for the variance estimation. See \bold{Details}.
#' \describe{
#'    \item{\code{boot}}{ nonparametric bootstrap approach}
#'    \item{\code{sand}}{ resampling-based sandwich estimator}
#' }
#' @param control a list of control parameters.
#'
#' @export
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711--730.
#' @references Wang, M.-C., Qin, J., and Chiang, C.-T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association}, \bold{96}(455): 1057--1065.
#' @references Ghosh, D. and Lin, D.Y. (2002). Marginal Regression Models for Recurrent and Terminal Events. \emph{Statistica Sinica}: 663--688.
#' @references Ghosh, D. and Lin, D.Y. (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877--885.
#' @references Huang, C.-Y. and Wang, M.-C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association}, \bold{99}(468): 1153--1165.
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017). Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association}, \bold{112}(518): 796--805.
#' @references Xu, G., Chiou, S.H.,Yan, J., Marr, K., and Huang, C.-Y. (2019). Generalized Scale-Change Models for Recurrent Event
#' Processes under Informative Censoring. \emph{Statistica Sinica}, \bold{30}: 1773--1795.
#'
#' @importFrom stats approxfun optim model.response
#' 
#' @seealso \code{\link{Recur}}, \code{\link{simGSC}}
#'
#' @example inst/examples/ex_reReg.R

reReg <- function(formula, data, subset,
                  model = "cox",
                  B = 0, se = c("boot", "sand"),
                  control = list()) {
    ## se = c("resampling", "bootstrap", "NULL"),
    ## se <- ifelse(is.null(se), "NULL", se)
    se <- match.arg(se)
    Call <- match.call()
    if (missing(formula)) stop("Argument 'formula' is required.")
    if (missing(data)) 
        data <- environment(formula)
    if (!missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (!is.logical(subIdx)) 
            stop("'subset' must be logical")
        subIdx <- subIdx & !is.na(subIdx)
        data <- data[subIdx, ]
    }    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$data <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    DF <- do.call(cbind, mf)
    DF <- as.data.frame(DF)
    obj <- model.response(mf)
    if (!is.Recur(obj)) stop("Response must be a `Recur` object")
    formula[[2]] <- NULL
    if (formula == ~ 1) DF$zero = 0 
    ctrl <- reReg.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    DF <- DF[order(DF$id, DF$time2), ]
    allModel <- apply(expand.grid(c("cox", "am", "gsc", "ar"),
                                   c("cox", "am", "gsc", "ar", ".")), 1, paste, collapse = "|")
    allModel <- c(allModel, "cox.LWYY", "cox.GL", "cox.HW", "am.GL", "am.XCHWY", "gsc.XCYH")
    model <- match.arg(model, c("cox", "am", "gsc", "ar", allModel))
    typeRec <- typeTem <- NULL
    if (grepl("|", model, fixed = TRUE)) {
        typeRec <- substring(model, 1, regexpr("[|]", model) - 1)
        typeTem <- substring(model, regexpr("[|]", model) + 1)
        model <- "general"
    }
    if (model %in% c("cox", "am", "gsc", "ar")) {
        typeRec <- model
        typeTem <- "."
        model <- "general"
    }
    ## Special cases:
    if (model == "cox.HW") {
        typeRec <- typeTem <- "cox"
        model <- "general"
    }
    if (model == "am.XCHWY") {
        typeRec <- typeTem <- "am"
        model <- "general"
    }
    if (model == "gsc.XCYH") {
        typeRec <- "gsc"
        typeTem <- "."
        model <- "general"        
    }
    if (model == "cox.LWYY") {
        typeRec <- "cox.LWYY"
        typeTem <- "."
    }
    if (model == "cox.GL") typeRec <- typeTem <- "cox.GL"
    if (model == "am.GL") typeRec <- typeTem <- "am.GL"
    if (length(unique(DF$time2[DF$event == 0])) == 1 & typeTem != ".") {
        typeTem <- "."
        cat("Only one unique censoring time is detected, terminal event model is not fitted.\n\n")
    }
    ## Temporary fix 
    if (typeRec != "gsc")  se <- "boot"
    engine.ctrl <- ctrl[names(ctrl) %in% names(attr(getClass(model), "slots"))]
    engine <- do.call("new", c(list(Class = model), engine.ctrl))
    engine@typeRec <- typeRec
    engine@typeTem <- typeTem
    if (se == "NULL" || B == 0)
        stdErr <- NULL
    else {
        stdErr.ctrl <- ctrl[names(ctrl) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class = se), stdErr.ctrl))
        stdErr@B <- B
    }
    ## initial values
    p <- ncol(DF) - ncol(mf[[1]])
    if (model == "general") {
        if (typeRec == "cox") {
            if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p + 1)
            if (length(engine@par1) == p) engine@par1 <- c(0, engine@par1)
            if (length(engine@par1) != (p + 1))
                stop("The length of initial value does not match with the number of covariates.")
            if (typeTem != ".") {
                engine@par3 <- engine@par2
                if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
                if (length(engine@par3) != p)
                    stop("The length of initial value does not match with the number of covariates.")
            }
        }
        if (typeRec == "gsc") {
            if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p)
            if (length(engine@par2) == 1) engine@par2 <- rep(engine@par2, p + 1)
            if (length(engine@par2) == p) engine@par2 <- c(0, engine@par2)
            if (length(engine@par1) != p | length(engine@par2) != (p + 1))
                stop("The length of initial value does not match with the number of covariates.")
            if (typeTem != ".") {
                if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
                if (length(engine@par4) == 1) engine@par4 <- rep(engine@par4, p)
                if (length(engine@par3) != p | length(engine@par4) != p)
                    stop("The length of initial value does not match with the number of covariates.")
            }
        }
        if (typeRec %in% c("ar", "am")) {
            if (length(engine@par1) == 1) engine@par1 <- rep(engine@par1, p)
            if (length(engine@par1) != p)
                stop("The length of initial value does not match with the number of covariates.")
            if (typeTem != ".") {
                engine@par3 <- engine@par2
                if (length(engine@par3) == 1) engine@par3 <- rep(engine@par3, p)
                if (length(engine@par3) != p)
                    stop("The length of initial value does not match with the number of covariates.")
            }
        }
    }
    engine@baseSE <- B > 0
    if (formula == ~1) {
        if (engine@baseSE) fit <- npFit(DF, B, typeTem)
        else fit <- npFit(DF, 0, typeTem)
        fit$typeRec <- "nonparametric"
        fit$typeTem <- typeTem
    } else {
        fit <- regFit(DF = DF, engine = engine, stdErr = stdErr)
    }    
    fit$DF <- DF
    fit$call <- Call
    fit$varNames <- names(DF)[-(1:6)]
    fit$se <- se
    if (engine@typeRec == "cox") fit$par1 <- fit$par1[-1]
    if (engine@typeRec == "gsc" & se != "boot") fit$par2 <- fit$par1 + fit$par2[-1]
    if (engine@typeRec == "gsc" & se == "boot") fit$par2 <- fit$par2[-1]
    if (se != "NULL" & se != "boot" & engine@typeRec == "gsc") fit$par2.se <- fit$par2.se[-1]   
    if (se != "NULL" & engine@typeRec == "cox") fit$par1.se <- fit$par1.se[-1]
    fit <- fit[order(names(fit))]
    class(fit) <- "reReg"
    return(fit)
}

#' Equation wrapper
#'
#' @noRd
#' @importFrom BB spg
#' @importFrom rootSolve uniroot.all
#' @keywords internal
eqSolve <- function(par, fn, solver, ...) {
    if (length(fn(par, ...)) == 1) {
        tmp <- uniroot.all(Vectorize(fn), interval = c(par - 10, par + 10))
        out <- NULL
        out$par <- tmp[which.min(abs(tmp - par))]
        out$convergence <- 0 ## 1 * (abs(fn(out$par, ...)) < 1e-5)
        return(out)
    }
    if (solver == "dfsane") {
        out <- dfsane(par = par, fn = function(z) fn(z, ...), 
                      alertConvergence = FALSE, quiet = TRUE, control = list(trace = FALSE))
        if (max(abs(out$par)) > 10) solver <- "BBsolve"
    }
    if (solver == "BBsolve")
        out <- BBsolve(par = par, fn = fn, ..., quiet = TRUE)
    if (solver == "BBoptim")
        out <- BBoptim(par = par, fn = function(z) sum(fn(z, ...)^2),
                       quiet = TRUE, control = list(trace = FALSE))
    if (solver == "optim")
        out <- optim(par = par, fn = function(z) sum(fn(z, ...)^2),
                     control = list(trace = FALSE))
    return(out)
}

reReg.control <- function(eqType = c("logrank", "gehan"),
                          solver = c("BB::dfsane", "BB::BBsolve", "BB::BBoptim", "optim"),
                          tol = 1e-7,
                          init = list(alpha = 0, beta = 0, eta = 0, theta = 0),
                          boot.parallel = FALSE, boot.parCl = NULL) {
    if (is.null(boot.parCl)) boot.parCl <- parallel::detectCores() / 2L
    solver <- match.arg(solver)
    if (solver == "BB::dfsane") solver <- "dfsane"
    if (solver == "BB::BBsolve") solver <- "BBsolve"
    if (solver == "BB::BBoptim") solver <- "BBoptim"
    eqType <- match.arg(eqType)
    list(tol = tol, eqType = eqType, solver = solver,
         par1 = init$alpha, par2 = init$beta, par3 = init$eta, par4 = init$theta,
         parallel = boot.parallel, parCl = boot.parCl)
}

##############################################################################
## Background functions...
## Probably need to clean these up someday
##############################################################################

#' R function for equation 8 of Ghosh & Lin (2002);
#' Marginal regression models for recurrent and terminal events.
#' 
#' @keywords internal
#' @noRd
coxGLeq <- function(beta, X, Y, T, cl, wgt) {
    p <- ncol(X)
    res <- vector("double", p)
    xb <- exp(X %*% beta)
    .C("coxGL", as.double(T), as.double(Y), as.double(X), as.double(xb), as.double(wgt),
       as.integer(length(T)), as.integer(cl), as.integer(c(0, cumsum(cl)[-length(cl)])),
       as.integer(nrow(X)), as.integer(p),        
       out = double(p), PACKAGE = "reReg")$out       
}

## varOut <- function(dat, na.rm = TRUE) {
##     dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
##     dat <- dat[complete.cases(dat),]
##     var(dat, na.rm = na.rm)
## }
