\name{plot.reReg}
\alias{plot.reReg}
\title{Plot a reReg Object.}
\description{
  Displays graph of the estimated baseline cumulative rate function
  and harzard function.
}
\usage{
\method{plot}{reReg}(x, se = FALSE, B = 200, breaks = 1000, ...)
}
\arguments{
  \item{x}{ an object of class \code{reReg}, usually returned by
    the \code{reReg} function.}
  \item{se}{a logical value, if TRUE and \code{x} is estimated with
    joint scale-change model, standard error will be estimated via a
    efficient resampling-based procedure. Empirical pointwise 95\%
    confidence intervals will be displayed as well.}
  \item{B}{ a numeric value specifies the resampling number for variance
    estimation when \code{se == TRUE}.}
  \item{breaks}{a numeric value specifies the number of breakpoints.}
  \item{...}{for future methods}
}
%% \details{
%% }
\value{
  a list containing the coordinates of the breakpoints on each of the
  curves as well as the empirical pointwise confidence bounds.
  The output contains the following components:
  \item{x}{a vector giving the breakpoints used in graph.}
  \item{ly}{a vector giving the estimated baseline cumulative rate
    function evaluated at \code{x}}
  \item{lyU}{a vector giving the empirical upper bound for baseline
    cumulative rate function evaluated at \code{x}}
  \item{lyL}{a vector giving the empirical lower bound for baseline
    cumulative rate function evaluated at \code{x}}
  \item{hy}{a vector giving the estimated baseline cumulative hazard
    function evaluated at \code{x}}
  \item{hyU}{a vector giving the empirical upper bound for baseline
    cumulative hazard function evaluated at \code{x}}
  \item{hyL}{a vector giving the empirical lower bound for baseline
    cumulative hazard function evaluated at \code{x}}
}

\references{
  Wang, M.C., Qin, J., and Chiang, C.T. (2001). Analyzing Recurrent
  Event Data with Informative Censoring. \emph{Journal of the American
  Statistical Association} \bold{96}{455}: 1057--1065.
  
  Huang, C.Y. and Wang, M.C. (2004). Joint Modeling and Estimation
  for Recurrent Event Processes and Failure Time Data. \emph{Journal
    of the American Statistical Association} \bold{99}(468): 1153--1165.
}


\seealso{\code{\link{reReg}}}
\examples{
## readmission data
data(readmission)
fit <- reReg(reSurv(t.start, id, death, t.stop) ~ sex + dukes + chemo, data = readmission, B = 0)
plot(fit)

## simulation data
simDat <- function(n, a, b, tau) {
  z <- rgamma(n, 2, 0.2)
  X1 <- rbinom(n, 1, 0.5)
  X2 <- runif(n, 0 ,1)
  lambda <- z * (1/10) * exp(cbind(X1, X2) \%*\% a)
  mt <- rpois(n, tau * lambda)
  eventT <- lapply(sapply(mt, function(x) runif(x) * tau), sort)
  D <- exp(- cbind(X1, X2) \%*\% b) * sqrt(800 * (-log(runif(n)) / z))
  C <- (X1 == 1) * rexp(n, 0.1) + (X1 == 0) * rexp(n, ifelse(z == 0, 10^-10, z)^2/300)
  Y <- pmin(C, tau, D)
  Delta <- 1 * (D <= Y)
  for (i in 1:n) {
    eventT[[i]] <- c(eventT[[i]][eventT[[i]] < Y[i]], Y[i])
  }
  out <- data.frame(id = rep(1:n, unlist(lapply(eventT, length))),
                    T = unlist(eventT),
                    X1 = rep(X1, unlist(lapply(eventT, length))),
                    X2 = rep(X2, unlist(lapply(eventT, length))),
                    Delta = rep(Delta, unlist(lapply(eventT, length))))
  out           
}
set.seed(123)
mydat <- simDat(30, a = c(-1, 1), b = c(-1, 1), tau = 10)
(fit.sim <- reReg(reSurv(T, id, Delta) ~ X1 + X2, data = mydat, B = 0))
plot(fit.sim, se = TRUE, B = 10)
}
\keyword{reReg}
\keyword{plot.reReg}
