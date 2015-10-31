\name{reReg}
\alias{reReg}
\title{Fits semiparametric regression models for recurrent events and
  failure time}
\description{
  Fits a joint model for the recurrent event process and the failure
  time that allows the censoring time to be informative about the
  recurrent event process.
  
  Estimating procedures include a joint scale-change model,
  Huang-Wang's method (Huang and Wang, 2004), and
  Ghosh-Lin's method (Ghosh and Lin, 2003). The variance is estiamted
  with an efficient bootstrap approach.
}
\usage{
reReg(formula, data, subset, method = "jsc", B = 200, contrasts = NULL)
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a
    `~' operator, and the terms on the right. The response must be a
    recurrent event survival object as returned by function
    \code{reSurv}.
  }
  \item{data}{
    an optional data frame in which to interpret the variables
    occurring in the `formula'.
  }
  \item{subset}{ an optional vector specifying a subset of observations
    to be used in the fitting process.
  }
  \item{method}{
    a character string specifying the underlying model and methods to
    estimate the regression parameter. The following are permitted: `jsc':
    joint scale-change model. `HW': Huang-Wang's
    method. `GL': Ghosh-Lin method. Among these, `jsc' is recommended.
    }
  \item{B}{
    a numeric value specifies the resampling number for variance
    estimation. 
    When B = 0, only the point estimates will be displayed.
  }
  \item{contrasts}{ an optional list.
  }
}
%% \details{
%%   The underlying models are different for different methods. 
%%   When \code{method = "jsc"}, the rate function of the counting process
%%   takes the form
%%   }
\value{
  An object of S3 class \code{"reReg"} representing the fit, with the
  following components:
  \item{alpha}{a vector of point estimates associated with rate function, \code{alpha}}
  \item{beta}{ a vector of point estimates associated with hazard function, \code{beta}}
  \item{va}{ estimated covariance matrix for \code{alpha}}
  \item{vb}{ estimated covariance matrix for \code{beta}}
  \item{lambda}{ a vector giving the estimated baseline cumulative rate
    function evaluated at the observation time of recurrence}
  \item{haz}{ a vector giving the estimated baseline cumulative hazard
    function evaluated at the observation time of recurrence}
  \item{zHat}{ a vector giving the estimated frailty random variable at the observation time of recurrence}
}
\references{
  Ghosh, D. and D.Y. Lin (2003). Semiparametric Analysis of Recurrent
  Events Data in the Presence of Dependent Censoring. \emph{Biometrics},
  \bold{59}: 877 -- 885.

  Huang, C.Y. and Wang, M.C. (2004). Joint Modeling and Estimation
  for Recurrent Event Processes and Failure Time Data. \emph{Journal
    of the American Statistical Association} \bold{99}(468): 1153--1165.

  %% Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. Joint
  %% Scale-change Models for Recurrent Events and Failure Time. 
}

\seealso{
  \code{\link{reSurv}}, \code{\link{plot.reReg}}
}
\examples{
## readmission data
data(readmission)
(fit <- reReg(reSurv(t.start, id, death, t.stop) ~ sex + dukes + chemo,
  data = readmission, B = 0))
summary(fit)
(fit.HW <- reReg(reSurv(t.start, id, death, t.stop) ~ sex + dukes + chemo,
  data = readmission, method = "HW", B = 0))
summary(fit.HW)

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
                    Delta = rep(Delta, unlist(lapply(eventT, length)))
                    )
  out           
}
set.seed(123)
mydat <- simDat(50, a = c(-1, 1), b = c(-1, 1), tau = 10)
(fit.sim <- reReg(reSurv(T, id, Delta) ~ X1 + X2, data = mydat, B = 0))
summary(fit.sim)
}
\keyword{reReg}


