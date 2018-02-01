\name{plotRate}
\alias{plotRate}
\title{Plot the baseline rate function after fitting reReg.}
\description{
  Plot the baseline rate function for the recurrent process after fitting reReg.
}
\usage{
plotRate(x, control = list(), ...)
}
\arguments{
  \item{x}{ an object of class \code{reReg}, usually returned by
    the \code{reReg} function.}
  \item{control}{ A list of control parameters.}
  \item{...}{for future methods}
}
\seealso{\code{\link{reReg}}}
\examples{
## readmission data
data(readmission)
set.seed(123)
fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
                   data = subset(readmission, id < 50),
                   method = "am.XCHWY", se = "resampling", B = 20)
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative rate function
plotRate(fit)
## Plot with user-specified labels
plotRate(fit, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title")) 
}
