\name{plotHaz}
\alias{plotHaz}
\title{Plot the baseline function after fitting reReg.}
\description{
  Plot the baseline cumulative hazard function for terminal events after fitting reReg.
}
\usage{
plotHaz(x, control = list(), ...)
}
\arguments{
  \item{x}{ an object of class \code{reReg}, usually returned by
    the \code{reReg} function.}
  \item{control}{ A list of control parameters.}
  \item{...}{for future methods}
}
\details{
  The argument control defaults to a list with values \code{xlab =
    "Time", ylab = "Subject", title = "Recurrent event plot",
    terminal.name = "Terminal event"}.
  The function allows user-specified labels.
}
\seealso{\code{\link{reReg}}}
\examples{
## readmission data
data(readmission)
set.seed(123)
fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
             data = subset(readmission, id < 50),
             method = "am.XCHWY", se = "resampling", B = 20, plot.ci = TRUE)
## Plot both the baseline cumulative rate and hazard function
plot(fit)
## Plot baseline cumulative hazard function
plotHaz(fit)
## Plot with user-specified labels
plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))  
}
