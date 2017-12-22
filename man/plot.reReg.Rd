\name{plot.reReg}
\alias{plot.reReg}
\title{Plot a reReg Object.}
\description{
  Displays graph of the estimated baseline cumulative rate function
  and harzard function if available.
}
\usage{
\method{plot}{reReg}(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{reReg}, usually returned by
    the \code{reReg} function.}
  \item{...}{for future methods}
}
\seealso{\code{\link{reReg}}}
\examples{
## readmission data
data(readmission)
fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
             data = subset(readmission, id < 50),
             method = "am.XCHWY", se = "resampling", B = 20, plot.ci = TRUE)
plot(fit)
}
\keyword{reReg}
\keyword{plot.reReg}
