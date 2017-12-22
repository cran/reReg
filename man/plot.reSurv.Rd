\name{plot.reSurv}
\alias{plot.reSurv}
\title{Plot a reSurv Object.}
\description{
  Displays event plot for an reSurv object.
}
\usage{
\method{plot}{reSurv}(x, control = list(), ...)
}
\arguments{
  \item{x}{ an object of class \code{reSurv}, usually returned by
    the \code{reSurv} function.}
  \item{control}{ A list of control parameters.}
  \item{...}{for future methods}
}
\details{
  The argument control defaults to a list with values \code{xlab =
    "Time", ylab = "Subject", title = "Recurrent event plot",
    terminal.name = "Terminal event", recurrent.name =
    "Recurrent event", recurrent.types = NULL}.
  The function allows user-specified labels.
}
\seealso{\code{\link{reSurv}}}
\examples{
data(readmission)
reObj <- with(readmission, reSurv(t.stop, id, event, death))
## Default labels
plot(reObj)
## User specified labels
plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))

## Multiple recurrent event types
readmission$event2 <- ifelse(readmission$event == 0, 0, sample(1:3,
nrow(readmission), TRUE))
reObj2 <- with(readmission, reSurv(t.stop, id, event2, death))
plot(reObj2)
## Also allows for user specified names
plot(reObj2, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title",
             recurrent.types = c("Recurrent type A", "Recurrent type B", "Recurrent type C")))
}
\keyword{plot.reSurv}
