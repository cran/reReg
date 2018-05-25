\name{plot.reSurv}
\alias{plot.reSurv}
\title{Plot a reSurv Object.}
\description{
  Displays an event plot for an reSurv object.
}
\usage{
\method{plot}{reSurv}(x, data, order = TRUE, return.grob = FALSE, control = list(), ...)
}
\arguments{
  \item{x}{an object of class \code{reSurv}, usually returned by
    the \code{reSurv} function.}
  \item{data}{
    an optional data frame in which to interpret the variables
    occurring in the `formula'.
  }
  \item{order}{an optional logical value. If "TRUE", the event plot
    is sorted by the terminal times.}
  \item{return.grob}{an optional logical value. If "TRUE", a ggplot2 plot grob
    will be returned.}
  \item{control}{ a list of control parameters.}
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
reObj <- with(subset(readmission, id <= 10), reSurv(t.stop, event, death, id))
## Default labels
plot(reObj)
plot(reObj, order = FALSE)
## User specified labels
plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))

## With multiple hypothetical event types
set.seed(1)
reObj2 <- with(readmission, reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id))
plot(reObj2)
}
\keyword{plot.reSurv}
