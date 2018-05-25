\name{plotEvents}
\alias{plotEvents}
\title{Produce event plot}
\description{
  Displays an event plot with more flexible options.
}
\usage{
plotEvents(formula, data, order = TRUE, return.grob = FALSE, control = list(), ...) 
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
\value{
  An object of S3 class \code{"reReg"} representing the fit, with the
}
\seealso{
  \code{\link{reSurv}}, \code{\link{plot.reSurv}}
}
\examples{
## readmission data
data(readmission)
plotEvents(reSurv(t.stop, event, death, id) ~ 1, data = readmission)

## Separate plots by gender
plotEvents(reSurv(t.stop, event, death, id) ~ sex, data = readmission)

## Separate plots by gender and chemo type
plotEvents(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission)

## With multiple hypothetical event types
plotEvents(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~
sex + chemo, data = readmission)
}


