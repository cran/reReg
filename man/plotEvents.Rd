\name{plotEvents}
\alias{plotEvents}
\title{Produce event plot}
\description{
  Fits a joint model for the recurrent event process and the failure
  time that allows the censoring time to be informative about the
  recurrent event process. The rate function of the underlying process
  for the recurrent event process be specified as a Cox-type model, an
  accelerated mean model, an accelerated failure time model, or a
  generalized scale-change model.
}
\usage{
plotEvents(formula, data, control = list(), ...) 
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
  \item{control}{ A list of control parameters.}
  \item{...}{for future methods}
}
\details{
  The underlying models and assumptions are different for different methods.
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
with(readmission, plotEvents(reSurv(t.stop, id, event, death)))

## Separate plots by gender
with(readmission, plotEvents(reSurv(t.stop, id, event, death) ~ sex))

## Separate plots by gender and chemo type
with(readmission, plotEvents(reSurv(t.stop, id, event, death) ~ sex:chemo))

## Multiple recurrent event types
readmission$event2 <- ifelse(readmission$event == 0, 0, sample(1:3,
nrow(readmission), TRUE))
plotEvents(reSurv(t.stop, id, event2, death) ~ sex, data = readmission)
## plotEvents allows for user specified names
plotEvents(reSurv(t.stop, id, event2, death) ~ sex, data = readmission,
           control = list(xlab = "User xlab", ylab = "User ylab", title = "User title",
           recurrent.types = c("Recurrent type A", "Recurrent type B", "Recurrent type C")))
}


