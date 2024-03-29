% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotEvents}
\alias{plotEvents}
\title{Produce Event Plots}
\usage{
plotEvents(
  formula,
  data,
  result = c("increasing", "decreasing", "asis"),
  calendarTime = FALSE,
  control = list(),
  ...
)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a "~" operator,
and the predictors on the right.
The response must be a recurrent event survival object as returned by function \code{Recur()}.}

\item{data}{an optional data frame in which to interpret the variables occurring in
the "\code{formula}".}

\item{result}{an optional character string specifying whether the event plot is
sorted by the subjects' terminal time. The available options are
\describe{
  \item{\code{increasing}}{sort the terminal time from in ascending order (default).
This places longer terminal times on top. }
  \item{\code{decreasing}}{sort the terminal time from in descending order.
This places shorter terminal times on top. }
  \item{\code{none}}{present the event plots as is, without sorting by the terminal times.}
}}

\item{calendarTime}{an optional logical value indicating whether to plot in calendar time.
When \code{calendarTime = FALSE} (default), the event plot will have patient time on the x-axis.}

\item{control}{a list of control parameters. See \bold{Details}.}

\item{...}{graphical parameters to be passed to methods.
These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.}
}
\value{
A \code{ggplot} object.
}
\description{
Plot the event plot for an \code{Recur} object.
The usage of the function is similar to that of \code{plot.Recur()} but with more flexible options.
}
\details{
The argument \code{control} consists of options with argument defaults to a list with
the following values:
\describe{
  \item{xlab}{customizable x-label, default value is "Time".}
  \item{ylab}{customizable y-label, default value is "Subject" for event plot and
"Cumulative mean" for MCF plot.}
  \item{main}{customizable title, the default value is "Recurrent event plot"
when \code{mcf = FALSE} and
"Sample cumulative mean function plot" when \code{mcf = TRUE}.}
  \item{terminal.name}{customizable label for terminal event,
the default value is "Terminal event".}
  \item{recurrent.name}{customizable legend title for recurrent event,
the default value is "Recurrent events".}
  \item{recurrent.types}{customizable label for recurrent event type,
the default value is \code{NULL}.}
  \item{alpha}{between 0 and 1, controls the transparency of points.}
}
The \code{xlab}, \code{ylab} and \code{main} parameters can be specified
outside of the \code{control} list.
}
\examples{
data(simDat)
plotEvents(Recur(t.start \%to\% t.stop, id, event, status) ~ 1, data = simDat,
           xlab = "Time in days", ylab = "Subjects arranged by terminal time")

## Separate plots by x1
plotEvents(Recur(t.start \%to\% t.stop, id, event, status) ~ x1, data = simDat)

## For multiple recurrent events
simDat$x3 <- ifelse(simDat$x2 < 0, "x2 < 0", "x2 > 0")
simDat$event <- simDat$event * sample(1:3, nrow(simDat), TRUE)
plotEvents(Recur(t.start \%to\% t.stop, id, event, status) ~ x1 + x3, data = simDat)
}
\seealso{
\code{\link{Recur}}, \code{\link{plot.Recur}}
}
\keyword{Plots}
