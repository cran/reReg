\name{reSurv}
\alias{reSurv}
\alias{is.reSurv}
\title{Create a Survival Object for Recurrent Event Data}
\description{
  Create a recurrent event survival object, used as a response
  variable in reReg model formula.
}
\usage{
  reSurv(time1, id, event, status, time2 = NULL)
  is.reSurv(x) 
}
\arguments{
  \item{id}{Observation subject's ID.}
  \item{time1}{When 'time2' is provided, this vector is treated as the
    starting time for the gap time between two successive recurrent
    events. In the absence of 'time2', this is the observation time of
    recurrence on calendar time scale, in which, the time corresponds to
    the time since entry/inclusion in the study.}
  \item{time2}{An optional vector for ending time for the gap time
    between two successive recurrent events.}
  \item{event}{ A binary vector used as the recurrent event indicator.}
  \item{status}{A binary vector used as the status indicator for the terminal event.}
  \item{x}{An \code{reSurv} object.}
}
%% \details{
%% }
%% \value{
%%   \item{rsDF}{a data frame with variable "id", "Time", and "event".}
%%   In the case of \code{is.reSurv}, a logical value \code{TRUE} if
%%   \code{x} inherits from class \code{"reSurv"}, otherwise an \code{FALSE}.
%% }
\seealso{\code{\link{reReg}}}
\examples{
data(readmission)
head(readmission)
(with(readmission[1:10,], reSurv(t.stop, id, event, death)))
plot(with(readmission[1:10,], reSurv(t.stop, id, event, death)))
}
\keyword{reReg}
\keyword{plot.reReg}