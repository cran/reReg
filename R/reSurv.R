reSurv <- function(time1, id, event, status, time2 = NULL) {
    if (sum(time1 <= 0) > 0 & is.null(time2))
        stop("Observation time must be positive.")
    if (!is.null(time2) & any(time1 > time2)) 
        stop("Stop time must be > start time")
    if (!is.numeric(time1))
        stop("Time variable is not numeric")
    if (sum(time1 < 0) > 0)
        stop("Observation time must be positive.")
    if (sum(is.wholenumber(id)) < length(id))
        stop("ID must be integer values.")
    if (!is.numeric(time2) & !is.null(time2)) 
        stop("Time variable is not numeric")
    if (length(event) != length(time1))
        stop("Time and status are different lengths")
    if (length(status) != length(time1))
        stop("Time and status are different lengths")
    event2 <- NULL
    if (is.logical(event)) 
        event <- as.numeric(event)
    else if (is.numeric(event)) {
        event2 <- event
        event <- ifelse(event == 0, 0, 1)
        ## temp <- (event == 0 | event == 1)
        ## event <- ifelse(temp, event, NA)
        ## if (sum(is.na(event)) > 0) stop("Event must be 0 or 1)")
    }
    else stop("Event must be logical or numeric")
    if (is.logical(status)) 
        status <- as.numeric(status)
    else if (is.numeric(status)) {
        temp <- (status == 0 | status == 1)
        status <- ifelse(temp, status, NA)
        if (sum(is.na(status)) > 0) 
            stop("Status must be 0 or 1)")
    }
    else stop("Status must be logical or numeric")
    if (is.null(time2))
        tab <- data.frame(id = id, Time = time1, event = event, status = status, recType = event2)
    else tab <- data.frame(id = id, Time = unlist(lapply(split(time2 - time1, id), cumsum)),
                           event = event, status = status, recType = event2)
    rc <- list(reDF = tab)
    class(rc) <- "reSurv"
    invisible(rc)
}

is.reSurv <- function(x) inherits(x, "reSurv")
is.reReg <- function(x) inherits(x, "reReg")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
