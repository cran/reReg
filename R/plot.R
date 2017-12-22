## Default Event plot
plot.reSurv <- function(x, control = list(), ...) {
    recType <- event <- status <- id <- Time <- NULL
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    ## event <- x$reDF$event
    ## status <- x$reDF$status
    ## id <- x$reDF$id
    ## Time <- x$reDF$Time
    if (!is.reSurv(x)) stop("Response must be a reSurv class")
    n <- length(unique(x$reDF$id))
    ldat0 <- data.frame(id = rep(unique(x$reDF$id), 2),
                        Time = c(rep(0, n), aggregate(Time ~ id, max, data = x$reDF)$Time),
                        status = c(rep(0, n), aggregate(status ~ id, max, data = x$reDF)$status))
    ldat0 <- ldat0[order(ldat0$id),]
    rownames(ldat0) <- NULL
    ldat <- x$reDF
    ## table(ldat$event:ldat$status)
    ## if (control$plotCen) {
    ##     ldat$event <- as.factor(ldat$event)
    ##     ldat$status <- as.factor(ldat$status)
    ##     p <- ggplot(ldat0, aes(x = Time, y = id, group = id)) +
    ##         geom_line(color = "gray55", size = 1.05) +
    ##         geom_point(data = ldat, aes(x = Time, y = id, color = event:status, shape = event:status),
    ##                    size = 1.5, alpha = .7) +
    ##         scale_color_manual(values = c("red", "red", "blue"), name = "Event types:",
    ##                            breaks=c("0:0", "0:1", "1:0"),
    ##                            labels = c("Censored\n terminal events",
    ##                                       "Uncensored\n terminal event", "Recurrent event")) +
    ##         scale_shape_manual(values=c(2, 17, 16), name = "Event types:",
    ##                        breaks=c("0:0", "0:1", "1:0"),
    ##                        labels = c("Censored\n terminal events",
    ##                                   "Uncensored\n terminal event", "Recurrent event")) +
    ##         ggtitle(control$title) +
    ##         theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
    ##         labs(x = control$xlab, y = control$ylab)
    ## }
    ## if (!control$plotCen) {
    if (sum(ldat$event + ldat$status == 0) > 0)
        ldat <- ldat[-which(ldat$event + ldat$status == 0),]
    ldat$event <- as.factor(ldat$event)
    ldat$status <- as.factor(ldat$status)
    ldat$recType <- as.factor(ldat$recType)
    if (identical(ldat$event, ldat$recType)) {
        p <- ggplot(ldat0, aes(x = Time, y = id, group = id)) +
            geom_line(color = "gray55", size = 1.05) +
            geom_point(data = ldat, aes(x = Time, y = id, shape = event:status),
                       size = 1.5, alpha = .7) +
            scale_shape_manual(values=c(16, 4), name = "Event types:",
                               breaks=c("0:1", "1:0"),
                               labels = c(ctrl$terminal.name, ctrl$recurrent.name)) +
            ggtitle(ctrl$title) +
            labs(x = ctrl$xlab, y = ctrl$ylab) + theme_bw()
    } else {
        tmp <- table(with(ldat, recType:status))
        k <- length(unique(ldat$recType)) - 1
        if (is.null(ctrl$recurrent.types))
            ctrl$recurrent.types <- sapply(1:k, function(x) paste("Recurrent Type", x, collapse=""))
        if (length(ctrl$recurrent.types) != k)
            stop("User specified recurrent name must be the same length as the number of recurrent events")
        p <- ggplot(ldat0, aes(x = Time, y = id, group = id)) +
            geom_line(color = "gray55", size = 1.05) +
            geom_point(data = ldat, aes(x = Time, y = id, shape = recType:status),
                       size = 1.5, alpha = .7) +
            scale_shape_manual(values = c(16, 1:k), name = "Event types:",
                               ## breaks = names(tmp[tmp > 0]),
                               labels = c(ctrl$terminal.name, ctrl$recurrent.types)) +
            ggtitle(ctrl$title) +
            labs(x = ctrl$xlab, y = ctrl$ylab) + theme_bw()
    }
    p
}

##################################################################################
## More flexibleEvent plots
##################################################################################

plotEvents <- function(formula, data, control = list(), ...) {
    recType <- event <- status <- id <- Time <- NULL
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    if (class(formula) == "reSurv") return(plot(formula))
    if (missing(data)) obj <- eval(formula[[2]], parent.frame()) 
    if (!missing(data)) obj <- eval(formula[[2]], data) 
    if (!is.reSurv(obj)) stop("Response must be a reSurv object")
    if (formula[[3]] == 1) {
        return(plot(obj))
    } else {
        if (length(attr(terms(formula), "term.labels")) > 1)
            stop("The plotEvents can only take one covaraite")
        if (missing(data)) DF <- cbind(obj$reDF, eval(formula[[3]], parent.frame()))
        if (!missing(data)) DF <- cbind(obj$reDF, eval(formula[[3]], data))
        colnames(DF) <- c(names(obj$reDF), paste0(formula[3], collapse = ""))
    }
    n <- length(unique(DF$id))
    ldat0 <- data.frame(id = rep(unique(DF$id), 2),
                        Time = c(rep(0, n), aggregate(Time ~ id, max, data = DF)$Time),
                        status = c(rep(0, n), aggregate(status ~ id, max, data = DF)$status))
    ldat0 <- ldat0[order(ldat0$id),]
    rownames(ldat0) <- NULL
    ldat0$tmp <- DF[match(ldat0$id, DF$id), 6]
    ldat <- DF
    if (grepl(":", names(ldat)[6])) names(ldat)[6] <- gsub(":", "", names(ldat)[6])
    names(ldat0)[4] <- names(ldat)[6]
    ldat[,6] <- as.factor(ldat[,6])
    if (sum(ldat$event + ldat$status == 0) > 0)
        ldat <- ldat[-which(ldat$event + ldat$status == 0),]
    ldat$event <- as.factor(ldat$event)
    ldat$status <- as.factor(ldat$status)
    ldat$recType <- as.factor(ldat$recType)
    if (identical(ldat$event, ldat$recType)) {
        ## k <- length(unique(ldat$recType)) - 1
        p <- ggplot(ldat0, aes(x = Time, y = id, group = id)) +
            geom_line(color = "gray55", size = 1.05) +
            facet_grid(eval(as.symbol(names(ldat)[6])) ~ ., scales = "free_y", space = "free_y") +
            geom_point(data = ldat, aes(x = Time, y = id, shape = event:status),
                       size = 1.5, alpha = .7) +
            scale_shape_manual(values=c(16, 4), name = "Event types:",
                               ## breaks=c("0:1", "1:0"),
                               labels = c(ctrl$terminal.name, ctrl$recurrent.name)) +
            ggtitle(ctrl$title) +
            labs(x = ctrl$xlab, y = ctrl$ylab) + theme_bw()
    } else {
        ## tmp <- table(with(ldat, recType:status))
        k <- length(unique(ldat$recType)) - 1
        if (is.null(ctrl$recurrent.types))
            ctrl$recurrent.types <- sapply(1:k, function(x) paste("Recurrent Type", x, collapse=""))
        if (length(ctrl$recurrent.types) != k)
            stop("User specified recurrent name must be the same length as the number of recurrent events")
        p <- ggplot(ldat0, aes(x = Time, y = id, group = id)) +
            geom_line(color = "gray55", size = 1.05) +
            facet_grid(eval(as.symbol(names(ldat)[6])) ~ ., scales = "free_y", space = "free_y") +
            geom_point(data = ldat, aes(x = Time, y = id, shape = recType:status),
                       size = 1.5, alpha = .7) +
            scale_shape_manual(values = c(16, 1:k), name = "Event types:",
                               ## breaks = names(tmp[tmp > 0]),
                               labels = c(ctrl$terminal.name, ctrl$recurrent.types)) +
            ggtitle(ctrl$title) +
            labs(x = ctrl$xlab, y = ctrl$ylab) + theme_bw()
    }
    p    
}

plotEvents.control <- function(xlab = "Time", ylab = "Subject", title = "Recurrent event plot",
                               terminal.name = "Terminal event",
                               recurrent.name = "Recurrent event",
                               recurrent.types = NULL) {
    list(xlab = xlab, ylab = ylab, title = title, terminal.name = terminal.name,
         recurrent.name = recurrent.name, recurrent.types = recurrent.types)
}

##################################################################################

plot.reReg <- function(x, ...) {
    options(warn = -1)
    if (!is.reReg(x)) stop("Response must be a reReg class")
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH") 
    par(mfrow = c(2, 1))
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.ly <- max(ly, lyU, lyL, na.rm = TRUE)
    win.hy <- max(hy, hyU, hyL, na.rm = TRUE)
    plot(t0, ly, type = "s",  xlab = "", ylab = "", ylim = c(0, win.ly),
         main = "Baseline Cumulative Rate Function")
    if (any(!is.na(x$lamU))) {
        lines(t0, lyU, col = 2)
        lines(t0, lyL, col = 2)
    }
    title(ylab = expression(hat(Lambda)[0](t)), xlab = "Time", line = 2.2)
    ##     title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.ly, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    ## op <- par(ask=TRUE)
    ## plot(t, hy, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, win.hy),
    ##      main = "Baseline Cumulative Hazard Function")
    plot(t0, hy, type = "s",  xlab = "", ylab = "", ylim = c(0, win.hy),
         main = "Baseline Cumulative Hazard Function")
    if (any(!is.na(x$hazU))) {
        lines(t0, hyU, col = 2, lty = 2)
        lines(t0, hyL, col = 2, lty = 2)
    }
    title(ylab = expression(hat(H)[0](t)), xlab = "Time", line = 2.2)
    ## title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.hy, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    ## par(op)
    par(mfrow = c(1, 1))
    out <- c(x, list(t0 = t0, lam = ly, lamU = lyU, lamL = lyL, haz = hy, hazU = hyU, hazL = hyL))
    options(warn = 0)
    invisible(out)
}

plotRate <- function(x, control = list(), ...) {
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH") 
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$title == "Recurrent event plot") # default value
        ctrl$title <- "Baseline Cumulative Rate Function"
    if (!is.reReg(x)) stop("Response must be a reReg class")
    options(warn = -1)
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.ly <- max(ly, lyU, lyL, na.rm = TRUE)
    plot(t0, ly, type = "s",  xlab = "", ylab = "", ylim = c(0, win.ly),
         main = ctrl$title)
    if (any(!is.na(x$lamU))) {
        lines(t0, lyU, col = 2, lty = 2, "s")
        lines(t0, lyL, col = 2, lty = 2, "s")
    }
    if (ctrl$ylab == "Subject") #default value 
        title(ylab = expression(hat(Lambda)[0](t)), xlab = ctrl$xlab, line = 2.2)
    else title(ylab = ctrl$ylab, xlab = ctrl$xlab, line = 2.2)
    options(warn = 0)
}


plotHaz <- function(x, control = list(), ...) {
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH")
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$title == "Recurrent event plot") # default value
        ctrl$title <- "Baseline Cumulative Hazard Function"
    if (!is.reReg(x)) stop("Response must be a reReg class")
    options(warn = -1)
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.hy <- max(hy, hyU, hyL, na.rm = TRUE)
    plot(t0, hy, type = "s",  xlab = "", ylab = "", ylim = c(0, win.hy),
         main = ctrl$title)
    if (any(!is.na(x$hazU))) {
        lines(t0, hyU, col = 2, lty = 2, "s")
        lines(t0, hyL, col = 2, lty = 2, "s")
    }
    if (ctrl$ylab == "Subject") #default value 
        title(ylab = expression(hat(Lambda)[0](t)), xlab = ctrl$xlab, line = 2.2)
    else title(ylab = ctrl$ylab, xlab = ctrl$xlab, line = 2.2)
    options(warn = 0)
}
