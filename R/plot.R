globalVariables(c("Time", "Yi", "id", "recType", "status", "tij"))

## Default Event plot
plot.reSurv <- function(x, data, order = TRUE, return.grob = FALSE, control = list(), ...) {
    plotEvents(x, data, order = order, return.grob = return.grob, control = control)
}

plotEvents <- function(formula, data, order = TRUE, return.grob = FALSE, control = list(), ...) {
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    nX <- 0
    if (is.reSurv(formula)) {dat <- formula$reTb
    } else {
        ## if (length(attr(terms(formula), "term.labels")) > 1)
        ##     stop("The current vision can only handle one covaraite.")
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        dat <- obj$reTb
        nX <- attr(terms(formula), "term.labels")
        if (formula[[3]] != 1) {
            if (missing(data)) DF <- cbind(obj$reDF, sapply(nX, function(x) eval(as.symbol(paste(x)), parent.frame())))
            if (!missing(data)) DF <- cbind(obj$reDF, sapply(nX, function(x) eval(as.symbol(paste(x)), data)))
            ## colnames(DF) <- c(names(obj$reDF), nX)
            suppressMessages(dat <- left_join(obj$reTb, unique(select(DF, id, nX))))
        }
    }
    if (order) {
        if (nX == 0 || formula[[3]] == 1) dat <- dat %>% mutate(id = rank(Yi, ties.method = "first"))
        else dat <- dat %>% group_by_at(6:ncol(dat)) %>% mutate(id = rank(Yi, ties.method = "first")) 
    }
    sz <- 1 + 8 / (1 + exp(length(unique(dat$id)) / 30))
    k <- length(unique(unlist(dat$recType)))
    shp.val <- c(17, rep(19, k))
    clr.val <- c(alpha("red", ctrl$alpha), hcl(h = seq(120, 360, length.out = k), l = 60, alpha = ctrl$alpha))
    rec.lab <- paste("r", 1:k, sep = "")
    if (k == 1) {
        shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.name)
    } 
    if (k > 1 & is.null(ctrl$recurrent.type)) {
        shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))
    }
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.type)
        } else {
            cat('The length of "recurrent.type" mismatched, default names are used.\n')
            shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))            
        }
    }
    names(shp.val) <- names(clr.val) <- c("Yi", rec.lab)
    ## ggplot starts here
    gg <- ggplot(dat, aes(id, Yi)) +
        geom_bar(stat = "identity", fill = "gray75") +
        geom_point(data = dat %>% filter(status > 0),
                   aes(id, Yi, shape = "Yi", color = "Yi"), size = sz) +
        geom_point(data = dat %>% filter(!map_lgl(tij, is.null)) %>%
                       unnest(tij, recType), # %>% select(id, tij, recType),
                   aes(id, tij, shape = factor(recType, labels = rec.lab), 
                       color = factor(recType, labels = rec.lab)),
                   size = sz) + 
        ## position = position_jitter(w = 0.1, h = 0)) +
        scale_shape_manual(
            name = "", values = shp.val,
            labels = shp.lab,
            breaks = c("Yi", rec.lab)) +
        scale_color_manual(
            name = "", values = clr.val,
            labels = shp.lab,
            breaks = c("Yi", rec.lab)) +
        coord_flip() + 
        theme(axis.line.y = element_blank(),
              axis.title.y = element_text(vjust = 0),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "both")
    if (!return.grob) {
        gg + theme(panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.key = element_rect(fill = "white", color = "white")) +
            scale_x_continuous(expand = c(0, 1)) +
            ggtitle(ctrl$title) + labs(x = ctrl$ylab, y = ctrl$xlab) +
            guides(shape = guide_legend(override.aes = list(size = 2.7)))
    } else {return(ggplotGrob(gg))}
}

plotEvents.control <- function(xlab = "Time", ylab = "Subject", title = "Recurrent event plot",
                               terminal.name = "Terminal event",
                               recurrent.name = "Recurrent event",
                               recurrent.type = NULL, alpha = .7) {
    list(xlab = xlab, ylab = ylab, title = title, terminal.name = terminal.name,
         recurrent.name = recurrent.name, recurrent.type = recurrent.type,
         alpha = alpha)
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
