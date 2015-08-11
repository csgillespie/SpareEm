#' @export
legFun <-
function (a = -1, b = 1) 
{
    legList <- list(sqrt(1/2) * c(1), sqrt(3/2) * c(0, 1), sqrt(5/2) * 
        (1/2) * c(-1, 0, 3), sqrt(7/2) * (1/2) * c(0, -3, 0, 
        5), sqrt(9/2) * (1/8) * c(3, 0, -30, 0, 35), sqrt(11/2) * 
        (1/8) * c(0, 15, 0, -70, 0, 63), sqrt(13/2) * (1/16) * 
        c(-5, 0, 105, 0, -315, 0, 231), sqrt(15/2) * (1/16) * 
        c(0, -35, 0, 315, 0, -693, 0, 429), sqrt(17/2) * (1/128) * 
        c(35, 0, -1260, 0, 6930, 0, -12012, 0, 6435), sqrt(19/2) * 
        (1/128) * c(0, 315, 0, -4620, 0, 18018, 0, -25740, 0, 
        12155), sqrt(21/2) * (1/256) * c(0, -63, 0, 3465, 0, 
        -30030, 0, 90090, 0, -109395, 0), sqrt(23/2) * (1/256) * 
        c(46189, 0, -693, 0, 15015, 0, -90090, 0, 218790, 0, 
            -230945, 0), sqrt(25/2) * (1/1024) * c(88179, 0, 
        231, 0, -18018, 0, 225225, 0, -1021020, 0, 2078505, 0, 
        -1939938))
    p <- length(legList) - 1
    coef <- sapply(0:p, function(i) c(legList[[i + 1]], rep(0, 
        p - i)))
    coef <- coef * sqrt(2/(b - a))
    rm(list = "legList")
    function(X, terms) {
        X <- as.matrix(X)
        X <- 2 * (X - a)/(b - a) - 1
        if (any(X < -1 | X > 1)) 
            warning(sprintf("Value(s) in \"X\" outside the interval [%.2f, %.2f]", 
                a, b))
        terms <- as.matrix(terms)
        if (ncol(X) != ncol(terms)) 
            stop("\"X\" and \"terms\" should have the same number of columns")
        if (any(terms != round(terms) | terms < 0 | terms > p)) 
            stop("Inappropriate values in \"terms\"")
        tmp <- outer(X, 0:p, "^")
        dd <- dim(tmp)
        dim(tmp) <- c(prod(dd[1:2]), dd[3])
        for (i in p:0) tmp[, i + 1] <- tmp %*% coef[, i + 1]
        dim(tmp) <- dd
        robj <- sapply(1:nrow(terms), function(i) {
            foo <- sapply(1:dd[2], function(j) tmp[, j, 1 + terms[i, 
                j]])
            dim(foo) <- dd[1:2]
            apply(foo, 1, prod)
        })
        dim(robj) <- c(dd[1], nrow(terms))
        robj
    }
}
