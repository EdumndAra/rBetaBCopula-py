#######################################################################################
##
##   R code by Ivan Kojadinovic and Bingqing Yi Copyright (C) 2023
##
##   The R implementation below is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This R code is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.
##
##   If you use this code, please cite R, the relevant R packages and the research
##   articles in which the algorithms below are introduced.
##
#######################################################################################

library(extraDistr)

##' @title Samping from the empirical beta copula
##'        see Kiriliouk, Segers and Tsukahara (2021)
##' @param n sample size
##' @param x matrix of d-dimensional observations (or multivariate ranks)
##'          from which the empirical beta copula will be computed
##' @param ranks logical indicating whether x contains observations
##'              or corresponding multivariate ranks
##' @return a random sample of size n from the empirical beta copula
##' @author Ivan Kojadinovic and Bingqing Yi

rBinCopula <- function(n, x, ranks = FALSE) {

    ## Checks
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(n >= 1L)

    ## Data dimensions
    m <- nrow(x)
    d <- ncol(x)

    ## Multivariate ranks from x
    if (!ranks)
        x <- apply(x, 2, rank)

    ## Generate n realizations
    v <- matrix(NA, n, d)
    for (i in 1:n) {
        I <- sample(1:m, size = 1) # select a row
        for (j in 1:d) {
            r <- x[I, j]
            v[i,j] <- rbeta(1, r, m + 1 - r)
        }
    }
    v
}

pbetab <- function(x, m, u, rho, lower.tail) {
    if (rho <= 1) {
        warning("rho is <= 1; using binomial instead")
        return(pbinom(x, size = m, prob = u, lower.tail = lower.tail))
    }
    ifelse((u == 0) | (u == 1),
           pbinom(x, size = m, prob = u, lower.tail = lower.tail), # degenerate cases
           pbbinom(x, size = m, alpha = (rho - m) * u / (1 - rho),
                   beta = (rho - m) / (1 - rho) * (1 - u),
                   lower.tail = lower.tail)
           )
}

kbetab <- function(u, r, m, rho = 4)
    pbetab(x = r - 1, m = m, u = u, rho = rho, lower.tail = FALSE)

kbetab.inv <- function(u, r, m, rho = 4) {
    f <- function(v)
        pbetab(x = r - 1, m = m, u = v, rho = rho, lower.tail = FALSE) - u
    uniroot(f, interval = 0:1)$root
}

##' @title Samping from the "BetaB4" smooth empirical copula
##' @param n sample size
##' @param x matrix of d-dimensional observations (or multivariate ranks)
##'          from which the "BetaB4" smooth empirical copula will be computed
##' @param ranks logical indicating whether x contains observations
##'              or corresponding multivariate ranks
##' @param rho parameter of the margins of the smoothing distributions;
##'            set to 4 by default
##' @return a random sample of size n from the "BetaB4" smooth empirical copula
##' @author Ivan Kojadinovic and Bingqing Yi

rBetaBCopula <- function(n, x, ranks = FALSE, rho = 4) {

    ## Checks
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(n >= 1L)

    ## Data dimensions
    m <- nrow(x)
    d <- ncol(x)

    ## Multivariate ranks from x
    if (!ranks)
        x <- apply(x, 2, rank)

    ## Generate n realizations from the empirical
    ## beta copula first
    u <- rBinCopula(n = n, x = x, ranks = TRUE)

    ## Transform these realizations marginally
    for (i in 1:n) {
        I <- sample(1:m, size = 1) # select a row
        for (j in 1:d) {
            r <- x[I, j]
            u[i,j] <- kbetab.inv(u[i,j], r = r, m = m, rho = rho)
        }
    }
    u
}
