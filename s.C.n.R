#######################################################################################
##
##   R code by Ivan Kojadinovic and Bingqing Yi Copyright (C) 2022
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
##   articles in which the estimators below are introduced.
##
#######################################################################################

library(copula)
library(extraDistr)

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

pbeta2 <- function(x, m, u, rho, ...)
    pbeta(x, shape1 = (m - rho) / rho * u, shape2 = (m - rho) / rho * (1 - u), ...)

##' @title Some possibly data-adaptive smooth empirical copulas
##' @param u matrix of d-dimensional evaluation points (points correspond to rows)
##' @param r matrix of d-dimensional multivariate ranks obtained from the data
##' @param marg smoothing margins: scaled binomial, scaled beta-binomial or beta
##' @param rho smoothing parameter for the scaled beta-binomial and beta margins
##' @param cop the smoothing d-dimensional survival copula
##' @return values of the smooth empirical copula at the rows of u
##' @author Ivan Kojadinovic and Bingqing Yi

s.C.n <- function(u, r, marg = c("binomial", "betabinomial", "beta"), rho = 4,
                  cop = indepCopula(dim = ncol(r))) {

    if(any(u < 0, 1 < u))
        stop("'u' must be in [0,1].")

    stopifnot((d <- ncol(u)) == ncol(r) && d == dim(cop))

    marg <- match.arg(marg)
    n <- nrow(r)
    m <- nrow(u)

    v <- matrix(NA, m, d)
    ec <- numeric(m)

    for (i in 1:n) {
        for (j in 1:d)
            v[,j] <- switch(marg,
                            binomial = pbinom(r[i,j] - 1, size = n, prob = u[,j], lower.tail = FALSE),
                            betabinomial = pbetab(r[i,j] - 1, m = n, u = u[,j], rho = rho, lower.tail = FALSE),
                            beta = pbeta2((r[i,j] - 0.5)/n, m = n, u = u[,j], rho = rho, lower.tail = FALSE))

        ec <- ec + pCopula(v, copula = cop)
    }
    ec / n
}

## Some examples
n <- 100
d <- 2
m <- 10
tau <- 0.5
theta <- iTau(gumbelCopula(), tau = tau)

data.cop <- gumbelCopula(theta, dim = d) # data generating copula
x <- rCopula(n, copula = data.cop) # data
u <- matrix(runif(m * d), m, d) # random evalatuion points
r <- apply(x, 2, rank) # multivariate ranks corresponding to x
ebc <- empCopula(X = r/(n+1), smoothing = "beta") # the empirical beta copula of x

## Some data-adaptive smooth empirical copulas evaluated at u
s.C.n(u, r, marg = "binomial", cop = ebc)
s.C.n(u, r, marg = "betabinomial", cop = ebc)
s.C.n(u, r, marg = "beta", cop = ebc)

## Tests to see whether they have standard uniform margins
v <- runif(m) # random points where to evaluate the margins of the estimators
w1 <- cbind(v, 1) # evaluations points margin 1
w2 <- cbind(1, v) # evaluations points margin 2
if (d > 2)
    for (j in 3:d) {
        w1 <- cbind(w1, 1)
        w2 <- cbind(w2, 1)
    }


stopifnot(all.equal(s.C.n(u = w1, r = r, marg = "binomial"), v))
stopifnot(all.equal(s.C.n(u = w2, r = r, marg = "binomial"), v))

stopifnot(all.equal(s.C.n(u = w1, r = r, marg = "betabinomial"), v))
stopifnot(all.equal(s.C.n(u = w2, r = r, marg = "betabinomial"), v))

print(all.equal(s.C.n(u = w1, r = r, marg = "beta"), v))
print(all.equal(s.C.n(u = w2, r = r, marg = "beta"), v))
