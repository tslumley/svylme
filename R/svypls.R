##
## Based on lme4pureR:pls.R by Steve Walker, Doug Bates
## GPL-2
##
## aweights: precision weights as in lmer
## yweights: observation sampling weights
## uweights: sampling weights for the unit-variance random effects
##
## Only usable when there is a well-defined sampling weight for each u,
## ie, nested case. 

## library(lme4pureR)
## library(lme4)
## library(minqa)
## set.seed(1)
## n <- 1000
## x <- rnorm(n)
## z <- rnorm(n)
## X <- cbind(1, x)
## ZZ <- cbind(1, z)
## grp <- gl(n/5,5)
## RE <- mkRanefStructures(list(grp), list(ZZ))
## Z <- t(RE$Zt)
## y <- as.numeric(X%*%rnorm(ncol(X)) + Z%*%rnorm(ncol(Z)) + rnorm(n))
## m <- lmer.fit(y,X,ZZ,grp)
## m$par
## Lambdat <- RE$Lambdat
## Lambdat@x <- m$par[RE$Lind]
## cov2cor(crossprod(Lambdat)[1:2,1:2])
## lmer(y ~ x + (z|grp))


svyseqlme <- function(y, mmFE, mmRE, grp,
                     aweights, offset = numeric(n),
                     REML = TRUE, yweights, uweights){
    n<-length(y)
    if(missing(aweights)) aweights <- rep(1,length(y))
    initRE <- mkRanefStructures(grp, mmRE)
    devfun <- with(initRE, {
        pls(mmFE,y,Zt,Lambdat,
            thfun = function(theta) theta[Lind],
            aweights = aweights, offset = offset,
            REML = REML, yweights=yweights, uweights=uweights)})
    varcomp<-with(initRE, {
        bobyqa(initRE$theta, devfun,
               lower = lower, upper = upper)})
    list(devfun=devfun, varcomp=varcomp, beta=environment(devfun)$beta)
}


pls <- function(X,y,Zt,Lambdat,thfun,aweights,
                offset = numeric(n),REML = FALSE,
                yweights=rep(1,n),uweights=rep(1,n),...)
{
    # SW: how to test for sparse matrices, without specifying the specific class?
    stopifnot(is.matrix(X)) #  is.matrix(Zt), is.matrix(Lambdat))
    n <- length(y); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n,
              nrow(Lambdat) == q, ncol(Lambdat) == q, is.function(thfun))
                                        # calculate weighted products
    Wahalf <- if (missing(aweights)) Diagonal(n=n) else Diagonal(x=sqrt(as.numeric(aweights)))
    Wuhalf<- Diagonal(x=sqrt(as.numeric(uweights)))
    Whalf <- Diagonal(x=sqrt(as.numeric(yweights)))%*%Wahalf
    
    WX <- Whalf %*% X
    Wy <- Whalf %*% y
    ZtW <- Zt %*% Whalf
    XtWX <- crossprod(WX)
    XtWy <- crossprod(WX, Wy)
    ZtWX <- ZtW %*% WX
    ZtWy <- ZtW %*% Wy
    rm(WX,Wy)
    local({                             # mutable values stored in local environment
        b <- numeric(q)                 # conditional mode of random effects
        beta <- numeric(p)              # conditional estimate of fixed-effects
        cu <- numeric(q)                # intermediate solution
        DD <- XtWX                      # down-dated XtWX
        L <- Cholesky(tcrossprod(Lambdat %*% ZtW)+Wuhalf, LDL = FALSE)
        Lambdat <- Lambdat              # stored here b/c x slot will be updated
        mu <- numeric(n)                # conditional mean of response
        RZX <- matrix(0,nrow=q,ncol=p)  # intermediate matrix in solution
        u <- numeric(q)                 # conditional mode of spherical random effects
        s2hat <- numeric(1)
        function(theta) {
            Lambdat@x[] <<- thfun(theta)
            L <<- Cholesky(tcrossprod(Lambdat %*% ZtW)+Wuhalf, LDL = FALSE)
                                        #update(L, Lambdat %*% ZtW, mult = 1) ## does this need Wuhalf? TL
                                        # solve eqn. 30
            cu[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWy, system="P"),
                                     system="L"))
                                        # solve eqn. 31
            RZX[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWX, system="P"),
                                      system="L"))
            ## downdate XtWX and form Cholesky factor (eqn. 32)
            DD <<- as(XtWX - crossprod(RZX), "dpoMatrix")
            ## conditional estimate of fixed-effects coefficients (solve eqn. 33)
            beta[] <<- as.vector(solve(DD, XtWy - crossprod(RZX, cu)))
            ## conditional mode of the spherical random-effects coefficients (eqn. 34)
            u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta, system = "Lt"),
                                    system="Pt"))
            b[] <<- as.vector(crossprod(Lambdat,u))
                                        # conditional mean of the response
            mu[] <<- as.vector(crossprod(Zt,b) + X %*% beta + offset)
            wtres <- Whalf*(y-mu)       # weighted residuals
            pwrss <- sum(wtres^2) + sum((Wuhalf*u)^2) # penalized, weighted residual sum-of-squares
            fn <- as.numeric(length(mu))
            Nhat <- sum(yweights)
            s2hat <<- sum(wtres^2)/Nhat
            ld <- 2*determinant(L,logarithm=TRUE)$modulus # log determinant
            if (REML) {
                ld <- ld + determinant(DD,logarithm=TRUE)$modulus
                fn <- fn - length(beta)
            }
            attributes(ld) <- NULL
                                        # profiled deviance or REML criterion
            ld + fn*(1 + log(2*pi*pwrss) - log(Nhat))-determinant(Wuhalf,logarithm=TRUE)$modulus
            #ld + fn*(1 + log(2*pi*pwrss))
        }
    })
}
