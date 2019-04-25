## Based on lme4pureR:pirls.R by Steve Walker and Doug Bates (GPL-2)
##
##  returns a function for evaluating the GLMM Laplace approximated deviance
##
pirls <- function(X,y,Zt,Lambdat,thfun,theta,
                  weights,offset=numeric(n),
                  eta=numeric(n),family=binomial,
                  tol = 10^-6, npirls = 30,
                  nstephalf = 10,verbose=0L,
                  ...){
    # FIXME: bad default starting value for eta

    n <- nrow(X); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n,
              nrow(Lambdat) == q, ncol(Lambdat) == q, is.function(thfun))

    if (is.function(family)) family <- family() # ensure family is a list

    local({    
        nth <- length(theta)
        betaind <- -seq_len(nth) # indices to drop 1:nth
        linkinv <- family$linkinv
        variance <- family$variance
        muEta <- family$mu.eta
        aic <- family$aic
        sqDevResid <- family$dev.resid
        mu <- linkinv(eta)
        beta <- numeric(p)
        u <- numeric(q)
        L <- Cholesky(tcrossprod(Lambdat %*% Zt), perm=FALSE, LDL=FALSE, Imult=1)
                                        # create function for conducting PIRLS
            function(thetabeta) {
                                        # initialize
                Lambdat@x[] <<- thfun(thetabeta[-betaind])
                LtZt <- Lambdat %*% Zt
                beta[] <<- thetabeta[betaind]
                offb <- offset + X %*% beta
                updatemu <- function(uu) {
                    eta[] <<- offb + as.vector(crossprod(LtZt, uu))
                    mu[] <<- linkinv(eta)
                    sum(sqDevResid(y, mu, weights)) + sum(uu^2)
                }
                u[] <<- numeric(q)
                olducden <- updatemu(u)
                cvgd <- FALSE
                for(i in 1:npirls){
                                        # update w and muEta
                    Whalf <- Diagonal(x = sqrt(weights / variance(mu)))
                                        # update weighted design matrix
                    LtZtMWhalf <- LtZt %*% (Diagonal(x = muEta(eta)) %*% Whalf)
                                        # update Cholesky decomposition
                    L <- update(L, LtZtMWhalf, 1)
                                        # alternative (more explicit but slower)
                                        # Cholesky update
                    # L <- Cholesky(tcrossprod(LtZtMWhalf), perm=FALSE, LDL=FALSE, Imult=1)
                                        # update weighted residuals
                    wtres <- Whalf %*% (y - mu)
                                        # solve for the increment
                    delu <- as.vector(solve(L, LtZtMWhalf %*% wtres - u))
                    if (verbose > 0L) {
                        cat(sprintf("inc: %12.4g", delu[1]))
                        nprint <- min(5, length(delu))
                        for (j in 2:nprint) cat(sprintf(" %12.4g", delu[j]))
                        cat("\n")
                    }
                                        # update mu and eta and calculate
                                        # new unscaled conditional log density
                    ucden <- updatemu(u + delu)
                    if (verbose > 1L) {
                        cat(sprintf("%6.4f: %10.3f\n", 1, ucden))
                    }

                    if(abs((olducden - ucden) / ucden) < tol){
                        cvgd <- TRUE
                        break
                    }
                                        # step-halving
                    if(ucden > olducden){
                        for(j in 1:nstephalf){
                            ucden <- updatemu(u + (delu <- delu/2))
                            if (verbose > 1L) {
                                cat(sprintf("%6.4f: %10.3f\n", 1/2^j, ucden))
                            }
                            if(ucden <= olducden) break
                        }
                        if(ucden > olducden) stop("Step-halving failed")
                    }
                                        # set old unscaled conditional log density
                                        # to the new value
                    olducden <- ucden
                                        # update the conditional modes (take a step)
                    u[] <<- u + delu
                }
                if(!cvgd) stop("PIRLS failed to converge")

                                        # create Laplace approx to -2log(L)
                ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
                attributes(ldL2) <- NULL
                # FIXME: allow for quadrature approximations too
                Lm2ll <- aic(y,rep.int(1,n),mu,weights,NULL) + sum(u^2) + ldL2 #+ (q/2)*log(2*pi)

                if (verbose > 0L) {
                    cat(sprintf("%10.3f: %12.4g", Lm2ll, thetabeta[1]))
                    for (j in 2:length(thetabeta)) cat(sprintf(" %12.4g", thetabeta[j]))
                    cat("\n")
                }

                Lm2ll
            }
    })
}
