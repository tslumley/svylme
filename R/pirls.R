## Based on lme4pureR:pirls.R by Steve Walker and Doug Bates (GPL-2)
##
##  returns a function for evaluating the GLMM Laplace approximated deviance
##


svyseqglme<-function(formula, family, design, REML=FALSE, scale=c("sample_size","effective_sample_size","raw")){
    data<-model.frame(design)
    m0 <-lme4::glmer(formula, family=family,data=data)
   
    ## remove missing from design
    if (!is.null(naa<-attr(m0@frame,"na.action"))){
        design<-design[-naa,]
    }

    X<-m0@pp$X
    Zt<-m0@pp$Zt
    Lambdat<-m0@pp$Lambdat
    thfun<- with(m0@pp, function(theta) theta[Lind])
    y<-m0@resp$y
    offset<-m0@resp$offset
    n<-length(y)

    ## check nesting
    model_clusters <- m0@flist
    is.nested<-function(fine, coarse){
        if(!(length(fine)==length(coarse)))
            stop("vectors must be the same length")
        length(unique(fine)) == length(unique(paste(fine,coarse,sep=":")))
    }
    last.nested<-function(g){
        a<-which(apply(design$cluster, 2, is.nested, fine=g))
        if (length(a))
            max(a)
        else
            NA
        }
    u_depth<-lapply(model_clusters, last.nested)
    design_depth<-NCOL(design$cluster)

    ## Don't need to worry about whether u_depth==design_depth is ok
    ## (is ok if whole-cluster sampling, not if element sampling)
    ## because if it isn't, the lmer() call won't have converged.

    if (any(is.na(unlist(u_depth)))){
        print(u_depth)
        stop("Model clusters not nested in sampling units")
        }

    yweights<-weights(design,"sampling")

    scaling_method <- match.arg(scale)
    clweights<-scale_weights(design, method=scaling_method)
    uweights<-numeric(NROW(Zt))
    nlevels<-length(m0@flist)
    for(i in 1:nlevels){
        idx<-match(rownames(Zt), unique(m0@flist[[i]]))
        uweights[!is.na(idx)]<-clweights[[ u_depth[[i]] ]][na.omit(idx)]
    }
    ###
    aweights<-rep(1,n) ## for now

    ###
    devfun <- pirls(X,y,Zt,Lambdat, theta=m0@theta,beta=m0@beta,
            thfun = thfun, aweights=aweights,
            offset = offset, eta=m0@resp$eta, u=m0@u,
            family=family, yweights=yweights, uweights=uweights)
    varcomp<-bobyqa(c(m0@theta,m0@beta), devfun,
                    lower = c(m0@lower,rep(-Inf,length(m0@beta))), upper = Inf)

    ## sandwich estimator
    p<-NCOL(X)
    beta <-environment(devfun)$beta
    b<-environment(devfun)$b
    ##DD<-environment(devfun)$DD
    ##scores<-(X*as.vector(y-X%*%beta-crossprod(Zt,b)))%*%solve(DD)
    ##Vbeta<-vcov(svytotal(scores,design))
    Vbeta<-matrix(NA,p,p)

    ## put variance parameters in printable form
    qi<-sapply(m0@cnms,length)
    q<-sum(qi)
    VC<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i))))
    ###(need indicator for where thetas go in the matrix)
    ThInd<-which((VC==1) & lower.tri(VC,diag=TRUE))
    Th<-matrix(0,q,q)
    Th[ThInd]<-varcomp$par
    VC[]<-tcrossprod(Th)

    names(beta)<-colnames(X)
    dimnames(Vbeta)<-list(colnames(X),colnames(X))

    znames<-unlist(m0@cnms)
    names(znames)<-rep(names(m0@cnms),sapply(m0@cnms,length))
    
    rval<-list(devfun=devfun, varcomp=varcomp, beta=beta, call=sys.call(),
               s2=environment(devfun)$s2hat, znames=znames, Vbeta=Vbeta,
               VC=VC)
    class(rval)<-"svyseqlme"
    rval
}


pirls <- function(X,y,Zt,Lambdat,thfun,theta,beta,
                  aweights,offset=numeric(n),
                  eta,u,family=binomial,
                  tol = 10^-5, npirls = 50,
                  nstephalf = 10,verbose=1L,
                  yweights=rep(1,n),uweights=rep(1,n),
                  ...){

    n <- nrow(X); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n,
              nrow(Lambdat) == q, ncol(Lambdat) == q, is.function(thfun))

    if (is.function(family)) family <- family() # ensure family is a list

    local({
           Wuhalf<- Diagonal(x=sqrt(as.numeric(uweights)))
           Wyhalf <- Diagonal(x=sqrt(as.numeric(yweights)))
           W_u <- Diagonal(x=(as.numeric(uweights)))
           ZtW <- Zt %*% Wyhalf
           yweights<-yweights
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
        L <- Cholesky(tcrossprod(Lambdat %*% ZtW)+W_u, perm=FALSE, LDL=FALSE)
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
                    sum(sqDevResid(y, mu, aweights*yweights)) + sum(uweights*(uu^2))
                }
                u[] <<- numeric(q)
                olducden <- updatemu(u)
                cvgd <- FALSE
                for(i in 1:npirls){
                                        # update w and muEta
                    Whalf <- Diagonal(x = sqrt(aweights / variance(mu)))%*% Wyhalf
                                        # update weighted design matrix
                    LtZtMWhalf <- LtZt %*% (Diagonal(x = muEta(eta)) %*% Whalf )
                                        # alternative (more explicit but slower)
                                        # Cholesky update
                    L <- Cholesky(tcrossprod(LtZtMWhalf)+W_u, perm=FALSE, LDL=FALSE)
                                        # update weighted residuals
                    wtres <- Whalf %*% (y - mu)
                                        # solve for the increment
                    delu <- as.vector(solve(L, LtZtMWhalf %*% wtres - uweights*u))
                    if (verbose > 1L) {
                        #cat(sprintf("inc: %12.4g", delu[1]))
                        #nprint <- min(5, length(delu))
                        #for (j in 2:nprint) cat(sprintf(" %12.4g", delu[j]))
                        print(summary(delu))
                        #cat("\n")
                    }
                                        # update mu and eta and calculate
                                        # new unscaled conditional log density
                    ucden <- updatemu(u + delu)
                    if (verbose > 0L) {
                        cat(sprintf("%6.4f: %10.3f\n", i, ucden))
                    }

                    if(abs((olducden - ucden) / ucden) < tol){
                        cvgd <- TRUE
                        break
                    }
                                        # step-halving
                    if(ucden > olducden){
                        if (verbose>0L) cat(sprintf("%6.4f: %10.3f\n", 1, olducden))
                        for(j in 1:nstephalf){
                            ucden <- updatemu(u + (delu <- delu/2))
                            if (verbose > 0L) {
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
                if(!cvgd) warning("PIRLS failed to converge")

                                        # create Laplace approx to -2log(L)
                ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
                attributes(ldL2) <- NULL
                # FIXME: allow for quadrature approximations too
                Lm2ll <- aic(y,rep.int(1,n),mu,aweights*yweights,NULL) + sum(u^2*uweights) + ldL2 #+ (q/2)*log(2*pi)

                if (verbose > 0L) {
                    cat(sprintf("%10.3f: %12.4g", Lm2ll, thetabeta[1]))
                    for (j in 2:length(thetabeta)) cat(sprintf(" %12.4g", thetabeta[j]))
                    cat("\n")
                }

                Lm2ll
            }
    })
}
