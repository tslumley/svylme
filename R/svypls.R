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
##
## currently only for the two-stage case:
##    multistage needs better bookkeeping to match u to uweights


scale_weights<-function(design, method){
    m<-NCOL(design$allprob)
    if (method=="raw"){
        w<-vector("list",m-1)
        wi<-t(apply(design$allprob,1, cumprod))
        for(i in seq_len(m-1)){
            w[[i]]<-wi[!duplicated(design$cluster[,i]),i]
        }
        w
    } else if(method=="sample_size"){
        w<-vector("list",m-1)
        wi<-t(apply(1/design$allprob,1, cumprod))
         for(i in seq_len(m-1)){
            ns<-design$fpc$sampsize[,i+1]
            if (!is.null(design$fpc$popsize))
                Ns<-design$fpc$popsize[,i+1]
            else
                Ns<-ave(1/design$allprob[,i+1], design$cluster[,i], FUN=sum)
            w[[i]]<-(wi[,i]*Ns/ns)[!duplicated(design$cluster[,i])]
         }
        w
    }
    else stop("not implemented yet")
}


svyseqlme<-function(formula, design, REML=FALSE, scale=c("sample_size","effective_sample_size","gk","raw")){
    data<-model.frame(design)
    m0 <-lme4::lmer(formula, data, REML=REML)
   
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

    yweights<-weights(design,"sampling")

    scaling_method <- match.arg(scale)
    clweights<-scale_weights(design, method=scaling_method)
    uweights<-rep(clweights[[1]], each=length(m0@cnms[[1]]))  #FIXME only for two-stage
    
    
    devfun <- pls(X,y,Zt,Lambdat,
            thfun = thfun,
            aweights = rep(1,n), offset = offset,
            REML = REML, yweights=yweights, uweights=uweights)
    varcomp<-bobyqa(m0@pp$theta, devfun,
                    lower = m0@lower, upper = Inf)

    ## sandwich estimator
    warning("Sandwich variance estimator may still contain nuts")
    p<-NCOL(X)
    beta <-environment(devfun)$beta
    b<-environment(devfun)$b
    DD<-environment(devfun)$DD
    scores<-(X*as.vector(y-X%*%beta-crossprod(Zt,b)))%*%solve(DD)
    Vbeta<-vcov(svytotal(scores,design))

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
    
    rval<-list(devfun=devfun, varcomp=varcomp, beta=beta, call=sys.call(),
               s2=environment(devfun)$s2hat, znames=m0@cnms[[1]], Vbeta=Vbeta,
               VC=VC)
    class(rval)<-"svyseqlme"
    rval
}




pls <- function(X,y,Zt,Lambdat,thfun,aweights,
                offset = numeric(n),REML = FALSE,
                yweights=rep(1,n),uweights=rep(1,n),...)
{
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
    W_u <- Diagonal(x=(as.numeric(uweights)))
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
        L <- Cholesky(tcrossprod(Lambdat %*% ZtW)+W_u, LDL = FALSE)
        Lambdat <- Lambdat              # stored here b/c x slot will be updated
        mu <- numeric(n)                # conditional mean of response
        RZX <- matrix(0,nrow=q,ncol=p)  # intermediate matrix in solution
        u <- numeric(q)                 # conditional mode of spherical random effects
        s2hat <- numeric(1)

        W_y<- Diagonal(x=(as.numeric(yweights)*as.numeric(aweights)))
        WX <- Whalf %*% X               #used in sandwich estimator
        Wy <- Whalf %*% y
        ZtW <- Zt %*% Whalf

        function(theta) {
            Lambdat@x[] <<- thfun(theta)
            L <<- Cholesky(tcrossprod(Lambdat %*% ZtW)+W_u, LDL = FALSE)
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
            s2hat <<- sum(pwrss)/Nhat
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



 
print.svyseqlme<-function(x,digits=max(3L, getOption("digits") - 3L),...){
    cat("Linear mixed model fitted by sequential pseudolikelihood\n")
    cat("Formula: ")
    cat(paste(deparse(x$call$formula),collapse="\n"))
    cat("\nRandom effects:\n")
    theta<-x$varcomp$par
    s<-sqrt(as.vector(x$s2))
    stdev<- matrix(s*sqrt(diag(x$VC)),ncol=1)
    rownames(stdev)<-x$znames
    colnames(stdev)<-"Std.Dev."
    print(round(stdev,digits))
    cat("Residual:\t",round(s,digits))
    cat("\n Fixed effects:\n")
    coef<- cbind(beta=x$beta,SE=sqrt(diag(x$Vbeta)),t=x$beta/sqrt(diag(x$Vbeta)))
    coef<-cbind(coef,p=2*pnorm(-abs(coef[,3])))
    colnames(coef)<-c("beta","SE","t","p")
    print(round(coef,digits))
    cat("\n")
    invisible(x)
    }


coef.svyseqlme<-function(object,...,random=FALSE){
    if (random) {
        VC<-x$VC
        s2<-drop(object$s2)
        dimnames(L)<-list(object$znames,object$znames)
        list(s2=s2, varb=VC*s2)
    } else 
        drop(object$beta)
}

vcov.svyseqlme<-function(object,...){
    object$Vbeta
}
