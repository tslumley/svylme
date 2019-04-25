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
    } else if(method=="effective_sample_size"){ ##FIXME check
        w<-vector("list",m-1)
        wi<-t(apply(1/design$allprob,1, cumprod))
         for(i in seq_len(m-1)){
            ns<-design$fpc$sampsize[,i+1]
            Swsq<-ave(1/design$allprob[,i+1], design$cluster[,i], FUN=sum)
            Sws<-ave(1/design$allprob[,i+1], design$cluster[,i], FUN=sum)
            w[[i]]<-(wi[,i]*Swsq/Sws)[!duplicated(design$cluster[,i])]
         }
        w
    } else stop("not implemented yet")
}


svyseqlme<-function(formula, design, REML=FALSE, scale=c("sample_size","effective_sample_size","raw")){
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
    
    devfun <- pls(X,y,Zt,Lambdat,
            thfun = thfun,
            aweights = rep(1,n), offset = offset,
            REML = REML, yweights=yweights, uweights=uweights)
    varcomp<-bobyqa(m0@pp$theta, devfun,
                    lower = m0@lower, upper = Inf)

    ## sandwich estimator
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

    znames<-unlist(m0@cnms)
    names(znames)<-rep(names(m0@cnms),sapply(m0@cnms,length))
    
    rval<-list(devfun=devfun, varcomp=varcomp, beta=beta, call=sys.call(),
               s2=environment(devfun)$s2hat, znames=znames, Vbeta=Vbeta,
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
            ## profiled deviance or REML criterion
            ld + fn*(1 + log(2*pi*pwrss) - log(Nhat))-determinant(Wuhalf,logarithm=TRUE)$modulus
            #ld + fn*(1 + log(2*pi*pwrss))
        }
    })
}



 
print.svyseqlme<-function(x,digits=max(3L, getOption("digits") - 3L),...){
    cat("Linear mixed model fitted by sequential pseudolikelihood\n")
    cat("Call: ")
    cat(paste(deparse(x$call),collapse="\n"))
    cat("\nRandom effects:\n")
    theta<-x$varcomp$par
    s<-sqrt(as.vector(x$s2))
    stdev<- rbind(matrix(s*sqrt(diag(x$VC)),ncol=1),round(s, digits))
    rownames(stdev)<-c(paste(names(x$znames),x$znames), "Residual:")
    colnames(stdev)<-"Std.Dev."
    print(round(stdev, digits))

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
        VC<-object$VC
        s2<-drop(object$s2)
        dimnames(VC)<-list(object$znames,object$znames)
        list(s2=s2, varb=VC*s2)
    } else 
        drop(object$beta)
}

vcov.svyseqlme<-function(object,...){
    object$Vbeta
}
