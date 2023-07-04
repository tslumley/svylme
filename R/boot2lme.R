
boot2lme<-function(model, rdesign, verbose=FALSE){

    if(is.null(model$devfun)) stop("model must be fitted with return.devfun=TRUE")
    
    naa<-environment(model$devfun)$naa
    if (!is.null(naa)){
        if (length(environment(model$devfun)$y)+length(naa) == NROW(rdesign))
            rdesign<-rdesign[-naa,,drop=FALSE]
        if(verbose) warning(paste(length(naa),"observations dropped because of missing values"))
        }
    if (length(environment(model$devfun)$y) != NROW(rdesign)){
        stop("number of rows of design does not match model")
    }

    basewts<-weights(rdesign, "sampling")
    replicates<-weights(rdesign, "analysis")
    scale<-rdesign$scale
    rscales<-rdesign$rscales

     
    nrep<-ncol(replicates)
    pwt0<-if (model$method=="nested") get("pwts",environment(model$devfun)) else get("pwt",environment(model$devfun))
    if (is.null(rscales)) rscales<-rep(1,nrep)

    ii<-get("ii", environment(model$devfun))
    jj<-get("jj", environment(model$devfun))  
    repwt<-(replicates/basewts)[ii,]
    repwtj<-(replicates/basewts)[jj,]
    if ((model$method=="nested") && (any(abs((repwt-repwtj)/(1+repwt+repwtj))>1e-5)))
        warning("replicate weights vary within cluster")
    else {
        repwt<-sqrt(repwt*repwtj)  
    }
    
    theta0<-model$opt$par
    thetastar<-matrix(nrow=nrep,ncol=length(theta0))
    betastar<-matrix(nrow=nrep,ncol=length(model$beta))
    s2star<-numeric(nrep)

    D<-get("L",environment(model$devfun))
    Dstar<-array(0,c(nrep,NROW(D),NCOL(D)))
    
    if (verbose) pb<-txtProgressBar(min = 0, max = nrep, style = 3)

    for(i in 1:nrep){
        if (verbose) setTxtProgressBar(pb, i)
        if (model$method=="nested"){
        thetastar[i,]<-bobyqa(theta0, model$devfun,
                              lower = model$lower,
                              upper = rep(Inf, length(theta0)), pwt=repwt[,i]*pwt0)$par
        } else {
            thetastar[i,]<-bobyqa(theta0, model$devfun,
                              lower = model$lower,
                              upper = rep(Inf, length(theta0)), pwt_new=repwt[,i]*pwt0, subtract_margins=model$subtract.margins)$par
         }
        betastar[i,]<-get("beta",environment(model$devfun))
        s2star[i]<-get("s2",environment(model$devfun))
        Dstar[i,,]<-get("L",environment(model$devfun))
    }

    if(verbose) close(pb)
    
    rval<-list(theta=thetastar, beta=betastar, s2=s2star, D=Dstar,scale=scale, rscales=rscales, formula=model$formula)

    class(rval)<-"boot2lme"
    rval   
}

print.boot2lme<-function(x,...){
    cat("boot2lme:",length(x$s2),"replicates from", deparse(x$formula))
    invisible(x)
}

vcov.boot2lme<-function(object, parameter=c("beta","theta","s2","relSD","SD","relVar","fullVar"),...){
    parameter<-match.arg(parameter)

    nthetas<-NCOL(object$theta)

    if (nthetas==1){
           V<-switch(parameter,
                      beta=svrVar(object$beta, object$scale,object$rscales),
                      theta=svrVar(object$theta, object$scale,object$rscales),
                      s2=svrVar(object$s2, object$scale, object$rscales),
                      relSD=svrVar(sqrt((apply(object$D,1, diag))), object$scale, object$rscales), ##FIXME: dimension decay when there's just one random effect
                      SD=svrVar(sqrt((apply(object$D,1, diag))*object$s2), object$scale, object$rscales),
                      relVar=svrVar((apply(object$D,1,c)), object$scale, object$rscales),
                      fullVar=svrVar((apply(object$D,1,c))*object$s2, object$scale, object$rscales)
                      )

        } else {
            V<-switch(parameter,
                      beta=svrVar(object$beta, object$scale,object$rscales),
                      theta=svrVar(object$theta, object$scale,object$rscales),
                      s2=svrVar(object$s2, object$scale, object$rscales),
                      relSD=svrVar(sqrt(t(apply(object$D,1, diag))), object$scale, object$rscales), ##FIXME: dimension decay when there's just one random effect
                      SD=svrVar(sqrt(t(apply(object$D,1, diag))*object$s2), object$scale, object$rscales),
                      relVar=svrVar(t(apply(object$D,1,c)), object$scale, object$rscales),
                      fullVar=svrVar(t(apply(object$D,1,c))*object$s2, object$scale, object$rscales)
                      )
    }

    as.matrix(V)

    }
