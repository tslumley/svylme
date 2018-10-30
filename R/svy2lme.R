

boot2lme<-function(model, basewts, replicates, scale, rscales=NULL,  verbose=FALSE){

    nrep<-ncol(replicates)
    pwt0<-get("pwts",environment(model$devfun))
    if (is.null(rscales)) rscales<-rep(1,nrep)

    ii<-get("ii", environment(model$devfun))
    jj<-get("jj", environment(model$devfun))  
    repwt<-(replicates/basewts)[ii,]
    repwtj<-(replicates/basewts)[jj,]
    if (any(abs((repwt-repwtj)/(1+repwt+repwtj))>1e-5)) warning("replicate weights vary within cluster")

    theta0<-model$opt$par
    thetastar<-matrix(nrow=nrep,ncol=length(theta0))
    betastar<-matrix(nrow=nrep,ncol=length(model$beta))
    s2star<-numeric(nrep)

    D<-get("L",environment(model$devfun))
    Dstar<-array(0,c(nrep,NROW(D),NCOL(D)))
    
    if (verbose) pb<-txtProgressBar(min = 0, max = nrep, style = 3)

    for(i in 1:nrep){
        if (verbose) setTxtProgressBar(pb, i)
        thetastar[i,]<-bobyqa(theta0, model$devfun,
                              lower = model$lower,
                              upper = rep(Inf, length(theta0)), pwt=repwt[,i]*pwt0)$par
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
    cat("boot2lme:",length(x$s2star),"replicates from", deparse(x$formula))
    invisible(x)
}

vcov.boot2lme<-function(object, parameter=c("beta","theta","s2","relSD","SD","relVar","fullVar"),...){
    parameter<-match.arg(parameter)

    switch(parameter,
           beta=svrVar(object$beta, object$scale,object$rscales),
           theta=svrVar(object$theta, object$scale,object$rscales),
           s2=svrVar(object$s2, object$scale, object$rscales),
           relSD=svrVar(sqrt(t(apply(object$D,1, diag))), object$scale, object$rscales),
           SD=svrVar(sqrt(t(apply(object$D,1, diag))*object$s2), object$scale, object$rscales),
           relVar=svrVar(t(apply(object$D,1,c)), object$scale, object$rscales),
           fullVar=svrVar(t(apply(object$D,1,c))*object$s2, object$scale, object$rscales)
           )

    }

is_close<-function(a,b, tolerance=1e-5){
    all(abs((as.matrix(a)-as.matrix(b))/(as.matrix(a)+as.matrix(b)))<tolerance)
    }

pi_from_design<-function(design, ii,jj){
    if (NCOL(design$allprob)==1){
        ## No multistage weights
        if (NCOL(design$cluster)>1)
            stop("you need weights/probabilities for each stage of sampling")
        
        if (NCOL(design$cluster)==1 && !any(duplicated(design$cluster))){
            ## ok, element sampling
            if(is.null(design$fpc$popsize)) #with replacement
                return(list(full=design$prob[ii]*design$prob[jj],
                       first=design$prob[ii],
                       cond=rep(1,length(ii))))
            else if(is_close(as.vector(design$allprob),
                              as.vector(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
                # srs, possibly stratified
                n<-design$fpc$sampsize
                N<-design$fpc$popsize
                return(list(full= n[ii]*(n[jj]-1)/( N[ii]*(N[jj]-1)),
                            first=n[ii]/N[ii],
                            cond=rep(1,length(ii))))
            } else {
                ## Hajek high entropy: Brewer p153
                pi<-design$allprob
                denom<-ave(1-pi, design$strata,FUN=sum)
                samestrata<-(design$strata[ii]==design$strata[jj])
                return(list(full=pi[ii]*pi[jj]*(1- ifelse(samestrata, (1-pi[ii])*(1-pi[jj])/denom, 0)),
                            first=pi[ii],
                            cond=rep(1,length(ii))))
            }
        } else if (all(by(design$prob, design$cluster[,1], function(x) length(unique(x)))==1)) {
            ## possibly ok, sampling of whole clusters
            warning("assuming no subsampling within clusters because multi-stage weights were not given")
             if(is.null(design$fpc$popsize)) #with replacement
                 return(list(full=design$prob[ii],
                             first=design$prob[ii],
                             cond=rep(1, length(ii))))
            else if(is_close(as.vector(design$allprob),
                              as.vector(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
                # srs, possibly stratified
                n<-design$fpc$sampsize
                N<-design$fpc$popsize
                return(list(full= n[ii]/N[ii],
                            first=n[ii]/N[ii],
                            cond=rep(1,length(ii))))
            } else {
                ## Hajek high entropy: Brewer p153
                pi<-design$allprob
                denom<-ave(1-pi, design$strata,FUN=sum)
                samestrata<-(design$strata[ii,1]==design$strata[jj,1])
                return(list(full=pi[ii,1],
                            first=pi[ii,1],
                            cond=rep(1,length(ii))))
            }
        } else {
            ## not ok
            stop("you need weights/probabilities for each stage of sampling")
        }       
    }

    ## If we're here, we have multistage weights
    if (ncol(design$allprob)!=ncol(design$cluster)){
        ## ? can't happen
        stop("number of stages of sampling does not match number of stages of weights")
    }

    if(is.null(design$fpc$popsize)){ #with replacement
        last<-ncol(design$allprob)
        return(list(full=design$prob[ii]*design$allprob[jj,last],
                    first=apply(design$allprob[ii,-last, drop=FALSE], 1, prod),
                    cond=design$allprob[ii,last]*design$allprob[jj,last]))
    }
    if(all.equal(as.matrix(design$allprob), as.matrix(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
        ## multistage stratified random sampling
        last<-ncol(design$allprob)
        n<-design$fpc$sampsize
        N<-design$fpc$popsize
        pstages <- n[ii,]*(n[jj,]-1)/(N[ii,]*(N[jj,]-1))
        return(list(full=apply((n[ii,]/N[ii,])[,-last,drop=FALSE],1,prod)*pstages[,last],
                    first=apply((n[ii,]/N[ii,])[,-last,drop=FALSE],1,prod),
                    cond=pstages[,last]))
    }

    ## Hajek high entropy: Brewer p153
    first<-cpwt<-rep_len(1,length(ii))
    for (i in 1:ncol(allprob)){
        pi<-design$allprob[,i]
        denom<-ave(1-pi, design$strata[,i],FUN=sum)
        samestrata<-(design$strata[ii,i]==design$strata[jj,i])
        if (i==ncol(allprob))
            cpwt<-cpwt*pi[ii]*pi[jj]*(1- ifelse(samestrata, (1-pi[ii])*(1-pi[jj])/denom, 0))
        else
            first<-first*pi[ii]
    }
    return(list(full=first*cpwt, first= first, cond=cpwt))

}

svy2lme<-function(formula,design,sterr=TRUE, return.devfun=FALSE){

    data<-model.frame(design)
    
    ## Use unweighted model to get starting values and set up variables
    m0<-lme4::lmer(formula,data,REML=FALSE)

    ## Extract varables
    y<-m0@resp$y
    X<-m0@pp$X

    ## cluster indicator
    g<-m0@flist[[1]]
    ## number of clusters
    n1<-length(unique(g))
    ## number of variance parameters
    q<-nrow(m0@pp$Zt)/n1

    ## Z in lme4::lmer has separate columns for each cluster; restructure it
    Z<-crossprod(m0@pp$Zt, (outer(1:(n1*q),1:q,function(i,j) ((i-j) %% q)==0)*1))

    n<-NROW(X)

    ## all pairs within same cluster
    ij<-subset(expand.grid(i=1:n,j=1:n), (g[i] == g[j]) & (i<j))
    ## columns of indices for first and second observation in a pair
    ii<-ij[,1]
    jj<-ij[,2]
    npairs<-nrow(ij)
    
    p<-NCOL(X)
    
    ## starting values from the unweighted model
    s2<-m0@devcomp$cmp["sigmaML"]^2
    theta<-theta0<- m0@theta
    beta<-beta0<-lme4::fixef(m0)

    ## second-order weights
    allpwts<-pi_from_design(design,ii,jj)
    pwts<-1/allpwts$full
    pwt2<-1/allpwts$cond
    p1<-allpwts$first
    
    ## if (is.null(N2)){
    ##     ## using probabilities
    ##     pwt2 <- (1/p2[ii])*(1/p2[jj])
    ##     pwts<- (1/p1[ii])*pwt2  ## with replacement at stage 2
    ## } else {
    ##     ## using cluster size
    ##     n2<-ave(as.numeric(g), g, FUN=length)
    ##     pwt2<-N2[ii]*(N2[jj]-1)/(n2[ii]*(n2[ii]-1))
    ##     pwts<-(1/p1[ii])*pwt2  ## SRS without replacement at stage 2
    ## }

    ## variance matrix of random effects
    qi<-sapply(m0@cnms,length)
    L<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i))))
    ###(need indicator for where thetas go in the matrix)
    ThInd<-which((L==1) & lower.tri(L,diag=TRUE))
    
    ## profile pairwise deviance
    devfun<-function(theta,pwt){
        ## variance parameters: Cholesky square root of variance matrix
        Th<-matrix(0,q,q)
        Th[ThInd]<-theta
        L[]<<-tcrossprod(Th)

        ## v11 is a vector of (1,1) entries of the matrix var(Y)
        ## for each pair, similarly for the others
        v11<-(rowSums(Z[ii,,drop=FALSE]*( Z[ii,,drop=FALSE]%*%L))+1)
        v12<-rowSums(Z[ii,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))
        v22<-(rowSums(Z[jj,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))+1)
        ## explicit 2x2 determinants
        det<-v11*v22-v12*v12
        ## explicit 2x2 inverses
        inv11<- v22/det
        inv22<- v11/det
        inv12<- -v12/det

        ## X matrices for first and second element of each pair
        Xii<-X[ii,,drop=FALSE]
        Xjj<-X[jj,,drop=FALSE]

        ## X^TWX
        xtwx<- crossprod(Xii,pwt*inv11*Xii)+
            crossprod(Xjj,pwt*inv22*Xjj)+
            crossprod(Xii,pwt*inv12*Xjj)+
            crossprod(Xjj,pwt*inv12*Xii)

        ## X^WY
        xtwy<-crossprod(Xii,pwt*inv11*y[ii])+
            crossprod(Xjj,pwt*inv22*y[jj])+
            crossprod(Xii,pwt*inv12*y[jj])+
            crossprod(Xjj,pwt*inv12*y[ii])

        ## betahat at the given variance parameter values
        beta<<-solve(xtwx,xtwy)
        Xbeta<-X%*%beta

        ## two residuals per pair
        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]

        ## -2 times Gaussian log profile pairwise likelihood
        qf<-crossprod(r1,pwt*inv11*r1)+
            crossprod(r2,pwt*inv22*r2)+
            crossprod(r1,pwt*inv12*r2)+
            crossprod(r2,pwt*inv12*r1)

        Nhat<-sum(pwt)*2
        s2<<-qf/Nhat
        
        ##sum(log(det)*pwt) + qf 
        sum(log(det)*pwt) + Nhat*log(qf*2*pi/Nhat)
        
    }

    ## Standard errors of regression parameters
    Vbeta<-function(theta,pwt){
        ## setup exactly as in devfun
        Th<-matrix(0,q,q)
        Th[ThInd]<-theta
        L<<-tcrossprod(Th)
        
        v11<-(rowSums(Z[ii,,drop=FALSE]*( Z[ii,,drop=FALSE]%*%L))+1)
        v12<-rowSums(Z[ii,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))
        v22<-(rowSums(Z[jj,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))+1)
        det<-v11*v22-v12*v12
        inv11<- v22/det
        inv22<- v11/det
        inv12<- -v12/det
        
        Xii<-X[ii,,drop=FALSE]
        Xjj<-X[jj,,drop=FALSE]
        
        xtwx<- crossprod(Xii,pwt*inv11*Xii)+
            crossprod(Xjj,pwt*inv22*Xjj)+
            crossprod(Xii,pwt*inv12*Xjj)+
            crossprod(Xjj,pwt*inv12*Xii)
        
        Xbeta<-X%*%beta
        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]

        ## score for betas
        xwr<-Xii*pwt2*(inv11*r1)+
            Xjj*pwt2*(inv22*r2)+
            Xii*pwt2*(inv12*r2)+
            Xjj*pwt2*(inv12*r1)

        ## cluster weights
        p1g<-p1[!duplicated(g[ii])]

        if (is.null(design)){
            ## sandwich estimator
            J<-crossprod((1/p1g)*rowsum(xwr,g[ii],reorder=FALSE)*sqrt(n1/(n1-1)))
            G<-solve(xtwx)
            G%*%J%*%G
        } else {
            inffun<-(1/p1g)*rowsum(xwr,g[ii],reorder=FALSE)%*%solve(xtwx)
            PSUg<-design$cluster[,1][ii[!duplicated(g[ii])]]
            
            inffun<-rowsum(inffun, PSUg,reorder=FALSE)
            stratPSU<-design$strata[,1][ii[!duplicated(design$cluster[,1][ii])]] ##CHECKME?

            one<-rep(1,NROW(inffun))
            ni<-ave(one,stratPSU,FUN=NROW)
            centered<- inffun-ave(inffun, stratPSU, FUN=mean)
            crossprod(centered*sqrt(ni/(ni-1)))
        }
    }
    

    ## Powell's derivative-free quadratic optimiser
    fit<-bobyqa(theta0, devfun,
                lower = m0@lower,
                upper = rep(Inf, length(theta)), pwt=pwts)

    ## variance of betas, if wanted
    Vbeta<-if (sterr) Vbeta(fit$par,pwts) else matrix(NA,q,q)
    
    ## return all the things
    rval<-list(opt=fit,
               s2=s2,
               beta=beta,
               Vbeta=Vbeta,
               formula=formula,
               znames=do.call(c,m0@cnms),
               L=L)
    
    ## for resampling
    if(return.devfun) {
        rval$devfun<-devfun
        rval$lower<-m0@lower
        }
    
    class(rval)<-"svy2lme"
    rval
}
        

 
print.svy2lme<-function(x,digits=max(3L, getOption("digits") - 3L),...){
    cat("Linear mixed model fitted by pairwise likelihood\n")
    cat("Formula: ")
    cat(paste(deparse(x$formula),collapse="\n"))
    cat("\nRandom effects:\n")
    theta<-x$opt$par
    s<-sqrt(as.vector(x$s2))
    stdev<- matrix(s*sqrt(diag(x$L)),ncol=1)
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


coef.svy2lme<-function(object,...,random=FALSE){
    if (random) {
        L<-object$L
        s2<-drop(object$s2)
        dimnames(L)<-list(object$znames,object$znames)
        list(s2=s2, varb=L*s2)
    } else 
        object$beta
}

vcov.svy2lme<-function(object,...){
    object$Vbeta
}
