## like all.equal only specialised to matrices and returns logical
is_close<-function(a,b, tolerance=1e-5){
    all(abs((as.matrix(a)-as.matrix(b))/(as.matrix(a)+as.matrix(b)))<tolerance)
}


## safe devision
"%//%"<-function(e1,e2) ifelse(e1==0,0, e1/e2)



## For nested sampling: assumes PSUs are not correlated and so a
## pair can't have observations from two PSUs
## $first are marginal PSU weights
## $cond are conditional pairwise weights given that the PSU is selected
## $all are the full pairwise weights
## first and cond are needed for the sandwich estimator, not for the
##   pairwise likelihood itself
pi_from_design<-function(design, ii,jj){

    if (design$pps && !is.null(design$dcheck)){
        ## We have pairwise probabilities already. Or, at least, covariances
        Deltacheck<-design$dcheck[ii,jj]
        indep<-design$prob[ii]*design$prob[jj]
        pi_ij<-(Deltacheck+1)*indep

        last<-ncol(design$allprob)
        n<-design$fpc$sampsize
        N<-design$fpc$popsize

        ## But sandwich standard errors would require fourth-order probabilities
        return(list(full=pi_ij,
                    first=NULL,
                    cond=NULL))
        }
    
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
                return(list(full= n[ii]*(n[jj]-1)%//%( N[ii]*(N[jj]-1)),
                            first=n[ii]/N[ii],
                            cond=rep(1,length(ii))))
            } else {
                ## Hajek high entropy: based on Brewer p153, equation 9.14
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
            else if(is_close(as.vector(design$allprob[[1]]),
                              as.vector(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
                # srs, possibly stratified
                n<-design$fpc$sampsize
                N<-design$fpc$popsize
                return(list(full= (n[ii]/N[ii]),
                            first=n[ii]/N[ii],
                            cond=rep(1,length(ii))))
            } else {
                ## Hajek high entropy: based on Brewer p153, equation 9.14
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
        samestrata<-(design$strata[ii, ]==design$strata[jj, ])
        pstages <-(n[ii,]/N[ii,])*(samestrata*((n[jj,]-1)%//%(N[jj,]-1)) + (1-samestrata)*(n[jj,]/N[jj,]))  ##FIXME divide by  zero when N==1
        return(list(full=apply((n[ii,]/N[ii,])[,-last,drop=FALSE],1,prod)*pstages[,last],
                    first=apply((n[ii,]/N[ii,])[,-last,drop=FALSE],1,prod),
                    cond=pstages[,last]))
    }

    ## Hajek high entropy: Brewer p153
    first<-cpwt<-rep_len(1,length(ii))
    for (i in 1:ncol(design$allprob)){
        pi<-design$allprob[,i]
        denom<-ave(1-pi, design$strata[,i],FUN=sum)
        samestrata<-(design$strata[ii,i]==design$strata[jj,i])
        if (i==ncol(design$allprob))
            cpwt<-cpwt*pi[ii]*pi[jj]*(1- ifelse(samestrata, (1-pi[ii])*(1-pi[jj])/denom, 0))
        else
            first<-first*pi[ii]
    }
    return(list(full=first*cpwt, first= first, cond=cpwt))

}

getpairs<-function(gp, TOOBIG=1000){
    n<-length(gp)
    if (n < TOOBIG){
        ij<-outer(gp,gp,"==")
        ij<-ij & upper.tri(ij)
        return(which(ij,arr.ind=TRUE))
    } 

    ng<-ave(1:n, gp, FUN=length)
    j<-rep(1:n,ng)
    i<-numeric(length(j))
    for (g in unique(gp)){
        this<- which(gp[j]==g)
        i[this] <-which(gp==g)
        }
    data.frame(i=i[i<j],j=j[i<j])
    }

svy2lme_nested<-function(formula,design,sterr=TRUE, return.devfun=FALSE){

    data<-model.frame(design)
    
    ## Use unweighted model to get starting values and set up variables
    m0<-lme4::lmer(formula,data,REML=FALSE)

    ## remove missing from design
    if (!is.null(naa<-attr(m0@frame,"na.action"))){
        design<-design[-naa,]
    }


    ## Extract varables
    y<-m0@resp$y
    X<-m0@pp$X

    ## cluster indicator
    g<-m0@flist[[1]]
    if (sterr){
        ll<-NCOL(design$cluster)
        gpsu<-g[!duplicated(design$cluster[,ll-1])]
        if(any(duplicated(gpsu))){
            stop("model clusters must be nested in design clusters")
        }
    }
    
    ## number of clusters
    n1<-length(unique(g))
    ## number of random effects
    qis<-sapply(m0@cnms,length)
    q<-sum(qis)
    n<-NROW(X)

    ## Z in lme4::lmer has separate columns for each cluster; restructure it
    Z<-matrix(nrow=n, ncol=sum(qis))
    pos<-0
    npos<-0
    for(qi in qis){
        Z[,pos+(1:qi)]<-as.matrix(crossprod(m0@pp$Zt[npos+(1:(n1*qi)),,drop=FALSE], outer(1:(n1*qi),1:qi,function(i,j) ((i-j) %% qi)==0)*1))
        pos<-pos+qi
        npos<-npos+qi*n1
    }
    
    ##Z<-crossprod(m0@pp$Zt, (outer(1:(n1*q),1:q,function(i,j) ((i-j) %% q)==0)*1))


    ## all pairs within same cluster
    ## Conceptually:
    ## ij<-subset(expand.grid(i=1:n,j=1:n), (g[i] == g[j]) & (i<j))
    ## but needs to work when n^2 is too big to construct
    ij<-getpairs(g)
    
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

        ## There could be multiple clusters in the same PSU
        ## Sum the cluster influence functions over PSU, then crossprod
        
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
            
            inffunS<-rowsum(inffun, PSUg,reorder=FALSE)
            stratPSU<-design$strata[,1][ii[!duplicated(design$cluster[,1][ii])]] ##FIXME to allow single-PSU strata?

            one<-rep(1,NROW(inffunS))
            ni<-ave(one,stratPSU,FUN=NROW)
            centering<-apply(inffunS,2,function(x) ave(x, stratPSU, FUN=mean))
            centered<- inffunS-centering
            crossprod(centered*sqrt(ni/ifelse(ni==1,1,(ni-1))))
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
               L=L, method="nested")
    
    ## for resampling
    if(return.devfun) {
        rval$devfun<-devfun
        rval$lower<-m0@lower
        }
    
    class(rval)<-"svy2lme"
    rval
}
        

 
print.svy2lme<-function(x,digits=max(3L, getOption("digits") - 3L),...){
    cat("Linear mixed model fitted by pairwise pseudolikelihood\n")
    if(!is.null(x$call)){
        cat("Call: ")
        cat(paste(deparse(x$call),collapse="\n"))
    }else{
            cat("Formula: ")
            cat(paste(deparse(x$formula),collapse="\n"))
    }
    cat("\nRandom effects:\n")
    theta<-x$opt$par
    s<-sqrt(as.vector(x$s2))
    stdev<- matrix(s*sqrt(diag(x$L)),ncol=1)
    if (!is.null(names(x$zname))){
        vcnames<-paste(names(x$znames), x$znames, sep=":")
    } else {
        vcnames<-x$znames
    }
    rownames(stdev)<-vcnames
    colnames(stdev)<-"Std.Dev."
    print(round(stdev,digits))
    cat("Residual:\t",round(s,digits))
    cat("\n Fixed effects:\n")
    coef<- cbind(beta=x$beta,SE=sqrt(diag(x$Vbeta)),t=x$beta/sqrt(diag(x$Vbeta)))
    coef<-cbind(coef,p=2*pnorm(-abs(coef[,3])))
    colnames(coef)<-c("beta","SE","t","p")
    printCoefmat(coef,digits=digits,P.values=TRUE,has.Pvalue=TRUE, signif.stars=FALSE)
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
        drop(object$beta)
}

vcov.svy2lme<-function(object,...){
    object$Vbeta
}
