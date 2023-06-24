

getallpairs<-function(gps, TOOBIG=1000){
    n<-length(gps[1])
    if (n < TOOBIG){
        ijall<-Reduce("|", lapply(gps, function(gp) outer(gp,gp,"==")), FALSE)
        diag(ijall)<-FALSE
        return(which(ijall, arr.ind=TRUE)) ## return(which(ijall & upper.tri(ijall), arr.ind=TRUE))
    } 

    alli<-list(); allj<-list()
    for (gp in gps){
        ng<-ave(1:n, gp, FUN=length)
        j<-rep(1:n,ng)
        i<-numeric(length(j))
        for (g in unique(gp)){
            this<- which(gp[j]==g)
            i[this] <-which(gp==g)
        }
        alli[[gp]]<-i
        allj[[gp]]<-j
    }
    i<-do.call(c,alli)
    j<-do.call(c,allj)
    rval<-data.frame(i=i[i<j],j=j[i<j])
    unique(rval)
}


svy2lme<-function(formula, design, sterr=TRUE, return.devfun=FALSE, method=c("general","nested"), all.pairs=FALSE, subtract.margins=FALSE){

    method<-match.arg(method)
    if(method=="nested"){
        if(all.pairs) stop("all.pairs=TRUE not allowed for method='nested'")
        return(svy2lme_nested(formula,design, sterr=sterr, return.devfun=return.devfun))
    }
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

    ## cluster indicators
    gs<-m0@flist
    
    ## number of clusters 
    n1s<-sapply(gs, function(gi) length(unique(gi)))
    
    ## number of random effects
    qis<-sapply(m0@cnms,length)
    q<-sum(qis)
    n<-NROW(X)
    
    Z<-t(m0@pp$Zt)
    
    ## need PSUs as well as clusters now
    psu<-design$cluster[[1]]
    
    if (all.pairs && !subtract.margins){
        ## unavoidably going to be big
        ij<-subset(expand.grid(i=1:n,j=1:n),i!=j)
    } else{
        ## all pairs within same cluster
        ## Conceptually, the union of 
        ## ij<-subset(expand.grid(i=1:n,j=1:n), (g[i] == g[j]) & (i<j))
        ## but needs to work when n^2 is too big to construct
        ij<-getallpairs(gs)
    }
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
    allpwts<-svylme:::all_pi_from_design(design,ii,jj)
    pwt<-1/allpwts$full
  
    ## variance matrix of random effects
    qi<-sapply(m0@cnms,length)
    L<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i)))) 
    ###(need indicator for where thetas go in the matrix)
    ThInd<-which((L==1) & lower.tri(L,diag=TRUE))
    Lambda<- lme4::getME(m0, "Lambda")
    Zt<-lme4::getME(m0,"Zt")
    
    ## profile pairwise deviance
    ## a whole heap of stuff is being passed by lexical scope
    devfun<-function(theta, pwt_new=NULL, subtract_margins=FALSE){
        if (!is.null(pwt_new)) pwt<-pwt_new  ##resampling
        
        ## variance parameters: Cholesky square root of variance matrix
        Lind<-lme4::getME(m0, "Lind")
        Lambda@x<- theta[Lind]
        ## Full (sparse) vcov(Y)
        Xi<-tcrossprod(crossprod(Zt, Lambda)) + Diagonal(n)
        D<-diag(Xi)
        
        ## v11 is a vector of (1,1) entries of the matrix var(Y)
        ## for each pair, similarly for the others
        v11<-D[ii]
        v22<-D[jj]
        v12<-Xi[cbind(ii,jj)]
        

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

        ## all pairs by subtraction
        ## nb: some observations may not be in *any* correlated pairs
        if (subtract_margins){
            v_margin <- D
            pw_uni<-weights(design)
            xtwx_margin<-crossprod(X,pw_uni*X/v_margin)
            xtwy_margin<-crossprod(X,pw_uni*y/v_margin)
            xtwx_ind<- crossprod(Xii,pwt*Xii/v11) + crossprod(Xjj,pwt*Xjj/v22)
            xtwy_ind<-crossprod(Xii,pwt*y[ii]/v11) + crossprod(Xjj,pwt*y[jj]/v22)     
            N<-sum(pw_uni)  ## population number of observations
            xtwx<-xtwx-xtwx_ind+2*(N-1)*xtwx_margin
            xtwy<-xtwy-xtwy_ind+2*(N-1)*xtwy_margin
        }

        ## betahat at the given variance parameter values
        beta<<-solve(xtwx,xtwy)
        Xbeta<-X%*%beta

        ## two residuals per pair
        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]

        Nhat<-sum(pwt)*2 ## population number of *correlated* pairs 

        ## -2 times Gaussian log profile pairwise likelihood
        qf<-crossprod(r1,pwt*inv11*r1)+
            crossprod(r2,pwt*inv22*r2)+
            crossprod(r1,pwt*inv12*r2)+
            crossprod(r2,pwt*inv12*r1)

        logdet<-sum(log(det)*pwt)
        
        ## all pairs by subtraction
        if (subtract_margins){
            qf_margin<-crossprod(r,pw_uni*r/v_margin)
            qf_ind<-crossprod(r1,pwt*r1/v11)+crossprod(r2,pwt*r2/v22)
            qf<-qf-qf_ind+(N-1)*qf_margin
            
            logdet_margin<-sum(log(v_margin)*pw_uni)
            logdet_ind<-sum(log(v11*v22)*pwt)
            logdet<- logdet-logdet_ind+(N-1)*logdet_margin

            Nhat<-N*(N-1)  ## population number of pairs
        } 
        s2<<-qf/Nhat
        
        logdet + Nhat*log(qf*2*pi/Nhat)
        
    }

    ## Standard errors of regression parameters
    ##
    ## If beta = (X^TWX)^{-1}(XTWY)
    ## the middle of the sandwich is the sum over design-correlated pairs
    ## of X^TW(Y-mu)^T(Y-mu)WX
    ##
    ## off-diag W is just off-diag Xi[ij]^{-1}/pi_{ij}, ie, inv12/pi_ij
    ## diag W is sum of diag Xi[ij]^{-1}/pi_{ij} for all pairs with i in them
    ## ie, sum_j(inv11/pi_ij) but being careful about indices
    ##
    ## The nested version was simpler because pairs were always in the same PSU
    
    Vbeta<-function(theta, subtract_margins=FALSE){
        ## setup exactly as in devfun
        ## variance parameters: Cholesky square root of variance matrix
        Lind<-lme4::getME(m0, "Lind")
        Lambda@x<- theta[Lind]
        ## Full (sparse) vcov(Y)
        Xi<-tcrossprod(crossprod(Zt, Lambda)) + Diagonal(n)
        D<-diag(Xi)
        
        ## v11 is a vector of (1,1) entries of the matrix var(Y)
        ## for each pair, similarly for the others
        v11<-D[ii]
        v22<-D[jj]
        v12<-Xi[cbind(ii,jj)]
        
        det<-v11*v22-v12*v12
        inv11<- v22/det
        inv22<- v11/det
        inv12<- -v12/det
        
        Xii<-X[ii,,drop=FALSE]
        Xjj<-X[jj,,drop=FALSE]

        ## sensitivity matrix ('bread')
        xtwx<- crossprod(Xii,pwt*inv11*Xii)+
            crossprod(Xjj,pwt*inv22*Xjj)+
            crossprod(Xii,pwt*inv12*Xjj)+
            crossprod(Xjj,pwt*inv12*Xii)

        if (subtract_margins){
            v_margin <- D
            pw_uni<-weights(design)
            xtwx_margin<-crossprod(X,pw_uni*X/v_margin)
            xtwx_ind<- crossprod(Xii,pwt*Xii/v11) + crossprod(Xjj,pwt*Xjj/v22)
            N<-sum(pw_uni)  ## population number of observations
             xtwx<-xtwx-xtwx_ind+2*(N-1)*xtwx_margin
        }
        
        Xbeta<-X%*%beta
        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]

        ## try making W explicitly
        W<-Matrix(0, n,n)
        W[cbind(ii,jj)]<-inv12*pwt
        idx<-which((1:n) %in% ii)
        W[cbind(idx,idx)]<-rowsum(inv11*pwt,ii,reorder=TRUE)
        if (subtract_margins){
            n_uncorr<-rep(n-1,n)
            n_uncorr[idx]<-n_uncorr[idx]-rowsum(rep(1,length(jj)),ii,reorder=TRUE)
            W[cbind(1:n,1:n)]<-W[cbind(1:n,1:n)]+pw_uni*(1/v_margin)*n_uncorr
        }
        xtwx<-crossprod(X, W%*%X)
        xwr<-X*(W%*%r)
        Delta<-survey:::Dcheck_multi(design$cluster, design$strata, design$allprob)
        xtwxinv<-solve(xtwx)
        V<-xtwxinv%*%crossprod(xwr, Delta%*%xwr)%*%xtwxinv
        return(V)
        
        ## ##  variability matrix ('cheese')
        ## ## score for betas 
        ## xwr<-Xii*pwt*(inv11*r1)+
        ##     Xjj*pwt*(inv22*r2)+
        ##     Xii*pwt*(inv12*r2)+
        ##     Xjj*pwt*(inv12*r1)

        ## if (subtract_margins){
        ##    stop("not done")
        ## }
        ## ## The grouping variables here are PSUs (not model clusters)
        
        ## if (is.null(design)){
        ##   stop("standard errors need a design argument")
        ## }
        ## inf<-rowsum( xwr%*%solve(xtwx), ii, reorder=FALSE)
        ## ##inf<- xwr%*%solve(xtwx)
        ## designi<-design[(1:n) %in% ii,]
        ## Dcheck<-survey:::Dcheck_multi(designi$cluster, designi$strata, designi$allprob)
        ## crossprod(inf, Dcheck %*% inf)
        
        ## ## stratPSU<-design$strata[,1][ii[!duplicated(psu[ii])]] ##FIXME to allow single-PSU strata?
        
        ## ## one<-rep(1,NROW(inffun))
        ## ## ni<-ave(one,stratPSU,FUN=NROW)
        ## ## centering<-apply(inffun,2,function(x) ave(x, stratPSU, FUN=mean))
        ## ## centered<- inffun-centering
        ## ## V <- crossprod(centered*sqrt(ni/ifelse(ni==1,1,(ni-1))))
        ## ## V
    }
    
    
    
    ## Powell's derivative-free quadratic optimiser
    fit<-minqa::bobyqa(theta0, devfun,
                lower = m0@lower,
                upper = rep(Inf, length(theta)), 
                subtract_margins=all.pairs && subtract.margins)

    ## variance of betas, if wanted
    Vb<-if (sterr ) Vbeta(fit$par,subtract_margins=all.pairs && subtract.margins) else matrix(NA,q,q)

    ## variance components
    Th<-matrix(0,q,q)
    Th[ThInd]<-fit$par
    L<-tcrossprod(Th)
    ## return all the things
    rval<-list(opt=fit,
               s2=s2,
               beta=beta,
               Vbeta=Vb,
               formula=formula,
               znames=do.call(c,m0@cnms),
               L=L, all.pairs=all.pairs,
               subtract.margins=subtract.margins, method="general")
    
    ## for resampling
    if(return.devfun) {
        rval$devfun<-devfun
        rval$lower<-m0@lower
        }
    
    class(rval)<-"svy2lme"
    rval
}


## pairwise probabilities: does *not* assume nesting
## probably does assume no more than two-stage sampling FIXME
##
## we only use $full, not the other components.
##
all_pi_from_design<-function(design, ii,jj){

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
            ## ok, element sampling, can't be same PSU
            if(is.null(design$fpc$popsize)) #with replacement
                return(list(full=design$prob[ii]*design$prob[jj],
                            first=design$prob[ii],
                            cond=rep(1,length(ii))))
            else if(is_close(as.vector(design$allprob),
                             as.vector(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
                ## srs, possibly stratified
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
            ## possibly ok, sampling of whole PSUs
            warning("assuming no subsampling within PSUs because multi-stage weights were not given")
            
            samePSU<-design$cluster[ii,1]==design$cluster[jj,1]

            if(is.null(design$fpc$popsize)){ #with replacement
                 return(list(full=ifelse(samePSU, design$prob[ii], design$prob[ii]*design$prob[jj]),
                             first=design$prob[ii],
                             cond=rep(1, length(ii))))
            } else if(is_close(as.vector(design$allprob[[1]]),
                              as.vector(design$fpc$sampsize/design$fpc$popsize),tolerance=1e-4)){
                # srs, possibly stratified
                n<-design$fpc$sampsize
                N<-design$fpc$popsize
                return(list(full= ifelse(samePSU, (n[ii]/N[ii]),(n[ii]/N[ii])*(n[jj]/N[jj])),
                            first=n[ii]/N[ii],
                            cond=rep(1,length(ii))))
            } else {
                ## Hajek high entropy: based on Brewer p153, equation 9.14
                pi<-design$allprob
                denom<-ave(1-pi, design$strata,FUN=sum)
                samestrata<-(design$strata[ii,1]==design$strata[jj,1])
                return(list(full=ifelse(samePSU, pi[ii,1], pi[ii,1]*pi[jj,1]*(1- ifelse(samestrata, (1-pi[ii,1])*(1-pi[jj,1])/denom, 0))),
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
    samePSU<-design$cluster[ii,1]==design$cluster[jj,1]

    if(is.null(design$fpc$popsize)){ #with replacement
        last<-ncol(design$allprob)
        return(list(full=ifelse(samePSU, design$prob[ii]*design$allprob[jj,last],design$prob[ii]*design$prob[jj]),
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
        return(list(full=ifelse(samePSU, apply((n[ii,]/N[ii,])[,-last,drop=FALSE],1,prod)*pstages[,last],design$prob[ii]*design$prob[jj]),
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
    return(list(full=ifelse(samePSU, first*cpwt,design$prob[ii]*design$prob[jj]), first= first, cond=cpwt))

}
