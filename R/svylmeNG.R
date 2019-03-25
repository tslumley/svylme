

getallpairs<-function(gps, TOOBIG=1000){
    n<-length(gps[1])
    if (n < TOOBIG){
        ijall<-Reduce("|", lapply(gps, function(gp) outer(gp,gp,"==")), FALSE)
        return(which(ijall & upper.tri(ijall), arr.ind=TRUE))
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


svy2lmeNG<-function(formula,design,sterr=TRUE, return.devfun=FALSE){

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
    
    ## Z in lme4::lmer has separate columns for each cluster; restructure it
    Z<-matrix(nrow=n, ncol=sum(qis))
    pos<-0
    npos<-0
    for(i in 1:length(qis)){
        qi<-qis[i]
        n1i<-n1s[i]
        Z[,pos+(1:qi)]<-as.matrix(crossprod(m0@pp$Zt[npos+(1:(n1i*qi)),,drop=FALSE], outer(1:(n1i*qi),1:qi,function(i,j) ((i-j) %% qi)==0)*1))
        pos<-pos+qi
        npos<-npos+qi*n1i
    }
    


    ## all pairs within same cluster
    ## Conceptually, the union of 
    ## ij<-subset(expand.grid(i=1:n,j=1:n), (g[i] == g[j]) & (i<j))
    ## but needs to work when n^2 is too big to construct
    ij<-getallpairs(gs)
    
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
    allpwts<-svylme:::pi_from_design(design,ii,jj)
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
        ##FIXME
        p1g<-p1[!duplicated(g[ii])]

        if (is.null(design)){
          stop("standard errors need a design argument")
        } 
        inffun<-(1/p1g)*rowsum(xwr,g[ii],reorder=FALSE)%*%solve(xtwx)
        PSUg<-design$cluster[,1][ii[!duplicated(g[ii])]]
        
        inffun<-rowsum(inffun, PSUg,reorder=FALSE)
        stratPSU<-design$strata[,1][ii[!duplicated(design$cluster[,1][ii])]] ##FIXME to allow single-PSU strata?
        
        one<-rep(1,NROW(inffun))
        ni<-ave(one,stratPSU,FUN=NROW)
        centering<-apply(inffun,2,function(x) ave(x, stratPSU, FUN=mean))
        centered<- inffun-centering
        phase2 <- crossprod(centered*sqrt(ni/ifelse(ni==1,1,(ni-1))))
        phase1 <- solve(xtwx)*s2
        phase1+phase2
    }
    
    
    
    ## Powell's derivative-free quadratic optimiser
    fit<-minqa::bobyqa(theta0, devfun,
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
