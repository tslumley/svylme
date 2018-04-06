
svy2lme<-function(formula,data, p1,p2,N2=NULL,sterr=TRUE){

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
    if (is.null(N2)){
        ## using probabilities
        pwt2 <- (1/p2[ii])*(1/p2[jj])
        pwt<- (1/p1[ii])*pwt2  ## with replacement at stage 2
    } else {
        ## using cluster size
        n2<-ave(as.numeric(g), g, FUN=length)
        pwt2<-N2[ii]*(N2[jj]-1)/(n2[ii]*(n2[ii]-1))
        pwt<-(1/p1[ii])*pwt2  ## SRS without replacement at stage 2
    }

    ## variance matrix of random effects
    qi<-sapply(m0@cnms,length)
    L<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i))))
    ###(need indicator for where thetas go in the matrix)
    ThInd<-which((L==1) & lower.tri(L,diag=TRUE))
    
    ## profile pairwise deviance
    devfun<-function(theta){
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

        s2<<-qf/sum(pwt)/2
        
        ##sum(log(det)*pwt) + qf 
        sum(log(det)*pwt) + 2*sum(pwt)*log(qf*2*pi/(2*sum(pwt)))
        
    }

    ## Standard errors of regression parameters
    Vbeta<-function(theta){
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
        p1g<-p1[ii][!duplicated(g[ii])]

        ## sandwich estimator
        J<-crossprod((1/p1g)*rowsum(xwr,g[ii],reorder=FALSE)*(n1/(n1-1)))
        G<-solve(xtwx)
        G%*%J%*%G
        }


    ## Powell's derivative-free quadratic optimiser
    fit<-bobyqa(theta0, devfun,
                lower = m0@lower,
                upper = rep(Inf, length(theta)))

    ## variance of betas, if wanted
    Vbeta<-if (sterr) Vbeta(fit$par) else matrix(NA,q,q)
    
    ## return all the things
    rval<-list(opt=fit,
               s2=s2,
               beta=beta,
               Vbeta=Vbeta,
               formula=formula,
               znames=do.call(c,m0@cnms),
               L=L)
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
        s2<-object$s2
        dimnames(L)<-list(object$znames,object$znames)
        list(s2=s2, varb=L*s2)
    } else 
        object$beta
}

vcov.svy2lme<-function(object,...){
    object$Vbeta
}
