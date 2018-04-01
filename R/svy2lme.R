
svy2lme<-function(formula,data, p1,p2,N2=NULL){
    
    m0<-lme4::lmer(formula,data,REML=FALSE)
    y<-m0@resp$y
    X<-m0@pp$X
    g<-m0@flist[[1]]
    n1<-length(unique(g))
    q<-nrow(m0@pp$Zt)/n1
    Z<-crossprod(m0@pp$Zt, (outer(1:(n1*q),1:q,function(i,j) ((i-j) %% q)==0)*1))

    matrix(diag(q), ncol = q, nrow = n1*q)

    n<-NROW(X)
    ij<-subset(expand.grid(i=1:n,j=1:n), (g[i]==g[j]) &  (i !=j))
    ii<-ij[,1]
    jj<-ij[,2]
    p<-NCOL(X)
    

    theta<-theta0<-c(2*log(m0@devcomp$cmp["sigmaML"]), m0@theta)
    beta<-beta0<-fixef(m0)

    if (is.null(N2)){
        pwt<- (1/p1[ii])*(1/p2[ii])*(1/p2[jj])  ## with replacement at stage 2
    } else {
        n2<-ave(as.numeric(g), g, FUN=length)
        pwt<-(1/p1[ii])*N2[ii]*(N2[jj]-1)/(n2[ii]*(n2[ii]-1))  ## SRS without replacement at stage 2
    }

    m<-nrow(ij)/n1
    L<-matrix(0,q,q)

    devfun<-function(theta){
        s2<-exp(theta[1])
        Th<-matrix(0,q,q)
        Th[lower.tri(Th,diag=TRUE)]<-theta[-1]
        L<<-tcrossprod(Th)

        v11<-(rowSums(Z[ii,,drop=FALSE]*( Z[ii,,drop=FALSE]%*%L))+1)*s2
        v12<-rowSums(Z[ii,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))*s2
        v22<-(rowSums(Z[jj,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))+1)*s2
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
        
        xtwy<-crossprod(Xii,pwt*inv11*y[ii])+
            crossprod(Xjj,pwt*inv22*y[jj])+
            crossprod(Xii,pwt*inv12*y[jj])+
            crossprod(Xjj,pwt*inv12*y[ii])
        
        beta<<-solve(xtwx,xtwy)
        Xbeta<-X%*%beta

        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]

        qf<-crossprod(r1,pwt*inv11*r1)+
            crossprod(r2,pwt*inv22*r2)+
            crossprod(r1,pwt*inv12*r2)+
            crossprod(r2,pwt*inv12*r1)
        
        (sum(log(det)*pwt) + qf)
        
    }

    Vbeta<-function(theta){
        s2<-exp(theta[1])
        Th<-matrix(0,q,q)
        Th[lower.tri(Th,diag=TRUE)]<-theta[-1]
        L<<-tcrossprod(Th)
        
        v11<-(rowSums(Z[ii,,drop=FALSE]*( Z[ii,,drop=FALSE]%*%L))+1)*s2
        v12<-rowSums(Z[ii,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))*s2
        v22<-(rowSums(Z[jj,,drop=FALSE]*(Z[jj,,drop=FALSE]%*%L))+1)*s2
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

        ##FIXME
        
        
        xwr<-Xii*pwt*(inv11*r1)+
            Xjj*pwt*(inv22*r2)+
            Xii*pwt*(inv12*r2)+
            Xjj*pwt*(inv12*r1)

        J<-crossprod(rowsum(xwr,g[ii]))
        G<-solve(xtwx)
        G%*%J%*%G
        }
    
    fit<-bobyqa(theta0, devfun,
                lower = c(-Inf,m0@lower),
                upper = rep(Inf, length(theta0)))

    
    
    rval<-list(opt=fit, beta=beta,Vbeta=Vbeta(fit$par), formula=formula,znames=m0@cnms[[1]],L=L)
    class(rval)<-"svy2lme"
    rval
    
}
        

 
print.svy2lme<-function(x,digits=max(3L, getOption("digits") - 3L),...){
    cat("Linear mixed model fitted by pairwise likelihood\n")
    cat("Formula: ")
    cat(deparse(x$formula),"\n")
    cat("Random effects:\n")
    theta<-x$opt$par[-1]
    s<-x$opt$par[1]
    stdev<- matrix(s*sqrt(diag(x$L)),ncol=1)
    rownames(stdev)<-x$znames
    colnames(stdev)<-"Std.Dev."
    print(stdev)
    cat("\n Fixed effects")
    coef<- cbind(beta=x$beta,SE=sqrt(diag(x$Vbeta)),t=x$beta/sqrt(diag(x$Vbeta)))
    coef<-cbind(coef,p=2*pnorm(-abs(coef[,"t"])))
    print(coef)
    cat("\n")
    }
