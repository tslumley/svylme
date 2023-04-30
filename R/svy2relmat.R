svy2relmer<-function(formula, design, sterr=TRUE, return.devfun=FALSE, relmat=NULL){

    data<-model.frame(design)

    formula_copy <-formula
    formula$relmat<-NULL

    ## Use unweighted model to get starting values and set up variables
    m0<-relmatLmer_naive(formula,data,relmat=relmat)

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
    

    ## all pairs within same cluster
    ## needs to be all correlated (in the model) pairs
    Lambda<- lme4::getME(m0, "Lambda")
    Zt<-lme4::getME(m0,"Zt")
    Xi<-tcrossprod(crossprod(Zt, Lambda)) + Diagonal(n)
    ij<-subset(expand.grid(i=1:n,j=1:n),i<j)
    ij<-ij[Xi[as.matrix(ij)]!=0,]
    
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
    if (sterr){
        if (is.null(allpwts$cond))
            stop("Can't get sandwich standard errors for this design")
        pwt2<-1/allpwts$cond
        p1<-allpwts$first
    }
    
    ## variance matrix of random effects
    qi<-sapply(m0@cnms,length)
    L<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i))))  ##FIXME: no, it's a lot more complicated
    ###(need indicator for where thetas go in the matrix)
    ThInd<-which((L==1) & lower.tri(L,diag=TRUE))
    Lambda<- lme4::getME(m0, "Lambda")
    Zt<-lme4::getME(m0,"Zt")
    
    ## profile pairwise deviance
    devfun<-function(theta,pwt){
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

       
        ## The grouping variables here are PSUs, not clusters
        pw1<-1/p1
        
        if (is.null(design)){
          stop("standard errors need a design argument")
        }
        inffun<-rowsum( (xwr*pw1)%*%solve(xtwx), psu[ii], reorder=FALSE)
        
        stratPSU<-design$strata[,1][ii[!duplicated(psu[ii])]] ##FIXME to allow single-PSU strata?
        
        one<-rep(1,NROW(inffun))
        ni<-ave(one,stratPSU,FUN=NROW)
        centering<-apply(inffun,2,function(x) ave(x, stratPSU, FUN=mean))
        centered<- inffun-centering
        V <- crossprod(centered*sqrt(ni/ifelse(ni==1,1,(ni-1))))
        V
    }
    
    
    
    ## Powell's derivative-free quadratic optimiser
    fit<-minqa::bobyqa(theta0, devfun,
                lower = m0@lower,
                upper = rep(Inf, length(theta)), pwt=pwts)

    ## variance of betas, if wanted
    Vbeta<-if (sterr) Vbeta(fit$par,pwts) else matrix(NA,q,q)

    ## variance components
    Th<-matrix(0,q,q)
    Th[ThInd]<-fit$par
    L<-tcrossprod(Th)
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
    
    class(rval)<-c("svy2lme","svy2relmer")
    rval


}


## From lme4qtl (github.com/variani/lme4qtl), GPL3
relmatLmer_naive <- function(formula, data = NULL, 
  start = NULL,
  relmat = list()
)
{
  # formula
  control <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    mc <- mcout <- match.call()

  # lmod
  lmod <- lme4::lFormula(formula, data, control = control)
  
  #-------------------------------
  # start of relmatLmer-specific code
  #-------------------------------
  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  relnms <- names(relmat)
  relfac <- relmat
  flist <- lmod$reTrms[["flist"]]   ## list of factors

  ind <- (relnms %in% names(flist))
  if(any(ind)) {
    ## random-effects design matrix components
    ##Ztlist <- lmod$reTrms[["Ztlist"]]
    
    asgn <- attr(flist, "assign")
    for(i in seq_along(relnms[ind])) {
      
      relmati <- relnms[ind][i]
      if(!(relmati %in% names(flist))) {
        stop("a relationship matrix must be (", relmati, ")",
          " associated with only one random effects term (", paste(names(flist), collapse = ", "), ")")
      }
      tn <- which(relmati == names(flist))
      fn <- names(flist)[tn]
      
      zn <- lmod$fr[, fn]
      zn<-as.factor(zn)
      zn.unique <- levels(zn)

      stopifnot(!is.null(rownames(relmat[[fn]])))
      rn <- rownames(relmat[[fn]])

      stopifnot(all(zn.unique %in% rn))
      
      # compute a relative factor R: K = R'R
      # See lme4qtl:::relfac
      K <- Matrix::Matrix(relmat[[fn]][zn.unique, zn.unique], sparse = TRUE)
      R <- Matrix::chol(K)
      relfac[[fn]] <- R

        pi <- length(lmod$reTrms$cnms[[i]])
        Zi_t <- lmod$reTrms$Ztlist[[i]] 
      Zi_t <- kronecker(R, diag(1, pi)) %*% Zi_t # t(Z*)

      # put the new t(Z*) back into the appropriate slot `Ztlist`
      lmod$reTrms$Ztlist[[i]] <- Zi_t

    }
  }
  lmod$reTrms[["Zt"]] <- do.call(rbind, lmod$reTrms$Ztlist)
    ##-------------------------------
    ## end of relmatLmer-specific code
    ##-------------------------------
    
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  
  devfun <- do.call(mkLmerDevfun, c(lmod,list(start = start)))#, verbose = verbose, control = control)))
  
    opt <-optimizeLmer(devfun,start=start)
        
    rval<-lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
    rval@optinfo$relfac<-list(relfac=relfac)
    rval
}
