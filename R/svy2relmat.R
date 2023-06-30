## FIXME: we're getting underestimation of variance components and non-zero score in models with relmat
## Also (or because) the estimates change slightly if you reorder the data. 

svy2relmer<-function(formula, design, sterr=TRUE, return.devfun=FALSE,
                     relmat=NULL,all.pairs=FALSE, subtract.margins=FALSE){

    data<-model.frame(design)

    formula_copy <-formula
    formula$relmat<-NULL

    ## Use unweighted model to get starting values and set up variables
    m0<-relmatLmer_naive(formula,data,relmat=relmat)
    ## now have to worry about ordering

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
        ## needs to be all correlated (in the model) pairs
        Lambda<- lme4::getME(m0, "Lambda")
        Zt<-lme4::getME(m0,"Zt")
        Xi<-tcrossprod(crossprod(Zt, Lambda)) + Diagonal(n)
        ij<-subset(expand.grid(i=1:n,j=1:n),i!=j)
        ij<-ij[Xi[as.matrix(ij)]!=0,]
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
    ##
    ## having this be a copy of the one in svy2lmeNG looks bad
    ## but it's to allow reference to big objects by lexical scope
    devfun<-function(theta,  subtract_margins=FALSE){
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


        ## assign to enclosing environment for resampling
        Th<-matrix(0,q,q)
        Th[ThInd]<-theta
        L<-tcrossprod(Th)


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

        Nhat<-sum(pwt) ## population number of correlated pairs 

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
            qf<-qf-qf_ind+2*(N-1)*qf_margin
            
            logdet_margin<-sum(log(v_margin)*pw_uni)
            logdet_ind<-sum(log(v11*v22)*pwt)
            logdet<- logdet-logdet_ind+2*(N-1)*logdet_margin

            Nhat<-N*(N-1)  ## population number of pairs
        } 
        s2<<-qf/(2*Nhat)
        
        logdet +2*Nhat*log(qf*2*pi/Nhat)
        
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


        
        Xbeta<-X%*%beta
        r<-y-Xbeta
        r1<-r[ii]
        r2<-r[jj]
        ## all pairs by subtraction
        ## nb: some observations may not be in *any* correlated pairs
        if (subtract_margins){
            v_margin <- D
            pw_uni<-weights(design)
            N<-sum(pw_uni)  ## population number of observations
        }
        
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
    
    }
    if (any(zero<-(theta0==m0@lower))){
        theta0[zero]<-0.5  ## relative variance, so 0.5 should be safe, but should see what lmer does
    }
     
    ## Powell's derivative-free quadratic optimiser
    fit<-minqa::bobyqa(theta0, devfun,
                lower = m0@lower,
                upper = rep(Inf, length(theta)), 
                subtract_margins=all.pairs && subtract.margins)

    ## variance of betas, if wanted
    Vbeta<-if (sterr) Vbeta(fit$par,subtract_margins=all.pairs && subtract.margins) else matrix(NA,q,q)

    ## variance components
    Th<-matrix(0,q,q)
    Th[ThInd]<-fit$par
    L<-tcrossprod(Th)

    ## names: get the relmat names into the output if possible
    znames<-do.call(c,m0@cnms)
    if (any(names(znames) %in% names(sys.call()$relmat))){
        tn<-names(znames)[names(znames) %in% names(sys.call()$relmat)]
        for(tni in tn){
            names(znames)[names(znames) %in% tni]<-deparse(sys.call()$relmat[[tni]])
        }
    }
    
    ## return all the things
    rval<-list(opt=fit,
               s2=s2,
               beta=beta,
               Vbeta=Vbeta,
               formula=formula,
               znames=znames,
               L=L,call=sys.call(),
               all.pairs=all.pairs,
               subtract.margins=subtract.margins)
    
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
  relmat =NULL
)
{
  ## lme4 formula
  control <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    mc <- mcout <- match.call()

  ## lme4 data setup
  lmod <- lme4::lFormula(formula, data, control = control)

    if (is.null(relmat)){
        warning("No relmat terms found")
        } else {
            ##-------------------------------
            ## start of relmatLmer-specific code
            ##-------------------------------
            if (!is.list(relmat)) stop("relmat must be a list")
            if (length(names(relmat)) != length(relmat)) stop("relmat terms must have names")

            relnms <- names(relmat)
            relfac <- relmat
            flist <- lmod$reTrms[["flist"]]   ## list of factors
            fnmns <- names(flist)
            
            ind <- (relnms %in% names(flist))
            if (any(!ind)) warning("some relmat terms are not used")
            
            if(any(ind)) {   
                asgn <- attr(flist, "assign")
                if (any(duplicated(asgn))) stop("a relmat term can have only one random-effect term")
                for(i in seq_along(fnmns)) {
                    fn <- fnmns[i]
                    if (!(fn %in% relnms)) next  ## not relmat
                    
                    ##tn <- which(relmati == names(flist))
                    ##fn <- names(flist)[tn]
                    
                    zn <- lmod$fr[, fn]
                    zn<-as.factor(zn)
                    zn.unique <- levels(zn)
                    
                    if(is.null(rownames(relmat[[fn]]))) stop("relmat matrices must have dimnames")
                    rn <- rownames(relmat[[fn]])
                    
                    if(!all(zn.unique %in% rn)) stop("relmat dimnames do not match factor levels")
                    
                    ## compute a relative factor R: K = R'R
                    ## See lme4qtl:::relfac
                    K <- Matrix::Matrix(relmat[[fn]][zn.unique, zn.unique], sparse = TRUE)
                    R <- Matrix::chol(K)
                    relfac[[fn]] <- R
                    
                    pi <- length(lmod$reTrms$cnms[[i]])
                    Zi_t <- lmod$reTrms$Ztlist[[i]] 
                    Zi_t <- kronecker(R, diag(1, pi)) %*% Zi_t ## t(Z*)
                    
                    ## put the new t(Z*) back into the appropriate slot `Ztlist`
                    lmod$reTrms$Ztlist[[i]] <- Zi_t
                    
                }
            }
            lmod$reTrms[["Zt"]] <- do.call(rbind, lmod$reTrms$Ztlist)
        }
    mcout$formula <- lmod$formula
    lmod$formula <- NULL
    
    devfun <- do.call(mkLmerDevfun, c(lmod,list(start = start)))
    
    opt <-optimizeLmer(devfun,start=start)
    
    rval<-lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
    rval@optinfo$relfac<-list(relfac=relfac)
    rval
}
