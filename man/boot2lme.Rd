\name{boot2lme}
\alias{boot2lme}
\alias{vcov.boot2lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resampling variances for svy2lme
}
\description{
  Computes variance estimates for the weighted-pairwise-likelihood
  linear mixed models fitted by \code{\link{svy2lme}} using replicate
  weights.  The replicate weights for a pair are obtained by dividing by
  the sampling weight and then multiplying by the replicate
  weight. There will be a warning if the ratio of replicate weight to
  sampling weight differs for observations in the same pair.  
}
\usage{
boot2lme(model, rdesign,  verbose = FALSE)
\method{vcov}{boot2lme}(object,
   parameter=c("beta", "theta","s2", "relSD" ,"SD","relVar","fullVar"),
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A model returned by \code{svy2lme} with the \code{devfun=TRUE} option
}
  \item{rdesign}{
  replicate-weights design corresponding to the design used to fit the model, see example
}
 
  \item{verbose}{
print progess information?
}
\item{object}{returned by \code{boot2lme}}
\item{\dots}{for method compatibility}
\item{parameter}{Variance matrix for: regression parameters, relative variance
  parameters on Cholesky square root scale, residual variance, relative
  standard errors of random effects, standard errors of random effects,
  entire relative variance matrix, entire variance matrix}
}
\details{
The variance is estimated from the replicates \code{thetastar} and original point estimate \code{theta} as \code{scale*sum(rscales* (thetastar-theta)^2)}. 
}
\value{
 For \code{boot2lme}, an object of class \code{boot2lme} with components
  \item{theta}{replicates of variance parameters (on Cholesky square
    root scale)}
  \item{beta}{replicates of regression parameters}
  \item{D}{replicates of relative variance matrix}
  \item{scale,rscales}{from the input}
  \item{formula}{model formula from the input}

  For the \code{vcov} method, a variance matrix.
}



\examples{

data(api, package="survey")

# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

m0<-svy2lme(api00~(1|dnum)+ell+mobility, design=dclus2,return.devfun=TRUE)
jkdes<-as.svrepdesign(dclus2, type="mrb")
jkvar<-boot2lme(m0,jkdes)

SE(jkvar, "beta")
SE(jkvar, "SD")
SE(jkvar,"s2")


m1<-svy2lme(api00~(1|dnum)+ell+mobility,
design=dclus2,return.devfun=TRUE, all.pairs=TRUE, subtract.margins=TRUE)
jk1var<-boot2lme(m1,jkdes)

SE(jk1var, "beta")
SE(jk1var, "SD")


\donttest{
##takes a few minutes
data(pisa)

pisa$w_condstuwt<-with(pisa, w_fstuwt/wnrschbw)
pisa$id_student<-1:nrow(pisa)

dpisa<-survey::svydesign(id=~id_school+id_student, weight=~wnrschbw+w_condstuwt, data=pisa)

m<-svy2lme(isei~(1+female|id_school)+female+high_school+college+one_for+both_for+test_lang,
	design=dpisa, return.devfun=TRUE,method="nested")

bpisa<-as.svrepdesign(dpisa, type="bootstrap", replicates=100)

b<-boot2lme(m, bpisa, verbose=TRUE)
str(b)

vcov(b,"beta")
vcov(b,"s2")

## SE() inherits the parameter= argument
SE(b,"beta")
SE(b,"theta")
SE(b,"SD")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
\keyword{survey}% __ONLY ONE__ keyword per line
