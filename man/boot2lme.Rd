\name{boot2lme}
\alias{boot2lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resampling variances for svy2lme
}
\description{
  Computes variance estimates for the weighted-pairwise-likelihood linear mixed models fitted by \code{\link{svy2lme}} using replicate weights.  The replicate weights for a pair are obtained by dividing by the sampling weight and then multiplying by the replicate weight. There will be a warning if the ratio of replicate weight to sampling weight differs for observations in the same pair. 
}
\usage{
boot2lme(model, basewts, replicates, scale, rscales = NULL,
   verbose = FALSE)
\method{vcov}{boot2lme}<-function(object,
   parameter=c("beta", "theta","s2", "relSD" ,"SD","relVar","fullVar"),
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A model returned by \code{svy2lme} with the \code{devfun=TRUE} option
}
  \item{basewts}{
vector of sampling weights
}
  \item{replicates}{
matrix of replicate weights
}
  \item{scale}{
see Details below
}
  \item{rscales}{
see Details below
}
  \item{verbose}{
print progess information
}
\item{object}{returned by \code{boot2lme}}
\item{parameter}{Variance matrix for: regression parameters, relative variance
  parameters on Cholesky square root scale, residual variance, relative
  standard errors of random effects, standard errors of random effects,
  entire relative variance matrix, entire variance matrix}
}
\details{
The variance is estimated from the replicates \code{thetastar} and original point estimate \code{theta} as \code{scale*sum(rscales* (thetastar-theta)^2)}. See the example for how to extract all the necessary information from a replicate-weights survey design object.

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
\donttest{
##takes a few minutes
data(pisa)

dpisa<-survey::svydesign(id=~id_school, weight=~w_fstuwt, data=pisa)

m<-svy2lme(isei~(1+female|id_school)+female+high_school+college+one_for+both_for+test_lang,
	data=pisa, p1=1/pisa$wnrschbw,p2=pisa$wnrschbw/pisa$w_fstuwt, design=dpisa,
	return.devfun=TRUE)

bpisa<-as.svrepdesign(dpisa, type="bootstrap", replicates=100)

b<-boot2lme(m, weights(bpisa,"sampling"), weights(bpisa,"analysis"),
	       bpisa$scale, bpisa$rscales,verbose=TRUE)
str(b)

vcov(b,"beta")
vcov(b,"SD")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
\keyword{survey}% __ONLY ONE__ keyword per line
