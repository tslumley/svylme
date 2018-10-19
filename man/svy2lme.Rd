\name{svy2lme}
\alias{svy2lme}
\alias{coef.svy2lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed models by pairwise likelihood
}
\description{
Fits linear mixed models by profile pairwise composite likelihood. At
the moment, it can only handle two-level models, where the model
clusters are the same as the sampling clusters at the second-last stage
of sampling. 
}
\usage{
svy2lme(formula, design, sterr=TRUE,  return.devfun=FALSE)
\method{coef}{svy2lme}(object,...,random=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    Model formula as in the \code{lme4} package, specifying a single
    grouping variable (cluster) for the random effects
  }
  \item{design}{
    A survey design object produced by \code{survey::svydesign}. The
    pairwise weights will be computed from this design, which must have
    separate probabilities or weights for each stage of sampling.
  }
  \item{sterr}{
    Estimate standard errors for fixed effects? Set to \code{FALSE} for
    greater speed when using resampling to get standard errors.
  }

  \item{return.devfun}{If \code{TRUE}, return the deviance function as a
    component of the object. This will increase the memory use
    substantially, but allows for bootstrapping}
    \item{object}{\code{svy2lme} object}
    \item{random}{if \code{TRUE}, the variance components rather than
    the fixed effects}
}
\details{
  The population pairwise likelihood would be the sum of the
    loglikelihoods for a pair of observations, taken over all pairs of
    observations from the same cluster.  This is estimated by taking a
    weighted sum over pairs in the sample, with the weights being the
    reciprocals of pairwise sampling probabilities. The advantage over
    standard weighted pseudolikelihoods is that there is no
    large-cluster assumption needed and no rescaling of weights. The
    disadvantage is some loss of efficiency (to an extent that isn't
    well understood), and simpler point estimation.
  
    The code uses \code{lme4::lmer} to parse the formula and produce
    starting values, profiles out the fixed effects and residual
    variance, and then uses \code{minqa::bobyqa} to maximise the
    resulting profile deviance.

    As with \code{lme4::lmer}, the default is to estimate the
    correlations of the random effects, since there is typically  no
    reason to assume these are zero. You can force two random effects to
    be independent by entering them in separate terms, eg
    \code{(1|g)+(-1+x|g)} in the model formula asks for a random slope
    and random intercept forced to be zero.
    
    The internal parametrisation of the variance components uses the
    Cholesky decomposition of the relative variance matrix (the variance
    matrix divided by the residual variance), as in
    \code{lme4::lmer}. The component \code{object$s2} contains the
    estimated residual variance and the component \code{object$opt$par}
    contains the elements of the Cholesky factor in column-major order,
    omitting any elements that are forced to be zero by requiring
    indepedent random effects. 
    
    Standard errors of the fixed effects are currently estimated using a
    "with replacement" approximation, valid when the sampling fraction
    is small and superpopulation (model, process) inference is
    intended. Tthe influence functions are added up within
    cluster, centered within strata, the residuals added up within strata, and then
    the crossproduct is taken within each stratum. The stratum variance
    is scaled by \code{ni/(ni-1)} where \code{ni} is the number of PSUs
    in the stratum, and then added up across strata. When the sampling
    and model structure are the same, this is the estimator of Yi et al,
    but it also allows for there to be sampling stages before the stages
    that are modelled.
  
}
\value{
  \code{svy2lme} returns an object with \code{print}, \code{coef}, and
  \code{vcov} methods.

  The \code{coef} method with \code{random=TRUE} returns a two-element
  list: the first element is the estimated residual variance, the second
  is the matrix of estimated variances and covariances of the random effects.
  
}
\references{
  J.N.K. Rao, François Verret and Mike A. Hidiroglou "A weighted composite likelihood approach to inference for two-level models from survey data" Survey Methodology, December 2013  Vol. 39, No. 2, pp. 263-282

  Grace Y. Yi , J. N. K. Rao and Haocheng Li "A WEIGHTED COMPOSITE LIKELIHOOD APPROACH FOR ANALYSIS OF SURVEY DATA UNDER TWO-LEVEL MODELS" Statistica Sinica Statistica Sinica
  Vol. 26, No. 2 (April 2016), pp. 569-587
}
\author{
  Thomas Lumley
}
\examples{

data(api, package="survey")

# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

svy2lme(api00~(1|dnum)+ell+mobility+api99, design=dclus1)
svy2lme(api00~(1|dnum)+ell+mobility+api99, design=dclus2)

lme4::lmer(api00~(1|dnum)+ell+mobility+api99, data=apipop)

## Simulated

set.seed(2000-2-29)

df<-data.frame(x=rnorm(1000*20),g=rep(1:1000,each=20), t=rep(1:20,1000), id=1:20000)
df$u<-with(df, rnorm(1000)[g])

df$y<-with(df, x+u+rnorm(1000,s=2))

## oversample extreme `u` to bias random-intercept variance
pg<-exp(abs(df$u/2)-2.2)[df$t==1]

in1<-rbinom(1000,1,pg)==1
in2<-rep(1:5, length(in1))

sdf<-subset(df, (g \%in\% (1:1000)[in1]) & (t \%in\% in2))

p1<-rep(pg[in1],each=5)
N2<-rep(20,nrow(sdf))

## Population values
lme4::lmer(y~x+(1|g), data=df)

## Naive estimator: higher intercept variance
lme4::lmer(y~x+(1|g), data=sdf)

##pairwise estimator
sdf$w1<-1/p1
sdf$w2<-20/5

design<-survey::svydesign(id=~g+id, data=sdf, weights=~w1+w2)
pair<-svy2lme(y~x+(1|g),design=design)
pair


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression }% use one of  RShowDoc("KEYWORDS")
\keyword{survey }% __ONLY ONE__ keyword per line
