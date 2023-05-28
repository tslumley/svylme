\name{svy2relmer}
\alias{svy2relmer}
\encoding{utf8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed models with correlated random effects
}
\description{
Fits linear mixed models by profile pairwise composite
likelihood. Allows for correlated random effects, eg, for genetic
relatedness (QTL) models
}
\usage{
svy2relmer(formula, design, sterr=TRUE,  return.devfun=FALSE, relmat=NULL,
 all.pairs=FALSE, subtract.margins=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    Model formula as in the \code{lme4} package, or with terms like
    \code{(1|id)} for correlated random effects together with the
    \code{relmat} argument.
    
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
    substantially, but allows for bootstrapping.}
  \item{relmat){
     Specifies a list of relatedness matrices that corresponds to one or
more random-effect groupings (eg \code{(1|id)} in the \code{formula}
    together with \code{relmat=list(id=Phi)} implies a covariance matrix
    of \code{Phi} for the random effects before scaling. See Details and
    the vignettes. 
    }
    \item{all.pairs}Use all pairs rather than just correlated pairs?}
    \item{subtract.margins}{If \code{TRUE} and \code{all.pairs=TRUE},
    compute with all pairs by the faster algorithm involving subtraction
    from the marginal likelihood}
}
\details{

    This function is very similar to \code{\link{svy2lme}} and only the
    differences are described here.

    Formula parsing and starting values use code based on the
    \code{lme4qtl} package.

    In \code{svy2lme} and \code{lme4::lmer}, the model is based on
    independent standard Normal random effects that are transformed to
    give random coefficients that might be correlated within observation
    but are either identical or independent between observations. In
    this function, the basic random effects in a term are multiplied by a square
    root of the \code{relmat} matrix for that term, giving basic random
    effects whose covariance between observations proportional to the
    \code{relmat} matrix. For example, in a quantitative trait locus
    model in genetics, the matrix would be a genetic relatedness matrix.


    The \code{return.devfun=TRUE} option is useful if you want to
    examine objects that aren't returned as part of the output. For
    example, \code{get("ij", environment(object$devfun))} is the set of
    pairs used in computation. 
  
}
\value{
  \code{svy2relmer} returns an object with \code{print}, \code{coef}, and
  \code{vcov} methods.

  
}
\references{
Ziyatdinov, A., Vázquez-Santiago, M., Brunel, H. et al. lme4qtl: linear mixed models with flexible covariance structure for genetic studies of related individuals. BMC Bioinformatics 19, 68 (2018). https://doi.org/10.1186/s12859-018-2057-x
}
\author{
  Thomas Lumley
}
\examples{
data(milk_subset)
herd_des<- svydesign(id = ~herd + id, prob = ~one + one2, data = milk_subset)

svy2lme(sdMilk ~ lact + log(dim) + (1|herd), design=herd_des, method="general")

svy2relmer(sdMilk ~ lact + log(dim) + (1|id) + (1|herd), design=herd_des,
    relmat = list(id = A_gen))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression }% use one of  RShowDoc("KEYWORDS")
\keyword{survey }% __ONLY ONE__ keyword per line
