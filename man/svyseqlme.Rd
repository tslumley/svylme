\name{svyseqlme}
\alias{svyseqlme}
\alias{coef.svyseqlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed models by sequential pseudolikelihood
}
\description{
Sequential pseudolikelihood is my name for the estimation approach introduced by Pfeffermann et al (1998) and expanded in a sequence of papers by Rabe-Hesketh and Skrondal. It is the procedure behind complex-survey estimation in the \code{gllamm} package of Stata code. The construction here is slightly different from that of the previous authors, since it is based on the penalised least squares computations behind the \code{lme4} package, but the results are the same for linear mixed models in designs without strata. 

}
\usage{
svyseqlme(formula, design, REML = FALSE,
scale = c("sample_size", "effective_sample_size", "gk", "raw"))
\method{coef}{svy2lme}(object,...,random=FALSE)
    }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A model formula, as in \code{lmer}
}
  \item{design}{
A survey design object produced by \code{survey::svydesign}, containing weights for each stage of sampling.
}
  \item{REML}{
Should (weighted pseudo)\code{REML} estimation be used instead of (weighted pseudo) likelihood?
}
\item{scale}{How to scale weights}
 \item{\dots}{for method compatibility}
    \item{random}{if \code{TRUE}, the variance components rather than
    the fixed effects}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Object of class \code{svyseqlme}
}
\references{
Pfeffermann, D., Skinner, C. J., Holmes, D. J., Goldstein, H., & Rasbash, J. (1998). Weighting for unequal selection probabilities in multilevel models. Journal of the Royal Statistical Society: series B (statistical methodology), 60(1), 23-40.

Rabe-Hesketh, S. and Skrondal, A. (2006). Multilevel modeling of complex survey data. Journal of the Royal Statistical Society, Series A 169 (4), 805-827
}



\seealso{
\code{\link{svy2lme}} for pairwise  pseudolikelihood modelleing
}
\examples{

## simulated example
## 
set.seed(2000-2-29)

df<-data.frame(x=rnorm(1000*20),g=rep(1:1000,each=20), t=rep(1:20,1000), id=1:20000)
df$u<-with(df, rnorm(1000)[g])

df$y<-with(df, x+u+rnorm(1000,s=2))

## oversample extreme `u` to bias random-intercept variance
pg<-exp(abs(df$u/2)-2.2)[df$t==1]

in1<-rbinom(1000,1,pg)==1
in2<-rep(1:5, length(in1))

sdf<-subset(df, (g %in% (1:1000)[in1]) & (t %in% in2))

p1<-rep(pg[in1],each=5)
N2<-rep(20,nrow(sdf))

sdf$w1<-1/p1
sdf$w2<-20/5

design<-survey::svydesign(id=~g+id, data=sdf, weights=~w1+w2)
svyseqlme(y~x+(1|g),design=design)

## from Stata
##            Variance  (sqrt() this for comparison)
## (Intercept)   1.034899
## Residual:	 3.728113
##                beta     SE       t     p
## (Intercept) -0.0756 0.0988 -0.7655 0.444
## x            1.0470 0.0795 13.1722 0.000


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
