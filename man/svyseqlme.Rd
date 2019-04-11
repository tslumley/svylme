\name{svyseqlme}
\alias{svyseqlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear mixed models by sequential pseudolikelihood
}
\description{
Sequential pseudolikelihood is my name for the estimation approach introduced by Pfeffermann et al (1998) and expanded in a sequence of papers by Rabe-Hesketh and Skrondal. It is the procedure behind complex-survey estimation in the \code{gllamm} package of Stata code. The construction here is slightly different from that of the previous authors, since it is based on the penalised least squares computations behind the \code{lme4} package, but the results are essentially identical in the case of linear mixed models. 

}
\usage{
svyseqlme(formula, data, REML = FALSE, aweights, yweights, uweights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A model formula, as in \code{lmer}
}
  \item{data}{
A data frame
}
  \item{REML}{
Should (weighted pseudo)\code{REML} estimation be used instead of (weighted pseudo) likelihood?
}
  \item{aweights}{
Precision weights for the observations (ie, inversely proportional to variances)
}
  \item{yweights}{
Sampling weights for the observations: the ultimate weights at the observation level
}
  \item{uweights}{
Sampling weights for the vector of spherical random effects. These are the cluster sampling weights for the clusters that the effects appear in. (Needs documentation as to ordering for multistage sampling!)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Pfeffermann, D., Skinner, C. J., Holmes, D. J., Goldstein, H., & Rasbash, J. (1998). Weighting for unequal selection probabilities in multilevel models. Journal of the Royal Statistical Society: series B (statistical methodology), 60(1), 23-40.

Rabe-Hesketh, S. and Skrondal, A. (2006). Multilevel modeling of complex survey data. Journal of the Royal Statistical Society, Series A 169 (4), 805-827
}



\seealso{
\code{\link{svy2lme}} for pairwise  pseudolikelihood modelleing
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
