\name{milk_subset}
\alias{milk_subset}
\docType{data}
\title{
  Milk production (subset)
}
\description{
A subset of a dataset from the \code{pedigreemm} package, created as an
example for the \code{lme4qtl} package.  The original data had records
of the milk production of 3397 lactations from first through fifty
parity Holsteins. These were 1,359 cows, daughters of 38 sires in 57
herds. The data was downloaded from the USDA internet site. All
lactation records represent cows with at least 100 days in milk, with an
average of 347 days. Milk yield ranged from 4,065 to 19,345 kg estimated
for 305 days, averaging 11,636 kg. There were 1,314, 1,006, 640, 334 and
103 records were from first thorough fifth lactation animals.  The
subset is of cows from 3 sires. 

  
}
\usage{data("milk_subset")}
\format{
  A data frame with 316 observations on the following 13 variables.
  \describe{
    \item{\code{id}}{numeric identifier of cow}
    \item{\code{lact}}{number of lactation for which production is measured}
    \item{\code{herd}}{a factor indicating the herd}
    \item{\code{sire}}{a factor indicating the sire}
    \item{\code{dim}}{number of days in milk for that lactation}
    \item{\code{milk}}{milk production estimated at 305 days}
    \item{\code{fat}}{fat production estimated at 305 days}
    \item{\code{prot}}{protein production estimated at 305 days}
    \item{\code{scs}}{the somatic cell score}
    \item{\code{sdMilk}}{\code{milk} scaled by cow-specific
      standard deviation}
    \item{\code{herd_id}}{a character vector indicating the herd}
    \item{\code{one}}{a vector of 1s for convenience in weighting}
    \item{\code{one2}}{another vector of 1s for convenience in weighting}
  }
}
\details{
This data example gives noticeably different results for full likelihood
and pairwise likelihood because the model is misspecified.  The best
fitting linear model for the large herd 89 is different, and that herd
gets relatively more weight in the pairwise analysis (because it has
more pairs).

}
\source{
  Constructed at \url{https://github.com/variani/lme4qtl/blob/master/vignettes/pedigreemm.Rmd}
}
\references{
   2010. A.I. Vazquez, D.M. Bates, G.J.M. Rosa, D. Gianola and K.A. Weigel.
   Technical Note: An R package for fitting generalized linear mixed models
   in animal breeding. Journal of Animal Science, 88:497-504.
 }
\examples{
data(milk_subset)
herd_des<- svydesign(id = ~herd + id, prob = ~one + one2, data = milk_subset)
lm(sdMilk ~ lact + log(dim),data=milk_subset,subset=herd==89)
lm(sdMilk ~ lact + log(dim),data=milk_subset,subset=herd!=89)
svy2lme(sdMilk ~ lact + log(dim) + (1|herd), design=herd_des,method="nested")
svy2lme(sdMilk ~ lact + log(dim) + (1|herd), design=herd_des,method="general")

## pairwise result is closer to herd 89 than to remainder
lmer(sdMilk ~ lact + log(dim) + (1|herd), data=milk_subset)
svy2relmer(sdMilk ~ lact + log(dim) + (1|id) + (1|herd), design=herd_des,
    relmat = list(id = A_gen))
}
\keyword{datasets}
