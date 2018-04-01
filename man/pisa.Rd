\name{pisa}
\alias{pisa}
\docType{data}
\title{
Data from the PISA international school survey
}
\description{
Data from the PISA survey of schools, obtained from Stata, who obtained it from Rabe-Hesketh \& Skrondal. 
}
\usage{data("pisa")}
\format{
  A data frame with 2069 observations on the following 11 variables.
  \describe{
    \item{\code{female}}{1 for female}
    \item{\code{isei}}{socioeconomic index}
    \item{\code{w_fstuwt}}{student sampling weight (total)}
    \item{\code{wnrschbw}}{school sampling weight}
    \item{\code{high_school}}{1 if highest level of parents' education is high school}
    \item{\code{college}}{1 if highest level of parents' education is college/uni}
    \item{\code{one_for}}{1 if one parent is foreign-born}
    \item{\code{both_for}}{1 if both parents are foreign-born}
    \item{\code{test_lang}}{1 if the test language is spoken at home}
    \item{\code{pass_read}}{1 if the student passed a reading proficiency test}
    \item{\code{id_school}}{school (sampling unit) identifier}
  }
}

\source{
Data downloaded from \url{http://www.stata-press.com/data/r15/pisa2000.dta}
}
\references{
Rabe-Hesketh, S., and A. Skrondal.  2006.  Multilevel modelling of complex survey data.Journal of the Royal StatisticalSociety,  Series  A169:  805–827
}
\examples{
data(pisa)

## This model doesn't make a lot of sense, but it's the one in the
## Stata documentation because the outcome variable is numeric.
svy2lme(isei~(1|id_school)+female+high_school+college+one_for+both_for+test_lang,
	data=pisa, p1=1/pisa$wnrschbw,p2=pisa$wnrschbw/pisa$w_fstuwt)

}
\keyword{datasets}
