\name{nzmaths}
\alias{nzmaths}
\docType{data}
\title{
Maths Performance Data from the PISA 2012 survey in New Zealand
}
\description{
Data on maths performance, gender, some problem-solving variables and some school resource variables.
}
\usage{data("nzmaths")}
\format{
  A data frame with 4291 observations on the following 26 variables.
  \describe{
    \item{\code{SCHOOLID}}{School ID}
    \item{\code{CNT}}{Country id: a factor with levels \code{New Zealand}}
    \item{\code{STRATUM}}{a factor with levels \code{NZL0101} \code{NZL0102} \code{NZL0202} \code{NZL0203}}
    \item{\code{OECD}}{Is the country in the OECD?}
    \item{\code{STIDSTD}}{Student ID}
    \item{\code{ST04Q01}}{Gender: a factor with levels \code{Female} \code{Male}}
    \item{\code{ST14Q02}}{Mother has university qualifications \code{No} \code{Yes}}
    \item{\code{ST18Q02}}{Father has university qualifications \code{No} \code{Yes}}
    \item{\code{MATHEFF}}{Mathematics Self-Efficacy: numeric vector}
    \item{\code{OPENPS}}{Mathematics Self-Efficacy:  numeric vector}
    \item{\code{PV1MATH},\code{PV2MATH},\code{PV3MATH},\code{PV4MATH},\code{PV5MATH} }{'Plausible values' (multiple imputations) for maths performance}
    \item{\code{W_FSTUWT}}{Design weight for student}
    \item{\code{SC35Q02}}{Proportion of maths teachers with professional development in maths in past year}
    \item{\code{PCGIRLS}}{Proportion of girls at the school}
    \item{\code{PROPMA5A}}{Proportion of maths teachers with ISCED 5A (math major)}
    \item{\code{ABGMATH}}{Does the school group maths students: a factor with levels \code{No ability grouping between any classes} \code{One of these forms of ability grouping between classes for s} \code{One of these forms of ability grouping for all classes}}
    \item{\code{SMRATIO}}{Number of students per maths teacher}
    \item{\code{W_FSCHWT}}{Design weight for school}
    \item{\code{condwt}}{Design weight for student given school}
      }
}

\source{
A subset extracted from the \code{PISA2012lite} R package, \url{https://github.com/pbiecek/PISA2012lite}
}
\references{
OECD (2013) PISA 2012 Assessment and Analytical Framework: Mathematics, Reading, Science, Problem Solving and Financial Literacy. OECD Publishing.
}
\examples{
data(nzmaths)
des<-svydesign(id=~SCHOOLID+STIDSTD, strata=~STRATUM, nest=TRUE,
	weights=~W_FSCHWT+condwt, data=nzmaths)

oo<-options(survey.lonely.psu="remove")

m1<-svy2lme(PV1MATH~ (1+ ST04Q01 |SCHOOLID)+ST04Q01*(PCGIRLS+SMRATIO)+MATHEFF+OPENPS, design=des)
m1

m3<-svyseqlme(PV1MATH~ (1+ ST04Q01 |SCHOOLID)+ST04Q01*(PCGIRLS+SMRATIO)+MATHEFF+OPENPS, 
	design=des,scale="sample_size")
m3

options(oo)
}
\keyword{datasets}
