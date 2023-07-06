
data(nzmaths)

nzmaths$cSTRATUM<- nzmaths$STRATUM
nzmaths$cSTRATUM[nzmaths$cSTRATUM=="NZL0102"]<-"NZL0202"


des<-svydesign(id=~SCHOOLID+STIDSTD, strata=~cSTRATUM, nest=TRUE,
	weights=~W_FSCHWT+condwt, data=nzmaths)

des<-update(des, centPCGIRLS=PCGIRLS-0.5)
jkdes<-as.svrepdesign(des)

m1<-svy2lme(PV1MATH~ (1+ ST04Q01 |SCHOOLID)+ST04Q01*(centPCGIRLS+SMRATIO)+MATHEFF+OPENPS, design=des, return.devfun=TRUE)
m2<-svy2lme(PV1MATH~ (1+ ST04Q01+MATHEFF+OPENPS |SCHOOLID)+ST04Q01*centPCGIRLS+MATHEFF+OPENPS, design=des, return.devfun=TRUE)

m1var<-boot2lme(m1,jkdes,verbose=TRUE)
m2var<-boot2lme(m2,jkdes,verbose=TRUE)
