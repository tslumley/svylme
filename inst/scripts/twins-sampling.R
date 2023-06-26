data("twinbmi",package="mets")
library(svylme)
library(Matrix)
I_twin<-with(twinbmi, Matrix(outer(1:nrow(twinbmi),1:nrow(twinbmi),function(i,j) (id[i]==id[j]) & (i!=j))))
I_mz<-with(twinbmi, Matrix(outer(1:nrow(twinbmi),1:nrow(twinbmi),function(i,j) (id[i]==id[j]) & (zyg[i]=="MZ") & (i!=j))))

n<-nrow(I_twin)
Phi_env<-I_twin+Diagonal(n)
Phi_add<-I_twin/2+I_mz/2+Diagonal(n)
Phi_dom<-I_twin/4+I_mz*3/4+Diagonal(n)

dimnames(Phi_env)<-list(twinbmi$id,twinbmi$id)
dimnames(Phi_add)<-list(twinbmi$id,twinbmi$id)
dimnames(Phi_dom)<-list(twinbmi$id,twinbmi$id)


twinbmi$id2<-twinbmi$id
twinbmi$id3<-twinbmi$id

## sampling

## whole twins

twinbmi$dbmi<-with(twinbmi, ave(bmi,id, FUN=function(v) if (length(v)>1) abs(diff(v)) else 0))

dup<-duplicated(twinbmi$id)
uid<-twinbmi$id[!dup]
udbmi<-twinbmi$dbmi[!dup]
twinbmi$strata<-cut(twinbmi$dbmi, quantile(udbmi,(1:5)/5), include.lowest=TRUE) 
nsample<-c(50,50,50,400)
names(nsample)<-levels(twinbmi$strata)
insample<- twinbmi$id %in% uid[stratsample(twinbmi$strata[!dup], nsample)]
twinbmi$fpc<-1383
des<-svydesign(id=~id,data=twinbmi[insample,], strata=~strata, fpc=~fpc)



## environment
a<-svy2lme(bmi ~ age+gender+(1|id), design=des)
b<-lme4::lmer(bmi ~ age+gender+(1|id),data=twinbmi[insample,])
##svy2lme(bmi ~ age+gender+(1|id), design=des,all.pairs=TRUE,subtract.margins=TRUE)

## environment plus additive genetic
d<-svy2relmer(bmi ~ age+gender+(1|id)+(1|id2), design=des,relmat=list(id=Phi_env,id2=Phi_add))
e<-lme4qtl::relmatLmer(bmi ~ age+gender+(1|id)+(1|id2), data=twinbmi[insample,],relmat=list(id=Phi_env,id2=Phi_add))

list(a=c(coef(a),coef(a, random=TRUE)),
     a=c(coef(a),coef(a, random=TRUE)),
     d=c(coef(d),coef(d, random=TRUE)))

## individuals

twinbmi$keep<-rbinom(nrow(twinbmi), 1, .5)
twinsubsample <- subset(twinbmi[insample,], keep==1)
twinsubsample$fpc2<-2
iid<-1:nrow(twinsubsample)
des2<-svydesign(id=~id+iid,data=twinsubsample, strata=~strata, fpc=~fpc+fpc2)

## environment
svy2lme(bmi ~ age+gender+(1|id), design=des2)
lme4::lmer(bmi ~ age+gender+(1|id),data=twinsubsample)
svy2lme(bmi ~ age+gender+(1|id), design=des2,all.pairs=TRUE,subtract.margins=TRUE)

## environment plus additive genetic
svy2relmer(bmi ~ age+gender+(1|id)+(1|id2), design=des2,relmat=list(id=Phi_env,id2=Phi_add))
lme4qtl::relmatLmer(bmi ~ age+gender+(1|id)+(1|id2), data=twinsubsample,relmat=list(id=Phi_env,id2=Phi_add))
svy2relmer(bmi ~ age+gender+(1|id)+(1|id2), design=des2,relmat=list(id=Phi_env,id2=Phi_add),all.pairs=TRUE,subtract.margins=TRUE)

