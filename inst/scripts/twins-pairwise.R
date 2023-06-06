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

des<-svydesign(id=~id,data=twinbmi)

## environment
svy2lme(bmi ~ age+gender+(1|id), design=des)
svy2relmer(bmi ~ age+gender+(1|id), design=des,relmat=list(id=Phi_env))
lme4::lmer(bmi ~ age+gender+(1|id),data=twinbmi)
lme4qtl::relmatLmer(bmi ~ age+gender+(1|id),data=twinbmi,relmat=list(id=Phi_env))
svy2lme(bmi ~ age+gender+(1|id), design=des,all.pairs=TRUE,subtract.margins=TRUE)

## environment plus additive genetic
svy2relmer(bmi ~ age+gender+(1|id)+(1|id2), design=des,relmat=list(id=Phi_env,id2=Phi_add))
lme4qtl::relmatLmer(bmi ~ age+gender+(1|id)+(1|id2), data=twinbmi,relmat=list(id=Phi_env,id2=Phi_add))


## environment plus additive and dominant genetic
svy2relmer(bmi ~ age+gender+(1|id)+(1|id2)+(1|id3), design=des,relmat=list(id=Phi_env,id2=Phi_add,id3=Phi_dom))
lme4qtl::relmatLmer(bmi ~ age+gender+(1|id)+(1|id2)+(1|id3), data=twinbmi,relmat=list(id=Phi_env,id2=Phi_add,id3=Phi_dom))
