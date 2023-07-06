library(pedigreemm)

data(milk)

milk <- within(milk, {
  id <- as.character(id)
  sdMilk <- milk / sd(milk)
})
system.time(
m0<-pedigreemm(sdMilk~lact+log(dim)+(1|id)+(1|herd),data=milk, pedigree=list(id=pedCowsR), REML=FALSE)
)

A_gen <- getA(pedCowsR)
ind <- rownames(A_gen) %in% milk$id
A_gen <- A_gen[ind, ind]


library(lme4qtl)
library(svylme)
library(sampling)

simMilk<-function(theta,model, n){
	Lambda<- getME(model, "Lambda")
    Zt<-getME(model,"Zt")
    Lind<-getME(model, "Lind")
    Lambda@x<- theta[Lind]
    s2<-model@devcomp$cmp["sigmaML"]
    m<-nrow(Zt)
  	u<-matrix(rnorm(m*n,0,1),ncol=n)
	U<-crossprod(Zt,Lambda)%*%u*sqrt(s2)
	Y<-drop(getME(model,"X")%*%model@beta)+U+matrix(rnorm(nrow(U)*n,0,s=sqrt(s2)),ncol=n)
	Y
}


set.seed(2023-6-29)
sim_milk<-simMilk(m0@optinfo$val, m0,2)
milk$simMilk<-sim_milk[,1]

herds<-aggregate(milk$herd,list(milk$herd),length)
herds$p<-herds[,2]*10/sum(herds[,2])

cfsvy<-function(model)  c(coef(model), unlist(coef(model,random=TRUE))[c(2,5,1)])

cflmer<-function(model){
    a<-VarCorr(model)
    c(fixef(model), as.vector(unlist(a[1:2])), attr(a,"sc")^2)
    }

Pi2<-UPtillepi2(herds$p)
dimnames(Pi2)<-list(herds[,1], herds[,1])



one.sim<-function(){
    sampled_herds<-as.logical(UPtille(herds$p))
    submilk<-subset(milk, herd %in% herds[sampled_herds,1])
    submilk$herd<-as.character(submilk$herd)
    
    p<-herds$p[sampled_herds]
    names(p)<-herds[sampled_herds,1]
    submilk$p<-p[submilk$herd]
    
    
    PI2_sub<-Pi2[sampled_herds,sampled_herds][submilk$herd,submilk$herd]
    
    
    
    sub_milk_des<-svydesign(id=~herd,data=submilk, prob=~p,pps=ppsmat(PI2_sub))
    
    
    m1a<-relmatLmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),data=submilk, relmat=list(id=A_gen),REML=FALSE)
    m2a<-svy2relmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),design=sub_milk_des, relmat=list(id=A_gen),return.devfun=TRUE)
    m3a<-svy2relmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),design=sub_milk_des, relmat=list(id=A_gen),all.pairs=TRUE, subtract.margins=TRUE)
    
    

    m1b<-relmatLmer(simMilk~lact+log(dim)+(1|id)+(1|herd),data=submilk, relmat=list(id=A_gen),REML=FALSE)
    m2b<-svy2relmer(simMilk~lact+log(dim)+(1|id)+(1|herd),design=sub_milk_des, relmat=list(id=A_gen),return.devfun=TRUE)
    m3b<-svy2relmer(simMilk~lact+log(dim)+(1|id)+(1|herd),design=sub_milk_des, relmat=list(id=A_gen),all.pairs=TRUE, subtract.margins=TRUE)
    
    rval<-c( cflmer(m1a),
            cfsvy(m2a),
            cfsvy(m3a),
            cflmer(m1b),
            cfsvy(m2b),
            cfsvy(m3b)
            )
}

results<-replicate(1000, tryCatch(one.sim(), error=function(e) rep(NA,36)))
save(results,file="~/milk-sampling.rda")
