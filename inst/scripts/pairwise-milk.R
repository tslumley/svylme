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
system.time(
m1<-relmatLmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),data=milk, relmat=list(id=A_gen))
)



library(svylme)
milk_des<-svydesign(id=~1,data=milk)
system.time(
m2<-svy2relmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),design=milk_des, relmat=list(id=A_gen))
)

system.time(
m3<-svy2relmer(sdMilk~lact+log(dim)+(1|id)+(1|herd),design=milk_des, relmat=list(id=A_gen),all.pairs=TRUE, subtract.margins=TRUE)
)


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


set.seed(2023-6-7)
sim_milk<-simMilk(m0@optinfo$val, m0,2)
milk$simMilk<-sim_milk[,1]
sim_milk_des<-svydesign(id=~1,data=milk)

m1a<-relmatLmer(simMilk~lact+log(dim)+(1|id)+(1|herd),data=milk, relmat=list(id=A_gen),REML=FALSE)
m2a<-svy2relmer(simMilk~lact+log(dim)+(1|id)+(1|herd),design=sim_milk_des, relmat=list(id=A_gen),return.devfun=TRUE)
m3a<-svy2relmer(simMilk~lact+log(dim)+(1|id)+(1|herd),design=sim_milk_des, relmat=list(id=A_gen),all.pairs=TRUE, subtract.margins=TRUE)




m0
m1
m2
m3

m1a
m2a
m3a