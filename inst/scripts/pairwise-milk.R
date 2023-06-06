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

set.seed(2023-6-1)
sim_milk<-simulate(m0)
milk$simMilk<-sim_milk$sim_1
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