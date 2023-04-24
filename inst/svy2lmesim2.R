
N1=400 
N2=400 
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long
overlap=ceiling(N2*3/4)


cflmer<-function(model){
    a<-VarCorr(model)
    c(fixef(model), as.vector(a[[1]]), attr(a,"sc")^2)
    }
cfsvy<-function(model){
    a<-coef(model, random=TRUE)
    c(coef(model), a$varb,a$s2)
    }



model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}


population<-model_cluster(population,overlap)


rr<-replicate(1000, {
    population$x<- population$long %% 40
    population$z<-rnorm(400*400)
    population$u<-sort(rnorm(400))[population$cluster]
    population$y<- with(population, x+z + u+rnorm(400*400))
    
    
    lmer(y~z+x+(1|cluster), data=population)
    
    
    population$strata<-(population$long-1) %/% 40
    population$uid<-1:nrow(population)
    
    stratsize<- c(20,5,4,3,2,2,3,4,5,20)
    names(stratsize)<-unique(population$strata)
    sstrat<-stratsample(population$strata[!duplicated(population$PSU)], stratsize)
    
    stage1psu<- population$PSU[!duplicated(population$PSU)][sstrat]
    stage1<- subset(population, PSU %in% stage1psu)
    
    
    stratsize2<-rep(c(20,8,20),c(1,66,1))
    names(stratsize2)<-unique(stage1$PSU)
    stage2<-stage1[stratsample(stage1$PSU, stratsize2),]
    
    
    stage2$fpc1<-400/10
    stage2$fpc2<-400
    des<-svydesign(id=~PSU+uid, fpc=~fpc1+fpc2, strata=~strata,data=stage2)
    
    
    c(
        cfsvy(svy2lme(y~x+z+(1|cluster),design=des)),
        cflmer(lmer(y~x+z+(1|cluster),population)),
        cflmer(lmer(y~x+z+(1|cluster),stage2))
    )
})


matrix(apply(rr, 1, median),byrow=TRUE,nrow=3)
matrix(apply(rr, 1, mad),byrow=TRUE,nrow=3)
