library(lme4)
library(svylme)

set.seed(2023-6-20)

N1=400 
N2=400 
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long
overlap=ceiling(N2*1/2)


cflmer<-function(model){
    a<-VarCorr(model)
    c(fixef(model), as.vector(unlist(a[1:2])), attr(a,"sc")^2, SE(model))
    }
cfsvy<-function(model){
    a<-coef(model, random=TRUE)
    c(coef(model), diag(a$varb),a$s2,SE(model))
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


f<-function(overlap,REPS=1000){
    
    population<-model_cluster(population,overlap)
    population$x<- population$long %% 40
    population$z<-rnorm(400*400)
    population$u<-sort(rnorm(400))[population$cluster]
    population$y<- with(population, x+z + u+rnorm(400*400))    
    
    population$strata<-(population$long-1) %/% 40
    population$uid<-1:nrow(population)
    
    true<-cflmer(lmer(y~x+z+(1|cluster), population))
    
    rr<-replicate(REPS, {
        
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
            cfsvy(svy2lme(y~x+z+(1|cluster), design=des)),
            ##cflmer(lmer(y~x+z+(1|cluster), population)),
            cflmer(lmer(y~x+z+(1|cluster), stage2))
        )
    })

    list(
        overlap=overlap/N2,
        true=true,
        median=matrix(apply(rr, 1, median),byrow=TRUE,nrow=2),
        mad=matrix(apply(rr, 1, mad),byrow=TRUE,nrow=2)
    )
}


results<- lapply(c(0.1,0.25,0.5,0.75,0.9,1)*N2, f)