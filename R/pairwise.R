
with_rep<-function(p1,p2){
    rval <- list(p1=p1, p2=p2)
    class(rval)<-"ppair"
    rval
    }

srswr<-function(p1, n2){
    rval <- list(p1=p1, p2=p2)
    class(rval)<-"ppair"
    rval
    }

overton<-function(p,n){
    (n-1)*outer(p,p,"*")/(n-outer(p,p,"+")/2)
}

HR<-function(p,psquared){
    (n-1)*outer(p,p,"*")/(n - outer(p,p,,"+") +psquared/n)
}

HR1<-function(p){
    psquaredhat <- sum(p)
    HR(p,psquaredhat)
    }

highentropy<-function(p, psquared=sum(p),n=length(p)){
    outer(p,p)*(1-outer(1-p,1-p)/(n-psquared))
        
    }
