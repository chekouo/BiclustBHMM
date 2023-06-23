F1<-function(TrueZ, TrueKappajk,probZmean,probkappajk,EstK){
  n=length(TrueZ);
  p=nrow(TrueKappajk);K=ncol(TrueKappajk)
  EstKappajk=matrix(0,p,EstK)
  EstZ=apply(probZmean,1, which.max)
  #EstZ=apply(1:n, function(x) which.max(probZmean[x,]))
  Nbr=sapply(1:EstK, function(x) sum(EstZ==x));
  Ks=setdiff(1:EstK,which(Nbr==0))
  
  NKr=list()
  for (k in 1:EstK){
    EstKappajk[,k]=sapply(1:p, function(x) which.max(probkappajk[x,k,]))
    NKr[[k]] =setdiff(1:3,which(sapply(1:3, function(x) sum(EstKappajk[,k]==x))==0));
  }
  F1=array(0, dim=c(length(Ks),3,K,3))
  for (k1 in 1:length(Ks)){
    for (l1 in 1:length(NKr[[Ks[k1]]])){
      for (k in 1:K){
        for (l in 1:3){
          rc=sum((TrueKappajk[,k]==l) & (EstKappajk[,Ks[k1]]==NKr[[Ks[k1]]][l1]))*sum((TrueZ==k)&(EstZ==Ks[k1]))
          NB=sum(TrueKappajk[,k]==l)*sum(TrueZ==k)
          NAA=sum(EstKappajk[,Ks[k1]]==NKr[[Ks[k1]]][l1])*sum(EstZ==Ks[k1])
          recall=rc/NB
          precision=rc/NAA
          #  if (NAA==0){F1[k1,l1,k,l]=0; } else {
          F1[k1,NKr[[Ks[k1]]][l1],k,l]=2/(1/recall+1/precision)
          # }
        }
      }
      #  maxF1=max(F11)
      #   F1final=F1final+maxF1
    }
  }
  
  0.5*(mean(apply(F1,c(1,2),function(x) max(x,na.rm = T)))+mean(apply(F1,c(3,4),function(x) max(x,na.rm = T))))
}

