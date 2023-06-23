GenDataHMM=function(OrderFeatured=FALSE,K=K,MeanOverExpre=1:K,MeanUnderExpre=-c(1:K),p=p,seed=1,n=n,
                    SigmaOE=rep(1,K),SigmaUE=rep(1,K))
  
  
{
  
  Zk=rep(0,n)
  
  for (i in 1:n){
    for (k in 1:K){
      if ((1+n*(k-1)/K<=i)&&(i<=n*k/K)){
        Zk[i]=k;
      }
    }
  }
  ### Create kappa, feature membership
  KappaKnown=OrderFeatured
  set.seed(1)
  kappajk=matrix(0,p,K)
  if (KappaKnown==TRUE){
    kappajk=matrix(3,p,K);
    kappajk[1:100,]=1;
    kappajk[101:200,]=2;
  } else {
    x1=sample(p,100)
    kappajk[x1,]=1;
    x2=sample(c(1:p)[-x1],100);
    kappajk[x2,]=2;
  }
  Y=matrix(0,n,p)
  set.seed(seed);
  for (j in 1:p){
    for (i in 1:n){
      for (k in 1:K){
        if (Zk[i]==k){
          if ((kappajk [j,k]==1)){
            Y[i,j]=rnorm(1,mean=MeanOverExpre[k],sd =SigmaOE[k])
          } else if ((kappajk [j,k]==2)){
            Y[i,j]=rnorm(1,mean=MeanUnderExpre[k],sd=SigmaUE[k])
          } else {
            Y[i,j]=rnorm(1,mean=0,sd=1)
          }
        }
        
      }
    }
    # K is the number of clusters
  }
  list(Y=Y,kappajk=kappajk,Zk=Zk)
}
  


  
