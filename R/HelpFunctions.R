## Hyperparameter setting
#alpha=rep(1,K)
#delt=rep(1,3) ## hyper for transition probbabilities
#### MCMC algorithm
MCMCfunction<-function(Method="HMM",Data=Data,K=K,t=0.01,alpha=0,rho=.1,gam=.1,alphaP=rep(1,K),
                       delt=rep(1,3),sigma20Mu=1000,
                       Nsample=1000, burnin=300,seed=1){
  set.seed(seed)
  ### Initialization
  TransMat=matrix(rep(c(1/3,1/3,1/3),3),3,byrow = T)
  n=nrow(Data)
  p=ncol(Data)
  Zk=rep(0,n)
  dat=Data
  for (i in 1:n){
    Zk[i]=sample(1:K,size=1)
  }
  kappajk=matrix(0,p,K)
  for (k in 1:K){
    for (j in 1:p){
      kappajk [j,k]=sample(1:3,size=1,prob=c(1,1,1))
    }
  }
  #kappajk=Dat$kappajk
  #Zk=Dat$Zk
  Muk=matrix(0,K,3)
  sigma2=array(1,dim = c(p,K,2))
  sigma2j3=rep(1,p);
  Mu=array(0,dim=c(p,K,3))
  sigma2Muk=matrix(1,K,2)
  probZMean=matrix(0,n,K);
  probkappajk=array(0, dim=c(p,K,3));
  TransMatMean=matrix(0,3,3)
  MuMean=matrix(0,K,3)
  SigmaMuMean=matrix(0,K,2);
  logpostsample=NULL;
  loglik=NULL
  for (tt in 1:(Nsample+burnin)){
    if ((tt)%%10==1){
      #print(paste("Number of mcmc iterations is ", tt))
      cat(paste("\rmcmc iteration", tt))
    }
    
 
    
    if (Method!="HMM"){
      Nr=rep(0,3)
      for (s in 1:3){
        Nr[s]=sum(kappajk==s)
      }
      ma=rdirichlet(1,as.vector(delt+Nr))
      TransMat=matrix(c(rep(ma,3)),3,byrow = T)
    }
    else {
      for (r in 1:3){
        Nr=rep(0,3)
        for (s in 1:3){
          for (j in 1:(p-1)){
            Nr[s]=Nr[s]+sum((kappajk[j,]==r) & (kappajk[j+1,]==s))
          }
        }
        TransMat[r,]=rdirichlet(1,as.vector(delt+Nr))
      }
    }
    if (tt>burnin) { TransMatMean=TransMatMean+TransMat/Nsample}
    
    
    ### Sample MU and Sigma
    
    for (j in 1:p){
      for (l in 1:3){
        for (k in 1:K){
          Nkl=0;Skl=0;
          #if (kappajk[j,k]==l){
          Nkl=(kappajk[j,k]==l)*sum(Zk==k)
          Skl=(kappajk[j,k]==l)*sum((Zk==k)*dat[,j]);
          
          if (l<=2){
            Mean=(Muk[k,l]*sigma2[j,k,l]/sigma2Muk[k,l]+Skl)/(sigma2[j,k,l]/sigma2Muk[k,l]+Nkl)
            SD2=1/(Nkl/sigma2[j,k,l]+1/sigma2Muk[k,l])
            # } else if (l==3){
            #Mean=(Muk[k,l]*sigma2j3[j]/sigma23+Skl)/(sigma2j3[j]/sigma23+Nkl)
            #   SD2=1/(Nkl/sigma2j3[j]+1/sigma23)
            #}
            
            Mu[j,k,l]=  rnorm(1, mean =Mean , sd = sqrt(SD2))
          }
          
          #Mu[k,l]=  rtruncnorm(1, a=t, b=Inf, mean =Mean , sd = sqrt(SD2))
          if (l<=2){
            S2kl=0;
            #if (kappajk[j,k]==l){
            S2kl=(kappajk[j,k]==l)*sum((Zk==k)*(dat[,j]-Mu[j,k,l])^2)
            #}
            sigma2[j,k,l]=rinvgamma(1, rho+Nkl/2, gam+S2kl/2)
          }
          
        }
        #for (i in 1:n){
        #N3=N3+sum(kappajk[j,Zk[i]]==3)
        # S32=S32+sum((kappajk[j,Zk[i]]==3)*(dat[i,j]-Mu[j,Zk[i],3])^2)
        #}
      }
      #if (l==3){
      N3=0;S32=0;
      N3=sum(kappajk[j,Zk]==3)
      S32=sum((kappajk[j,Zk]==3)*(dat[,j])^2)
      sigma2j3[j]=rinvgamma(1, rho+N3/2, gam+S32/2)
      #}
      
    }
    
    
    for (l in 1:2){
      for (k in 1:K){
        Nkl=0;Skl=0;
        Nkl=sum(kappajk[,k]==l)
        Skl=sum((kappajk[,k]==l)*Mu[,k,l])
        Mean=(alpha*sigma2Muk[k,l]/sigma20Mu+Skl)/(sigma2Muk[k,l]/sigma20Mu+Nkl)
        SD2=1/(Nkl/sigma2Muk[k,l]+1/sigma20Mu)
        if (l==1){
          Muk[k,l]=  rtruncnorm(1, a=t, b=Inf, mean =Mean , sd = sqrt(SD2))
        } else if (l==2){
          Muk[k,l]=  rtruncnorm(1, a=-Inf, b=-t, mean =Mean , sd = sqrt(SD2))
        }
        S2kl=sum((kappajk[,k]==l)*(Mu[,k,l]-Muk[k,l])^2)
        sigma2Muk[k,l]=rinvgamma(1, rho+Nkl/2, gam+S2kl/2)
        
      }
    }
    #N3=sum(kappajk==3)
    
    if (tt>burnin) {MuMean=MuMean+Muk/Nsample;SigmaMuMean=SigmaMuMean+sigma2Muk/Nsample}
    ### Sample sig2 variance of the noise
    
    
    ### Sample latent kappa and Z
    ## Sample kappa
    for (k in 1:K){
      for (j in 1:p){
        logprobkappa=rep(0,3);
        for (l in 1:3){
          MN=Mu[j,k,l];
          if (l<=2){
            SD=sqrt(sigma2[j,k,l]);
          } else  if (l==3){SD=sqrt(sigma2j3[j]);}
          if (j==1){
            logprobkappa[l]=log(TransMat[l,kappajk[2,k]])+
              sum((Zk==k)*dnorm(dat[,j],mean =MN ,sd=SD,log = T))
            
          } else if (j<=p-1) {
            logprobkappa[l]=log(TransMat[kappajk[j-1,k],l])+log(TransMat[l,kappajk[j+1,k]])+
              sum((Zk==k)*dnorm(dat[,j],mean = MN,sd=SD,log = T))
          } else if (j==p){
            logprobkappa[l]=log(TransMat[kappajk[j-1,k],l])+
              sum((Zk==k)*dnorm(dat[,j],mean = MN,sd=SD,log = T))
          }
        }
        kappajk[j,k]=sample(1:3,size=1,prob=exp(logprobkappa-max(logprobkappa)))
        if (tt>burnin){
          for (l in 1:3){
            probkappajk[j,k,l]= probkappajk[j,k,l]+(kappajk[j,k]==l)
          }
        }
      }
    }
    
    #### Sample cluster probabilities
    Omeg=rdirichlet(1,as.vector(alphaP+sapply(1:K,function(x) sum(Zk==x))))
    ## Sample Z
    
    for (i in 1:n){
      logprobZ=rep(0,K);
      for (k in 1:K){
        logprobZ[k]=sum((kappajk[,k]==1)*dnorm(dat[i,],mean = Mu[,k,1],sd=sqrt(sigma2[,k,1]),log = T))+
          sum((kappajk[,k]==2)*dnorm(dat[i,],mean = Mu[,k,2],sd=sqrt(sigma2[,k,2]),log=T))+
          sum((kappajk[,k]==3)*dnorm(dat[i,],mean = 0,sd=sqrt(sigma2j3),log=T))+log(Omeg[k])
        
      }
      Zk[i]=sample(1:K,size=1,prob=exp(logprobZ-max(logprobZ)))
      if (tt>burnin) {
        for (k in 1:K){
          probZMean[i,k]=probZMean[i,k]+(Zk[i]==k);
        }
      }
    }
    #print(Zk)
    LogP=logpost(dat,Mu,Muk,sigma2Muk,sigma2,sigma2j3,Zk,kappajk,TransMat,Omeg,alphaP,t,alpha,rho,gam,delt,K);
    
    logpostsample[tt]=LogP$logpost
    loglik[tt]=LogP$loglik
  }
  cat("\n")
  ttmax=which.max(logpostsample)
  
  list(probZMean=probZMean/(Nsample),probkappajk=probkappajk/Nsample,TransMatMean=TransMatMean,MuMean=MuMean,
       logpostsample=logpostsample,DIC=-4*mean(loglik)+2*loglik[ttmax])
}

logpost<-function(dat,Mu,Muk,sigma2Muk,sigma2,sigma2j3,Zk,kappajk,TransMat,Omeg,alphaP,t,alpha,rho,gam,delt,K){
  loglik=0;logpost=0;
  n=nrow(dat);p=ncol(dat);
  
  for (i in 1:n){
    for (k in 1:K){
      if (Zk[i]==k){
        loglik=loglik+sum((kappajk[,k]==1)*dnorm(dat[i,],mean = Mu[,k,1],sd=sqrt(sigma2[,k,1]),log = T))+
          sum((kappajk[,k]==2)*dnorm(dat[i,],mean = Mu[,k,2],sd=sqrt(sigma2[,k,2]),log=T))+
          sum((kappajk[,k]==3)*dnorm(dat[i,],mean = Mu[,k,3],sd=sqrt(sigma2j3),log=T))
      }
    }
  }
  #logprior for Zk
  for (k in 1:K){
    logpost=logpost+sum(Zk==k)*log(Omeg[k])
    for (j in 1:p){
      for (l in 1:3){
        if (j==1){
          logpost=logpost+(kappajk[1,k]==l)*log(1/3);
        } else {
          logpost=logpost+(kappajk[j,k]==l)*log(TransMat[kappajk[j-1,k],l])
        }
      }
    }
  }
  logpost=logpost+ddirichlet(Omeg,as.vector(alphaP))
  
  for (r in 1:3){
    logpost=logpost+ddirichlet(TransMat[r,],as.vector(delt))
  }
  
  for (k in 1:K){
    for (l in 1:2){
      if (l==1){
        logpost=logpost+dtruncnorm(Muk[k,l], a=t, b=Inf, mean =alpha , sd = sqrt(sigma2Muk[k,l]))
      } else if (l==2){
        logpost=logpost+dtruncnorm(Muk[k,l], a=-Inf, b=-t, mean =alpha , sd = sqrt(sigma2Muk[k,l]))
      }
      logpost=logpost+sum(dnorm(Mu[,k,l],mean = Muk[k,l],sd=sqrt(sigma2Muk[k,l]),log = T))
      
      logpost=logpost+sum(dinvgamma(sigma2[,k,l], rho, gam))
      logpost=logpost+dinvgamma(sigma2Muk[k,l], rho, gam)
    }}
  logpost=logpost+sum(dinvgamma(sigma2j3, rho, gam))
  list(loglik=loglik,logpost=logpost+loglik)
}


    
  
          
    