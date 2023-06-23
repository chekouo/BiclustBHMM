BiclustBHMM<-function(method="HMMBi-C",Data=Data,K=K,TruncValue=.2,Mu0=0,sigma20Mu=1000,hyperErrVar=c(1,1),alphaP=rep(1,K),delt=rep(0.5,3),
                      Nsample=3000, burnin=1000,seed=1){
# We have in total 4 methods as described in the main manuscript
  #1. HMMBi-C 
 # 2. HMMBi-NoC 
 # 3. NoHMMBi-C 
  #4. NoHMMBi-NoC
  
library(truncnorm) ## generate truncated normal distributions
library(mc2d) ## generate dirichlet
library(invgamma)
  
  if (!method%in%c("HMMBi-C","HMMBi-NoC","NoHMMBi-C","NoHMMBi-NoC")){
    stop("The method should be one of these: HMMBi-C,HMMBi-NoC,NoHMMBi-C,NoHMMBi-NoC")
  }
  rho=hyperErrVar[1]
  gam=hyperErrVar[2]
  alpha=Mu0
if (method=="HMMBi-C"){
  Res=MCMCfunction(Method = "HMM",Data=Data,K=K,t=TruncValue,alpha=alpha,rho=rho,gam=gam,alphaP=alphaP,delt=delt,sigma20Mu=sigma20Mu,
                   Nsample=Nsample, burnin=burnin,seed=1)
}
if (method=="HMMBi-NoC"){
  Res=MCMCfunction(Method = "HMM",Data=Data,K=K,t=-Inf,alpha=alpha,rho=rho,gam=gam,alphaP=alphaP,delt=delt,sigma20Mu=sigma20Mu,
                   Nsample=Nsample, burnin=burnin,seed=1) }
if (method=="NoHMMBi-C"){
  Res=MCMCfunction(Method = "NoHMM",Data=Data,K=K,t=TruncValue,alpha=alpha,rho=rho,gam=gam,alphaP=alphaP,delt=delt,sigma20Mu=sigma20Mu,
                   Nsample=Nsample, burnin=burnin,seed=1) }
if (method=="NoHMMBi-NoC"){
  Res=MCMCfunction(Method = "NoHMM",Data=Data,K=K,t=-Inf,alpha=alpha,rho=rho,gam=gam,alphaP=alphaP,delt=delt,sigma20Mu=sigma20Mu,
                   Nsample=Nsample, burnin=burnin,seed=1) }

return(Res)

}