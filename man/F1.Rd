\name{F1}
\alias{F1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute criteria to evaluate biclustering performamce based on the F1-measure
}
\description{
	This function computes  F1-measure (F1) as defined in the manuscript. 
}

\usage{
F1(TrueZ, TrueKappajk,probZmean,probkappajk,EstK)
}

\arguments{
  \item{TrueZ}{ True membership values of subjects/samples of dimension \eqn{n}, sample size. Possible values are \eqn{1,...,K}.}
\item{TrueKappajk}{ Binary vector with 1 as positive and 0 as negative.}
\item{probZmean}{ (Posterior or Estimated) probability of cluster membership for each sample/subject.}
\item{probkappajk}{(Posterior or Estimated) probability of group membership for each feature within each sample cluster. It's an array of dimension \eqn{p\times K\times 3}.}
\item{EstK}{(Estimated) number of cluster }
}

\details{
The function returns F1 as compute in the manuscript
}
\references{
Thierry Chekouo and Himadri Mukherjee (2023), \emph{A Bayesian Hierarchical Hidden Markov Model for clustering and gene selection: Application to Kidney cancer gene expression data}, \emph{submitted}.
}

\seealso{
\code{\link{GenDataHMM}}
\code{\link{BiclustBHMM}}
}

\examples{
### We run Setting 1 in the manuscript
library(BiclustBHMM);
s2=1 # Variance of the error
K=2 # number of clusters
n=100 # number of subjects
p=1000 # number of features
Dat=GenDataHMM(OrderFeatured=TRUE,K=K,p=p,seed=1,n=n,SigmaOE=rep(s2,K),
SigmaUE=rep(s2,K),MeanOverExpre=c(2:(K+1))/2,MeanUnderExpre=-c(2:(K+1))/2)

dat=Dat$Y

Result=BiclustBHMM(method="HMMBi-C",Data=dat,K=K,TruncValue=.2,Mu0=0,
sigma20Mu=1000,hyperErrVar=c(1,1), alphaP=rep(1,K),delt=rep(0.5,3),Nsample=300, 
burnin=100,seed=1)
str(Result)
F11=F1(TrueZ=Dat$Zk,TrueKappajk=Dat$kappajk,
              probZmean=Result$probZMean,probkappajk=Result$probkappajk,EstK=K)
F11
}
