\name{BiclustBHMM}
\alias{BiclustBHMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	An MCMC algorithm to perform a Bayesian model-based biclustering using Hidden Markov Models 
}
\description{
The method accounts for prior dependence between features, and performs a clustering of samples. Within each cluster, features are classified in 3 groups: over-expressed, under-expressed and irrelevant features. Hence, if K is the number of clusters, we have in total 3*K biclusters. The algorithm is implemented on the four methods: i) HMMBi-C: this method determines biclusters by assuming an order structure between features through  a hidden Markov structure on the data, and imposing a constraint on the mean parameters \eqn{\mu_{kl}} (positive, negative and zero depending on the feature state) in order to have a better biological interpretation of obtained clusters; ii)  HMMBi-NoC: this method assumes an order structure between features through  a hidden Markov structure on the data, but does not impose a constraint on the mean parameters \eqn{\mu_{kl}}; iii) NoHMMBi-C: this method does not assume an order structure between features, but does impose a constraint on the mean parameters \eqn{\mu_{kl}}; and iv) NoHMMBi-NoC: this method does not assume an order structure between features, and  does not impose a constraint on the mean parameters \eqn{\mu_{kl}}.  The algorithm computes mainly (i) the posterior probabilities of inclusions for each subject to be in each cluster and ii) For each cluster, the posterior probabilities of inclusions for each feature to be in a group (or state). 
}

\usage{
BiclustBHMM(method="HMMBi-C",Data=Data,K=K,TruncValue=.2,Mu0=0,sigma20Mu=1000,
hyperErrVar=c(1,1),alphaP=rep(1,K),delt=rep(0.5,3),Nsample=3000, burnin=1000,seed=1)
}

\arguments{
	\item{method}{ It's one of the four methods: "HMMBi-C","HMMBi-NoC","NoHMMBi-C" or "NoHMMBi-NoC".  
	}
	\item{Data}{ An expression data matrix of dimension \eqn{n\times p} where \eqn{n} is the number of subjects, and \eqn{p} is the number of features. Features are ordered with respect to some criteria. For instance, we sorted features using similarities computed Gene Ontology.}
	\item{K}{ It's the number of subject clusters.}
	\item{TruncValue}{ It's the minimum value imposed on the mean parameters to identify either under-expressed or over-expressed features.}
	\item{Mu0}{ Prior mean of \eqn{\mu_{kl}}, parameter means of biclusters. }
	\item{sigma20Mu}{ Prior variance of \eqn{\mu_{kl}}, parameter means of biclusters. }
	\item{hyperErrVar }{ Shape and rate parameters of the inverse prior distribution of both the error variance and the variance of \eqn{\mu_{kl}}. They are denoted \eqn{\alpha_0} and \eqn{\beta_0} in the manuscript.}
	\item{alphaP}{ Scale parameters of the dirichlet prior distributions of proba. of inclusion of subjects. They are denoted \eqn{\alpha_1}, ..., \eqn{\alpha_K} in the manuscript.}
	\item{delt}{ Scale parameters of the dirichlet prior distributions of transition prob. of features for each cluster. They are denoted \eqn{\delta_1}, \eqn{\delta_3} and \eqn{\delta_3} in the manuscript.}
\item{Nsample}{ Total number of MCMC draws. It must be larger than burnin.}
\item{burnin}{ 	Number of draws to discard for burn-in.} 
\item{seed}{ Set a seed number to generate distributions in the MCMC algorithm.}
}
\details{
	The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}
\value{
\item{probZMean}{ Posterior probability of cluster membership for each sample/subject}
\item{probkappajk}{ Posterior probability of group membership for each feature within each sample cluster. It's an array of dimension \eqn{p\times K\times 3}. }
\item{TransMatMean}{ Overall (posterior) transition matrix between the three states (or groups) of featurues.}
\item{MuMean}{ Estimated posterior mean of cluster means. The last column is zero, and the first column should be means of over-expressed biclusters, the second cloumn is the mean of the under-expressed bicluster. The dimension is \eqn{K\times 3}..}
\item{logpostsample}{ Provide the log-posterior for each mcmc sample.}
\item{DIC}{ Provide the Deviance information criteria (DIC)}
}

\references{
Thierry Chekouo and Himadri Mukherjee (2023), \emph{A Bayesian Hierarchical Hidden Markov Model for clustering and gene selection: Application to Kidney cancer gene expression data}, \emph{submitted}.
}

\seealso{
\code{\link{GenDataHMM}}
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
Result=BiclustBHMM(method="NoHMMBi-NoC",Data=dat,K=K,TruncValue=.2,Mu0=0,
sigma20Mu=1000,hyperErrVar=c(1,1), alphaP=rep(1,K),delt=rep(0.5,3),Nsample=300, 
burnin=100,seed=1)
str(Result)
F11=F1(TrueZ=Dat$Zk,TrueKappajk=Dat$kappajk,
probZmean=Result$probZMean,probkappajk=Result$probkappajk,EstK=K)

F11

\dontrun{
# this is a long running example on the real data application
data(SimilarityBP)
HC=hclust(as.dist(1-SimilarityBP), method = "average", members = NULL)
OrderEntrez=HC$order
data(mRNAExpression)
dat=mRNAExpression[OrderEntrez,]
K=2
Result=BiclustBHMM(method="HMMBi-C",Data=dat,K=K,TruncValue=.2,Mu0=0,
sigma20Mu=1000,hyperErrVar=c(1,1), alphaP=rep(1,K),delt=rep(0.5,3),Nsample=300, 
burnin=100,seed=1)
}

}

