\name{GenDataHMM}
\alias{GenDataHMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of simulated data as explained in the reference manuscript.	      
}
\description{ 
	This function generates data described in the manuscript.
}

\usage{
GenDataHMM(OrderFeatured=FALSE,K=K,MeanOverExpre=1:K,MeanUnderExpre=-c(1:K),p=p,seed=1,n=n,
                    SigmaOE=rep(1,K),SigmaUE=rep(1,K))
  
}

\arguments{
	        \item{OrderFeatured}{ If TRUE, then consecutive features belong to the same group. Otherwise, features in any group are chosen randomly without a predefined order.}
		\item{K}{ Number of sample clusters}
		\item{MeanOverExpre}{ (positive) Means of over-expressed biclusters.}
			\item{MeanUnderExpre}{ (negative) Means of under-expressed biclusters..}
		\item{p}{ Number of features.}
		\item{seed}{ Seed to generate random numbers.}
				\item{n}{ Number of subjects.}
		\item{SigmaOE}{ Variance to generate over-expressed features.}
				\item{SigmaUE}{ Variance to generate under-expressed feature.}
}

\details{
The function will generate data as explained in the manuscript. To see the results, use the ``$" operator. 
} 
\value{
\item{Y}{ A expression matrix of dimension \eqn{n\times p}}
\item{kappajk}{ A matrix of feature membership of dimension \eqn{p\times K}. Values of the matrix are 1, 2 or 3 that correspond respectively to over-expressed, under-expressed and irrelevant features. }
\item{TransMatMean}{ Overall (posterior) transition matrix between the three states (or groups) of featurues.}
\item{Zk}{ A vector of length \eqn{n} of subject membership. Values are \eqn{1,..., K}.}
}

\references{
Thierry Chekouo and Himadri Mukherjee (2023), \emph{A Bayesian Hierarchical Hidden Markov Model for clustering and gene selection: Application to Kidney cancer gene expression data}, \emph{submitted}.
}
\seealso{
\code{\link{BiclustBHMM}}
}

\examples{
library(BiclustBHMM);
s2=1 # Variance of the error
K=2 # number of clusters
n=100 # number of subjects
p=1000 # number of features
Dat=GenDataHMM(OrderFeatured=TRUE,K=K,p=p,seed=1,n=n,SigmaOE=rep(s2,K),
SigmaUE=rep(s2,K),MeanOverExpre=c(2:(K+1))/2,MeanUnderExpre=-c(2:(K+1))/2)

str(Dat)
}
