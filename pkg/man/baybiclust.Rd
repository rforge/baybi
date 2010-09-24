\name{baybiclust}
\alias{baybiclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Biclustering Function with Importance}
\description{Imports data matrix and transformed parameters and gives bitrees}
\usage{
baybiclust(x, x.id, transformed.par, row.labels = rownames(x), col.labels = colnames(x))
}
\arguments{
  \item{x}{data matrix}
  \item{x.id}{an integer vector declaring replications of a clustering subject}
  \item{transformed.par}{transformed model parameters, a vector of length 5}
  \item{row.labels}{row labels }
  \item{col.labels}{column labels}
}
\details{There is no error control! Segfault may happen easily.}
\value{
  \item{rowtree}{dendrogram tree of rows}
  \item{coltree}{dendrogram tree of columns}
  \item{data}{reordered data \code{x} according to dendrogram orders}
}
\author{Vahid Partovi Nia}
\seealso{loglikelihood in bclust package}
\examples{
library(baybi)
library(bclust)
library(labeltodendro)
data(julien)

#unreplicated example

y=julien
y.select= !(apply(y==0,1,sum)>0)
y=log(y[y.select,])
x=(y-median(as.vector(y)))/sd(as.vector(y))
x.id<-1:nrow(y)
mc.julien<-meancss(x,x.id)
stn.ratio=1
     optimfunc<-function(theta)
     {
     -loglikelihood(x.mean=mc.julien$mean,x.css=mc.julien$css,
     repno=mc.julien$repno,transformed.par=c(theta[1],-100,theta[1]+log(stn.ratio),theta[2:3]),
     var.select=FALSE)#compute - log likelihood
     }
     
theta<-optim(rep(0,3),optimfunc,method="BFGS")$par
transpar<-c(theta[1],-100,theta[1]+log(stn.ratio),theta[2:3])

bi=baybiclust(x=x[1:20,],transformed.par=transpar)

baybiplot(bi,image.width=1,vertdendrogram.size=1.5)

plot(bi$rowtree$clust.number,bi$rowtree$logposterior,type="b")
abline(h=max(bi$rowtree$logposterior))
viplot(bi$rowtree$imp.unordered,xlab=bi$rowtree$labels.unordered,xlab.mar=20)


plot(bi$coltree$clust.number,bi$coltree$logposterior,type="b")
abline(h=max(bi$coltree$logposterior))
viplot(bi$coltree$imp.unordered,xlab=bi$coltree$labels.unordered)

# replicated examples

library(bclust)
data(gaelle)
transpar=c(-1.8390106,-0.9863396,1.6291749,0.0829399,-3.33)
gaelle.id=c(rep(1,3),rep(2:5,each=4))
gaelle.lab<-c("ColWT","d172","d263","isa2","sex4")

bi=baybiclust(x=gaelle[1:19,1:20],x.id=gaelle.id,
row.labels=gaelle.lab,transformed.par=transpar)

baybiplot(bi)

}
