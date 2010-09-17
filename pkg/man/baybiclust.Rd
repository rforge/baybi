\name{baybiclust}
\alias{baybiclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Biclustering Function with Importance}
\description{Imports data matrix and transformed parameters and gives bitrees}
\usage{
baybiclust(x, transformed.par, row.labels = rownames(x), col.labels = colnames(x))
}
\arguments{
  \item{x}{data matrix}
  \item{transformed.par}{vector of length 5}
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
library(bclust)
library(labeltodendro)
data(julien)
y=julien
y.select= !(apply(y==0,1,sum)>0)
y=log(y[y.select,])
x=(y-median(as.vector(y)))/sd(as.vector(y))

x=(y-median(as.vector(y)))/sd(as.vector(y))
x.id<-rep(1:nrow(y))
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

transpar=c( -0.69,-100.0000000,-0.69,0.1157195,6.5472181)
bi=baybiclust(julien[1:20,],transformed.par=transpar)

plot(bi$rowtree)
abline(h=bi$rowtree$cut)
plot(bi$rowtree$clust.number,bi$rowtree$logposterior,type="b")


plot(bi$coltree)
abline(h=bi$coltree$cut)
plot(bi$coltree$clust.number,bi$coltree$logposterior,type="b")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
