\name{baybiplot}
\alias{baybiplot}
\title{plots a baybimp object}
\description{ this is a generic function with a lot of flexible options}
\usage{
baybiplot(x, xlab = x$coltree$labels.unordered, 
ylab = x$rowtree$labels.unordered, xlab.cex = 1, ylab.cex = 1, 
horizdendrogram.lwd = 1, horizdendrogram.size = 2, vertdendrogram.lwd = 1, 
vertdendrogram.size = 2, xlab.mar = 6, ylab.mar = 3, 
image.col = cm.colors(12), horizbar.col = rev(c(heat.colors(5)[-4], "white")), 
vertbar.col = rev(c(heat.colors(5)[-4], "white")), horizbar.size = 0.25, 
vertbar.size = 0.25, horizteeth.size = 0.25, vertteeth.size = 0.25, 
image.width = 5, image.height = 3)
}
\arguments{
  \item{x}{a baybimp object}
  \item{xlab}{row labels appearing in x axis}
  \item{ylab}{column labels appearing in y axis}
  \item{xlab.cex}{cex for xlab}
  \item{ylab.cex}{cex for ylab}
  \item{horizdendrogram.lwd}{lwd for the horizontal dendrogram}
  \item{horizdendrogram.size}{size of the horizontal dendrogram}
  \item{vertdendrogram.lwd}{lwd for the vertical dendrogram}
  \item{vertdendrogram.size}{size of the vertical dendrogram}
  \item{xlab.mar}{margins for xlab}
  \item{ylab.mar}{margins for ylab}
  \item{image.col}{colors for the image plot}
  \item{horizbar.col}{colors for the horizontal bar of importances}
  \item{vertbar.col}{colors for the vertical bar of importances}
  \item{horizbar.size}{size for the horizontal bar of importances}
  \item{vertbar.size}{size for the vertical bar of importances}
  \item{horizteeth.size}{size for the horizontal teeth plot}
  \item{vertteeth.size}{size for the vertical teeth plot}
  \item{image.width}{width of the image plot}
  \item{image.height}{height of the image plot}
}
\details{see also ditplot and dptplot in bclust package}
\author{Vahid Partovi Nia}
\seealso{baybimp}
\examples{
library(bclust)
data(gaelle)
transpar=c(-1.8390106,-0.9863396,1.6291749,0.0829399,-3.33)
bi=baybiclust(x=gaelle[1:10,1:20],transformed.par=transpar)
baybiplot(bi)
  }

