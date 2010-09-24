baybiclust=function(x,x.id=1:nrow(x),transformed.par,row.labels=rownames(x),col.labels=colnames(x))
{
nocrossing.order<-function(mergemat)
{
	myorder<-matrix(NA,nrow(mergemat),nrow(mergemat)+1)
			na.add <- function(vec)
	{
		return(c(vec,rep(NA,nrow(mergemat)+1-length(vec))))
	}

	myorder[1,]<-na.add(mergemat[1,])
			i<-2
			while (i<=nrow(mergemat))
	{
		helpvec1<-mergemat[i,][1]
				helpvec2<-mergemat[i,][2]
				if (helpvec1>0)
		{
			helpvec1<-na.exclude(myorder[helpvec1,])
		}
		if (helpvec2>0)
		{
			helpvec2<-na.exclude(myorder[helpvec2,])
		}
		myorder[i,]<-na.add(c(helpvec1,helpvec2))
				i<-i+1
	}
	return(-myorder[i-1,])
}
typenomaker<-function(label)
{
helpvec<-c()
for (i in 1:max(label))
      {
      helpvec[i]<-sum(label==i)
      }
return(helpvec)
}
id.order<-order(relabel(x.id))
repno<-as.vector(table(x.id[id.order]))
baybi=.C("RbaybiclustG",PACKAGE="baybi",as.double(x), as.integer(nrow(x)), 
as.integer(ncol(x)), as.integer(repno), as.integer(length(repno)), 
as.double(transformed.par),rowmerge=as.integer(rep(0,(length(repno)-1)*2)), 
colmerge=as.integer(rep(0,(ncol(x)-1)*2)), 
rowheight=as.double(rep(0,(length(repno)-1))),colheight=as.double(rep(0,ncol(x)-1)))
rowbtree=NULL
rowbtree$repno.unordered=repno
rowbtree$merge=matrix(baybi$rowmerge,ncol=2,byrow=TRUE)
rowbtree$order=nocrossing.order(rowbtree$merge)
rowbtree$height=baybi$rowheight
rowbtree$logposterior<- -rowbtree$height
rowminheight<-min(rowbtree$height)
rowminindex<-which(rowbtree$height==rowminheight)
rowincrement<-diff(rowbtree$height,lag=1)
rowbtree$height<-diffinv(abs(rowincrement),lag=1)
rowbtree$clust.number<-((length(repno)-1):1)
rowbtree$cut<-rowbtree$height[rowminindex]
oldClass(rowbtree)="hclust"
rowbtree$optim.alloc<-cutree(rowbtree,h=rowbtree$cut)
rowbtree$optim.clustno<-max(rowbtree$optim.alloc)
rowbtree$labels.unordered<-row.labels
rowbtree$labels.ordered<-row.labels[rowbtree$order]

colbtree=NULL
colbtree$merge=matrix(baybi$colmerge,ncol=2,byrow=TRUE)
colbtree$order=nocrossing.order(colbtree$merge)
colbtree$height=baybi$colheight
colbtree$logposterior<- -colbtree$height
colminheight<-min(colbtree$height)
colminindex<-which(colbtree$height==colminheight)
colincrement<-diff(colbtree$height,lag=1)
colbtree$height<-diffinv(abs(colincrement),lag=1)
colbtree$clust.number<-((ncol(x)-1):1)
colbtree$cut<-colbtree$height[colminindex]
colbtree$optim.alloc<-cutree(colbtree,h=colbtree$cut)
colbtree$optim.clustno<-max(colbtree$optim.alloc)
oldClass(colbtree)="hclust"
colbtree$labels.unordered=col.labels
colbtree$labels.ordered=col.labels[colbtree$order]
roworder=rowbtree$order
colorder=colbtree$order
y=x[,colorder]
newy=c()
for (i in roworder)
{
newy=rbind(newy,y[which(rep(1:length(repno),repno)==i),])
}
y=newy
repno=repno[roworder]
rowbtree$repno.ordered=repno

rownames(y)=row.labels[rep(rowbtree$order,repno)]
coltypeno=rep(1,ncol(y))
rowtypeno=typenomaker(relabel(rowbtree$optim.alloc[rowbtree$order]))
lvsrow1 <- matrix(.C("Rlogmarg1dataG",PACKAGE="baybi", y=as.double(y), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)),repno=as.integer(repno), 
nrepno=as.integer(length(repno)),rowtypeno=as.integer(rowtypeno), 
nrowtypeno=as.integer(length(rowtypeno)), coltypeno=as.integer(coltypeno), ncoltypeno=as.integer(length(coltypeno)), 
theta=as.double(transformed.par),result=as.double(rep(0,length(rowtypeno)*length(coltypeno))))$result,
ncol=length(coltypeno))
rowtypeno=length(repno)
lvsrow0 <- matrix(.C("Rlogmarg1dataG",PACKAGE="baybi", y=as.double(y), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)),repno=as.integer(repno), 
nrepno=as.integer(length(repno)),rowtypeno=as.integer(rowtypeno), 
nrowtypeno=as.integer(length(rowtypeno)), coltypeno=as.integer(coltypeno), ncoltypeno=as.integer(length(coltypeno)), 
theta=as.double(transformed.par),result=as.double(rep(0,length(rowtypeno)*length(coltypeno))))$result,ncol=length(coltypeno))
colbtree$imp.ordered=as.vector(apply(lvsrow1,2,sum)-lvsrow0)
colbtree$imp.unordered=colbtree$imp.ordered[order(colorder)]

rowtypeno=rep(1,length(repno))
coltypeno=typenomaker(relabel(colbtree$optim.alloc[colbtree$order]))

lvscol1 <- matrix(.C("Rlogmarg1dataG",PACKAGE="baybi", y=as.double(y), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)),repno=as.integer(repno), 
nrepno=as.integer(length(repno)),rowtypeno=as.integer(rowtypeno), 
nrowtypeno=as.integer(length(rowtypeno)), coltypeno=as.integer(coltypeno), ncoltypeno=as.integer(length(coltypeno)), 
theta=as.double(transformed.par),result=as.double(rep(0,length(rowtypeno)*length(coltypeno))))$result,ncol=length(coltypeno))

coltypeno=ncol(y)
lvscol0 <- matrix(.C("Rlogmarg1dataG",PACKAGE="baybi", y=as.double(y), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)),repno=as.integer(repno), 
nrepno=as.integer(length(repno)),rowtypeno=as.integer(rowtypeno), 
nrowtypeno=as.integer(length(rowtypeno)), coltypeno=as.integer(coltypeno), ncoltypeno=as.integer(length(coltypeno)), 
theta=as.double(transformed.par),result=as.double(rep(0,length(rowtypeno)*length(coltypeno))))$result,nrow=1)
rowbtree$imp.ordered=as.vector(apply(lvscol1,1,sum)-lvscol0)
rowbtree$imp.unordered=rowbtree$imp.ordered[order(roworder)]
baybiobj=list(rowtree=rowbtree,coltree=colbtree,data.ordered=y,data.unordered=x)
class(baybiobj)="baybimp"
return(baybiobj)
}






baybiplot<-function(x,
xlab=x$coltree$labels.unordered,
ylab=x$rowtree$labels.unordered,
xlab.cex=1,
ylab.cex=1,
horizdendrogram.lwd=1,
horizdendrogram.size=2,
vertdendrogram.lwd=1,
vertdendrogram.size=2,
xlab.mar=6,
ylab.mar=3,
image.col=cm.colors(12),
horizbar.col=rev(c(heat.colors(5)[-4],"white")),
vertbar.col=rev(c(heat.colors(5)[-4],"white")),
horizbar.size=0.25,
vertbar.size=0.25,
horizteeth.size=0.25,
vertteeth.size=0.25,
image.width=5,
image.height=3)
{
######### Defining some functions
teethplotvh<-function(x,teeth.space=0.25,teeth.lwd=1,horiz=FALSE)
{
breakpoints<-function(label)
{
result<-c()
j<-1
for (i in 1:(length(label)-1))
    {
    if (!(label[i]==label[i+1])) {result[j]<-i;j<-j+1}
    }
return(result)
}
typenomaker<-function(label)
{
helpvec<-c()
for (i in 1:max(label))
      {
      helpvec[i]<-sum(label==i)
      }
return(helpvec)
}



if (horiz) {
label.dr<-x$optim.alloc[order.dendrogram(as.dendrogram(x))]
if ((teeth.space>0.25)| (teeth.space<0)) stop("teeth.space is wrongly adjusted")
if (teeth.lwd<0) stop("teeth.lwd is wrongly adjusted")
label.dr<-relabel(label.dr)
breakpoints.dr<-breakpoints(label.dr)
image(1:3, 1:length(label.dr), matrix(0,3,length(label.dr)), xlim = c(0, 
        1.2 ), ylim = c(0.5,length(label.dr) + 0.5), axes=FALSE,xlab = "", 
        ylab = "", col="white") #this makes plots comparable in axes

points(c(0,1),c(0.5+(0.5-teeth.space),0.5+(0.5-teeth.space)),
type="l",lwd=teeth.lwd)
points(c(0,1),c(length(label.dr)+teeth.space,
length(label.dr)+teeth.space),type="l",lwd=teeth.lwd)

if (max(label.dr)>1){ 
for (i in 1:length(breakpoints.dr))
{
points(c(0,1),c(breakpoints.dr[i]+0.5+teeth.space,breakpoints.dr[i]+0.5+teeth.space),
type="l",lwd=teeth.lwd)
points(c(0,1),c(breakpoints.dr[i]+0.5-teeth.space,breakpoints.dr[i]+0.5-teeth.space),
type="l",lwd=teeth.lwd)
}}

typeno.dr<-typenomaker(label.dr)
points(c(1,1),c(sum(typeno.dr)+teeth.space,
sum(typeno.dr)-typeno.dr[length(typeno.dr)]+0.5+teeth.space),type="l",lwd=teeth.lwd)

points(c(1,1),c(0.5+(0.5-teeth.space),
typeno.dr[1]+0.5-teeth.space),type="l",lwd=teeth.lwd)

if (length(typeno.dr)>2)
 {
for (i in 2:(length(typeno.dr)-1))
  {
  points(c(1,1),c(0.5+breakpoints.dr[i-1]+teeth.space,
  breakpoints.dr[i]+0.5-teeth.space),type="l",lwd=teeth.lwd)
  }
 }
} else {




label.dr<-x$optim.alloc[order.dendrogram(as.dendrogram(x))]
if ((teeth.space>0.25)| (teeth.space<0)) stop("teeth.space is wrongly adjusted")
if (teeth.lwd<0) stop("teeth.lwd is wrongly adjusted")
label.dr<-relabel(label.dr)
breakpoints.dr<-breakpoints(label.dr)
image(1:length(label.dr),1:3, t(matrix(0,3,length(label.dr))), ylim = c(0, 
        1.2 ), xlim = c(0.5,length(label.dr) + 0.5), axes=FALSE,xlab = "", 
        ylab = "", col="white") #this makes plots comparable in axes

points(c(0.5+(0.5-teeth.space),0.5+(0.5-teeth.space)),c(0,1),
type="l",lwd=teeth.lwd)
points(c(length(label.dr)+teeth.space,
length(label.dr)+teeth.space), c(0,1),type="l",lwd=teeth.lwd)

if (max(label.dr)>1){ 
for (i in 1:length(breakpoints.dr))
{
points(c(breakpoints.dr[i]+0.5+teeth.space,breakpoints.dr[i]+0.5+teeth.space),c(0,1),
type="l",lwd=teeth.lwd)
points(c(breakpoints.dr[i]+0.5-teeth.space,breakpoints.dr[i]+0.5-teeth.space),c(0,1),
type="l",lwd=teeth.lwd)
} }

typeno.dr<-typenomaker(label.dr)
points(c(sum(typeno.dr)+teeth.space,
sum(typeno.dr)-typeno.dr[length(typeno.dr)]+0.5+teeth.space),c(0,0),type="l",lwd=teeth.lwd)

points(c(0.5+(0.5-teeth.space),
typeno.dr[1]+0.5-teeth.space),c(0,0),type="l",lwd=teeth.lwd)

if (length(typeno.dr)>2)
 {
for (i in 2:(length(typeno.dr)-1))
  {
  points(c(0.5+breakpoints.dr[i-1]+teeth.space,
  breakpoints.dr[i]+0.5-teeth.space),c(0,0),type="l",lwd=teeth.lwd)
  }
 }
}
}

relabel<-function(current.label)
{
help.vec<-rep(NA,length(current.label))
j<-1
labelvalue<-current.label[1]
help.vec [which(current.label==labelvalue)]<-j
current.label[which(current.label==labelvalue)]<-0

for (i in 2:length(current.label))
        {
        if ((current.label[i]!=labelvalue) & (current.label[i]!=0))
                {
                labelvalue<-current.label[i]
                j<-j+1
                help.vec [which(current.label==labelvalue)]<-j
                current.label[which(current.label==labelvalue)]<-0
                } 
        }
return(help.vec)
}

gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}
logbfcoder<-function(x)
{
result<-c()
for (i in 1:length(x))
   {
           if (x[i]<= 0 ) {result[i]<- 0}
                    else{if (x[i]< 1 ) {result[i]<- 1}
                          else{if (x[i]< 3 ) {result[i]<- 2}
                                else{if (x[i]< 5 ) {result[i]<- 3}
                                      else{result[i]<- 4}
                                    }
                               }     
                        }
                   
     }
return(result)
}

cutplot.dendrogramh = function(x, h, cluscol=NULL, leaflab= "none", lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    
    # Not nice, but necessary
    pn  = stats:::plotNode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]],axes=FALSE, horiz=TRUE,leaflab="none", yaxs="i",xaxs="i",...)
    
    x = x[[2]]
    K = length(x)
    if (is.null(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=TRUE)
        x1 = x2 + 1
   }
abline(v=h,col="gray",lwd=lwd,lty=2)   
}

#cutplot.dendrogramh(as.dendrogram(x$rowtree), h=x$rowtree$cut)


cutplot.dendrogramv = function(x, h, cluscol=NULL, leaflab= "none", lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    
    # Not nice, but necessary
    pn  = stats:::plotNode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]],axes=FALSE, horiz=FALSE,leaflab="none", yaxs="i",xaxs="i",...)
    
    x = x[[2]]
    K = length(x)
    if (is.null(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=FALSE)
        x1 = x2 + 1
   }
abline(h=h,col="gray",lwd=lwd,lty=2)   
}



logbfcoder<-function(x)
{
result<-c()
for (i in 1:length(x))
   {
           if (x[i]<= 0 ) {result[i]<- 0}
                    else{if (x[i]< 1 ) {result[i]<- 1}
                          else{if (x[i]< 3 ) {result[i]<- 2}
                                else{if (x[i]< 5 ) {result[i]<- 3}
                                      else{result[i]<- 4}
                                    }
                               }     
                        }
                   
     }
return(result)
}

######### Defining some functions Finished


x.data<-(matsum(x$data.unordered,x$rowtree$repno.unordered)/x$rowtree$repno.unordered)

#cutplot.dendrogramv(as.dendrogram(x$coltree), h=x$coltree$cut)


layout(matrix(c(0,0,4,0,0,0,7,0,2,6,1,3,0,0,5,0),4,4,byrow=TRUE),
c(horizdendrogram.size,vertbar.size,image.width,horizteeth.size), 
c(vertdendrogram.size,horizbar.size,image.height,vertteeth.size), respect=TRUE) # defines the space of plot

#layout.show(lo)

if (is.null(xlab)) {image.bmargin<-0} else {image.bmargin<-xlab.mar} #bottom margin for xlabels
if (is.null(ylab)) {image.rmargin<-0.2} else {image.rmargin<-ylab.mar} #right margin for ylabels

            xrow.dendro <- as.dendrogram(x$rowtree)
            xcol.dendro <- as.dendrogram(x$coltree)
           rowInd <- order.dendrogram(xrow.dendro)
           colInd <- order.dendrogram(xcol.dendro)
            x.data<-x.data[rowInd,colInd] 
            ylab<-ylab[rowInd]
            xlab<-xlab[colInd]
# plot 1
par(mar=c(image.bmargin,0,0,image.rmargin))
image(1:ncol(x.data), 1:nrow(x.data), t(x.data), axes = FALSE, xlim = c(0.5, 
       ncol(x.data) + 0.5), ylim = c(0.5, nrow(x.data) + 0.5), xlab = "", 
       ylab = "", col=image.col)
#writes labels
if(!is.null(xlab)) {axis(1, 1:ncol(x.data), las = 2, line = -0.5, tick = 0, 
            labels = xlab, cex.axis = xlab.cex)}
if(!is.null(ylab)){axis(4, 1:nrow(x.data), las = 2, line = -0.5, tick = 0, 
        labels = ylab, cex.axis =ylab.cex)} #r.cex)

par(mar=c(image.bmargin,0.5,0,0))
# plot 2
cutplot.dendrogramh(xrow.dendro, h=x$rowtree$cut, lwd=horizdendrogram.lwd)




#plot 3
    par(mar=c(image.bmargin,0,0,0))

    teethplotvh(x$rowtree,horiz=TRUE)


# plot 4
par(mar=c(0,0,0,image.rmargin))
cutplot.dendrogramv(xcol.dendro, h=x$coltree$cut, lwd=vertdendrogram.lwd)
teethplotvh(x$coltree,horiz=FALSE)

horizimp=logbfcoder(x$rowtree$imp.ordered)

par(mar=c(image.bmargin,0,0,0.5))
if(sum(horizimp>0)>0){
image(1,1:length(rowInd),matrix(horizimp,nrow=1), 
axes = FALSE, xlab = "",ylab ="",col=horizbar.col)} else {
cat("no variable is important \n")
image(1,1:length(rowInd),matrix(horizimp,nrow=1), 
axes = FALSE, xlab = "",ylab ="",col="white")
} 

vertimp=logbfcoder(x$coltree$imp.ordered)

par(mar=c(0.5,0,0,image.rmargin))
if(sum(vertimp>0)>0){
image(1:ncol(x.data),1,matrix(vertimp,ncol=1), 
axes = FALSE, xlab = "",ylab ="",col=vertbar.col)} else {
cat("no subject is important \n")
image(1:ncol(x.data),1,matrix(vertimp,ncol=1), 
axes = FALSE, xlab = "",ylab ="",col="white")
}
layout(1)
}




