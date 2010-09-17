baybiclust=function(x,transformed.par,row.labels=rownames(x),col.labels=colnames(x))
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
repno=rep(1,nrow(x))
baybi=.C("RbaybiclustG",PACKAGE="baybi",as.double(x), as.integer(nrow(x)), 
as.integer(ncol(x)), as.integer(repno), as.integer(length(repno)), 
as.double(transformed.par),rowmerge=as.integer(rep(0,(length(repno)-1)*2)), 
colmerge=as.integer(rep(0,(ncol(x)-1)*2)), 
rowheight=as.double(rep(0,(length(repno)-1))),colheight=as.double(rep(0,ncol(x)-1)))
rowbtree=NULL
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
rowbtree$labels<-row.labels
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
colbtree$labels=col.labels
roworder=rowbtree$order
colorder=colbtree$order
y=x[,colorder]
newy=c()
for (i in roworder)
{
newy=rbind(newy,y[which(rep(1:length(repno),repno)==i),])
}
y=newy
rownames(y)=row.labels[rowbtree$order]
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
colbtree$imp=as.vector(apply(lvsrow1,2,sum)-lvsrow0)

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
rowbtree$imp=as.vector(apply(lvscol1,1,sum)-lvscol0)
return(list(rowtree=rowbtree,coltree=colbtree,data=y))
}

