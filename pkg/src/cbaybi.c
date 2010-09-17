# include <R.h>
# include <math.h>
# include <Rmath.h>
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include <R_ext/Utils.h>
     
 void R_CheckUserInterrupt(void);



void logmarg0tmG (double sigma2, double tau2eta, double mu, double mean[], double css[], int n,int r[], int nr, double  *result) /* here r works like repno */
{
	int i=0,j;
	while (i<n)
		{
		for (j=0; j<nr;j++)
			{
			double R=r[j]; /*avoids problem of division of an integer by an integer like R/2 */
			result[i] = - R/2 * log(2 * M_PI) - (R-1)/2 * log (sigma2)- log(R*tau2eta+sigma2)/2  - 
			css[i]/(2*sigma2)-pow(mean[i]-mu,2)/(2*(tau2eta+sigma2/R));
			i++;
			}
		}
}

void logmarg1tmG (double sigma2, double tau2eta, double tau2theta, double mu, double mean[], double css[], int n, int r[], int nr, double  *result)
{
logmarg0tmG(sigma2, tau2eta+tau2theta, mu, mean, css, n, r, nr, result);
}

void logmargtmG (double theta[], double mean[], double css[], int n, int r[], int nr, double *result)
{
double sigma2,tau2eta,tau2theta,mu,p,suml=0,l0[n],l1[n];
/*l0=(double *) R_alloc(n,sizeof(double));
l1=(double *) R_alloc(n,sizeof(double));*/
int i;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]);
mu=theta[3];
p=1/(1+exp(-theta[4]));
logmarg0tmG(sigma2, tau2eta, mu, mean, css, n, r, nr, l0);
logmarg1tmG(sigma2, tau2eta,tau2theta, mu, mean, css, n, r, nr, l1);
	for (i=0;i<n;i++)
	{
	if ( (l1[i]-l0[i])>0) 
	suml+=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml+=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=suml;
}

void RlogmargtmG (double *theta, double *mean, double *css, int *n, int *r, int *nr, double *result)
{
logmargtmG (theta, mean, css, *n, r, *nr,result);
}


void cssofmatrix_oneside (double X[],int nX, int N[], int nN, double *meanXN, double *cssXN)
	{
	int i=0,j,k,l=0;
	double sumx=0, sumx2=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				sumx2 += pow(X[i],2);
				i++;
				}
			meanXN[l] = sumx / N[j];
			cssXN[l]  = sumx2 - pow(sumx,2) / N[j];
			l++;
			sumx  = 0;
			sumx2= 0;
			}
		}
	}

/* matsum is like cssofmatrix, but with the result in integer instead of double and does not give the corrected sum of squares*/
void matsum_oneside (double X[],int nX, int N[], int nN, double *result)
	{
	int i=0,j,k,l=0;
	double sumx=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				i++;
				}
			result[l] = sumx;
			l++;
			sumx  = 0;
			}
		}
	}

void integermatsum_oneside (int X[],int nX, int N[], int nN, int *result)
	{
	int i=0,j,k,l=0;
	int sumx=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				i++;
				}
			result[l] = sumx;
			l++;
			sumx  = 0;
			}
		}
	}


	void intmatsum_oneside (int X[],int nX, int N[], int nN, int *result)
	{
	int i=0,j,k,l=0;
	double sumx=0;
	while (i<nX)
		{
		for (j=0;j<nN;j++)
			{
			for (k=0;k<N[j];k++)
				{
				sumx   += X[i];
				i++;
				}
			result[l] = sumx;
			l++;
			sumx  = 0;
			}
		}
	}


void cssofmatrix_twoside (double X[],int nrowX, int ncolX, int Nrow[], int nNrow, int Ncol[], int nNcol,double meanX[],double cssX[])
	{
	int krow=0,kcol=0,i=0,j=0,k=0,cumNrow=0,cumNcol=0;
	double sumx=0, sumx2=0;
	cumNcol=0;
	cumNrow=0;
	for (krow=0;krow<nNrow;krow++)
		{
		for (kcol=0;kcol<nNcol;kcol++)
			{
			for (i=cumNrow;i<(cumNrow+Nrow[krow]);i++)
				{
				for (j=cumNcol;j<(cumNcol+Ncol[kcol]);j++)
					{
					sumx=sumx+X[(i+j* nrowX)];
					sumx2=sumx2+X[(i+j*nrowX)]*X[(i+j*nrowX)];
					/*printf("index is %d \n", (1+));*/
					}
				}
			meanX[k]=sumx/ (Ncol[kcol]*Nrow[krow]);
			cssX[k]=sumx2-sumx*sumx/(Ncol[kcol]*Nrow[krow]);
			k++;
			sumx=0;
			sumx2=0;
			cumNcol=cumNcol+Ncol[kcol];
			}
		cumNcol=0;
		cumNrow=cumNrow+Nrow[krow];
		}
	}

void Rcssofmatrix_twoside (double *X,int *nrowX, int *ncolX, int *Nrow, int *nNrow, int *Ncol, int *nNcol,double *meanX,double *cssX)
	{
	cssofmatrix_twoside ( X, *nrowX,  *ncolX, Nrow, *nNrow, Ncol, *nNcol,meanX,cssX);
	}
	

void transpose(double y[],int nrowy,int ncoly,double *result)
{
int i,j,l=0;
	for (i=0;i<nrowy;i++)
	{
		for (j=0;j<ncoly;j++)
		{
		result[i+j*nrowy]=y[l];
		l++;
		}
	}
}
	
void matsum_twoside(double X[],int nrowX, int ncolX, int Nrow[], int nNrow, int Ncol[], int nNcol,
	double sumX[])
	{
	int krow=0,kcol=0,i=0,j=0,k=0,cumNrow=0,cumNcol=0;
	double sumx=0,tsum[nNcol*nNrow];
	/*tsum=(double *) R_alloc(nNcol*nNrow,sizeof(double));*/
	cumNcol=0;
	cumNrow=0;
		for (krow=0;krow<nNrow;krow++)
		{
		for (kcol=0;kcol<nNcol;kcol++)
			{
			for (i=cumNrow;i<(cumNrow+Nrow[krow]);i++)
				{
				for (j=cumNcol;j<(cumNcol+Ncol[kcol]);j++)
					{
					sumx=sumx+X[(i+j* nrowX)];
					}
				}
			tsum[k]=sumx;
			k++;
			sumx=0;
			cumNcol=cumNcol+Ncol[kcol];
			}
		cumNcol=0;
		cumNrow=cumNrow+Nrow[krow];
		}
	transpose(tsum,nNrow,nNcol,sumX);
	}
	
void Rmatsum_twoside(double *X,int *nrowX, int *ncolX, int *Nrow, int *nNrow, int *Ncol, int *nNcol,
	double *sumX)
	{
        matsum_twoside(X,*nrowX, *ncolX, Nrow, *nNrow, Ncol,  *nNcol,sumX);
	}


void Rcssofmatrix_oneside (double *X,int *nX, int *N, int *nN, double *meanXN, double *cssXN)
	{
	cssofmatrix_oneside (X,*nX, N, *nN, meanXN, cssXN);
	}

	
void logmarg0dataG (double y[], int nrowy, int ncoly, int repno[] , int nrepno, 
	 int rowtypeno[], int nrowtypeno, int coltypeno[], int ncoltypeno, 
	double theta[],double  *result)
{
double sigma2,tau2eta, mu, ybar[nrepno*ncoly],css[nrepno*ncoly], l0[nrepno*ncoly];
	/*l0=(double *) R_alloc((nrepno*ncoly), sizeof(double));
	ybar=(double *) R_alloc((nrepno*ncoly), sizeof(double));
	css=(double *) R_alloc((nrepno*ncoly), sizeof(double));*/
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
mu=theta[3];
cssofmatrix_oneside (y, (nrowy*ncoly), repno, nrepno, ybar, css);
logmarg0tmG(sigma2, tau2eta, mu, ybar, css, (nrepno*ncoly ) , repno , nrepno, l0);
matsum_twoside (l0, nrepno,ncoly , rowtypeno, nrowtypeno, 
	coltypeno, ncoltypeno, result);
}

void Rlogmarg0dataG (double *y, int *nrowy, int *ncoly, int *repno , int *nrepno, 
	 int *rowtypeno, int *nrowtypeno, int *coltypeno, int *ncoltypeno, 
	double *theta,double  *result)
{
logmarg0dataG (y, *nrowy,  *ncoly, repno ,  *nrepno, 
	  rowtypeno,  *nrowtypeno,  coltypeno,  *ncoltypeno, 
	 theta,  result);
}


void covmaker(double sigma2, double tau2eta, double tau2theta, int rawr[], int nrawr, int ncol,double *Cov)
{
int i, j, k,rsum=0,counter=0,rdone,rremain,ddone,dremain,r[nrawr*ncol],nr=nrawr*ncol;
double offone=tau2theta,offtwo=tau2theta+tau2eta,diag=tau2theta+tau2eta+sigma2;
/*r=(int *) R_alloc(nrawr*ncol, sizeof(int));*/
k=0;
for (j=0;j<ncol;j++)
	{
	for (i=0;i<nrawr;i++)
		{
		r[k]=rawr[i];
		k++;
		}
	}
for (i=0;i<nr;i++){rsum+=r[i];}
rdone=0;
rremain=rsum-r[0];
for (i=0;i<nr;i++)
	{
	ddone=0;
	dremain=r[i];
	for (j=0;j<r[i];j++)
		{
		if (rdone>0)
			{
			for (k=0;k<(rdone);k++)
				{
				Cov[counter]=offone;
				counter++;
				}
			}	
		if (ddone>0) 
			{
			for (k=0;k<(ddone);k++)
				{
				Cov[counter]=offtwo;
				counter++;
				}
			}
		Cov[counter]=diag;
		counter++;
		dremain--;
		ddone++;
		if (dremain>0)
			{
			for (k=0;k<(dremain);k++)
				{
				Cov[counter]=offtwo;
				counter++;
				}
			}
		if(rremain>0)
			{
			for (k=0;k<(rremain);k++)
				{
				Cov[counter]=offone;
				counter++;
				}
			}
		}
rdone=rdone+r[i];
if (i<nr-1) {rremain=rremain-r[i+1];}
	}
}

void dmvnorm(double X[], int nX, double Mu[], double Sigma[],double *logdensity)
{

char myL='L',myN='N';
int info=0,N=nX,KD=nX,LD=KD+1,one=1,LDA=nX,LDB=nX,i=0;

/*for (i=0;i<nX*nX;i++){Rprintf("cov is %f i is %d \n ",Sigma[i],(i+1));}*/

double Chol[nX*nX],Centered[nX];
/*Chol=(double *) R_alloc(nX*nX, sizeof(double));/*used as Cholesky factor of Sigma*/
/*Centered=(double *) R_alloc(nX, sizeof(double)); /*used to centered density values*/




double d2=0;/*gives Mahalanobis distance*/
double sumlogeigen=0; /*gives sum of log eigen values of Cholesky factor, useful as to be traslated as log determinant of Sigma*/
for (i=0;i<(nX*nX);i++)
	{
	Chol[i]=Sigma[i];
	}
for (i=0;i<nX;i++)
	{
	Centered[i]=X[i]-Mu[i];
	}
	

F77_CALL(dpbtrf)(&myL, &N, &KD , Chol, &LD,&info);
/*Now Sigma is B, a lower triangular matrix that Sigma =B'B */

F77_CALL(dtrtrs)(&myL, &myN, &myN, &N, &one, Chol, &LDA,  Centered ,&LDB, &info);
/* Now Centred is inverse of Chol times (X-Mu)*/

	for (i=0;i<nX;i++)
	{
	d2 += Centered[i]*Centered[i];
	sumlogeigen += log(Chol[i*nX+i]);
	}
logdensity[0]=-0.5*(nX*log(2*M_PI)+d2)-sumlogeigen;
}



void logmarg1dataG(double y[], int nrowy,int ncoly, int repno[] , 
int nrepno, int rowtypeno[], int nrowtypeno,int coltypeno[], int ncoltypeno, 
double theta[],double  *result)
{
int maxdim=imax2(nrowy,ncoly),rtsum[nrepno],subrepno[maxdim*nrepno];
/*double sigma2,tau2eta, tau2theta, mu, *logdensity, 
	*submean,*subcov, *suby,*tresult;
submean=(double *) R_alloc(maxdim*maxdim, sizeof(double));
subcov=(double *) R_alloc(maxdim*maxdim, sizeof(double));
subrepno=(int *) R_alloc(maxdim*nrepno, sizeof(int));
suby=(double *) R_alloc(maxdim*maxdim, sizeof(double));
subcov=(double *) R_alloc(maxdim*maxdim*maxdim*maxdim,sizeof(double));
rtsum=(int *) R_alloc(nrepno, sizeof(int));
logdensity=(double *) R_alloc(1, sizeof(double));
tresult=(double *) R_alloc(nrowtypeno*ncoltypeno, sizeof(double));
*/
	
double sigma2,tau2eta, tau2theta, mu, logdensity[1], 
	submean[maxdim*maxdim],subcov[maxdim*maxdim*maxdim*maxdim], 
	suby[maxdim*maxdim],tresult[nrowtypeno*ncoltypeno];
/*rtsum=(int *) R_alloc(nrepno, sizeof(int));*/


int krow,kcol,i=0,j=0,k=0,l=0,cumrtrow=0,cumtcol=0,cumtrow=0;
sigma2=exp(theta[0]);
tau2eta=exp(theta[1]);
tau2theta=exp(theta[2]); 
mu=theta[3];
intmatsum_oneside(repno,nrepno,rowtypeno,nrowtypeno,rtsum);
/*Rprintf("nrowtypeno is %d, ncoltypeno is %d \n ",nrowtypeno,ncoltypeno);*/
for (krow=0;krow<nrowtypeno;krow++)
		{
		for (kcol=0;kcol<ncoltypeno;kcol++)
			{
				for (j=cumtcol;j<(cumtcol+coltypeno[kcol]);j++)
					{
						for (i=cumrtrow;i<(cumrtrow+rtsum[krow]);i++)
							{
							suby[k]=y[(i+j* nrowy)];
							/*Rprintf("krow is %d, kcol is %d, y index is %d \n ",krow+1,kcol+1,(i+j* nrowy)+1);*/
							k++;
						}
					}
				cumtcol=cumtcol+coltypeno[kcol];
				for (j=0;j<coltypeno[kcol];j++)
					{
						for (i=cumtrow;i<(cumtrow+rowtypeno[krow]);i++)
							{
							subrepno[j*rowtypeno[krow]+(int)(i-cumtrow)]=repno[i];
							}
					}
				for (i=0;i<k;i++){submean[i]=mu;}

/*				for (i=0;i<k;i++){Rprintf("suby is %f i is %d \n ",suby[i],(i+1));}*/

				covmaker(sigma2,tau2eta,tau2theta,subrepno,rowtypeno[krow],coltypeno[kcol],subcov);
				dmvnorm(suby,k,submean,subcov,logdensity);
				/*printf("krow is %d kcol is %d index is %d log denstiy is %f \n",krow+1,kcol+1,l,*logdensity);*/
				tresult[l]=*logdensity;l++;
				k=0;
			}
		cumtcol=0;
		cumrtrow=cumrtrow+rtsum[krow];
		cumtrow=cumtrow+rowtypeno[krow];
		}
transpose(tresult,nrowtypeno,ncoltypeno,result);
}



void Rlogmarg1dataG(double *y, int *nrowy,int *ncoly, int *repno , 
int *nrepno, int *rowtypeno, int *nrowtypeno,int *coltypeno, int *ncoltypeno, 
double *theta,double  *result)
{
	logmarg1dataG(y,  *nrowy, *ncoly,  repno , 
        *nrepno, rowtypeno, *nrowtypeno, coltypeno, *ncoltypeno, 
       theta,result);
}

void Rcovmaker(double *sigma2, double *tau2eta, double *tau2theta, int *r, int *nr, int *ncol, double *Cov)
{
	covmaker(*sigma2, *tau2eta, *tau2theta, r, *nr, *ncol, Cov);
}

void Rdmvnorm(double *X, int *nX, double *Mu, double *Sigma,double *logdensity)
{
	dmvnorm(X, *nX, Mu, Sigma,logdensity);
}


void sumlogmargdataG (double y[], int nrowy, int ncoly, int repno[] , int nrepno, int rowtypeno[], int nrowtypeno, int coltypeno[], int ncoltypeno, double theta[],double  *result)
{
double l0[nrowtypeno*ncoltypeno], l1[nrowtypeno*ncoltypeno],p, suml=0 ;
int i;
/*l0=(double *) R_alloc(nrowtypeno*ncoltypeno, sizeof(double));
l1=(double *) R_alloc(nrowtypeno*ncoltypeno, sizeof(double));*/

p=1/(1+exp(-theta[4]));
logmarg0dataG (y, nrowy, ncoly, repno , nrepno, rowtypeno, nrowtypeno,coltypeno, ncoltypeno, theta, l0);
logmarg1dataG (y, nrowy, ncoly, repno , nrepno, rowtypeno, nrowtypeno,coltypeno, ncoltypeno, theta, l1);
for (i=0; i<(nrowtypeno*ncoltypeno); i++)
	{
	if ( (l1[i]-l0[i])>0) /* this "if" calculates marginal likelihood without overflow error, more precisely*/
	suml+=l1[i] + log(p + (1-p) * exp(l0[i]-l1[i]));
	else 
	suml+=l0[i] + log((1-p) + p * exp(l1[i]-l0[i])); 
	}
*result=suml;
}


void RsumlogmargdataG (double *y, int *nrowy, int *ncoly, int *repno , int *nrepno, 
int *rowtypeno, int *nrowtypeno, int *coltypeno, int *ncoltypeno, double *theta,double  *result)
{
	sumlogmargdataG (y, *nrowy, *ncoly, repno , *nrepno, 
	rowtypeno, *nrowtypeno, coltypeno, *ncoltypeno, theta,result);
}



void eliminaterowsofmat(double y[],  int nrowy, int ncoly, int irow1,int nrowpick1, 
int irow2, int nrowpick2,int startfill, double *result)
{
int i,j,l=0;
	for (i=1; i<=nrowy;i++)
	{
	for (j=1 ; j<=ncoly; j++ )
		{
		if( ((i>=irow1) && (i< (irow1+nrowpick1))) || ((i>=irow2) && (i< (irow2+nrowpick2))) ) 
			{break;}
			else {
				result[l+startfill]=y[((i-1)+nrowy*(j-1) )];
				l++;
			}
		}
	}
}




void integereliminaterowsofmat(int y[],  int nrowy, int ncoly, int irow1,int nrowpick1, 
int irow2, int nrowpick2,int startfill, int *result)
{
int i,j,l=0;
	for (i=1; i<=nrowy;i++)
	{
	for (j=1 ; j<=ncoly; j++ )
		{
		if( ((i>=irow1) && (i< (irow1+nrowpick1))) || ((i>=irow2) && (i< (irow2+nrowpick2))) ) 
			{break;}
			else {
				result[l+startfill]=y[((i-1)+nrowy*(j-1) )];
				l++;
			}
		}
	}
}



void eliminatecolsofmat(double y[],  int nrowy, int ncoly, int jcol1,int ncolpick1, 
int jcol2, int ncolpick2,int startfill, double *result)
{
int i,j,l=0;
for (j=1 ; j<=ncoly; j++ )
	{
		{
		for (i=1; i<=nrowy;i++)
		if( !(((j>=jcol1) && (j< (jcol1+ncolpick1))) || ((j>=jcol2) && (j< (jcol2+ncolpick2)))) ) 
			{
				result[l+startfill]=y[((i-1)+nrowy*(j-1) )];
				l++;
			} 
		}
	}
}

void Reliminatecolsofmat(double *y,  int *nrowy, int *ncoly, int *jcol1,int *ncolpick1, 
int *jcol2, int *ncolpick2,int *startfill, double *result)
{
	eliminatecolsofmat(y,  *nrowy, *ncoly, *jcol1,*ncolpick1,*jcol2, *ncolpick2,*startfill, result);
}


void ldirich(int N[], int nN,double *result)
{
int i;
double sumN=0, sumlagmmaN1=0 ;
	for (i=0;i<nN;i++)
	{
		sumN += N[i];
		sumlagmmaN1 += lgammafn((double)(N[i]+1));
	}

*result= lgammafn((double)nN)+sumlagmmaN1-log((double) sumN)-lgammafn((double) (sumN+nN));
}

void pickuprowsofmat(double y[],  int nrowy, int ncoly, int irowstart,int nrowpick, int startfill, double *result)
{
int i,j,l=0;
	for (j=1; j<=nrowpick;j++)
	{
	for (i=0 ; i<ncoly; i++ )
		{
		result[l+startfill]=y[i*nrowy+j+irowstart-2];
		l++;
		}
	}
}

void integerpickuprowsofmat(int y[],  int nrowy, int ncoly, int irowstart,int nrowpick, int startfill, int *result)
{
int i,j,l=0;
	for (j=1; j<=nrowpick;j++)
	{
	for (i=0 ; i<ncoly; i++ )
		{
		result[l+startfill]=y[i*nrowy+j+irowstart-2];
		l++;
		}
	}
}


void pickupcolsofmat(double y[], int nrowy, int ncoly, int jcolstart,int ncolpick, int startfill, double *result)
{
int j,l=0;
		for (j=0; j<(nrowy*ncolpick);j++)
		{
		result[l+startfill]=y[(jcolstart-1)*nrowy+j];
		l++;
		}
	
}


void Rpickupcolsofmat(double *y, int *nrowy, int *ncoly, int *jcolstart,int *ncolpick, int *startfill, double *result)
{
	pickupcolsofmat(y, *nrowy, *ncoly, *jcolstart,*ncolpick, *startfill, result);
}


void yrowreorder(double y[], int nrowy, int ncoly, int repno[], int nrepno, int typeno[], int ntypeno, int noi, int  noj, double *result)
{
	/* for noi=1 or noj=1 it won't work*/
double tresult[nrowy*ncoly], irowstart1=1,irowstart2=1;
int i,j, startfill=0,repnosum[ntypeno];
	

/*repnosum=(int *) R_alloc(ntypeno, sizeof(int));
tresult=(double *) R_alloc(nrowy*ncoly, sizeof(double));*/

integermatsum_oneside (repno, nrepno, typeno, ntypeno, repnosum);

for (i=0;i<(noi-1);i++)
	{
	irowstart1 += repnosum[i];
	}
pickuprowsofmat(y,  nrowy, ncoly, (int)irowstart1, (int) repnosum[(noi-1)]  , startfill, tresult); 
/*it picks up the first part related to noi  and put it in the first line of "result" matrix*/

for (i=0;i<(noj-1);i++)
	{
	irowstart2 += repnosum[i];
	}
startfill=(repnosum[(noi-1)] ) * ncoly;
/*it picks up the second part related to noj  and put it in "result" immidiately after first part taken in the previous step*/
pickuprowsofmat(y,  nrowy, ncoly, (int)irowstart2, (int) repnosum[(noj-1)]  , startfill, tresult);

startfill=(repnosum[(noi-1)]+repnosum[(noj-1)] ) * ncoly;

/*it eiminates the values taken in the previous 2 steps and out it at last*/
eliminaterowsofmat(y,nrowy,ncoly, irowstart1, repnosum[(noi-1)] , 
irowstart2, repnosum[(noj-1)]  ,startfill, tresult);

/*we have arranged by rows but we read the matrix by columns so we need to transpose the result matrix*/
transpose(tresult, nrowy, ncoly, result);
}



void ycolreorder(double y[], int nrowy, int ncoly, int typeno[], int ntypeno, int noi, int  noj, double *result)
{
	/* for noi=1 or noj=1 it won't work*/
double tresult[nrowy*ncoly], icolstart1=1,icolstart2=1;
int i,j, startfill=0,*repnosum;
/*tresult=(double *) R_alloc(nrowy*ncoly, sizeof(double));*/

for (i=0;i<(noi-1);i++)
	{
	icolstart1 += typeno[i];
	}

for (i=0;i<(noj-1);i++)
	{
	icolstart2 += typeno[i];
	}

/*it picks up the first part related to noi  and put it in the first line of "result" matrix*/
pickupcolsofmat(y,  nrowy, ncoly, icolstart1, typeno[(noi-1)]  , startfill, result); 
startfill=nrowy*typeno[(noi-1)];
/*it picks up the second part related to noj  and put it in "result" immidiately after first part taken in the previous step*/
pickupcolsofmat(y,  nrowy, ncoly, icolstart2, typeno[(noj-1)]  , startfill, result); 
startfill=nrowy*(typeno[(noi-1)]+typeno[(noj-1)]);


/*it eiminates the values taken in the previous 2 steps and out it at last*/
eliminatecolsofmat(y,nrowy,ncoly, icolstart1, typeno[(noi-1)] , 
icolstart2, typeno[(noj-1)]  ,startfill, result);
}

void Rycolreorder(double *y, int *nrowy, int *ncoly, int *typeno, int *ntypeno, int *noi, int  *noj, double *result)
{
	ycolreorder(y, *nrowy, *ncoly, typeno, *ntypeno, *noi, *noj, result);
}

void repnoreorder(int repno[], int nrepno, int typeno[], int ntypeno, int noi, int noj, int  *result)
{
int irowstart1=1, irowstart2=1;
int i;
	
for (i=0;i<(noi-1);i++)
	{
	irowstart1 += typeno[i];
	}
integerpickuprowsofmat(repno,  nrepno, 1, irowstart1 , typeno[(noi-1)]  , 0, result); 

for (i=0;i<(noj-1);i++)
	{
	irowstart2 += typeno[i];
	}

integerpickuprowsofmat(repno,  nrepno, 1, irowstart2 , typeno[(noj-1)] , typeno[(noi-1)] , result); 
integereliminaterowsofmat(repno,nrepno,1,irowstart1 , typeno[(noi-1)], irowstart2, typeno[(noj-1)] ,typeno[(noi-1)]+typeno[(noj-1)]  , result);
}



void distmatrixG(double y[], int nrowy, int ncoly, int repno[], int nrepno, 
int rowtypeno[], int nrowtypeno,
int coltypeno[], int ncoltypeno, double theta[], double *minvalue,int  *minindex,int *isrow)
{
int i,j,minrowindex1,minrowindex2,mincolindex1,mincolindex2;
double newy[nrowy*ncoly], initlogmarg[1], initrowdirich[1],initcoldirich[1], distmat,minrowvalue,mincolvalue;
int newrepno[nrepno], newrowtypeno[nrowtypeno],newcoltypeno[ncoltypeno];
/*newy=(double *) R_alloc(nrowy*ncoly, sizeof(double));

newrepno=(int *) R_alloc(nrepno, sizeof(int));
newrowtypeno=(int *) R_alloc(nrowtypeno, sizeof(int));
newcoltypeno=(int *) R_alloc(ncoltypeno, sizeof(int));
*/
if (nrowtypeno>1)
{
newrowtypeno[0]=rowtypeno[0]+rowtypeno[1];
integereliminaterowsofmat(rowtypeno,  nrowtypeno, 1 , 1, 1, 2, 1, 1, newrowtypeno);
/*col remains fix*/
sumlogmargdataG (y, nrowy, ncoly, repno, nrepno, newrowtypeno, (nrowtypeno-1),coltypeno, ncoltypeno, theta, initlogmarg);
ldirich(newrowtypeno,(nrowtypeno-1),initrowdirich);
ldirich(coltypeno,ncoltypeno,initcoldirich); 
minrowvalue=-initlogmarg[0]/*-initrowdirich[0]-initcoldirich[0]*/;
minrowindex1=1; minrowindex2=2;

	for (i=1; i<=(nrowtypeno-1);i++)
		{
		for (j=(i+1); j<=nrowtypeno; j++)
			{
				R_CheckUserInterrupt();
				yrowreorder(y, nrowy, ncoly, repno, nrepno, rowtypeno, nrowtypeno, i, j, newy);
				repnoreorder(repno, nrepno, rowtypeno, nrowtypeno, i, j, newrepno);
				newrowtypeno[0]=rowtypeno[(i-1)]+rowtypeno[(j-1)];
				integereliminaterowsofmat(rowtypeno,  nrowtypeno, 1 , i, 1, j, 1, 1, newrowtypeno);
				ldirich(newrowtypeno,(nrowtypeno-1),initrowdirich);
				sumlogmargdataG (newy, nrowy, ncoly, newrepno, nrepno, 
				newrowtypeno, (nrowtypeno-1),coltypeno, ncoltypeno, theta, initlogmarg);
				distmat=-initlogmarg[0]/*-initrowdirich[0]-initcoldirich[0]*/;
				/*Rprintf("row: i=%d, j=%d, dist=%f \n",i,j,distmat);*/
				if (distmat<minrowvalue)
					{
					minrowvalue=distmat;
					minrowindex1=i; minrowindex2=j;
					}
			}			
		}
}

if (ncoltypeno>1)
{
newcoltypeno[0]=coltypeno[0]+coltypeno[1];
integereliminaterowsofmat(coltypeno,  ncoltypeno, 1 , 1, 1, 2, 1, 1, newcoltypeno);
/*row remaines fixed*/
sumlogmargdataG (y, nrowy, ncoly, repno, nrepno, rowtypeno, nrowtypeno,newcoltypeno, (ncoltypeno-1), theta, initlogmarg);
ldirich(rowtypeno,nrowtypeno,initrowdirich);
ldirich(newcoltypeno,(ncoltypeno-1),initcoldirich); 
mincolvalue=-initlogmarg[0]/*-initrowdirich[0]-initcoldirich[0]*/;
mincolindex1=1; mincolindex2=2;
	for (i=1; i<=(ncoltypeno-1);i++)
		{
			for (j=(i+1); j<=ncoltypeno; j++)
				{
				R_CheckUserInterrupt();
				ycolreorder(y, nrowy, ncoly, coltypeno, ncoltypeno, i, j, newy);
				newcoltypeno[0]=coltypeno[(i-1)]+coltypeno[(j-1)];
				integereliminaterowsofmat(coltypeno,  ncoltypeno, 1 , i, 1, j, 1, 1, newcoltypeno);
				sumlogmargdataG (newy, nrowy, ncoly, repno, nrepno, 
				rowtypeno, nrowtypeno,newcoltypeno, (ncoltypeno-1), theta, initlogmarg);
				distmat=-initlogmarg[0]/*-initrowdirich[0]-initcoldirich[0]*/;
				/*Rprintf("col: i=%d, j=%d,dist=%f\n",i,j,distmat);*/
				if (distmat<=mincolvalue)
					{
					mincolvalue=distmat;
					mincolindex1=i; mincolindex2=j ;
 					}
				}
		}
}

			
	/*distance calculations finishes here*/
	/*Rprintf("row is %f, col is %f\n",minrowvalue,mincolvalue);*/
        
		if ((ncoltypeno>1)&(minrowvalue<=mincolvalue))
		{
		*isrow=1;
		*minvalue=minrowvalue;
		minindex[0]=minrowindex1;
		minindex[1]=minrowindex2;
		}
	 
	
		if ((nrowtypeno>1)&(mincolvalue<=minrowvalue))
		{
		*isrow=0;
		*minvalue=mincolvalue;
		minindex[0]=mincolindex1;
		minindex[1]=mincolindex2;
		}
	
	if (ncoltypeno<=1)
	{
	*isrow=1;
	*minvalue=minrowvalue;
	minindex[0]=1;
	minindex[1]=2;
	}
	if (nrowtypeno<=1)
	{
		*isrow=0;
		*minvalue=mincolvalue;
		minindex[0]=1;
		minindex[1]=2;
	} 
	
	
}



void RdistmatrixG(double *y, int *nrowy, int *ncoly, int *repno, int *nrepno, 
int *rowtypeno, int *nrowtypeno,
int *coltypeno, int *ncoltypeno, double *theta, double *minvalue,int  *minrowindex,int *isrow)
{
	distmatrixG(y, *nrowy, *ncoly, repno, *nrepno, rowtypeno, *nrowtypeno,
 coltypeno, *ncoltypeno, theta, minvalue,minrowindex,isrow);
}



void baybiclustG(double y[], int nrowy, int ncoly, int repno[], int nrepno, double theta[],int  *rowmerge, int *colmerge, double *rowheight,double *colheight)
{
double yclust[nrowy*ncoly], minvalue[1], yresult[nrowy*ncoly];
int rowtypelabel[nrepno], rowtypelabelresult[nrepno], coltypelabel[ncoly], coltypelabelresult[ncoly],rowtypenoclust[nrepno],coltypenoclust[ncoly], repnoresult[nrepno],
	repnoclust[nrepno], rowtypenoresult[nrepno],coltypenoresult[ncoly],minindex[2],isrow[1];
/*
rowtypelabel=(int *) R_alloc(nrepno, sizeof(int));
rowtypelabelresult=(int *) R_alloc(nrepno, sizeof(int));
coltypelabel=(int *) R_alloc(ncoly, sizeof(int));
coltypelabelresult=(int *) R_alloc(ncoly, sizeof(int));
yclust=(double *) R_alloc(nrowy*ncoly, sizeof(double));
repnoclust=(int *) R_alloc(nrepno, sizeof(int));
rowtypenoclust=(int *) R_alloc(nrepno, sizeof(int));
coltypenoclust=(int *) R_alloc(ncoly, sizeof(int));
yresult=(double *) R_alloc(nrowy*ncoly, sizeof(double));
repnoresult=(int *) R_alloc(nrepno, sizeof(int));
rowtypenoresult=(int *) R_alloc(nrepno, sizeof(int));
coltypenoresult=(int *) R_alloc(ncoly, sizeof(int));
*/



int i,irow,icol,j, nrepnoclust=nrepno, nrowtypenoclust=nrepno,ncoltypenoclust=ncoly;
/*initialization*/
for 	(i=0; i<nrepno;i++)
	{
		rowtypelabel[i]=-(i+1);
		repnoclust[i]=repno[i];
		rowtypenoclust [i]=1;
	}
for 	(i=0; i<ncoly;i++)
	{
		coltypelabel[i]=-(i+1);
		coltypenoclust [i]=1;
	}

for ( i=0;i<(nrowy*ncoly);i++)
	{
	yclust[i]=y[i];
	}
irow=1;
icol=1;
while ( (nrowtypenoclust>1) || (ncoltypenoclust>1) )
	{
		distmatrixG(yclust, nrowy,  ncoly, repnoclust, nrepnoclust, rowtypenoclust, nrowtypenoclust, coltypenoclust, ncoltypenoclust, theta, minvalue, minindex,isrow);
		if (*isrow==1) 
			{
				/*Rprintf("dist=%f ",*minvalue);*/
				rowmerge[2*irow-2]=rowtypelabel[(minindex[0]-1)];
				rowmerge[2*irow-1]=rowtypelabel[(minindex[1]-1)];
				rowtypelabelresult[0]=irow;
				integereliminaterowsofmat(rowtypelabel,  nrowtypenoclust, 1 , (minindex[0]), 1, (minindex[1]), 1, 1, rowtypelabelresult);
				for (j=0;j<nrowtypenoclust;j++) {rowtypelabel[j]=rowtypelabelresult[j];}
				rowheight[irow-1]=minvalue[0];
				yrowreorder(yclust, nrowy, ncoly, repnoclust, nrepnoclust, rowtypenoclust, nrowtypenoclust, (minindex[0]), (minindex[1]), yresult);
				for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
				repnoreorder(repnoclust, nrepnoclust, rowtypenoclust, nrowtypenoclust, (minindex[0]), (minindex[1]),repnoresult);
				for (j=0;j<nrepnoclust;j++){repnoclust[j]=repnoresult[j];}
				rowtypenoresult[0]=rowtypenoclust[(minindex[0]-1)]+rowtypenoclust[(minindex[1]-1)];
				integereliminaterowsofmat(rowtypenoclust,  nrowtypenoclust, 1 , (minindex[0]), 1, (minindex[1]), 1, 1, rowtypenoresult);
				for (j=0;j<nrowtypenoclust;j++){rowtypenoclust[j]=rowtypenoresult[j];}
				nrowtypenoclust--;
				Rprintf("nrow=%d, ncol=%d\n",nrowtypenoclust,ncoltypenoclust);
				irow++;
			} else
			{
				/*Rprintf("dist=%f ",*minvalue);*/
				colmerge[2*icol-2]=coltypelabel[(minindex[0]-1)];
				colmerge[2*icol-1]=coltypelabel[(minindex[1]-1)];
				coltypelabelresult[0]=icol;
				integereliminaterowsofmat(coltypelabel,  ncoltypenoclust, 1 , (minindex[0]), 1, (minindex[1]), 1, 1, coltypelabelresult);
				for (j=0;j<ncoltypenoclust;j++) {coltypelabel[j]=coltypelabelresult[j];}

				colheight[icol-1]=minvalue[0];
				ycolreorder(yclust, nrowy, ncoly, coltypenoclust, ncoltypenoclust, (minindex[0]), (minindex[1]), yresult);
				for (j=0;j<(nrowy*ncoly);j++){yclust[j]=yresult[j];}
				coltypenoresult[0]=coltypenoclust[(minindex[0]-1)]+coltypenoclust[(minindex[1]-1)];
				integereliminaterowsofmat(coltypenoclust,  ncoltypenoclust, 1 , (minindex[0]), 1, (minindex[1]), 1, 1, coltypenoresult);
				for (j=0;j<ncoltypenoclust;j++){coltypenoclust[j]=coltypenoresult[j];}
				ncoltypenoclust--;
				Rprintf("nrow=%d, ncol=%d\n",nrowtypenoclust,ncoltypenoclust);
				icol++;
			}
	}
}

void RbaybiclustG(double *y, int *nrowy, int *ncoly, int *repno, int *nrepno, double *theta,int *rowmerge, int *colmerge, double *rowheight,double *colheight)
{
	baybiclustG(y, *nrowy, *ncoly, repno, *nrepno, theta,rowmerge, colmerge, rowheight,colheight);
}


int maxofvector(int vector[], int nvector)
	{
		int result;
		int i;
		result=vector[0];
		for (i=1;i<nvector;i++)
		{
		if (result<vector[i]) 
			{
				result=vector[i];
			}
		}
	return result;
	}
	
void findindexoflabel (int label[],int nlabel,int desiredlabel,int *result, int *nresult)
	{
		int i,j=0;
		for (i=0;i<nlabel;i++)
		{
			if (label[i]==desiredlabel)
			{
				result[j]=i+1;
				j++;
			}
		}
		*nresult=j;
	}

void yarrangebyrowlabel(double y[], int nrowy, int ncoly, int repno[], 
	int nrepno, int label[], int nlabel,double *ypickedup,int *repnopickedup, int *typenopickedup, int *ntypenopickedup)
	{
	int nindexofdesiredlabel[1],maxlab=maxofvector(label,nlabel);
	int indexofdesiredlabel[nlabel];
	double typickedup[ncoly*nrowy]; 
	int repnostartfill=0,desiredlabel=1,j,i=1,ystartfill=0,irowstarty,k;

		
		for (i=1;i<=maxlab;i++)
		{
			desiredlabel=i;
			findindexoflabel(label, nlabel,desiredlabel,indexofdesiredlabel,nindexofdesiredlabel);
			for  (j=0;j<*nindexofdesiredlabel;j++)
			{
			k=0;irowstarty=1;
			while (k<(indexofdesiredlabel[j])-1)
				{
				irowstarty+=repno[k];
				k++;
				}
			repnopickedup[repnostartfill]=repno[(indexofdesiredlabel[j]-1)];
			repnostartfill+=1;
			pickuprowsofmat(y,  nrowy, ncoly, irowstarty, (repno[indexofdesiredlabel[j]-1]) , ystartfill, typickedup); 
/*			printf("startfrom %d pick up %d start writing from %d \n",irowstarty,asint(repno[asint(indexofdesiredlabel[j])-1]),ystartfill);*/
			ystartfill+= ((repno[(indexofdesiredlabel[j])-1])*ncoly);
		}
	        typenopickedup[i-1]=*nindexofdesiredlabel;
		}
transpose(typickedup,nrowy,ncoly,ypickedup);
	*ntypenopickedup=maxlab;
	}
	
void Ryarrangebyrowlabel (double *y, int *nrowy, int *ncoly, int *repno, 
	int *nrepno, int *label, int *nlabel,double *ypickedup, int *repnopickedup, int *typenopickedup,int *ntypenopickedup)
	{
       yarrangebyrowlabel(y, *nrowy, *ncoly, repno, *nrepno, label, *nlabel,ypickedup, repnopickedup,typenopickedup,ntypenopickedup);
	}
	

	
	/*
void yarrangebytwolabels(double y[], int nrowy, int ncoly, int repno[], 
	int nrepno, double rowlabel[], int nrowlabel, double collabel[], int ncollabel, double *ypickedup,
       int *repnopickedup, int *rowtypenopickedup, int *nrowtypenopickedup,int *coltypenopickedup, int *ncoltypenopickedup)
	{
	int nindexofdesiredlabel[1],maxlab=maxofvector(label,nlabel);
	double indexofdesiredlabel[nlabel], typickedup[ncoly*nrowy]; 
	int repnostartfill=0,desiredlabel=1,j,i=1,ystartfill=0,irowstarty,k;

		
		for (i=1;i<=maxlab;i++)
		{
			desiredlabel=i;
			findindexoflabel(label, nlabel,desiredlabel,indexofdesiredlabel,nindexofdesiredlabel);
			for  (j=0;j<*nindexofdesiredlabel;j++)
			{
			k=0;irowstarty=1;
			while (k<(asint(indexofdesiredlabel[j])-1))
				{
				irowstarty+=repno[k];
				k++;
				}
			repnopickedup[repnostartfill]=repno[asint(indexofdesiredlabel[j]-1)];
			repnostartfill+=1;
			pickuprowsofmat(y,  nrowy, ncoly, irowstarty, asint(repno[asint(indexofdesiredlabel[j])-1]) , ystartfill, typickedup); 
			printf("startfrom %d pick up %d start writing from %d \n",irowstarty,asint(repno[asint(indexofdesiredlabel[j])-1]),ystartfill);
			ystartfill+= (asint(repno[asint(indexofdesiredlabel[j])-1])*ncoly);
		}
	        typenopickedup[i-1]=*nindexofdesiredlabel;
		}
transpose(typickedup,nrowy,ncoly,ypickedup);
	*ntypenopickedup=maxlab;
	}

*/



