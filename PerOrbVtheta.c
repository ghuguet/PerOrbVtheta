//compile: g++ -o program PerOrbVtheta.c rk78.c -lm -lgsl -lgslcblas
//program output: A b suc z0 z1 op01 opn T1...Top01
//if the numbers obtained to compute the periodic orbit work fine suc(success)=1, otherwise suc=0 and results are not reliable.
//z0={0,1}, z1={0,1} depending whether it is present or not
//op01 integer that indicates number of periodic orbits with period >=2 which hit S0 and S1
//opn integer that indicates number of periodic orbits with period >=2 which hit at least one Sn with n>=2
//T1... Top01 period of the periodic orbits counted in op01


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>  
#include <stdlib.h> 
#include "rk78.h" 

using std::cout;
using std::endl;
using namespace std;

double A;
double twopi=8*atan(1.);
double T=0.5;
double dT;
double a=0.08;
double b;
double c=0.53;
double Iapp;

#pragma omp threadprivate(A,Iapp,b)
void vfield(double t, double *x, int ndim, double *dx);
int partition( double * a, int * index, int l, int r);
void quickSort( double * a, int * index, int l, int r);
void nper (double *xa, double *xf, int ifin, int *scount,std::ofstream& outfile);
int stroboscopic(double* x);
void sortsequence(int q,int *seq,int *perm);
int ismaximin(int q, int *X);
int invstroboscopic(double* x);
double rotnum(int q, int *X);

int main(int argc, char * argv[])
{

int cycles;
int number;
char name[50];
FILE *output;
int i,j;
double v, vmax;
int n;
double x0,x1;
int nite=100;
int ifin=0;
double *xf,*xa;
int newo;
double xap,tap;
int nxap;
double opn;
double bini,bfin;
double Aini=2.2;
double Afin=11;
int l;
double pas=0.001;
double pasA=1;
int Nb,NA;
int spike;
double x1max=15.;
  
  if(argc != 4){
    printf("usage: %s bini bfin Nb\n", argv[0]);
    abort();
  }

  sscanf(argv[1], "%lg", &bini);
  sscanf(argv[2], "%lg", &bfin);
 
  dT=0.5*T;
  
  sscanf(argv[3],"%d",&Nb);
  NA=2000;
  pasA=(Afin-Aini)/double(NA);

//#pragma omp parallel for private(j,i,ifin,x1,nxap,newo,xap,opn,x0,spike)
for (l=0;l<Nb;l++){
  b=bini+double(l)*pas;
  ofstream outfile;
  std::ostringstream s;
  s<<"output_T0d5_"<<"_b"<<b<<".tna";
  outfile.precision(20);
  outfile.open(s.str().c_str());	
    //#pragma omp parallel for private(i,ifin,x1,nxap,newo,xap,opn,x0,spike)
    for  (j=0;j<NA+1;j++){
      double *xa=new double[200];
      double *xf=new double[200];
      double *thetaf=new double[200];
      double *thetaa=new double[200];
      double *x=new double[2];
      int *scount=new int[200];
      A=Aini+double(j)*pasA;
      ifin=0;
      for (i=0;i<200;i++){
       xf[i]=0;
       thetaf[i]=0;
       xa[i]=0;
       thetaa[i]=0;
       scount[i]=0;
      } 
      x1=0;
      while (x1<x1max){
      x0=0;
      while (x0<x1){
       x[0]=x0;
       x[1]=x1;

    /*integrate*/
    nxap=0;
       spike=stroboscopic(x);
       //if (spike==1 || spike==0){
       if (spike==0){
	 for (i=1;i<nite;i++){
	   spike=stroboscopic(x);
	   if (i==nite-2 && nxap==0){
	     xap=x[0];
             tap=x[1];
	     nxap=1;
	   }
	 }
	 newo=0;
	 for (i=0;i<ifin && newo==0;i++){
	  if (((x[0]-xf[i])*(x[0]-xf[i])+(x[1]-thetaf[i])*(x[1]-thetaf[i]))<1.0e-6){newo=1;}
	 }
	 
	 if (newo==0){
	   xf[ifin]=x[0];
	   thetaf[ifin]=x[1];
	   xa[ifin]=xap;
	   thetaa[ifin]=tap;
	   scount[ifin]=spike;
	   ifin+=1;
	 }
       }
       x0+=0.1;
      }
      x1+=0.1;
      }  

//to add points 

      int *index1=new int[200];
      int *index2=new int[200];
      double *xa_aux=new double[200];
      double *xf_aux=new double[200];
      double *xtmp=new double[2];
      int suc;
      int nitew=0;
      double dist;
	suc=0;
	spike=0;
	while (suc==0 && nitew<30 && spike!=-1){

	  for (i=0;i<ifin;i++){
	    index1[i]=i;	
	    index2[i]=i;
	  }
  
	  for (i=0;i<ifin;i++){
	   xa_aux[i]=xa[i];
	   xf_aux[i]=xf[i];
	  }

	  quickSort(xa_aux, index1,0,ifin-1);
	  quickSort(xf_aux, index2,0,ifin-1);
 
 	  suc=1;
	  for (i=0;i<ifin && suc==1;i++){
	   if ((xa_aux[i]-xf_aux[i])>1.0e-3){ 
		suc=0;
	        x[0]=xf_aux[i];
		x[1]=thetaf[index2[i]];
	        xa[ifin]=x[0];
		thetaa[ifin]=x[1];
	        spike=stroboscopic(x);
	        xf[ifin]=x[0]; 
	        thetaf[ifin]=x[1];
                scount[ifin]=spike;
	        ifin++;
	        }

	  if ((xf_aux[i]-xa_aux[i])>1.0e-3){
		suc=0;
	        x[0]=xa_aux[i];
		x[1]=thetaa[index1[i]];
	        xf[ifin]=x[0];
	        thetaf[ifin]=x[1];
		xtmp[0]=x[0];
		xtmp[1]=x[1];
	        spike=invstroboscopic(x);
		xa[ifin]=x[0]; 
        	thetaa[ifin]=x[1];
                scount[ifin]=spike;
		stroboscopic(x);
		dist=pow(x[0]-xtmp[0],2.);
		dist+=pow(x[1]-xtmp[1],2.);
		if (dist>1e-6){
		  spike=-1;
		}
        	ifin++;
        	}
   
	   }
	   if (spike==-1){
	     suc=0;
	   }
  nitew++;
}

      if (spike!=-1){
      	nper(xa,xf,ifin,scount,outfile);
      }
      else{
	outfile<<A<<" "<<b<<" "<<0<<endl;
      }


      delete [] xa;
      delete [] xf;
      delete [] thetaf;
      delete [] thetaa;
      delete [] x;
      delete [] scount; 
      delete [] index1;
      delete [] index2;
      delete [] xa_aux;
      delete [] xf_aux;
      delete [] xtmp;


    }
  outfile.close();
}
  
}

void nper (double *xa, double *xf, int ifin, int * scount,std::ofstream& outfile){
  int i,j,k;
  int * index1, *index2;
  int *pos1, *pos2;
  int *map;
  int *vper, *full;
  double *vs;
  double *rot;
  double *rot2;
  int *arrow;
  int nper;
  int sc;
  int z0=0,z1=0,op01=0,opn=0;
  int count, count2;
  int *Tper,*seq,*imx;
  int *tper, *ixper;
  int snhit;
  int suc;
    
  index1=(int *)calloc(ifin, sizeof(int));
  index2=(int *)calloc(ifin, sizeof(int));
  pos1=(int *)calloc(ifin, sizeof(int));
  pos2=(int *)calloc(ifin, sizeof(int));
  map=(int *)calloc(ifin, sizeof(int));
  arrow=(int *)calloc(ifin, sizeof(int));
  vs=(double *)calloc(ifin, sizeof(double));
  vper=(int *)calloc(ifin, sizeof(int));
  full=(int *)calloc(ifin, sizeof(int));
  seq=(int *)calloc(ifin, sizeof(int));
  imx=(int *)calloc(ifin, sizeof(int));
  tper=(int *)calloc(ifin, sizeof(int));
  ixper=(int *)calloc(ifin, sizeof(int));
  rot=(double *)calloc(ifin, sizeof(double));
  rot2=(double *)calloc(ifin, sizeof(double));



  for (i=0;i<ifin;i++){
    index1[i]=i;
    index2[i]=i;
  }
  
  quickSort(xa, index1,0,ifin-1);
  quickSort(xf, index2,0,ifin-1);

  suc=1;
  for (i=0;i<ifin;i++){
   if ((xa[i]-xf[i])>1.0e-3){suc=0;}
  }
  

  for (i=0;i<ifin;i++){
   pos1[index1[i]]=i;
   pos2[index2[i]]=i;
  }

  for (i=0;i<ifin;i++){
   map[pos1[i]]=pos2[i];
   arrow[i]=scount[index1[i]];
  }
  
  count=0;
  for (i=0;i<ifin;i++){
   if (full[i]==1) continue;
   full[i]=1;
   j=map[i];
   snhit=0;
   if (arrow[i]>1){snhit=1;}
   sc=arrow[i]; //counts the number of spikes of the periodic orbit
   nper=0;
   seq[nper]=arrow[i];
   nper++;
   while (j !=i){
    full[j]=1;
    k=j;
    j=map[k];
    if (arrow[k]>1){snhit=1;}
    sc+=arrow[k];
    seq[nper]=arrow[k];
    nper++;
   }
   vper[count]=nper; //counts the period of the periodic orbit
   if (snhit==1){sc=-1;} 
   if (snhit==0 && nper>1){
    sc=ismaximin(nper,seq); // the vector seq has indexes from 0 to nper-1, thus length nper
   }
   rot[count]=rotnum(nper,seq);
   for(k=0;k<nper;k++){seq[k]=0;}
   vs[count]=sc; 
//the vs vector saves different information: if period is=1 tells you if it is z0 or z1, while if the period is greater than 1 and hits only s0 and s1 returns whether it is maximin (0 or 1), if an orbit hits sn, n>2 then -1
   count++;
  }

  z0=0;
  z1=0;
  count2=0;
  op01=0;
  opn=0;
  for (i=0;i<count;i++){
   if (vper[i]==1 && vs[i]==0){z0=1;}
   if (vper[i]==1 && vs[i]==1){z1=1;}
   if (vper[i]>1 && vs[i]!=-1){op01++; tper[count2]=vper[i]; ixper[count2]=vs[i];rot2[count2]=rot[i]; count2++;}
   if (vper[i]>1 && vs[i]==-1){opn++;}
  }  

  outfile<<A<<" "<<b<<" "<<suc<<" "<<z0<<" "<<z1<<" "<<op01<<" "<<opn;
  for(i=0;i<count2;i++){outfile<<" "<<tper[i]<<" "<<ixper[i]<<" "<<rot2[i]<<" ";}
  outfile<<endl;


  free(index1);
  free(index2);
  free(pos1);
  free(pos2);
  free(map);
  free(vper);
  free(vs);
  free(arrow);
  free(full);
  free(seq);
  free(imx);
  free(tper);
  free(ixper);
  free(rot);
  free(rot2);
} 

int partition( double * a, int * index, int l, int r) {
   int i, j, k;
   double pivot,aux;

   pivot = a[l];
	i = l; j = r+1;

	while( 1)
	{
		do ++i; while( a[i] <= pivot && i <= r );
		do --j; while( a[j] > pivot );
		if( i >= j ) break;
		aux = a[i]; a[i] = a[j]; a[j] = aux;
		k=index[i]; index[i]=index[j]; index[j]=k;
	}
	aux = a[l]; a[l] = a[j]; a[j] = aux;
	k=index[l]; index[l]=index[j]; index[j]=k;
   return j;
}


void quickSort( double * a, int * index, int l, int r)
{
	int j;

	if( l < r )
	{
   	// divide and conquer
        j = partition( a, index, l, r);
	   quickSort( a, index, l, j-1);
	   quickSort( a, index, j+1, r);
   }

}


void vfield(double t, double *x, int ndim, double *dx){

double sgm;
double Vr=0.1, tau=2.;  
double ti;

dx[0]=-x[0]+Vr+Iapp;
sgm=a+exp(b*(x[0]-c));
dx[1]=-(x[1]-sgm)/tau;
}

int invstroboscopic(double* x){
  double t=0.; 
  double *dx=new double[2];
  int i;
  double hini=T/100;
  double h=hini;
  double hmin=1.0E-4;
  double hmax=1;
  double hmax2=0.01;//This is for the Newton
  double tol=1.0e-13;
  int ndim=2;
  int out=0;
  int scount=0;
  double Newtol=1e-10;
  int maxiter=100;
  double hprev;
  
 ini_rk78(ndim);
 Iapp=0;
 
 h=-hini;
 while (t>-(T-dT) && out==0){
  if(x[0]>x[1]){
    end_rk78(ndim);
    delete[] dx;
    return -1;
  }
  rk78(&t,x,&h,tol,hmin,hmax,ndim,vfield);
 }

 h=-(t+T-dT);
 rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);

 Iapp=A;
 h=-hini;
 out=0;
 while (t>-T  && out==0){
  rk78(&t,x,&h,tol,hmin,hmax,ndim,vfield);
    if (x[0]<0){
      hprev=h;
      h=hini;
      while (x[0]<0){
	rk78(&t,x,&h,tol,hmin,hmax2,ndim,vfield);
      }
      i=0;
      while (fabs(x[0])>Newtol && i<maxiter){
       vfield(t,x,ndim,dx);
       h=-x[0]/dx[0];
       rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);
       i++;
      }
      if (i>=maxiter){
       printf("# Newton fails computing invstroboscopic\n");
       exit(1);
      }
      else{
	x[0]=0;
      }
      if (t<-T){
	out=1;
	h=-(t+T);
      }
      else{
	h=hprev;
	x[1]=x[1]-0.3;
	x[0]=x[1];
	scount++;
	if (x[1]<0 || x[0]<0){
	  end_rk78(ndim);
	  delete[] dx;
	  return -1;
	}
      }
    }
 }
   h=-(t+T);
   rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);
 end_rk78(ndim);
 delete[] dx;
 return scount;
}


int stroboscopic(double* x){
  double t=0.; 
  double *dx=new double[2];
  int i;
  double hini=0.05;
  double h=hini;
  double hmin=1.0E-4;
  double hmax=1;
  double hmax2=0.01;//This is for the Newton
  double tol=1.0e-13;
  double **aux;
  int ndim=2;
  int out=0;
  int scount=0;
  
 ini_rk78(ndim);
 Iapp=A;
 
 while (t<dT && out==0){
  rk78(&t,x,&h,tol,hmin,hmax,ndim,vfield);
  if (x[0] >= x[1]){
   //We go backwards for a better seed
     while (x[0]>x[1]){
      h=-hini;
      rk78(&t,x,&h,tol,hmin,hmax2,ndim,vfield);
     }
     while (fabs(x[0]-x[1])>1.0e-10){
      vfield(t,x,ndim,dx);
      h=-(x[0]-x[1])/(dx[0]-dx[1]);
      rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);  
      }
     if (t>dT){
       out=1;
       h=-(t-dT);
       rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);  
     }
     else{
     h=hini; 
     x[0]=0.;
     x[1]+=0.3;
     scount++;
     }
   }
 }

 if(out==0){
 h=-(t-dT);
 rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);
 }

 Iapp=0;
 h=hini;
 out=0;
 while (t<T && out==0){
  rk78(&t,x,&h,tol,hmin,hmax,ndim,vfield);
   if (x[0] >= x[1]){
     while (fabs(x[0]-x[1])>1.0e-10){
      vfield(t,x,ndim,dx);
      h=-(x[0]-x[1])/(dx[0]-dx[1]);
      rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);  
      }
     if (t>T){
       out=1;
       h=-(t-T);
       rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);  
     }
     else{
     h=hini; 
     x[0]=0.;
     x[1]+=0.3;
     scount++;
     }
   }
 }
 if (out==0){
 h=-(t-T);
 rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);
 }
 end_rk78(ndim);
 delete[] dx;
 return scount;
}

int ismaximin(int q, int *X){
  //This may fail of X is all 0's or 1's
  //X: sequence of 0's and 1's
  //q: length of the sequence
  int i,j,inc,newinc;
  int **siterates=new int*[q];//Here we store the iterates of X by the shift
  for (i=0;i<q;i++){
    siterates[i]=new int[q];
  }
  int *binseq=new int[q];
  int *perm=new int[q];
  
  binseq[0]=0;
  for (j=0;j<q;j++){
    siterates[0][j]=X[j];
    binseq[0]=binseq[0]+pow(2,q-1-j)*siterates[0][j];
  }
  for (i=1;i<q;i++){
    binseq[i]=0;
    for (j=0;j<q;j++){
      siterates[i][j]=siterates[i-1][(j+1)%q];
      binseq[i]=binseq[i]+ pow(2,q-1-j)*siterates[i][j];
    }
  }
  sortsequence(q,binseq,perm);
  inc=perm[1]-perm[0];
  if (inc<0){
    inc+=q;
  }
  newinc=inc;
  i=2;
  while (newinc==inc && i<q){
    newinc=perm[i]-perm[i-1];
    if (newinc<0){
      newinc+=q;
    }
    i++;
  }
  for (i=0;i<q;i++){
    delete[] siterates[i];
  }
  delete[] siterates;
  delete[] binseq;
  delete[] perm;

  if (i==q && newinc==inc){
    return 1;
  }
  else{
    return 0;
  }
}
void sortsequence(int q,int *seq,int *perm){
  //sorts a sequence of numbers, seq, of length q
  int i,j,tmp;
  for (i=0;i<q;i++){
    perm[i]=i;
  }
  for (i=0;i<q;i++){
    for (j=1;j<q-i;j++){
      if (seq[j-1]>seq[j]){
	tmp=seq[j-1];
	seq[j-1]=seq[j];
	seq[j]=tmp;
	tmp=perm[j-1];
	perm[j-1]=perm[j];
	perm[j]=tmp;
      }
    }
  }
}

double rotnum(int q, int *X){
  int i,sum;
  sum=0;
  for (i=0;i<q;i++){
    sum+=X[i];
  }
  return double(sum)/double(q);
}
