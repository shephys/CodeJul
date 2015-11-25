#include <rsf.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_cblas.h>

#include "funaux.h"

#define Id(v,i,j) v[(i)*(Nz)+(j)] 
/*^*/

#ifndef _funaux_h
#define _funaux_h

typedef struct msr MSR;
/*^*/

typedef struct csr CSR;
/*^*/

typedef struct coo COO;
/*^*/

struct msr{
	float *AA;
	int *JA;
};
/*^*/

struct csr{
	float *AA;
	int *JA;
	int *IA;
};
/*^*/

struct coo{
	float *AA;
	int *JR;
	int *JC;
};
/*^*/
#endif

void sparserand (float *A, float densy, int Nx, int Nz, int a, int b)
/*< Sparse Matrix >*/
{
  int i,j,n,x,y,N;
  int seed;
  seed = time(NULL);
  init_genrand((unsigned long) seed);

  N= Nx*Nz;
  
  n= (int) N*densy; 
  
  
  
  if(densy<1.0 && n>=Nx)
  {
    
    for(i=0;i<Nx;i++)
    {
      Id(A,i,i) = genrand_real1()*(b-a+1.0)+a;
    }
    
    j=Nz;
    
    while(j<n)
    {
    
      x =  genrand_real1()*Nx;
      y =  genrand_real1()*Nz;
      
      if(!Id(A,x,y))
      {
        Id(A,x,y)=genrand_real1()*(b-a+1.0)+a;
        j++;
      }
     
    }
  }

}


int formatMSR(float *A, int Nx, int Nz, MSR *sparse)
/*< Sparse Matrix in MSR>*/
{
	int i,j;
	int cont=0;

	for(i=0;i<Nx;i++)
	{
			
			sparse->AA[i]= Id(A,i,i);
	}	
 	

	sparse->AA[Nx]=-1;
	cont=Nx+1;
	for(i=0; i<Nz ; i++)
	{
		sparse->JA[i]=cont;
		for(j=0; j<Nx; j++)
		{
			if((Id(A,i,j)!=0)&&(i!=j))
			{
				sparse->AA[cont]=Id(A,i,j);
				sparse->JA[cont]=j;
				cont++;
			}
		}
  }

	sparse->JA[Nx]=cont;  
	return(cont);
}

void mulMsrVec(MSR sparse, float *V, float *R, int Nx, int nz)
/*< Sparse Matrix MSR Vector>*/
{
	int i,j;
	for(i=0; i < Nx; i++)
	{
		R[i]= V[i]*sparse.AA[i];
	}
	j=0;
	for(i=Nx+1; i < nz; i++)
	{
    while(sparse.JA[j+1]==i) j++;
		R[j] += (V[sparse.JA[i]] * sparse.AA[i]);
	}
}


int formatCSR(float *A, int Nx, int Nz, CSR *sparse)
/*< Sparse Matrix CSR>*/
{
  int i,j;
  int cont=0;
  for(i=0; i<Nx; i++)
  {
    sparse->IA[i]=cont;
    for(j=0; j<Nx; j++)
    {
      if(Id(A,i,j)!=0)
      {
        sparse->AA[cont]=Id(A,i,j);
        sparse->JA[cont]=j;
        cont++;
      }
    }
  }
  sparse->IA[Nx]=cont;
  return cont;
}

void mulCsrVec(CSR sparse, float *V, float *R, int Nx, int nz)
/*< Sparse Matrix CSR Vector>*/
{
  int i,j;
  
  for(i=0; i < Nx; i++)
  {
    R[i]=0;
  }
  j=0;
  
  for(i=0; i < nz; i++)
  {
    while(sparse.IA[j+1]==i)j++;
    R[j] += (V[sparse.JA[i]] * sparse.AA[i]);
  }
}

int formatCOO(float *A, int Nx, int Nz, COO *sparse)
/*< Sparse Matrix COO>*/
{
  int i,j;
  int cont=0;
  
  for(i=0; i < Nx; i++)
  {
    for(j=0; j < Nz; j++)
    {
      if(Id(A,i,j)!=0)
      {
        sparse->AA[cont]=Id(A,i,j);
        sparse->JR[cont]=i;
        sparse->JC[cont]=j;
        cont++;
      }
    }
  }
 return cont;
}

void mulCooVec(COO sparse, float *V, float *R, int Nx, int nz)
/*< Sparse Matrix COO Vector>*/
{
  int i;
  
  for(i=0; i < Nx; i++)
  {
    R[i]=0;
  }
   
  for(i=0; i < nz; i++)
  {
    R[sparse.JR[i]] += (V[sparse.JC[i]] * sparse.AA[i]);
  }
}



void mulAVec(float *A, float *V, float *R, int Nx, int Nz)
/*<Matrix Vector>*/
{
  int i, j;
  
  for(i=0;i<Nx;i++)
	{
		for(j=0;j<Nz;j++) R[i] +=Id(A,i,j)*V[j];
	}
}





