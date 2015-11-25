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

typedef struct msr *MSR;
/*^*/
#endif

struct msr{
	float *AA;
	int *JA;
};


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


int formatMSR(float *A, int Nx, int Nz, MSR sparse)
/*< Sparse Matrix in MSR>*/
{
	int i,j;
	int cont=0;
	//sparse->AA = sf_floatalloc(Nx*Nz);
	//sparse->JA = sf_intalloc(Nx*Nz);

	for(i=0;i<Nx;i++)
	{
			sf_warning("%f",Id(A,i,i));
			sparse->AA[i]= Id(A,i,i);
			//sparse->AA[i]=1.0;
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





