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



