#include <rsf.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_cblas.h>

#include "funaux.h"

int main (int argc, char* argv[])
{
    /* init RSF */  
    sf_init (argc,argv);
    
    sf_file in, out;
    in  = sf_input("in");    
    out = sf_output ("out");
    
		int i; 
    int Nx, Nz;
	  float densy; 
    float *A=NULL;    


    if (!sf_histint(in,"n1",&Nz))   sf_error("no n1");  // n1 
    if (!sf_histint(in,"n2",&Nx))   sf_error("no n2");  // n2 
		if (!sf_getfloat("den",&densy))  densy=0.1;         // Density 
       
    A = sf_floatalloc(Nx*Nz);
    for(i=0;i<Nx*Nz;i++) A[i]= 0.0;
    
    
    sparserand(A, densy,Nx,Nz,10,50);   
    sf_floatwrite(A,Nx*Nz,out); 
   
    free(A);

return 0;
}

