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
    
		int i,j; 
    int Nx, Nz;
	  float densy; 
    float *A=NULL, *AT=NULL;    
		MSR sparse=NULL;		
		//float *sparseAA=NULL;
		//int *sparseJA=NULL;


    if (!sf_histint(in,"n1",&Nz))   sf_error("no n1");  // n1 
    if (!sf_histint(in,"n2",&Nx))   sf_error("no n2");  // n2 
		if (!sf_getfloat("den",&densy))  densy=0.1;         // Density 
       
    A        = sf_floatalloc(Nx*Nz);
 		AT       = sf_floatalloc(Nx*Nz);
		//sparseAA = sf_floatalloc(Nx*Nz);
		//sparseJA = sf_intalloc(Nx*Nz);

    for(i=0;i<Nx*Nz;i++) A[i]= 0.0;
    
    
    sparserand(A, densy,Nx,Nz,10,5);
		
		for(i=0;i<Nx;i++)
		{
			for(j=0;j<Nz;j++) Id(AT,i,j)=Id(A,j,i);
		}
		   
		int nz; 
		nz = formatMSR(AT, Nx, Nz, sparse);

		sf_warning("los no cero %d",nz); 

		//sf_warning("%d",sparse);
		//for(i=0;i<nz;i++) sf_warning("%d \t %f", &sparse.JA[i], sparse.AA[i]);
		sf_floatwrite(A,Nx*Nz,out); 
   
    free(A);
		free(AT);
		//free(sparseAA);
		//free(sparseJA);

return 0;
}

