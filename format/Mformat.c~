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
    
		
		int i,j,k; 
    int Nx, Nz;
	  float densy; 
    float *A=NULL, *AT=NULL, *V=NULL, *R=NULL, *R2=NULL, *R3=NULL, *R4=NULL;    
	  MSR sparseMSR;	
	  CSR sparseCSR;	
	  COO sparseCOO;
	
    if (!sf_getint("Nx",&Nx))  Nx=500; 
    if (!sf_getint("Nz",&Nz))  Nz=500; 
		if (!sf_getfloat("den",&densy))  densy=0.1;         // Density 
       
    A            = sf_floatalloc(Nx*Nz);
    V            = sf_floatalloc(Nx);
    
    R            = sf_floatalloc(Nx);
    R2           = sf_floatalloc(Nx);
    R3           = sf_floatalloc(Nx);
    R4           = sf_floatalloc(Nx);
        
 		AT           = sf_floatalloc(Nx*Nz);
 		sparseMSR.AA = sf_floatalloc(Nx*Nz);
		sparseMSR.JA = sf_intalloc(Nx*Nz);
		
		sparseCSR.AA = sf_floatalloc(Nx*Nz);
		sparseCSR.IA = sf_intalloc(Nx*Nz);
		sparseCSR.JA = sf_intalloc(Nx*Nz);
		
		
		sparseCOO.AA = sf_floatalloc(Nx*Nz);
		sparseCOO.JR = sf_intalloc(Nx*Nz);
		sparseCOO.JC = sf_intalloc(Nx*Nz);
		
  
    for(i=0;i<Nx*Nz;i++) A[i]= 0.0;
    int nz, nz2, nz3; 
		struct timespec ts1, ts2, ts3, ts4, ts5; 
		for(i=0;i<Nx;i++) V[i]=1.0;
  
    for(k=0;k<9;k++)
    {  
    
      densy=0.1+k*0.1;
    
      sparserand(A,densy,Nx,Nz,10,5);
      				
		  for(i=0;i<Nx;i++)
		  {
			  for(j=0;j<Nz;j++) Id(AT,i,j)=Id(A,j,i);
		  }                              
    
			
		  nz  = formatMSR(AT, Nx, Nz, &sparseMSR);		
		  nz2 = formatCSR(AT, Nx, Nz, &sparseCSR);	
		  nz3 = formatCOO(AT, Nx, Nz, &sparseCOO); 
		
		  clock_gettime( CLOCK_REALTIME, &ts1 );
		
		  mulAVec(AT, V, R, Nx, Nz);
		  clock_gettime( CLOCK_REALTIME, &ts2 );
		  
		  mulCooVec(sparseCOO, V, R2, Nx, nz);
		  clock_gettime( CLOCK_REALTIME, &ts3 );
		
		  mulCsrVec(sparseCSR, V, R3, Nx, nz2);
		  clock_gettime( CLOCK_REALTIME, &ts4 );
		  
		  mulMsrVec(sparseMSR, V, R4, Nx, nz3);
		  clock_gettime( CLOCK_REALTIME, &ts5);		 
		
		  float seg  =0.0;
		  float seg2 =0.0;
		  float seg3 =0.0;
		  float seg4 =0.0;
		
		  //seg= (double)(sec2 - sec1)/CLOCKS_PER_SEC;
		  seg= (float) ( 1.0*(1.0*ts2.tv_nsec - ts1.tv_nsec*1.0)*1e-9 + 1.0*ts2.tv_sec - 1.0*ts1.tv_sec );
		  seg2= (float) ( 1.0*(1.0*ts3.tv_nsec - ts2.tv_nsec*1.0)*1e-9 + 1.0*ts3.tv_sec - 1.0*ts2.tv_sec );
		
		  seg3= (float) ( 1.0*(1.0*ts4.tv_nsec - ts3.tv_nsec*1.0)*1e-9 + 1.0*ts4.tv_sec - 1.0*ts3.tv_sec );
		
		  seg4= (float) ( 1.0*(1.0*ts5.tv_nsec - ts4.tv_nsec*1.0)*1e-9 + 1.0*ts5.tv_sec - 1.0*ts4.tv_sec );
		
		  sf_warning("%.0f %e %e %e %e", densy*100,seg,seg2,seg3,seg4);
		  
		  nz=0; nz2=0; nz3=0;
		  for(i=0;i<Nx*Nz;i++) A[i]= 0.0;
		}
		
		sf_floatwrite(A,Nx,out); 
   
      
		  //free(sparseAA);
		  //free(sparseJA);

return 0;
}

