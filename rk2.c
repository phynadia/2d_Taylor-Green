#include "ns2d.h"
#include <math.h>
#include<fftw3.h>
void rk2(double d,fftw_complex *uk1, fftw_complex *uk2,fftw_complex *Nkx1,fftw_complex *Nky2){
 
    	for(i=0; i<N; i++){  
    	   if(i <= N/2)
	       k1  = i;
	    else   
	       k1 = i - N;
	  for(j=0;j<(N/2+1);j++){
	     k2 = j;
	     ksqr = k1*k1+k2*k2;
	     m = j + (N/2+1)*i;
	     
	     uk1[m][0] = uk1[m][0]*(1-d*ksqr+d*d*ksqr*ksqr/2.0)+Nkx1[m][0]*(1-d*ksqr/2.0)*dt; 
	     uk1[m][1] = uk1[m][1]*(1-d*ksqr+d*d*ksqr*ksqr/2.0)+Nkx1[m][1]*(1-d*ksqr/2.0)*dt; 
	     
	     uk2[m][0] = uk2[m][0]*(1-d*ksqr+d*d*ksqr*ksqr/2.0)+Nky2[m][0]*(1-d*ksqr/2.0)*dt;
	     uk2[m][1] = uk2[m][1]*(1-d*ksqr+d*d*ksqr*ksqr/2.0)+Nky2[m][1]*(1-d*ksqr/2.0)*dt;  
	     
	}
     }
}
void euler(double d,fftw_complex *uk1, fftw_complex *uk2,fftw_complex *Nkx1,fftw_complex *Nky2){

        for(i=0; i<N; i++){  
    	   if(i <= N/2)
	       k1  = i;
	    else   
	       k1 = i - N;
	  for(j=0;j<(N/2+1);j++){
	     k2 = j;
	     ksqr = k1*k1+k2*k2;
	     m = j + (N/2+1)*i;
	     
	     uk1[m][0] = (uk1[m][0]*(1-d*ksqr)+Nkx1[m][0]*dt); 
	     uk1[m][1] = (uk1[m][1]*(1-d*ksqr)+Nkx1[m][1]*dt); 
	     
	     uk2[m][0] = (uk2[m][0]*(1-d*ksqr)+Nky2[m][0]*dt);
	     uk2[m][1] = (uk2[m][1]*(1-d*ksqr)+Nky2[m][1]*dt);
	     
	  }
	}     
}	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
    	
