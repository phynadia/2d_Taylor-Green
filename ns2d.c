/*....................................................................................................
.  									                             .
.  C PROGRAM FOR NUMERICAL SOLUTION OF 2D NAVIER-STOKE EQUATION WITH TAYLOR-GREEN INITIAL CONDITIONS .
.                                      BY PSEUDO-SPECTRAL METHOD.                                    .                                    
.  We have used FFTW library for Fast Fourier Transform and Rk-2 for time marching.                  .
.												     .
. 									NADIA BIHARI PADHAN          .
.....................................................................................................*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include "ns2d.h"
#include "file_open.c"


int main(int argc, char *argv[]){
	N=256;
	tmax = 0.1;
        nu = 1.6e-3;
        dt = 0.001;

	system ("mkdir -p data");
	scale = 1.0/(double)(N*N);
        d = nu*dt;
        
        if (argc > 1)
        N = atoi(argv[1]);
        if (argc > 2)
        tmax = atof(argv[2]);
        
        allocate();
        file_open();
//..........Creating plan for FFT.................................................... 
        fftw_plan p1,p2,p3,p4,p5,p6,p7;	
	p1 = fftw_plan_dft_r2c_2d(N,N,u1,uk1,FFTW_ESTIMATE);
	p2 = fftw_plan_dft_r2c_2d(N,N,u2,uk2,FFTW_ESTIMATE);
	p3 = fftw_plan_dft_r2c_2d(N,N,N1,Nk1,FFTW_ESTIMATE);
	p4 = fftw_plan_dft_r2c_2d(N,N,N2,Nk2,FFTW_ESTIMATE);
	p5 = fftw_plan_dft_c2r_2d(N,N,uk1,u1,FFTW_ESTIMATE);
	p6 = fftw_plan_dft_c2r_2d(N,N,uk2,u2,FFTW_ESTIMATE);
	p7 = fftw_plan_dft_c2r_2d(N,N,omegak3,omegaz,FFTW_ESTIMATE);
//..........Setting the range of x and y.............................................	
	for(i=0;i<N;i++){
	    j = i - N/2 ;
	    x[i]=2*M_PI*j/(double)N;
	    y[i]=2*M_PI*j/(double)N;
	   
	 }
//..........Using Taylor-Green Initial conditions....................................	 
	for(i=0;i<N;i++){ 
	   for(j=0;j<N;j++){ 
	    m = j+N*i;
	    u1[m]= sin(x[i])*cos(y[j]); 
	    u2[m]= -cos(x[i])*sin(y[j]); 	    
	   } 
	}
//..........Taking fourier transform of u1 and u2 and save into uk1 and uk2...........	
	fftw_execute(p1);
	fftw_execute(p2);
//....................................................................................	         
      for(t=0.0;t<=tmax;t = t+dt){
	        fftw_execute(p1);
		fftw_execute(p2);
	      
	    for(i=0;i<N;i++){
	        if(i <= N/2)
	           k1  = i;
	        else   
	           k1 = i - N;
	      for(j=0;j<(N/2+1);j++){ 
	          k2 = j;
		  m = j+(N/2+1)*i;
//................Calculating omega in k-space.......................................		
	          omegak3[m][0] = -(k1*uk2[m][1]-k2*uk1[m][1]);
	          omegak3[m][1] =   k1*uk2[m][0]-k2*uk1[m][0]; 
//...................................................................................	        
	      }
	  } 
//..............Taking omega into real space.........................................	  
	    fftw_execute(p7);
//...................................................................................	    
	   for(i=0;i<N;i++){ 
	      for(j=0;j<N;j++){
	         m = j+N*i;
	         
	         omegaz[m] = scale*omegaz[m];
//...............Evaluating non-linear terms in real space..........................	         
	         N1[m] =  u2[m]*omegaz[m];
	         N2[m] = -u1[m]*omegaz[m];
	        // printf("%g\n",N2[m]);
//..................................................................................	     
	     }
	  }
//..............Taking FFT of N1 and N2 and storing in Nk1 and Nk2..................	        
	  fftw_execute(p3);
	  fftw_execute(p4);  
//.................................................................................	  
	  for(i=0;i<N;i++){ 
	     if(i <= N/2)
	       k1  = i;
	     else   
	       k1 = i - N;
	     for(j=0;j<(N/2+1);j++){ 
	         k2 = j;
		 m = j+(N/2+1)*i; 
		 ksqr = k1*k1+k2*k2;       
	         if(ksqr > 1e-5){
	           p11 = 1.0 - k1*k1/(double)ksqr;
	           p22 = 1.0 - k2*k2/(double)ksqr;
	           p12 = -k1*k2/(double)ksqr;
	           }else{
	                 p11 = 1.0/2.0;
	                 p22 = p11;
	                 p12 = -1.0/2.0;
	                }       
	        Nkx1[m][0] = (Nk1[m][0]*p11 + Nk2[m][0]*p12);
	        Nkx1[m][1] = (Nk1[m][1]*p11 + Nk2[m][1]*p12);
	        Nky2[m][0] = (Nk2[m][0]*p22 + Nk1[m][0]*p12);
	        Nky2[m][1] = (Nk2[m][1]*p22 + Nk1[m][1]*p12);
	       // printf("%g",Nkx1[m][0]);
	        
	    }
	  } 
	     
	   rk2(d,uk1,uk2,Nkx1,Nky2);  //Time marching
	   
//.........IFFT of uk1 and uk2.................................................	   
	   fftw_execute(p5);
	   fftw_execute(p6);
//.............................................................................	   
	   E = 0.0;
	   E1 = 0.0;
	   for(i=0;i<N;i++){ 
	    for(j=0;j<N;j++){
	      m = j+N*i;
	      
	      u1exact[m] =  sin(x[i])*cos(y[j])*exp(-2*nu*(t+dt));
	      u2exact[m] = -cos(x[i])*sin(y[j])*exp(-2*nu*(t+dt));  
	      u1[m]=scale*u1[m];
	      u2[m]=scale*u2[m];
	      E += (u1[m]*u1[m]+u2[m]*u2[m])*scale;
	      E1 += u1exact[m]*u1exact[m]+u2exact[m]*u2exact[m];
	      
	   }
	 } 
	  fprintf(f1,"%g \t  %g \t %g\n",t,E,E1); 
      }   //.......End of the time loop........................................
//.............................................................................	
	  for(i=0;i<N;i++){ 
	    for(j=0;j<N;j++){
	     m = j+N*i;
	     
	     err1[m] = fabs(u1[m]-u1exact[m]);
	     err2[m] = fabs(u2[m]-u2exact[m]);
	     
	     fprintf(f2,"%g ",u1[m]);
	     fprintf(f3,"%g ",u2[m]);
	     fprintf(f4,"%g ",err1[m]);
	     fprintf(f5,"%g ",err2[m]);
	    } 
	     fprintf(f2,"\n");
	     fprintf(f3,"\n");
	     fprintf(f4,"\n");
	     fprintf(f5,"\n");
        }
	return 0;
}	
//..............................................................................	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	       
