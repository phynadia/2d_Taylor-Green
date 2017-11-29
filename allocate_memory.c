#include "ns2d.h"
#include<fftw3.h>
void allocate(){
        x = fftw_malloc(N*sizeof(double));
        y = fftw_malloc(N*sizeof(double));
	u1 = fftw_malloc(N*N*sizeof(double));
	u2 = fftw_malloc(N*N*sizeof(double));
	u1exact = fftw_malloc(N*N*sizeof(double));
	u2exact = fftw_malloc(N*N*sizeof(double));
	err1 = fftw_malloc(N*N*sizeof(double));
	err2 = fftw_malloc(N*N*sizeof(double));
	N1 = fftw_malloc(N*N*sizeof(double));
	N2 = fftw_malloc(N*N*sizeof(double));
	omegaz = fftw_malloc(N*N*sizeof(double));
	
	uk1 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	uk2 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	omegak3 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	Nk1 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	Nk2 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	Nkx1 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	Nky2 = fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
	
}	
