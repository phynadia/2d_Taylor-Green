#ifndef _NS2D_H
#define _NS2D_H
#include<fftw3.h>

void allocate();
void rk2(double d,fftw_complex *uk1,fftw_complex *uk2,fftw_complex *Nkx1,fftw_complex *Nky2);
void euler(double d,fftw_complex *uk1,fftw_complex *uk2,fftw_complex *Nkx1,fftw_complex *Nky2);
void file_open();

FILE *f1,*f2,*f3,*f4,*f5;
int N,i,j,m,k1,k2,ksqr;
double t;
double p11,p22,p12;
double scale;
double E1,E;
double tmax,nu,dt;
double d;
double *x,*y,*u1,*err1,*err2,*u2,*N1,*N2,*u1exact,*u2exact,*omegaz;
fftw_complex *uk1,*uk2,*Nk1,*Nk2,*Nkx1,*Nky2,*omegak3;

#endif
