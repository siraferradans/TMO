#include "Imagen.h"
#include <stdio.h>
#include <fftw.h>
#include <rfftw.h>
#include <math.h>

fftw_complex * Imagen::fft(rfftwnd_plan p)
{
  int fil=dim[0];
  int col=dim[1];
  int M=p->n[0];  
  int N=p->n[1];  
  int escala=fil*col;
  fftw_complex * F=new fftw_complex[M*N];
  rfftwnd_one_real_to_complex(p, this->datos, F);
  return F;
}

Imagen & invfft(fftw_complex * F,rfftwnd_plan pinv, int fil, int col )
{

  int M=pinv->n[0];  
  int N=pinv->n[1];  
  int escala=fil*col;
  float invescala=1.0/escala;
  Imagen * res= new Imagen(fil,col,0.0);
  //  fftw_real d[fil*col];
  //  fftwnd_one_complex_to_real(pinv, F, &d[0]);
  rfftwnd_one_complex_to_real(pinv, F, res->datos);
  
  //  for(int i=0; i<escala; i++)
  //    res->datos[i]=d[i];
  
  (*res)*=invescala;
  return(*res);
}


fftw_complex * producto(fftw_complex * A, fftw_complex * B, int fil, int col)
{
  int largo=fil*(col/2+1);
  fftw_complex * C= new fftw_complex[largo];
  for (int i = 0; i < largo; i++)
      {
	C[i].re = (A[i].re * B[i].re - A[i].im * B[i].im);
	C[i].im = (A[i].re * B[i].im + A[i].im * B[i].re);
      }

  return(C);
}

fftw_complex * producto(fftw_complex * A, Imagen & k)
{
  int fil=k.fils();
  int col=k.cols();
  int largo=fil*(col/2+1);
  fftw_complex * C= new fftw_complex[largo];
  for (int i = 0; i < fil; i++)
    for (int j = 0; j < col/2+1; j++) 
      {
	int ij = i*(col/2+1) + j;   
	C[ij].re = A[ij].re * k(i,j) ;
	C[ij].im = A[ij].im * k(i,j);
	//	C[ij].re = A[ij].re * k(i,col/2+j-1) ;
	//	C[ij].im = A[ij].im * k(i,col/2+j-1);
      }

  return(C);
}


void Imagen::fftshift()
{
  int fil=dim[0];
  int col=dim[1];
  double tmp;

  for(int i=0;i<fil/2;i++)
    for(int j=0;j<col/2;j++)
      {
	tmp=(*this)(i,j);
	(*this)(i,j)=(*this)(i+fil/2,j+col/2);
	(*this)(i+fil/2,j+col/2)=tmp;
	tmp=(*this)(i+fil/2,j);
	(*this)(i+fil/2,j)=(*this)(i,j+col/2);
	(*this)(i,j+col/2)=tmp;
      }
}
	

Imagen & modulo(fftw_complex * F, int fil, int col )
{

  Imagen * res= new Imagen(fil,col,0.0);
  int largo=fil*col;
  double re, im;
  for(int i=0; i< largo; i++)
    {
      re=F[i].re;
      im=F[i].im; 
      (*res).datos[i]=sqrt(re*re+im*im)/largo;
    }
  return(*res);
}


Imagen & nucleo_gaussiano(int fil, int col, double sigma)
{
  Imagen * res= new Imagen(fil,col,0.0);
  double normaliza=1.0/(sqrt(2*M_PI)*sigma+1e-6);
  int i_;
  int mitfil=fil/2;
  int mitcol=col/2;
  for(int i=0;i<fil;i++)  
    for(int j=0;j<col;j++)
      (*res)(i,j)=normaliza*exp( -((i-mitfil)*(i-mitfil)+(j-mitcol)*(j-mitcol) )/(2*sigma*sigma) );

  return(*res);
}


Imagen & nucleo_exp(int fil, int col, double alfa)
{

  Imagen * res= new Imagen(fil,col,0.0);
  
  float mitfil=((float)fil)/2.0;
  float mitcol=((float)col)/2.0;
  
  float c1=mitfil*mitfil;
  float c2=mitcol*mitcol;
  double x;


  for(int i=0;i<fil;i++)  
    for(int j=0;j<col;j++)
      {
	//	x=sqrt( (i-mitfil)*(i-mitfil)+(j-mitcol)*(j-mitcol) );
	x=sqrt( (i-mitfil)*(i-mitfil)/c1+(j-mitcol)*(j-mitcol)/c2 );
	(*res)(i,j)=0.5*(exp(-alfa*x)+exp(-alfa*alfa*x*x));
      }

  return(*res);
}


Imagen & nucleo_lineal(int fil, int col, double alfa)
{

  Imagen * res= new Imagen(fil,col,0.0);
  
  float mitfil=((float)fil)/2.0;
  float mitcol=((float)col)/2.0;
  
  float c1=mitfil*mitfil;
  float c2=mitcol*mitcol;
  double d;

  //  c2*=0.5;

  for(int i=0;i<fil;i++)  
    for(int j=0;j<col;j++)
      {
	//d=sqrt( (i-mitfil)*(i-mitfil)+(j-mitcol)*(j-mitcol) );
	d=sqrt( (i-mitfil)*(i-mitfil)/c1+(j-mitcol)*(j-mitcol)/c2 );
	//	d*=alfa;
	(*res)(i,j)=1.0/(alfa+d);
      }

  return(*res);
}



/***
fftw_complex * Imagen2complex(Imagen k)
{
  int fil=k.fils();
  int col=k.cols();
  int largo=fil*(col/2+1);
  fftw_complex * C= new fftw_complex[largo];
  for (int i = 0; i < largo; i++)
      {
	C[i].re = (A[i].re * B[i].re - A[i].im * B[i].im);
	C[i].im = (A[i].re * B[i].im + A[i].im * B[i].re);
      }

  return(C);
}
***/
