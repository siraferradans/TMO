#ifndef IMAGEN_H
#define IMAGEN_H
#include <math.h>
#include <fftw.h>
#include <rfftw.h>

inline int signo(double a)
{ 
  if(a!=0.0)
    if(a>0)
      return(1);
    else
      return(-1);
  else 
    return(0);
} ;


inline double max(double a, double b)
{
  if(a<=b)
    return(b);
  else
    return(a);
};


inline double min(double a, double b)
{
  if(a>=b)
    return(b);
  else
    return(a);
};


inline float minmod(float x,float y)
{
  return(signo(x)*max(0,min(fabs(x),y*signo(x))));
};


inline double mediana(float a, float b, float c)
{
  if(a <= b)
    if(c<=a)
      return(a);
    else
      if(c<=b)
	return(c);
      else
	return(b);
  else // b < a
    if(c<=b)
      return(b);
    else
      if(c<=a)
	return(c);
      else
	return(a);
}


class Imagen
{
 public:
  int dim[2];
  double * datos; 
  Imagen();
  Imagen(Imagen & im2);
  Imagen(int fil, int col);
  Imagen(int fil, int col, float val);
  ~Imagen();
  void limpia();
  void  operator*=(float dt);
  void  operator+=(float val);
  void  operator+=(Imagen & im2);
  void  operator-=(Imagen & im2);
  void  operator*=(Imagen & im2);
  void  operator/=(Imagen & im2);
  void  operator=(Imagen & im2);
  int fils(){return dim[0];}
  int cols(){return dim[1];} 
  int area(){return dim[0]*dim[1];}
  inline double & operator()(int i, int j)
    {
      return(datos[i*dim[1]+j]);
    };
  float maxval();
  float minval();
  float medval();

  void recorta(float M,float m);
  Imagen & desplaza(int dim, int paso);
  Imagen & Laplaciano(float dx,float dy);
  Imagen & KModGrad();
  Imagen & KModGrad(Imagen & bandera);
  Imagen & anisotropica(Imagen & g);
  Imagen & minmod2(Imagen & c);
  Imagen * K(int tipo);
  Imagen * K();
  Imagen & D2uDuDu(Imagen & bandera);
  Imagen & K(Imagen & bandera,int tipo);
  void inpaint(Imagen & bandera, float dt);
  void inpaint(Imagen & bandera, float dt, Imagen & Lap0);
  void D3u(Imagen & bandera, float dt, int metodo, int cons,int upw,Imagen & cx_ant, Imagen & cy_ant);
  Imagen & derivada1(Imagen & bandera, int dir, int dif=0);
  Imagen & derivada1(int dir, int dif=0);
  Imagen & derivada2(int dir);
  Imagen & derivada2_O3(int dir);
  Imagen & Laplaciano(Imagen & bandera);
  Imagen & ModGrad();
  Imagen * ModGrad(int dummy);
  void umbral(double nivel, double M, double m);
  void copia_cercano(Imagen & b,int dis);
  void potencia(double p);
  void rotar90();
  void rotar_90();

  void escala(float p, float q);
fftw_complex * fft(rfftwnd_plan p);
void fftshift();
 double MSE(Imagen & Im2, float escala);
 Imagen & minmod(Imagen & c);
};
void escribeHDR_color(char * archima, Imagen * imagen);
void leeHDR_color(Imagen * imagen, char * archima) ;
Imagen & lee(char * archima, int campo);
Imagen * lee_ptr(char * archima, int campo);
void lee_ptr(char * archima, Imagen * Ip);
int compara_dims(Imagen & im1, Imagen & im2);
void escribe(char * archima, Imagen & imR, Imagen & imG, Imagen & imB);
void escribe(char * archima, Imagen & imR);


Imagen * optica(Imagen * Im,Imagen & g,Imagen & Zbuf,int W);
Imagen * optica2(Imagen * Im,Imagen & g,Imagen & Zbuf,int ancho_max);
Imagen * optica3(Imagen * Im,Imagen & Zbuf,int ancho_max,float zmax,float zmin,float foco);

Imagen & motionblur(Imagen & I,Imagen & d1,Imagen &d2);
Imagen * rgb2rfp(Imagen & imR,Imagen & imG,Imagen & imB);
Imagen * rfp2rgb(Imagen & imr,Imagen & imf,Imagen &imp);
Imagen * rgb2YCbCr(Imagen & imR,Imagen & imG,Imagen & imB);
Imagen * YCbCr2rgb(Imagen & imY,Imagen & imCb,Imagen & imCr);

Imagen * distancia_entera(Imagen & b,int ancho=21);
Imagen *  extension(Imagen & U1, Imagen & B);

float average(float * m, int n);
float torben(float *m, int n);
int itorben(float *m, int n);


Imagen & invfft(fftw_complex * F,rfftwnd_plan pinv, int fil, int col );
Imagen & nucleo_gaussiano(int fil, int col, double sigma);
fftw_complex * producto(fftw_complex * A, fftw_complex * B, int fil, int col);
fftw_complex * producto(fftw_complex * A, Imagen & k);
Imagen & modulo(fftw_complex * F, int fil, int col );



#endif

