#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MINDOUBLE 1.0e-2

#include <Imagen.h>
#include <Magick++.h>
#include <iostream> 
using namespace std; 
using namespace Magick; 

#include <fftw.h>
#include <rfftw.h>
#include <time.h>


void uso()
{
  fprintf(stderr,"Use: CM input output rho inv_alfa \n rho: -1.0 -> 5.0 ; the bigger is rho the darker the final image \n inv_alfa:(1->100) parameter that controls the amount of detail in the final image, side of the weighting function for the second (variational) step.\n");
  exit(1);
};


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */


#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }
//#define ELEM_SWAP(a,b) { register float t=a;a=b;b=t; }

double quick_select(double arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP





main(int argc,char **argv){
  
  if(argc<5)
    uso();


  argv++;
  char * entrada_char=*argv++;
  char * salida_char=*argv++;
  float rho=atof(*argv++);
  float invalfa = atof(*argv++);
  int colores=3;

  Imagen Im[3],LMS[3];
 
  leeHDR_color(Im,entrada_char);

  for(int color=0;color<colores;color++)
    {
      if ( Im[color].minval() < 0 ) 
	Im[color].recorta(Im[color].maxval(), 0.0);

      LMS[color]=Im[color];
    }


  int nfil=Im[0].dim[0];
  int ncol=Im[0].dim[1];
  int largo=nfil*ncol;
double alfa=min(nfil,ncol)/invalfa;
  Imagen g;
  g=nucleo_gaussiano(nfil,ncol,alfa);
  g.escala(1.0,0.0);
        
  //  Des-comentar para grabar una imagen con el nucleo usado
  //  escribe("nucleo.gif",g,g,g);
  
  g.fftshift();
  double suma=0;
  for(int i=0;i<largo;i++)
    suma+=g.datos[i];
  g*=(1.0/suma); //ahora la integral del nucleo vale 1
  

  rfftwnd_plan p, pinv;
  p=rfftw2d_create_plan(nfil, ncol, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  pinv = rfftw2d_create_plan(nfil, ncol, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  int fil=nfil;
  int col=ncol;
  largo=fil*col;
  
  int KBR;

  float lambda;
  float dt=0.2;//1e-1;
  double umbral_diferencia=dt/20.0;

  fftw_complex * G=g.fft(p);
  
  int iteracion=0;
  char indice[3];
  double diferencia0=2000.0;
  double diferencia=1000.0;      

  //  fprintf(stderr,"\n se crea el nucleo \n");


  Imagen RGB[3];
  Imagen RGBorig[3];
  float med[3];
  float valmax[3];
  float valmin[3];
  
  //fprintf(stderr,"\n nfil %d ncol %d \n",nfil,ncol);

  Imagen M(3,3), iM(3,3);

  //stockman + sharpe, apendice, lambda=575nm
  M(0,0)=5.8411e-01 ;
  M(0,1)=4.1581e-01 ;
  M(0,2)= 8.0518e-05;

  //stockman + sharpe, apendice, lambda=540nm
  M(1,0)=  4.8020e-01   ;
  M(1,1)=  5.1763e-01 ;
  M(1,2)=  2.1673e-03;

  //stockman + sharpe, apendice, lambda=440nm
  M(2,0)=  4.5961e-02 ;
  M(2,1)=  7.0566e-02 ; 
  M(2,2)=  8.8347e-01 ; 


  // iM*M=eye(3)
  iM(0,0)= 5.0420172 ;  
  iM(0,1)=  -4.0514967; 
  iM(0,2)= 0.0094795;

  iM(1,0)= -4.6778641; 
  iM(1,1)= 5.6913998  ;
  iM(1,2)= -0.0135357;

  iM(2,0)=  0.1113314;
  iM(2,1)=  -0.2438159;
  iM(2,2)= 1.1324845;

  // DESCOMENTAR LO SIGUIENTE SI SE QUIERE USAR RGB EN VEZ DE LMS
   for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      {
  	M(i,j)=0;
  	iM(i,j)=0;
      }
  for(int i=0;i<3;i++)
    {
      M(i,i)=1;
      iM(i,i)=1;
    }

  Imagen aux;
  float medianLMS[3],meanLMS[3],mu[3];
  float minLMS[3],maxLMS[3];
  float rangof;


  int cInit=clock();  //empieza a contar el tiempo para la parte 1

  // PASAMOS DE RGB ORIGINAL AL ESPACIO LMS DE LOS CONOS
  for(int color=0;color<3;color++)
    for(int d=0;d<largo;d++)
      LMS[color].datos[d]=M(color,0)*Im[0].datos[d]+M(color,1)*Im[1].datos[d]+M(color,2)*Im[2].datos[d];


  for(int k=0;k<3;k++)
    {
      LMS[k]+=1e-6;
      aux=LMS[k];
      medianLMS[k]=quick_select(aux.datos,largo);
      meanLMS[k]=aux.medval();
      minLMS[k]=aux.minval();
      maxLMS[k]=aux.maxval();
      if ( k ==-1)
	mu[k]=medianLMS[k];
      else 
          mu[k]=pow(meanLMS[k],0.5)*pow(medianLMS[k],0.5);

    }

  // LA CONSTANTE DE SEMISATURACIÓN SIGMA DEPENDE DE LA ILUMINACIÓN
  // DEL BACKGROUND. DATOS DE LA TABLA1 DE VALETON+VAN NORREN. 
  // TOMAMOS COMO REFERENCIA EL CANAL CON ILUMINACIÓN MÁXIMA.
  // SI NO HACEMOS ESTE PASO (SI SIMPLEMENTE IGUALAMOS SIGMA
  // A MU ORIGINAL, O MU+RHO) LOS COLORES NO QUEDAN BIEN
  float muMax=max(max(mu[0],mu[1]),mu[2]);
  for(int k=0;k<3;k++)
  {
	//DESPLAZAMOS log(mu) POR IGUAL PARA LOS 3 CANALES: rho ES EL UNICO
	//PARAMETRO DE NUESTRO ALGORITMO
	float x=log10(muMax)-log10(mu[k]);
	float y;

	y=log10(mu[k])+4-rho;

	float z= - 0.37*y + 1.9;
    float w=z;
    mu[k]*=pow(10,w);
  }
  

  Imagen WF[3],NR[3],mez[3];
  // VALETON + VAN NORREN: 
  // RANGO = 4 ÓRDENES; r ES LA MITAD DEL RANGO
  // n=0.74 : EXPONENTE EN FÓRMULA DE NAKA-RUSHTON
  float r=2;
  float n=0.74;
  for(int k=0;k<3;k++)
    {
      //ENCONTRAMOS CONSTANTES DE WEBER-FECHNER A PARTIR DE NAKA-RUSHTON
      float sigma=mu[k];

      float logs=log10(sigma);
      float I0=sigma/pow(10,1.2);

      // WYSZECKI-STILES, PÁG. 530: FECHNER FRACTION
      float K_=100.0/1.85;
      if(k==2)
 	K_= 100.0/8.7; 
      float Ir=pow(10,logs+r);
      float mKlogc=pow(Ir,n)/(pow(Ir,n)+pow(sigma,n))-K_*log10(Ir+I0);

    
      //APLICAMOS WEBER-FECHNER
      WF[k]=LMS[k];
      for(int d=0;d<largo;d++)
      	WF[k].datos[d]=K_*log10(LMS[k].datos[d]+I0) +mKlogc;

        //CALCULAMOS NAKA-RUSHTON 
      NR[k]=LMS[k];
       float sigma_n=pow(mu[k],n);
      for(int d=0;d<largo;d++)
	NR[k].datos[d]=pow(LMS[k].datos[d],n)/(pow(LMS[k].datos[d],n)+sigma_n);
  
      //MEZCLAMOS W-F Y N-R
      mez[k]=LMS[k];
      for(int d=0;d<largo;d++)
      	{
      	  float x=log10(LMS[k].datos[d]);
	  float In= pow(LMS[k].datos[d],n);
	  float srn=pow(pow(10,logs+r),n);
	  // ANTES DE logs+r APLICAMOS W-F, DESPUÉS N-R
      	  if(x<=logs+r)
	    mez[k].datos[d]=K_*log10( LMS[k].datos[d] + I0)+mKlogc;      	
	  else
	    mez[k].datos[d]= In/(In+sigma_n);

      	}

      float minmez=mez[k].minval();
      mez[k]+= -minmez;
      float escalamez=1.0/(mez[k].maxval()+1e-12);
      mez[k]*=escalamez;

    }

  //VOLVEMOS DEL ESPACIO LMS AL ESPACIO RGB
  for(int color=0;color<3;color++)
    for(int d=0;d<largo;d++)
      Im[color].datos[d]=iM(color,0)*mez[0].datos[d]+iM(color,1)*mez[1].datos[d]+iM(color,2)*mez[2].datos[d];

    char* aux_char = new char[1000];
    sprintf(aux_char,(char*)"./%sAtt_rho=%f.png",salida_char,rho);
    escribe(aux_char,Im[0],Im[1],Im[2]);
 
    int cFin=clock();
  fprintf(stderr,"\n Step 1 done in: %f \n", (float)(cFin - cInit)/(float)CLOCKS_PER_SEC);

  for(int color=0;color<colores;color++)
    {
      RGB[color]=Im[color];
      RGBorig[color]=RGB[color];
      med[color]=RGB[color].medval();
    }
  
  lambda=1;
  KBR=0;
  
  int c1,c0=clock();
  while(diferencia>umbral_diferencia)
    {
      iteracion++;
      diferencia0=diferencia;
      diferencia=0.0;
      for(int color=0;color<colores;color++)
	{
	  Imagen u0=RGB[color];
	  Imagen RGB0=RGB[color];
	  
	  Imagen u(u0);
	  Imagen u2(u);
	  u2*=u;
	  Imagen u3(u2);
	  u3*=u;
	  Imagen u4(u3);
	  u4*=u;
	  Imagen u5(u4);
	  u5*=u;
	  Imagen u6(u5);
	  u6*=u;
	  Imagen u7(u6);
	  u7*=u;
	  
	  fftw_complex * U=u0.fft(p);
	  fftw_complex * U2=u2.fft(p);
	  fftw_complex * U3=u3.fft(p);
	  fftw_complex * U4=u4.fft(p);
	  fftw_complex * U5=u5.fft(p);
	  fftw_complex * U6=u6.fft(p);
	  fftw_complex * U7=u7.fft(p);
	  
	  fftw_complex * UG=producto(U,G,fil,col);
	  fftw_complex * U2G=producto(U2,G,fil,col);
	  fftw_complex * U3G=producto(U3,G,fil,col);
	  fftw_complex * U4G=producto(U4,G,fil,col);
	  fftw_complex * U5G=producto(U5,G,fil,col);
	  fftw_complex * U6G=producto(U6,G,fil,col);
	  fftw_complex * U7G=producto(U7,G,fil,col);
	  
	  
	  
	  Imagen & iu=invfft(UG,pinv,fil, col);
	  Imagen & iu2=invfft(U2G,pinv,fil, col);
	  Imagen & iu3=invfft(U3G,pinv,fil, col);
	  Imagen & iu4=invfft(U4G,pinv,fil, col);
	  Imagen & iu5=invfft(U5G,pinv,fil, col);
	  Imagen & iu6=invfft(U6G,pinv,fil, col);
	  Imagen & iu7=invfft(U7G,pinv,fil, col);
	  
	  
	  delete[] U;
	  delete[] U2;
	  delete[] U3;
	  delete[] U4;
	  delete[] U5;
	  delete[] U6;
	  delete[] U7;
	  
	  delete[] UG;
	  delete[] U2G;
	  delete[] U3G;
	  delete[] U4G;
	  delete[] U5G;
	  delete[] U6G;
	  delete[] U7G;
	  
	  for(int i=0;i<largo;i++)
	    {
	      float Ip=u0.datos[i];
	      float I=iu.datos[i];
	      float I2=iu2.datos[i];
	      float I3=iu3.datos[i];
	      float I4=iu4.datos[i];
	      float I5=iu5.datos[i];
	      float I6=iu6.datos[i];
	      float I7=iu7.datos[i];
	      float gr1=Ip-I;
	      float gr2=Ip*Ip-2*Ip*I+I2;
	      float gr3=Ip*Ip*Ip-3*Ip*Ip*I+3*Ip*I2-I3;
	      float gr4=Ip*Ip*Ip*Ip+4*Ip*Ip*I2+I4-4*Ip*I3+2*Ip*Ip*(I2-2*Ip*I);
	      float gr5=Ip*gr4-Ip*Ip*Ip*Ip*I-4*Ip*Ip*I3-I5+4*Ip*I4-2*Ip*Ip*(I3-2*Ip*I2);
	      float gr6=Ip*Ip*gr4-Ip*Ip*Ip*Ip*Ip*I-4*Ip*Ip*Ip*I3-Ip*I5+4*Ip*Ip*I4-2*Ip*Ip*Ip*(I3-2*Ip*I2)-(I*Ip*gr4-Ip*Ip*Ip*Ip*I2-4*Ip*Ip*I4-I6+4*Ip*I5-2*Ip*Ip*(I4-2*Ip*I3)    );
	      float gr7=Ip*Ip*(Ip*(  Ip*Ip*Ip*Ip+4*Ip*Ip*I2+I4-4*Ip*I3+2*Ip*Ip*(I2-2*Ip*I)     )-Ip*Ip*Ip*Ip*I-4*Ip*Ip*I3-I5+4*Ip*I4-2*Ip*Ip*(I3-2*Ip*I2))-2*Ip*(Ip*(Ip*Ip*Ip*Ip*I+4*Ip*Ip*I3+I5-4*Ip*I4+2*Ip*Ip*(I3-2*Ip*I2)) -Ip*Ip*Ip*Ip*I2-4*Ip*Ip*I4-I6+4*Ip*I5-2*Ip*Ip*(I4-2*Ip*I3))+Ip*(Ip*Ip*Ip*Ip*I2+4*Ip*Ip*I4+I6-4*Ip*I5+2*Ip*Ip*(I4-2*Ip*I3)) -Ip*Ip*Ip*Ip*I3-4*Ip*Ip*I5-I7+4*Ip*I6-2*Ip*Ip*(I5-2*Ip*I4);
	   	      u0.datos[i]=( -7.7456e+00)*gr7+ (3.1255e-16)*gr6+(1.5836e+01)*gr5+(-1.8371e-15)*gr4+(-1.1013e+01)*gr3+(4.4531e-16)*gr2+(3.7891e+00)*gr1+ 1.2391e-15 ;//Arctan, pendiente 10
	    }
	  
	  u0*=0.5;
	  float alpha=1;
	  u0+= alpha*med[color];
	  Imagen I0(RGBorig[color]);
	  I0*=lambda;
	  u0+=I0;
	  u0*=dt;

	  Imagen ut(RGB[color]);
	  ut*= (1-dt*(lambda+1));
	  ut+=u0;
	  RGB[color]=ut;
	  
	  RGB[color].recorta(1,0);
	  
	  diferencia+=RGB0.MSE(RGB[color],1.0);
	  
	  iu2.~Imagen();
	  iu3.~Imagen();
	  iu4.~Imagen();
	  iu5.~Imagen();
	  iu6.~Imagen();
	  iu7.~Imagen();
	  iu.~Imagen();
	  
	}
      c1=clock();
      
    }
      
      
  fprintf(stderr,"\n Complete execution done in: %f \n", (float)(c1 - cInit)/(float)CLOCKS_PER_SEC);
  
  float med2[3];
  for(int c=0;c<3;c++)
    {
      RGB[c].recorta(1,0);
      RGB[c].escala(1,0);
    }
  
  char entradatmp_char[256];
  strcpy(entradatmp_char,salida_char);
 
  escribe(salida_char,RGB[0],RGB[1],RGB[2]);
  printf("Saved output image as:%s\n",salida_char);
  return 0;
  
  
}


double Imagen::MSE(Imagen & Im2, float escala)
{
  int largo=dim[0]*dim[1];
  double res=0.0;
  double tmp;
  float v=1.0/escala;
  float v2=1.0/1.0;
  for(int i=0;i<largo;i++)
    {
      tmp=(datos[i])*v-(Im2.datos[i])*v2;
      tmp=fabs(tmp);
      res+=tmp;
    }
  return(res/largo);

}


void Imagen::escala(float maxv, float minv)
{

  int largo=dim[0]*dim[1];
  float M=this->maxval();
  float m=this->minval();

  double R;
  double eps=1e-6;
  double s=(maxv-minv)/(M-m);
  for(int i=0;i<largo;i++)
    {
      R=datos[i];
      datos[i]=minv+s*(R-m);
    }
}

