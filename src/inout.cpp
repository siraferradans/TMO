#include <Imagen.h>
#include <stdio.h>
#include <stdlib.h>

#include <cstring>

#include <Magick++.h>
#include <iostream> 
using namespace std; 
using namespace Magick; 

#define maxpaleta 65535.0
#define maxmipaleta 255.0
/* El tama#o maximo de paleta es 65536, especificado por QuantumLeap al compilar ImageMagick*/


void escribeHDR_color(char * archima, Imagen * imagen)
{
  /*
    Esta funcion escribe los numeros correctos de las HDR
    asume que la imagen tiene 3 planos de color
  */

  FILE* fid = fopen( archima, "w" );
  // read header

  fprintf(fid, "PF\n%d\t%d\n-1\n", imagen->cols(),imagen->fils() );

  float aux =0.0;
  int filas=imagen->fils();
  int columnas=imagen->cols();

	   
  for (int i = filas-1; i >= 0; i--)
      {
	for (int j = 0; j < columnas; j++)
	  {
	    
	    for (int color= 0; color < 3; color++)
	      {
		aux = (float)imagen[color](i,j);//.datos[i*filas+j];
	
		fwrite(&aux,sizeof(float),1,fid);	
	      }
	  }
      }

 fclose(fid);

}


void leeHDR_color(Imagen * imagen, char * archima) 
{
  /*
    Esta funcion coge los numeros correctos de las HDR, como doubles
  */

  FILE* fid = fopen( archima, "r" );
  // read header
  char* id=new char[4];
  int a= fscanf(fid, "%s\n", id);
  if ( a == EOF) cout<<"leeHDR_color:Error reading "<<archima<<"\n";
  
  int filas, columnas;
  float scale; 
  a=fscanf(fid, "%d %d\n", &columnas,&filas );
  if ( a == EOF) cout<<"leeHDR_color:Error reading "<<archima<<"\n";
  a=fscanf(fid, "%f\n", &scale );
  if ( a == EOF) cout<<"leeHDR_color:Error reading "<<archima<<"\n";
 
  scale = abs((int) scale );

  float aux=0.0;
 
  Imagen dummy(filas,columnas);
  if (strcmp(id,"PF")==0)
    {
      
      for (int color =0; color < 3; color++)
	imagen[color]=dummy;
      //       	imagen[color].Inicializa(filas,columnas,0.0);
      
      for (int i = 1; i <= filas; i++)
	{
	  for (int j = 0; j < columnas; j++)
	    {
	      for (int color= 0; color < 3; color++)
		{
		  a=fread(&aux,sizeof(float),1, fid);
		  
		  imagen[color].datos[(filas-i)*columnas+j] =(double)aux;
		  
		}
	      
	    }
	}
      
    }
  else if (strcmp(id,"Pf")==0)
    {
      cout<<"leeHDR_color:Imagen de un solo plano de color: NO TESTEADO!\n";
      //imagen->Inicializa(filas,columnas,0.0);
      for (int i = 1; i <= filas; i++)
	{
	  for (int j = 0; j < columnas; j++)
	    {
	      
	      a=fread(&aux,sizeof(float),1, fid);
	      imagen->datos[(filas-i)*columnas+j] =(double)aux;
	      
	    }
	}
      
    }
  else
    cout<<"leeHDR_color:Error! tipo de archivo (PF) no valido\n";
  
  fclose(fid);
  
}



Imagen  & lee(char * archima, int campo)
{
  Image image(archima);

  int alto=image.rows();
  int ancho=image.columns();  
  
  Imagen * im=new Imagen(alto,ancho,0.0);
  //  double c=maxmipaleta/maxpaleta;
  double c=1;

   for(int i=0;i<alto;i++)
    for(int j=0;j<ancho;j++)
      switch(campo)
	{
	case 0:
	  (*im)(i,j)=c*(image.pixelColor(j,i)).redQuantum();
	  break;
	case 1:
	  (*im)(i,j)=c*(image.pixelColor(j,i)).greenQuantum();
	  break;
	case 2:
	  (*im)(i,j)=c*(image.pixelColor(j,i)).blueQuantum();
	  break;
	default:
	  fprintf(stderr,"\n Error en inout.cpp al leer archivo \n");
	  break;
	} 
  
   return(*im);
};

int compara_dims(Imagen & im1, Imagen & im2)
{
  if(im1.fils()!=im2.fils() | im1.cols()!=im2.cols())
    return 0;
  return 1;
};


// void escribe(char * archima, Imagen & imR, Imagen & imG, Imagen & imB)
// {
//   int ancho=imR.cols();
//   int alto=imR.fils();
//   if(compara_dims(imR,imG)!=1 | compara_dims(imG,imB)!=1 )
//     {
//       fprintf(stderr,"\n Error al escribir en inout.cpp: imagenes R,G,B con distintas dimensiones. \n");
//       return;
//     }
//       fprintf(stderr,"\n Sin error al escribir en inout.cpp. \n");

//   Geometry geom(ancho,alto);
//       fprintf(stderr,"\n Sin error al escribir en inout.cpp. \n");
//   Image im(geom);
//   // float MR=imR.maxval();
//   // float MG=imG.maxval();
//   // float MB=imB.maxval();
//   // int maxuno=1;
//   // if(MR>1 || MG>1 || MB>1)
//   //   maxuno=0;

//       fprintf(stderr,"\n Sin error al escribir en inout.cpp. \n");


//   for(int i=0;i<alto;i++)
//     for(int j=0;j<ancho;j++)
//       {
// 	float R=imR(i,j);
// 	float G=imG(i,j);
// 	float B=imB(i,j);
// 	//	if(maxuno)
// 	if(1)
// 	  {
// 	    // R*=255;
// 	    // G*=255;
// 	    // B*=255;
// 	    R*=65535;
// 	    G*=65535;
// 	    B*=65535;
// 	  }
// 	Color col((int)R,(int)G,(int)B);
// 	//Color col(R,G,B);
//         im.pixelColor(j,i,col);
//       }

//   im.write(archima);
//   return;
// };


void escribe(char * archima, Imagen & imR, Imagen & imG, Imagen & imB)
{
  int ancho=imR.cols();
  int alto=imR.fils();
  if(compara_dims(imR,imG)!=1 | compara_dims(imG,imB)!=1 )
    {
      fprintf(stderr,"\n Error al escribir en inout.cpp: imagenes R,G,B con distintas dimensiones. \n");
      return;
    }

  Geometry geom(ancho,alto);
  Image im(geom,"white");
  float MR=imR.maxval();
  float MG=imG.maxval();
  float MB=imB.maxval();
  int maxuno=1;
  if(MR>1 || MG>1 || MB>1)
    maxuno=0;
  for(int i=0;i<alto;i++)
    for(int j=0;j<ancho;j++)
      {
	float R=imR(i,j);
	float G=imG(i,j);
	float B=imB(i,j);
	//	if(maxuno)
	if(1)
	  {
	    R*=65535;
	    G*=65535;
	    B*=65535;
	  }
	Color col((int)R,(int)G,(int)B);
        im.pixelColor(j,i,col);
      }

  im.write(archima);
  return;
};


void escribe(char * archima, Imagen & imR)
{
  int ancho=imR.cols();
  int alto=imR.fils();

  Geometry geom(ancho,alto);
  Image im(geom,"white");

//   float M=imR.maxval();
//   float m=imR.minval();
//   float c=maxpaleta;
//   if(m>=0.0 && m<=maxmipaleta && M>=0.0 && M<=maxmipaleta)
//     {
//       m=0.0;
//       M=1.0;
//       c=maxpaleta/maxmipaleta;
//     }

  for(int i=0;i<alto;i++)
    for(int j=0;j<ancho;j++)
      {
	float R=imR(i,j);
	//	R=c*(R-m)/(M-m);
	//R/=255.0;
	Color col((int)R,(int)R,(int)R);
        im.pixelColor(j,i,col);
      }

  im.write(archima);
  return;
}

