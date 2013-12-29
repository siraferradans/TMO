
#ifndef INOUT_H
#define INOUT_H

#include "Imagen.h"
#include "Connections.h"



//Funciones inout (inout.cpp)
void lee(Imagen * imagen, char * archima, int campo);
void escribe(char * archima, Imagen & imR, Imagen & imG, Imagen & imB);
bool lee_carpetaLDR( char * input_char, map<double, Imagen*>& Input_Images, char * format);

void escribe(char * archima, Imagen & imR);
void escribeTxt(char * archima, Imagen & imR); //se guardan valores en formato txt, para debugging


//HDR 
Imagen* readHDRfolder( char * folder_char, map<double,Imagen*>& Input_Images);
bool saveHDRfolder( char * folder_char, map<double,Imagen*>& Input_Images, Imagen* g);
void escribeHDR_color( char * archima, Imagen * imagen) ;
void escribeHDR_BW(char * archima, Imagen * imagen);
void leeHDR_color(Imagen * imagen, char * archima) ;


//Connections
void escribeTxt(char * archima, Connections imR);
void escribeTxtConnections(char * archima, Connections imR);
void escribeTxtExaminator(char * archima,char * archima2, Connections imR);


#endif
