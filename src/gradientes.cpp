#include <Imagen.h>
#include <stdio.h>
#include <math.h>


Imagen & Imagen::minmod(Imagen & c)
{
  int fil=dim[0];
  int col=dim[1];
  Imagen * res = new Imagen(c);
  double Ixat, Ixad, Iyat, Iyad;
  double Ixatm, Ixadm, Iyatm, Iyadm;
  double IxatM, IxadM, IyatM, IyadM;
  for(int i=1; i<fil-1; i++)
    for(int j=1; j<col-1;j++)
      {
	Ixat=(*this)(i,j)-(*this)(i-1,j);//backward derivative
	Ixad=(*this)(i+1,j)-(*this)(i,j);//forward derivative
	Iyat=(*this)(i,j)-(*this)(i,j-1);
	Iyad=(*this)(i,j+1)-(*this)(i,j);
	
	Ixatm=min(Ixat,0);
	IxatM=max(Ixat,0);
	Ixadm=min(Ixad,0);
	IxadM=max(Ixad,0);
	Iyatm=min(Iyat,0);
	IyatM=max(Iyat,0);
	Iyadm=min(Iyad,0);
	IyadM=max(Iyad,0);
	
	if (c(i,j)>0)
	  (*res)(i,j)*=sqrt( Ixatm*Ixatm + IxadM*IxadM + Iyatm*Iyatm + IyadM*IyadM );
	else
	  (*res)(i,j)*=sqrt( IxatM*IxatM + Ixadm*Ixadm + IyatM*IyatM + Iyadm*Iyadm );
      }
  
  return(*res);
}



Imagen & Imagen::derivada1(int dir, int dif)
{
  /* dir=0 para filas, dir=1 para columnas*/
  /* dif=0 para centrada, 1 forward, -1 backward */

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
 
  switch(dir)
    {
    case 0:
      switch(dif)
	{
	case 0:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      (*res)(i,j)=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
	  break;
	case 1:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      (*res)(i,j)=(*this)(i+1,j) - (*this)(i,j) ;
	  break;
	case 2:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      (*res)(i,j)=(*this)(i,j) - (*this)(i-1,j) ;
	  break;
	default:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      (*res)(i,j)=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
	  break;
	}

      for(int j=0;j<col;j++)
	(*res)(0,j)=( (*this)(1,j) - (*this)(0,j) );
      
      for(int j=0;j<col;j++)
	(*res)(fil-1,j)=( (*this)(fil-1,j) - (*this)(fil-2,j) );
      break;

    case 1:
      switch(dif)
	{
	case 0:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      (*res)(i,j)=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
	  break;
	case 1:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      (*res)(i,j)=(*this)(i,j+1) - (*this)(i,j) ;
	  break;
	case 2:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      (*res)(i,j)=(*this)(i,j) - (*this)(i,j-1) ;
	  break;
	default:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      (*res)(i,j)=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
	  break;
	}

      
      for(int i=0;i<fil;i++)
	(*res)(i,0)=( (*this)(i,1) - (*this)(i,0) );
      
      for(int i=0;i<fil;i++)
	(*res)(i,col-1)=( (*this)(i,col-1) - (*this)(i,col-2) );
      
      break;

    default:
      break;
    }


  return((*res));
}


Imagen & Imagen::derivada1(Imagen & bandera,int dir, int dif)
{
  /* dir=0 para filas, dir=1 para columnas*/
  /* dif=0 para centrada, 1 forward, -1 backward */

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
 
  switch(dir)
    {
    case 0:
      switch(dif)
	{
	case 0:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      if(bandera(i,j)>0)
		(*res)(i,j)=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
	  break;
	case 1:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=(*this)(i+1,j) - (*this)(i,j) ;
	  break;
	case 2:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=(*this)(i,j) - (*this)(i-1,j) ;
	  break;
	default:
	  for(int i=1;i<fil-1;i++)
	    for(int j=0;j<col;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
	  break;
	}

      for(int j=0;j<col;j++)
	(*res)(0,j)=( (*this)(1,j) - (*this)(0,j) );
      
      for(int j=0;j<col;j++)
	(*res)(fil-1,j)=( (*this)(fil-1,j) - (*this)(fil-2,j) );
      break;

    case 1:
      switch(dif)
	{
	case 0:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
	  break;
	case 1:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=(*this)(i,j+1) - (*this)(i,j) ;
	  break;
	case 2:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=(*this)(i,j) - (*this)(i,j-1) ;
	  break;
	default:
	  for(int i=0;i<fil;i++)
	    for(int j=1;j<col-1;j++)
	      if(bandera(i,j)>0)
	      (*res)(i,j)=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
	  break;
	}

      
      for(int i=0;i<fil;i++)
	(*res)(i,0)=( (*this)(i,1) - (*this)(i,0) );
      
      for(int i=0;i<fil;i++)
	(*res)(i,col-1)=( (*this)(i,col-1) - (*this)(i,col-2) );
      
      break;

    default:
      break;
    }


  return((*res));
}



Imagen & Imagen::derivada2(int dir)
{
  /* dir=0 para filas, dir=1 para columnas*/

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
 
  switch(dir)
    {
    case 0:
      for(int i=1;i<fil-1;i++)
	for(int j=0;j<col;j++)
	  (*res)(i,j)= (*this)(i+1,j) + (*this)(i-1,j) -2*(*this)(i,j) ;
      
      for(int j=0;j<col;j++)
	(*res)(0,j)=(*res)(1,j);
      
      for(int j=0;j<col;j++)
	(*res)(fil-1,j)=(*res)(fil-2,j);
      break;

    case 1:
      for(int i=0;i<fil;i++)
	for(int j=1;j<col-1;j++)
	  (*res)(i,j)= (*this)(i,j+1) + (*this)(i,j-1) -2*(*this)(i,j) ;
      
      for(int i=0;i<fil;i++)
	(*res)(i,0)=(*res)(i,1);
      
      for(int i=0;i<fil;i++)
	(*res)(i,col-1)=(*res)(i,col-2);
      
      break;

    default:
      break;
    }


  return((*res));
}



Imagen & Imagen::Laplaciano(float dx,float dy)
{

  int fil=dim[0];  
  int col=dim[1];
  float dx2=dx*dx;
  float dy2=dy*dy;

  Imagen * res= new Imagen(fil,col,0.0);
 
      for(int i=1;i<fil-1;i++)
	for(int j=1;j<col-1;j++)
	  (*res)(i,j)= ( (*this)(i+1,j) + (*this)(i-1,j) -2*(*this)(i,j))/dx2 +( (*this)(i,j+1) + (*this)(i,j-1) -2*(*this)(i,j))/dy2  ;
      
 
      for(int j=1;j<col-1;j++)
	{
	  (*res)(0,j)=((*this)(1,j) -1*(*this)(0,j))/dx2 + ( (*this)(0,j+1) + (*this)(0,j-1) -2*(*this)(0,j))/dy2 ;
	  (*res)(fil-1,j)=((*this)(fil-2,j) -1*(*this)(fil-1,j))/dx2 + ((*this)(fil-1,j+1) + (*this)(fil-1,j-1) -2*(*this)(fil-1,j))/dy2;
	}

      for(int i=1;i<fil-1;i++)
	{
	  (*res)(i,0)= ((*this)(i+1,0) + (*this)(i-1,0) -2*(*this)(i,0))/dx2 + ((*this)(i,1) - 1*(*this)(i,0))/dy2;
	  (*res)(i,col-1)= ((*this)(i+1,col-1) + (*this)(i-1,col-1) -2*(*this)(i,col-1))/dx2 + ( (*this)(i,col-2)-2*(*this)(i,col-1))/dy2 ;
	}     

      (*res)(0,0)=0.5*((*res)(0,1)+(*res)(1,0)); 
      (*res)(0,col-1)=0.5*((*res)(0,col-2)+(*res)(1,col-1)); 
      (*res)(fil-1,col-1)=0.5*((*res)(fil-1,col-2)+(*res)(fil-2,col-1)); 
      (*res)(fil-1,0)=0.5*((*res)(fil-1,1)+(*res)(fil-2,0)); 

  return((*res));
}


Imagen & Imagen::Laplaciano(Imagen & bandera)
{

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
 
      for(int i=1;i<fil-1;i++)
	for(int j=1;j<col-1;j++)
	  if(bandera(i,j)>0)
	  (*res)(i,j)= (*this)(i+1,j) + (*this)(i-1,j) -4*(*this)(i,j) + (*this)(i,j+1) + (*this)(i,j-1) ;
      
  return((*res));
}



Imagen & Imagen::minmod2(Imagen & c)
{
  int fil=dim[0];
  int col=dim[1];
  Imagen * res = new Imagen(c);
  double Ixat, Ixad, Iyat, Iyad;
  double Ixatm, Ixadm, Iyatm, Iyadm;
  double IxatM, IxadM, IyatM, IyadM;
  for(int i=1; i<fil-1; i++)
    for(int j=1; j<col-1;j++)
      {
	Ixat=(*this)(i,j)-(*this)(i-1,j);//backward derivative
	Ixad=(*this)(i+1,j)-(*this)(i,j);//forward derivative
	Iyat=(*this)(i,j)-(*this)(i,j-1);
	Iyad=(*this)(i,j+1)-(*this)(i,j);
	
	Ixatm=min(Ixat,0);
	IxatM=max(Ixat,0);
	Ixadm=min(Ixad,0);
	IxadM=max(Ixad,0);
	Iyatm=min(Iyat,0);
	IyatM=max(Iyat,0);
	Iyadm=min(Iyad,0);
	IyadM=max(Iyad,0);
	
	if (c(i,j)>0)
	  (*res)(i,j)*=sqrt( Ixatm*Ixatm + IxadM*IxadM + Iyatm*Iyatm + IyadM*IyadM );
	else
	  (*res)(i,j)*=sqrt( IxatM*IxatM + Ixadm*Ixadm + IyatM*IyatM + Iyadm*Iyadm );
      }
  
  return(*res);
}

Imagen & Imagen::KModGrad()
{
  int fil=dim[0];  
  int col=dim[1];
  double MINDOUBLE=1e-5;

  Imagen * res= new Imagen(fil,col,0.0);
  float fxx,fyy,fx,fy,fxy;

  for(int i=1;i<fil-1;i++)
    for(int j=1;j<col-1;j++)
      {
	fx=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
	fy=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
	fxx= (*this)(i+1,j) + (*this)(i-1,j) -2*(*this)(i,j) ;
	fyy= (*this)(i,j+1) + (*this)(i,j-1) -2*(*this)(i,j) ;
	fxy= ((*this)(i+1,j+1) - (*this)(i-1,j+1) - (*this)(i+1,j-1) + (*this)(i-1,j-1))/4.0;
		(*res)(i,j)=(fyy*fx*fx+fxx*fy*fy-2*fx*fy*fxy)/(fx*fx+fy*fy+MINDOUBLE);
      }


  // i=0 // 
    for(int j=1;j<col-1;j++)
      {
	fx=( (*this)(1,j) - (*this)(0,j) );
	fy=0.5*( (*this)(0,j+1) - (*this)(0,j-1) );
	fxx= (*this)(1,j) - (*this)(0,j) ;
	fyy= (*this)(0,j+1) + (*this)(0,j-1) -2*(*this)(0,j) ;
	fxy= ((*this)(1,j+1) - (*this)(0,j+1) - (*this)(1,j-1) + (*this)(0,j-1))/4.0;
		(*res)(0,j)=(fyy*fx*fx+fxx*fy*fy-2*fx*fy*fxy)/(fx*fx+fy*fy+MINDOUBLE);
      }

  // i=fil-1 // 
    for(int j=1;j<col-1;j++)
      {
	fx=( (*this)(fil-1,j) - (*this)(fil-2,j) );
	fy=0.5*( (*this)(fil-1,j+1) - (*this)(fil-1,j-1) );
	fxx= (*this)(fil-2,j) -(*this)(fil-1,j) ;
	fyy= (*this)(fil-1,j+1) + (*this)(fil-1,j-1) -2*(*this)(fil-1,j) ;
	fxy=( (*this)(fil-1,j+1) - (*this)(fil-2,j+1) - (*this)(fil-1,j-1) + (*this)(fil-2,j-1))/4.0;
		(*res)(fil-1,j)=(fyy*fx*fx+fxx*fy*fy-2*fx*fy*fxy)/(fx*fx+fy*fy+MINDOUBLE);
      }

    //j=0//
  for(int i=1;i<fil-1;i++)
      {
	fx=0.5*( (*this)(i+1,0) - (*this)(i-1,0) );
	fy=( (*this)(i,1) - (*this)(i,0) );
	fxx= (*this)(i+1,0) + (*this)(i-1,0) -2*(*this)(i,0) ;
	fyy= (*this)(i,1) -(*this)(i,0) ;
	fxy=( (*this)(i+1,1) - (*this)(i-1,1) - (*this)(i+1,0) + (*this)(i-1,0))/4.0;
		(*res)(i,0)=(fyy*fx*fx+fxx*fy*fy-2*fx*fy*fxy)/(fx*fx+fy*fy+MINDOUBLE);
      }

  // j=col-1//
  for(int i=1;i<fil-1;i++)
      {
	fx=0.5*( (*this)(i+1,col-1) - (*this)(i-1,col-1) );
	fy=( (*this)(i,col-1) - (*this)(i,col-2) );
	fxx= (*this)(i+1,col-1) + (*this)(i-1,col-1) -2*(*this)(i,col-1) ;
	fyy= (*this)(i,col-1)  -(*this)(i,col-2) ;
	fxy= ((*this)(i+1,col-1) - (*this)(i-1,col-1) - (*this)(i+1,col-2) + (*this)(i-1,col-2))/4.0;
		(*res)(i,col-1)=(fyy*fx*fx+fxx*fy*fy-2*fx*fy*fxy)/(fx*fx+fy*fy+MINDOUBLE);
      }
  
  (*res)(0,0)=0.5*((*res)(0,1)+(*res)(1,0)); 
  (*res)(0,col-1)=0.5*((*res)(0,col-2)+(*res)(1,col-1)); 
  (*res)(fil-1,col-1)=0.5*((*res)(fil-1,col-2)+(*res)(fil-2,col-1)); 
  (*res)(fil-1,0)=0.5*((*res)(fil-1,1)+(*res)(fil-2,0)); 


  return(*res);
}


Imagen & Imagen::ModGrad()
{

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
  Imagen ux(fil,col,0.0);
  Imagen uy(fil,col,0.0);
  double uxij, uyij;

      for(int i=1;i<fil-1;i++)
	for(int j=0;j<col;j++)
	  ux(i,j)=0.5*( (*this)(i+1,j) - (*this)(i-1,j) );
      
      for(int j=0;j<col;j++)
	{
	  ux(0,j)=( (*this)(1,j) - (*this)(0,j) );
	  ux(fil-1,j)=( (*this)(fil-1,j) - (*this)(fil-2,j) );
	}

      for(int i=0;i<fil;i++)
	for(int j=1;j<col-1;j++)
	  uy(i,j)=0.5*( (*this)(i,j+1) - (*this)(i,j-1) );
      
      for(int i=0;i<fil;i++)
	{
	  uy(i,0)=( (*this)(i,1) - (*this)(i,0) );
	  uy(i,col-1)=( (*this)(i,col-1) - (*this)(i,col-2) );
	}

      for(int i=0;i<fil;i++)
	for(int j=0;j<col;j++)
	  {
	    uxij=ux(i,j);
	    uyij=uy(i,j);
	    (*res)(i,j)=sqrt(uxij*uxij+uyij*uyij);
	  }


  return((*res));
}


Imagen & Imagen::derivada2_O3(int dir)
{
  /* O3 quiere decir que es una aproximacion de 3er orden */
  /* dir=0 para filas, dir=1 para columnas*/

  int fil=dim[0];  
  int col=dim[1];

  Imagen * res= new Imagen(fil,col,0.0);
 
  switch(dir)
    {
    case 0:
      for(int i=2;i<fil-2;i++)
	for(int j=0;j<col;j++)
	  (*res)(i,j)= (4.0/3.0)*((*this)(i+1,j)+(*this)(i-1,j)-2*(*this)(i,j)) -(1.0/12.0)*((*this)(i+2,j) -2*(*this)(i,j) + (*this)(i-2,j) ) ;
      
      for(int j=0;j<col;j++)
	{
	  (*res)(1,j)=(*this)(2,j) + (*this)(0,j) -2*(*this)(1,j) ;
	  (*res)(0,j)=(*res)(1,j);
	  (*res)(fil-2,j)=(*this)(fil-1,j) + (*this)(fil-3,j) -2*(*this)(fil-2,j) ;
	  (*res)(fil-1,j)=(*res)(fil-2,j);
	}
      break;

    case 1:
      for(int i=0;i<fil;i++)
	for(int j=2;j<col-2;j++)
	  (*res)(i,j)= (4.0/3.0)*((*this)(i,j+1)+(*this)(i,j-1)-2*(*this)(i,j)) -(1.0/12.0)*((*this)(i,j+2) -2*(*this)(i,j) + (*this)(i,j-2) ) ;
      
      for(int i=0;i<fil;i++)
	{
	  (*res)(i,1)=(*this)(i,2) + (*this)(i,0) -2*(*this)(i,1) ;
	  (*res)(i,0)=(*res)(i,1);
	  (*res)(i,col-2)=(*this)(i,col-1) + (*this)(i,col-3) -2*(*this)(i,col-2) ;
	  (*res)(i,col-1)=(*res)(i,col-2);
	}

      break;

    default:
      break;
    }


  return((*res));
}


Imagen & Imagen::K(Imagen & bandera, int tipo)
{
  int fil=dim[0];  
  int col=dim[1];

  Imagen vx(fil,col,0.0);
  Imagen vy(fil,col,0.0);
  Imagen * res= new Imagen(fil,col,0.0);
  double modgrad_u, uxij,uyij;
  double MINDOUBLE=1e-5;
  switch(tipo)
    {
    case 0:
      // 00
      for(int i=0;i<fil-1;i++)
	for(int j=0;j<col-1;j++)
	    {
	      uxij= (*this)(i+1,j) - (*this)(i,j) ;
	      uyij= (*this)(i,j+1) - (*this)(i,j) ;
	      modgrad_u=sqrt(uxij*uxij+uyij*uyij+MINDOUBLE);
	      vx(i,j)=uxij/modgrad_u;
	      vy(i,j)=uyij/modgrad_u;
	    }

      for(int i=1;i<fil;i++)
	for(int j=1;j<col;j++)
	  if(bandera(i,j)>0)
	    (*res)(i,j)=vx(i,j)-vx(i-1,j)+vy(i,j)-vy(i,j-1);
      break;
       
    case 1:
      // 01
      for(int i=0;i<fil-1;i++)
	for(int j=1;j<col;j++)
	  {
	    uxij= (*this)(i+1,j) - (*this)(i,j) ;
	    uyij= (*this)(i,j) - (*this)(i,j-1) ;
	    modgrad_u=sqrt(uxij*uxij+uyij*uyij+MINDOUBLE);
	    vx(i,j)=uxij/modgrad_u;
	    vy(i,j)=uyij/modgrad_u;
	  }
      for(int i=1;i<fil;i++)
	for(int j=0;j<col-1;j++)
	  if(bandera(i,j)>0)
	    (*res)(i,j)=vx(i,j)-vx(i-1,j)+vy(i,j+1)-vy(i,j);
      break;


    case 2:
      // 10
      for(int i=1;i<fil;i++)
	for(int j=0;j<col-1;j++)
	  {
	    uxij= (*this)(i,j) - (*this)(i-1,j) ;
	    uyij= (*this)(i,j+1) - (*this)(i,j) ;
	    modgrad_u=sqrt(uxij*uxij+uyij*uyij+MINDOUBLE);
	    vx(i,j)=uxij/modgrad_u;
	    vy(i,j)=uyij/modgrad_u;
	  }
      for(int i=0;i<fil-1;i++)
	for(int j=1;j<col;j++)
	  if(bandera(i,j)>0)
	    (*res)(i,j)=vx(i+1,j)-vx(i,j)+vy(i,j)-vy(i,j-1);
      break;
      
      
      
    case 3:
      // 11
      for(int i=1;i<fil;i++)
	for(int j=1;j<col;j++)
	  {
	    uxij= (*this)(i,j) - (*this)(i-1,j) ;
	    uyij= (*this)(i,j) - (*this)(i,j-1) ;
	    modgrad_u=sqrt(uxij*uxij+uyij*uyij+MINDOUBLE);
	    vx(i,j)=uxij/modgrad_u;
	    vy(i,j)=uyij/modgrad_u;
	  }
      for(int i=0;i<fil-1;i++)
	for(int j=0;j<col-1;j++)
	  if(bandera(i,j)>0)
	    (*res)(i,j)=vx(i+1,j)-vx(i,j)+vy(i,j+1)-vy(i,j);
      break;
      
    default:
      break;
    }


  return(*res);
};


