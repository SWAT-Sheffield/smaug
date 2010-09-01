#include "cudapars.h"
#include "paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "step.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////



__device__ __host__
int encode_dc2 (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int fencode_dc2 (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
}

__device__ __host__
real evalgrad_dc2(real fi, real fim1, real fip2, real fim2,struct params *p,int dir)
{
 //real valgrad_dc2;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx[0])))*(fi-fim1)-(1.0/(12.0*(p->dx[0])))*(fip2-fim2);
   //return((1.0/(2.0*(p->dx[0])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[0])))*(fi-fim1)):((1.0/(12.0*(p->dx[0])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(fi-fim1)-(1.0/(12.0*(p->dx[1])))*(fip2-fim2);
     // return((2.0/(1.0*(p->dx[1])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[1])))*(fi-fim1)):((1.0/(12.0*(p->dx[1])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }

 return -1;
}


__device__ __host__
real grad_dc2(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_dc2;

  if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx[0])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx[0])))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
//return((1.0/(2.0*(p->dx[0])))*(wmod[fencode_dc2(p,i+1,j,field)]-wmod[fencode_dc2(p,i-1,j,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_dc2(p,i+1,j,field)]-NVAR*wmod[fencode_dc2(p,i-1,j,field)]+wmod[fencode_dc2(p,i-2,j,field)]-wmod[fencode_dc2(p,i+2,j,field)])/6.0):wmod[fencode_dc2(p,i+1,j,field)]-wmod[fencode_dc2(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dx[1])))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
// return((1.0/(2.0*(p->dx[1])))*(wmod[fencode_dc2(p,i,j+1,field)]-wmod[fencode_dc2(p,i,j-1,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_dc2(p,i,j+1,field)]-NVAR*wmod[fencode_dc2(p,i,j-1,field)]+wmod[fencode_dc2(p,i,j-2,field)]-wmod[fencode_dc2(p,i,j+2,field)])/6.0):wmod[fencode_dc2(p,i,j+1,field)]-wmod[fencode_dc2(p,i,j-1,field)])/(2.0*(p->dx[1]))    );
}
 return 0;
}



__device__ __host__
real ddotcurrentb (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  //real ddc=0;

  real fi, fim1, fip2=0, fim2=0;
  real ddc1,ddc2;
  real ddcx,ddcy;

  switch(direction)
  {
	case 0:
	       fi=w[fencode_dc2(p,ix+1,iy,mom1)]*w[fencode_dc2(p,ix+1,iy,b1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,mom1)]*w[fencode_dc2(p,ix-1,iy,b1)]/w[fencode_dc2(p,ix-1,iy,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,mom1)]*w[fencode_dc2(p,ix+2,iy,b1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,mom1)]*w[fencode_dc2(p,ix-2,iy,b1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);
	       fi=w[fencode_dc2(p,ix,iy+1,mom1)]*w[fencode_dc2(p,ix,iy+1,b2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,mom1)]*w[fencode_dc2(p,ix,iy-1,b2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,mom1)]*w[fencode_dc2(p,ix,iy+2,b2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,mom1)]*w[fencode_dc2(p,ix,iy-2,b2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);
        break;
	case 1:
	       fi=w[fencode_dc2(p,ix+1,iy,mom2)]*w[fencode_dc2(p,ix+1,iy,b1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,mom2)]*w[fencode_dc2(p,ix-1,iy,b1)]/w[fencode_dc2(p,ix-1,iy,rho)];

               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,mom2)]*w[fencode_dc2(p,ix+2,iy,b1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,mom2)]*w[fencode_dc2(p,ix-2,iy,b1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);

	       fi=w[fencode_dc2(p,ix,iy+1,mom2)]*w[fencode_dc2(p,ix,iy+1,b2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,mom2)]*w[fencode_dc2(p,ix,iy-1,b2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,mom2)]*w[fencode_dc2(p,ix,iy+2,b2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,mom2)]*w[fencode_dc2(p,ix,iy-2,b2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);
        break;
	case 2:
	       fi=w[fencode_dc2(p,ix+1,iy,mom3)]*w[fencode_dc2(p,ix+1,iy,b1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,mom3)]*w[fencode_dc2(p,ix-1,iy,b1)]/w[fencode_dc2(p,ix-1,iy,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,mom3)]*w[fencode_dc2(p,ix+2,iy,b1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,mom3)]*w[fencode_dc2(p,ix-2,iy,b1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);
	       fi=w[fencode_dc2(p,ix,iy+1,mom3)]*w[fencode_dc2(p,ix,iy+1,b2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,mom3)]*w[fencode_dc2(p,ix,iy-1,b2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,mom3)]*w[fencode_dc2(p,ix,iy+2,b2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,mom3)]*w[fencode_dc2(p,ix,iy-2,b2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);

        break;
  }
  ddc1=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);



  switch(direction)
  {
	case 0:
	       fi=w[fencode_dc2(p,ix+1,iy,b1)]*w[fencode_dc2(p,ix+1,iy,mom1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,b1)]*w[fencode_dc2(p,ix-1,iy,mom1)]/w[fencode_dc2(p,ix-1,iy,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,b1)]*w[fencode_dc2(p,ix+2,iy,mom1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,b1)]*w[fencode_dc2(p,ix-2,iy,mom1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);
	       fi=w[fencode_dc2(p,ix,iy+1,b1)]*w[fencode_dc2(p,ix,iy+1,mom2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,b1)]*w[fencode_dc2(p,ix,iy-1,mom2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,b1)]*w[fencode_dc2(p,ix,iy+2,mom2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,b1)]*w[fencode_dc2(p,ix,iy-2,mom2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);
        break;
	case 1:
	       fi=w[fencode_dc2(p,ix+1,iy,b2)]*w[fencode_dc2(p,ix+1,iy,mom1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,b2)]*w[fencode_dc2(p,ix-1,iy,mom1)]/w[fencode_dc2(p,ix-1,iy,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,b2)]*w[fencode_dc2(p,ix+2,iy,mom1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,b2)]*w[fencode_dc2(p,ix-2,iy,mom1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);
	       fi=w[fencode_dc2(p,ix,iy+1,b2)]*w[fencode_dc2(p,ix,iy+1,mom2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,b2)]*w[fencode_dc2(p,ix,iy-1,mom2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,b2)]*w[fencode_dc2(p,ix,iy+2,mom2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,b2)]*w[fencode_dc2(p,ix,iy-2,mom2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);
        break;
	case 2:
	       fi=w[fencode_dc2(p,ix+1,iy,b3)]*w[fencode_dc2(p,ix+1,iy,mom1)]/w[fencode_dc2(p,ix+1,iy,rho)];
	       fim1=w[fencode_dc2(p,ix-1,iy,b3)]*w[fencode_dc2(p,ix-1,iy,mom1)]/w[fencode_dc2(p,ix-1,iy,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix+2,iy,b3)]*w[fencode_dc2(p,ix+2,iy,mom1)]/w[fencode_dc2(p,ix+2,iy,rho)];
	       fim2=w[fencode_dc2(p,ix-2,iy,b3)]*w[fencode_dc2(p,ix-2,iy,mom1)]/w[fencode_dc2(p,ix-2,iy,rho)];
               }
	       ddcx=evalgrad_dc2(fi,fim1,fip2,fim2,p,0);
	       fi=w[fencode_dc2(p,ix,iy+1,b3)]*w[fencode_dc2(p,ix,iy+1,mom2)]/w[fencode_dc2(p,ix,iy+1,rho)];
	       fim1=w[fencode_dc2(p,ix,iy-1,b3)]*w[fencode_dc2(p,ix,iy-1,mom2)]/w[fencode_dc2(p,ix,iy-1,rho)];
               if(p->sodifon)
               {
	       fip2=w[fencode_dc2(p,ix,iy+2,b3)]*w[fencode_dc2(p,ix,iy+2,mom2)]/w[fencode_dc2(p,ix,iy+2,rho)];
	       fim2=w[fencode_dc2(p,ix,iy-2,b3)]*w[fencode_dc2(p,ix,iy-2,mom2)]/w[fencode_dc2(p,ix,iy-2,rho)];
               }
	       ddcy=evalgrad_dc2(fi,fim1,fip2,fim2,p,1);
        break;
  }
  ddc2=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);
  return(ddc1-ddc2);
  //return -ddc2;

}

__device__ __host__
real ddotcurrentenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

 // real ddc=0;
  real dd1,dd2,dd3;
 
  real ddcx,ddcy;
  real fi, fim1,fip2=0, fim2=0;
  //real dpi, dpim1;//, dpip2=0, dpim2=0;


  //int field=energy;

  //fi=w[fencode_dc2(p,ix+1,iy,energy)]*w[fencode_dc2(p,ix+1,iy,mom1)]/w[fencode_dc2(p,ix,iy,rho)];
  //fim1=w[fencode_dc2(p,ix-1,iy,energy)]*w[fencode_dc2(p,ix-1,iy,mom1)]/w[fencode_dc2(p,ix-1,iy,rho)];
if(p->sodifon==1)
{
  fip2=w[fencode_dc2(p,ix+2,iy,energy)]*w[fencode_dc2(p,ix+2,iy,mom1)]/w[fencode_dc2(p,ix+2,iy,rho)];
  fim2=w[fencode_dc2(p,ix-2,iy,energy)]*w[fencode_dc2(p,ix-2,iy,mom1)]/w[fencode_dc2(p,ix-2,iy,rho)];
}
 // ddcx=evalgrad_dc2(fi,fim1,0,0,p,0);
  ddcx=evalgrad_dc2(w[fencode_dc2(p,ix+1,iy,energy)]*w[fencode_dc2(p,ix+1,iy,mom1)]/w[fencode_dc2(p,ix+1,iy,rho)],w[fencode_dc2(p,ix-1,iy,energy)]*w[fencode_dc2(p,ix-1,iy,mom1)]/w[fencode_dc2(p,ix-1,iy,rho)],fip2,fim2,p,0);

 // fi=w[fencode_dc2(p,ix,iy+1,energy)]*w[fencode_dc2(p,ix,iy+1,mom2)]/w[fencode_dc2(p,ix,iy+1,rho)];
 // fim1=w[fencode_dc2(p,ix,iy-1,energy)]*w[fencode_dc2(p,ix,iy-1,mom2)]/w[fencode_dc2(p,ix,iy-1,rho)];
if(p->sodifon==1)
{
  fip2=w[fencode_dc2(p,ix,iy+2,energy)]*w[fencode_dc2(p,ix,iy+2,mom2)]/w[fencode_dc2(p,ix,iy+2,rho)];
  fim2=w[fencode_dc2(p,ix,iy-2,energy)]*w[fencode_dc2(p,ix,iy-2,mom2)]/w[fencode_dc2(p,ix,iy-2,rho)];
}
  //ddcy=evalgrad_dc2(fi,fim1,0,0,p,1);
  ddcy=evalgrad_dc2(w[fencode_dc2(p,ix,iy+1,energy)]*w[fencode_dc2(p,ix,iy+1,mom2)]/w[fencode_dc2(p,ix,iy+1,rho)],w[fencode_dc2(p,ix,iy-1,energy)]*w[fencode_dc2(p,ix,iy-1,mom2)]/w[fencode_dc2(p,ix,iy-1,rho)],fip2,fim2,p,1);

  dd1=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);


 // dpi=(w[fencode_dc2(p,ix+1,iy,b1)]*w[fencode_dc2(p,ix+1,iy,mom1)]+w[fencode_dc2(p,ix+1,iy,b2)]*w[fencode_dc2(p,ix+1,iy,mom2)]+w[fencode_dc2(p,ix+1,iy,b3)]*w[fencode_dc2(p,ix+1,iy,mom3)])/w[fencode_dc2(p,ix+1,iy,rho)];
 // dpim1=(w[fencode_dc2(p,ix-1,iy,b1)]*w[fencode_dc2(p,ix-1,iy,mom1)]+w[fencode_dc2(p,ix-1,iy,b2)]*w[fencode_dc2(p,ix-1,iy,mom2)]+w[fencode_dc2(p,ix-1,iy,b3)]*w[fencode_dc2(p,ix-1,iy,mom3)])/w[fencode_dc2(p,ix-1,iy,rho)];
  //dpip2=(w[fencode_dc2(p,ix+2,iy,b1)]*w[fencode_dc2(p,ix+2,iy,mom1)]+w[fencode_dc2(p,ix+2,iy,b2)]*w[fencode_dc2(p,ix+2,iy,mom2)]+w[fencode_dc2(p,ix+2,iy,b3)]*w[fencode_dc2(p,ix+2,iy,mom3)])/w[fencode_dc2(p,ix+2,iy,rho)];
 // dpim2=(w[fencode_dc2(p,ix-2,iy,b1)]*w[fencode_dc2(p,ix-2,iy,mom1)]+w[fencode_dc2(p,ix-2,iy,b2)]*w[fencode_dc2(p,ix-2,iy,mom2)]+w[fencode_dc2(p,ix-2,iy,b3)]*w[fencode_dc2(p,ix-2,iy,mom3)])/w[fencode_dc2(p,ix-2,iy,rho)];

 // fi=dpi*w[fencode_dc2(p,ix+1,iy,b1)];
 // fim1=dpim1*w[fencode_dc2(p,ix-1,iy,b1)];
  //fip2=dpip2*w[fencode_dc2(p,ix+2,iy,b1)];
 // fim2=dpim2*w[fencode_dc2(p,ix-2,iy,b1)];
 // ddcx=evalgrad_dc2(fi,fim1,0,0,p,0);
 //  ddcx=evalgrad_dc2(((w[fencode_dc2(p,ix+1,iy,b1)]*w[fencode_dc2(p,ix+1,iy,mom1)]+w[fencode_dc2(p,ix+1,iy,b2)]*w[fencode_dc2(p,ix+1,iy,mom2)]+w[fencode_dc2(p,ix+1,iy,b3)]*w[fencode_dc2(p,ix+1,iy,mom3)])/w[fencode_dc2(p,ix+1,iy,rho)])*w[fencode_dc2(p,ix+1,iy,b1)],((w[fencode_dc2(p,ix-1,iy,b1)]*w[fencode_dc2(p,ix-1,iy,mom1)]+w[fencode_dc2(p,ix-1,iy,b2)]*w[fencode_dc2(p,ix-1,iy,mom2)]+w[fencode_dc2(p,ix-1,iy,b3)]*w[fencode_dc2(p,ix-1,iy,mom3)])/w[fencode_dc2(p,ix-1,iy,rho)])*w[fencode_dc2(p,ix-1,iy,b1)],0,0,p,0);

if(p->sodifon==1)
{
  fip2=wd[fencode_dc2(p,ix+2,iy,bdotv)]*w[fencode_dc2(p,ix+2,iy,b1)];
  fim2=wd[fencode_dc2(p,ix-2,iy,bdotv)]*w[fencode_dc2(p,ix-2,iy,b1)];
}

  ddcx=evalgrad_dc2(wd[fencode_dc2(p,ix+1,iy,bdotv)]*w[fencode_dc2(p,ix+1,iy,b1)],wd[fencode_dc2(p,ix-1,iy,bdotv)]*w[fencode_dc2(p,ix-1,iy,b1)],fip2,fim2,p,0);

 // dpi=(w[fencode_dc2(p,ix,iy+1,b1)]*w[fencode_dc2(p,ix,iy+1,mom1)]+w[fencode_dc2(p,ix,iy+1,b2)]*w[fencode_dc2(p,ix,iy+1,mom2)]+w[fencode_dc2(p,ix,iy+1,b3)]*w[fencode_dc2(p,ix,iy+1,mom3)])/w[fencode_dc2(p,ix,iy+1,rho)];
 // dpim1=(w[fencode_dc2(p,ix,iy-1,b1)]*w[fencode_dc2(p,ix,iy-1,mom1)]+w[fencode_dc2(p,ix,iy-1,b2)]*w[fencode_dc2(p,ix,iy-1,mom2)]+w[fencode_dc2(p,ix,iy-1,b3)]*w[fencode_dc2(p,ix,iy-1,mom3)])/w[fencode_dc2(p,ix,iy-1,rho)];  
  //dpip2=(w[fencode_dc2(p,ix,iy+2,b1)]*w[fencode_dc2(p,ix,iy+2,mom1)]+w[fencode_dc2(p,ix,iy+2,b2)]*w[fencode_dc2(p,ix,iy+2,mom2)]+w[fencode_dc2(p,ix,iy+2,b3)]*w[fencode_dc2(p,ix,iy+2,mom3)])/w[fencode_dc2(p,ix,iy+2,rho)];
  //dpim2=(w[fencode_dc2(p,ix,iy-2,b1)]*w[fencode_dc2(p,ix,iy-2,mom1)]+w[fencode_dc2(p,ix,iy-2,b2)]*w[fencode_dc2(p,ix,iy-2,mom2)]+w[fencode_dc2(p,ix,iy-2,b3)]*w[fencode_dc2(p,ix,iy-2,mom3)])/w[fencode_dc2(p,ix,iy-2,rho)];

 // fi=dpi*w[fencode_dc2(p,ix,iy+1,b2)];
 // fim1=dpim1*w[fencode_dc2(p,ix,iy-1,b2)];
if(p->sodifon==1)
{
  fip2=wd[fencode_dc2(p,ix,iy+2,bdotv)]*w[fencode_dc2(p,ix,iy+2,b2)];
  fim2=wd[fencode_dc2(p,ix,iy-2,bdotv)]*w[fencode_dc2(p,ix,iy-2,b2)];
}

//fi=w[fencode_dc2(p,ix,iy+1,b2)];
//  fim1=w[fencode_dc2(p,ix,iy-1,b2)];
  ddcy=evalgrad_dc2(wd[fencode_dc2(p,ix,iy+1,bdotv)]*w[fencode_dc2(p,ix,iy+1,b2)],wd[fencode_dc2(p,ix,iy-1,bdotv)]*w[fencode_dc2(p,ix,iy-1,b2)],fip2,fim2,p,1);
//ddcx=0;
//ddcy=evalgrad_dc2(((w[fencode_dc2(p,ix,iy+1,b1)]*w[fencode_dc2(p,ix,iy+1,mom1)]+w[fencode_dc2(p,ix,iy+1,b2)]*w[fencode_dc2(p,ix,iy+1,mom2)]+w[fencode_dc2(p,ix,iy+1,b3)]*w[fencode_dc2(p,ix,iy+1,mom3)])/w[fencode_dc2(p,ix,iy+1,rho)])*w[fencode_dc2(p,ix,iy+1,b2)],((w[fencode_dc2(p,ix,iy-1,b1)]*w[fencode_dc2(p,ix,iy-1,mom1)]+w[fencode_dc2(p,ix,iy-1,b2)]*w[fencode_dc2(p,ix,iy-1,mom2)]+w[fencode_dc2(p,ix,iy-1,b3)]*w[fencode_dc2(p,ix,iy-1,mom3)])/w[fencode_dc2(p,ix,iy-1,rho)])*w[fencode_dc2(p,ix,iy-1,b2)],0,0,p,1);

  dd2=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);



  ddcx=wd[fencode_dc2(p,ix,iy,pressuret)]*grad_dc2(w,p,ix,iy,mom1,0)/w[fencode_dc2(p,ix,iy,rho)];
  ddcy=wd[fencode_dc2(p,ix,iy,pressuret)]*grad_dc2(w,p,ix,iy,mom2,1)/w[fencode_dc2(p,ix,iy,rho)];


  dd3=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);

  
  return(dd1+dd2+dd3);
 //return dd1;
 // return ( ddc);
}


__device__ __host__
int derivcurrentb (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field, int direction) {

  int status=0;
        dw[fencode_dc2(p,ix,iy,field)]= -ddotcurrentb(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivcurrentenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=energy;
        dw[fencode_dc2(p,ix,iy,field)]= -ddotcurrentenergy(dw,wd,w,p,ix,iy);

  return ( status);
}

//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void derivcurrent2 (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case energy:
       derivcurrentenergy(dw,wd,w,p,ix,iy);
     break;
     case b1:
       derivcurrentb(dw,wd,w,p,ix,iy,field,0);
     break;
     case b2:
      derivcurrentb(dw,wd,w,p,ix,iy,field,1);
     break;
     case b3:
      derivcurrentb(dw,wd,w,p,ix,iy,field,2);
     break;
  }
  //return ( status);
}



__global__ void derivcurrent2_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
 // int index;
  //int ni=p->n[0];
  //int nj=p->n[1];
 // real dt=p->dt;
  //real dy=p->dx[1];
 // real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/(p->n[0]);
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*(p->n[0]));


  //if(i>(1+(p->sodifon==1)) && j >(1+(p->sodifon==1)) && i<((p->n[0])-1-(p->sodifon==1)) && j<((p->n[1])-1-(p->sodifon==1)))
if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
               /*for(int f=rho; f<=b3; f++)               
                  wmod[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               computebdotv(wmod,wd,p,i,j);*/
          
               for(int f=energy; f<=b3; f++)
               {              
                  derivcurrent2(dwn1+(NVAR*(p->n[0])*(p->n[1])*order),wd,wmod,p,i,j,f);
                 // dwn1[fencode_dc2(p,i,j,f)]=1.0;
                  //__syncthreads();
               }
               
               /*for(int f=rho; f<=b3; f++) 
                  wmod[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)]+0.5*dt*dwn1[fencode_dc2(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn2,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)]+0.5*dt*dwn2[fencode_dc2(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn3,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)]+dt*dwn3[fencode_dc2(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn4,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  {
                  wnew[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)]+(dt/6.0)*(
                     dwn1[fencode_dc2(p,i,j,f)]+2.0*dwn2[fencode_dc2(p,i,j,f)]
                         +2.0*dwn3[fencode_dc2(p,i,j,f)]+dwn4[fencode_dc2(p,i,j,f)]);
               }*/
              //  __syncthreads();
              /* for(int f=rho; f<=b3; f++)
                   wnew[fencode_dc2(p,i,j,f)]=w[fencode_dc2(p,i,j,f)]+dt*dwn1[fencode_dc2(p,i,j,f)];
               computej(wnew,wd,p,i,j);
               computepk(wnew,wd,p,i,j);
               computept(wnew,wd,p,i,j);*/ 


	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_dc2(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchroneous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}




int cuderivcurrent2(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
 //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
 // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
     derivcurrent2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}


