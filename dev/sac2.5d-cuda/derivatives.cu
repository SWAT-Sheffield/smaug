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
int encode (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}

__device__ __host__
float evalgrad(float fi, float fim1, float fip2, float fim2,struct params *p,int dir)
{
 //float valgrad;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx)))*(fi-fim1)-(1.0/(12.0*(p->dx)))*(fip2-fim2);
   return(1.0/(1.0*(p->dx)))*(fi-fim1);
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(fi-fim1)-(1.0/(12.0*(p->dy)))*(fip2-fim2);
      return(1.0/(1.0*(p->dy)))*(fi-fim1);
 }

 return -1;
}


__device__ __host__
float grad(float *wmod,struct params *p,int i,int j,int field,int dir)
{
 //float valgrad;

 if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx)))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
return(1.0/(1.0*(p->dx)))*(wmod[fencode(p,i+1,j,field)]-wmod[fencode(p,i-1,j,field)]);
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dy)))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
 return(1.0/(1.0*(p->dy)))*(wmod[fencode(p,i,j+1,field)]-wmod[fencode(p,i,j-1,field)]);

 }

 return -1;
}

__device__ __host__
void computej(float *wmod,float *wd,struct params *p,int i,int j)
{
 // int status=0;

 // float dbzdy, dbydz;
 // float dbzdx, dbxdz;
 // float dbydx, dbxdy;

 // dbzdy=grad(wmod,p,i,j,b3,1);
 // dbydz=0.0;
 // dbzdx=grad(wmod,p,i,j,b3,0);
//  dbxdz=0.0;
 // dbydx=grad(wmod,p,i,j,b2,0);
 // dbxdy=grad(wmod,p,i,j,b1,1);

  wd[fencode(p,i,j,0)]=(grad(wmod,p,i,j,b3,1))/(p->mu);
  wd[fencode(p,i,j,1)]=(grad(wmod,p,i,j,b3,0))/(p->mu);
  wd[fencode(p,i,j,2)]=(grad(wmod,p,i,j,b2,0)-grad(wmod,p,i,j,b1,1))/(p->mu);
 
  //return ( status);
}

__device__ __host__
void computebdotv(float *wmod,float *wd,struct params *p,int i,int j)
{
 // int status=0;
 //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
//  wd[fencode(p,i,j,4)]=  wd[fencode(p,i,j,3)]+0.5*(wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)]);

wd[fencode(p,i,j,bdotv)]=(wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,mom1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,mom2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,mom3)])/wmod[fencode(p,i,j,rho)];
 // return ( status);
}


__device__ __host__
void computepk(float *wmod,float *wd,struct params *p,int i,int j)
{
 // int status=0;
 //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
  wd[fencode(p,i,j,4)]=  wd[fencode(p,i,j,3)]+0.5*(wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)]);
 // return ( status);
}
__device__ __host__
void computept(float *wmod,float *wd,struct params *p,int i,int j)
{
  //int status=0;
  //float momsq=wmod[fencode(p,i,j,mom1)]*wmod[fencode(p,i,j,mom1)]+wmod[fencode(p,i,j,mom2)]*wmod[fencode(p,i,j,mom2)]+wmod[fencode(p,i,j,mom3)]*wmod[fencode(p,i,j,mom3)];
  //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
  wd[fencode(p,i,j,3)]=((p->gamma)-1)*(wmod[fencode(p,i,j,energy)]- 0.5*(wmod[fencode(p,i,j,mom1)]*wmod[fencode(p,i,j,mom1)]+wmod[fencode(p,i,j,mom2)]*wmod[fencode(p,i,j,mom2)]+wmod[fencode(p,i,j,mom3)]*wmod[fencode(p,i,j,mom3)])/wmod[fencode(p,i,j,rho)]-0.5*(wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)]) );
  //return ( status);
}

__device__ __host__
float sourcerho (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

 // float src=0;
 // int field=rho;
 
  return 0;
}

__device__ __host__
float sourcemom (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  //float src=0;
  switch(direction)
  {
	case 0:
         return(w[fencode(p,ix,iy,rho)]*(p->g1))-grad(wd,p,ix,iy,pressuret,0);
	break;
	case 1:
         return(w[fencode(p,ix,iy,rho)]*(p->g2))-grad(wd,p,ix,iy,pressuret,1);
	break;
	case 2:
         return(w[fencode(p,ix,iy,rho)]*(p->g3))-grad(wd,p,ix,iy,pressuret,2);
	break;
  }
  return 0;
}

__device__ __host__
float sourceb (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  //float src=0;
  switch(direction)
  {
	case 0:
         return(p->eta)*grad(wd,p,ix,iy,current3,1);
	break;
	case 1:
         return -(p->eta)*grad(wd,p,ix,iy,current3,0);
	break;
	case 2:
         return (p->eta)*(grad(wd,p,ix,iy,current2,0)-grad(wd,p,ix,iy,current1,1));
	break;
  }
  return 0;
}

__device__ __host__
float sourceenergy (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

 // float src=0;
  float srcg,srcb;
  int field=energy;
  float ddcx,ddcy;
  float fi,fim1;//fip2,fim2;
      srcg=(p->g1)*w[fencode(p,ix,iy,mom1)]+(p->g2)*w[fencode(p,ix,iy,mom2)]+(p->g3)*w[fencode(p,ix,iy,mom3)];

       fi=(w[fencode(p,ix+1,iy,b2)]*wd[fencode(p,ix+1,iy,current3)]-w[fencode(p,ix+1,iy,b3)]*wd[fencode(p,ix+1,iy,current2)]);
       fim1=(w[fencode(p,ix-1,iy,b2)]*wd[fencode(p,ix-1,iy,current3)]-w[fencode(p,ix-1,iy,b3)]*wd[fencode(p,ix-1,iy,current2)]);
      // fip2=(w[fencode(p,ix+2,iy,b2)]*wd[fencode(p,ix+2,iy,current3)]-w[fencode(p,ix+2,iy,b3)]*wd[fencode(p,ix+2,iy,current2)]);
     //  fim2=(w[fencode(p,ix-2,iy,b2)]*wd[fencode(p,ix-2,iy,current3)]-w[fencode(p,ix-2,iy,b3)]*wd[fencode(p,ix-2,iy,current2)]);
      // ddcx=evalgrad(fi,fim1,fip2,fim2,p,0);
      ddcx=evalgrad(fi,fim1,0,0,p,0);

       fi=(w[fencode(p,ix+1,iy,b3)]*wd[fencode(p,ix+1,iy,current1)]-w[fencode(p,ix+1,iy,b1)]*wd[fencode(p,ix+1,iy,current3)]);
       fim1=(w[fencode(p,ix,iy-1,b3)]*wd[fencode(p,ix,iy-1,current1)]-w[fencode(p,ix,iy-1,b1)]*wd[fencode(p,ix,iy-1,current3)]);
     //  fip2=(w[fencode(p,ix,iy+2,b3)]*wd[fencode(p,ix,iy+2,current1)]-w[fencode(p,ix,iy+2,b1)]*wd[fencode(p,ix,iy+2,current3)]);
     //  fim2=(w[fencode(p,ix,iy-2,b3)]*wd[fencode(p,ix,iy-2,current1)]-w[fencode(p,ix,iy-2,b1)]*wd[fencode(p,ix,iy-2,current3)]);
      // ddcx=evalgrad(fi,fim1,fip2,fim2,p,0);
      ddcy=evalgrad(fi,fim1,0,0,p,1);

      srcb=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);

 // src=srcg+srcb;
  return ( srcg+srcb);
}


__device__ __host__
float ddotcurrentrho (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

  float ddc=0;
//  int field=rho;

      ddc= grad(w,p,ix,iy,mom1,0)+grad(w,p,ix,iy,mom2,1);
  return ( isnan(ddc)?0:ddc);
}

__device__ __host__
float ddotcurrentmom (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  float ddc=0;
  float fi, fim1;
  //float  fip2=0, fim2=0;
  float ddc1,ddc2;
  float ddcx,ddcy;
   //     ddc= grad(w,p,ix,iy,mom1,0)+grad(w,p,ix,iy,mom2,1);
//evalgrad(float fi, float fim1, float fip2, float fim2,struct params *p,int dir)
  //fi=w(fencode(p,ix,iy,rho))
  //calculate momentum current

//w[fencode(p,ix,iy,rho)])=1;
//w[fencode(p,ix-1,iy,rho)])=1;
//w[fencode(p,ix+2,iy,rho)])=1;
//w[fencode(p,ix-2,iy,rho)])=1;
//w[fencode(p,ix,iy,rho)])=1;
//w[fencode(p,ix,iy-1,rho)])=1;
//w[fencode(p,ix,iy+2,rho)])=1;
//w[fencode(p,ix,iy-2,rho)])=1;

  switch(direction)
  {
    case 0:
       fi=(w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix+1,iy,rho)])*w[fencode(p,ix+1,iy,mom1)];
       fim1=(w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)])*w[fencode(p,ix-1,iy,mom1)];
    //   fip2=(w[fencode(p,ix+2,iy,mom1)]/w[fencode(p,ix+2,iy,rho)])*w[fencode(p,ix+2,iy,mom1)];
     //  fim2=(w[fencode(p,ix-2,iy,mom1)]/w[fencode(p,ix-2,iy,rho)])*w[fencode(p,ix-2,iy,mom1)];
      // ddcx=evalgrad(fi,fim1,fip2,fim2,p,0);
      ddcx=evalgrad(fi,fim1,0,0,p,0);
       //ddcx=fi-fim1;
       fi=(w[fencode(p,ix,iy+1,mom1)]/w[fencode(p,ix,iy+1,rho)])*w[fencode(p,ix,iy+1,mom2)];
       fim1=(w[fencode(p,ix,iy-1,mom1)]/w[fencode(p,ix,iy-1,rho)])*w[fencode(p,ix,iy-1,mom2)];
      // fip2=(w[fencode(p,ix,iy+2,mom1)]/w[fencode(p,ix,iy+2,rho)])*w[fencode(p,ix,iy+2,mom2)];
      // fim2=(w[fencode(p,ix,iy-2,mom1)]/w[fencode(p,ix,iy-2,rho)])*w[fencode(p,ix,iy-2,mom2)];
       //ddcy=fi;
       ddcy=evalgrad(fi,fim1,0,0,p,1);
       //ddcy=evalgrad(0,0,fip2,fim2,p,1);
    break;
    case 1:
       fi=(w[fencode(p,ix+1,iy,mom2)]/w[fencode(p,ix+1,iy,rho)])*w[fencode(p,ix+1,iy,mom1)];
       fim1=(w[fencode(p,ix-1,iy,mom2)]/w[fencode(p,ix-1,iy,rho)])*w[fencode(p,ix-1,iy,mom1)];
      // fip2=(w[fencode(p,ix+2,iy,mom2)]/w[fencode(p,ix+2,iy,rho)])*w[fencode(p,ix+2,iy,mom1)];
      // fim2=(w[fencode(p,ix-2,iy,mom2)]/w[fencode(p,ix-2,iy,rho)])*w[fencode(p,ix-2,iy,mom1)];
       ddcx=evalgrad(fi,fim1,0,0,p,0);
       fi=(w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)])*w[fencode(p,ix,iy+1,mom2)];
       fim1=(w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)])*w[fencode(p,ix,iy-1,mom2)];
      // fip2=(w[fencode(p,ix,iy+2,mom2)]/w[fencode(p,ix,iy+2,rho)])*w[fencode(p,ix,iy+2,mom2)];
      // fim2=(w[fencode(p,ix,iy-2,mom2)]/w[fencode(p,ix,iy-2,rho)])*w[fencode(p,ix,iy-2,mom2)];
       ddcy=evalgrad(fi,fim1,0,0,p,1);
    break;
    case 2:
       fi=(w[fencode(p,ix+1,iy,mom3)]/w[fencode(p,ix+1,iy,rho)])*w[fencode(p,ix+1,iy,mom1)];
       fim1=(w[fencode(p,ix-1,iy,mom3)]/w[fencode(p,ix-1,iy,rho)])*w[fencode(p,ix-1,iy,mom1)];
      // fip2=(w[fencode(p,ix+2,iy,mom3)]/w[fencode(p,ix+2,iy,rho)])*w[fencode(p,ix+2,iy,mom1)];
     //  fim2=(w[fencode(p,ix-2,iy,mom3)]/w[fencode(p,ix-2,iy,rho)])*w[fencode(p,ix-2,iy,mom1)];
       ddcx=evalgrad(fi,fim1,0,0,p,0);
       fi=(w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)])*w[fencode(p,ix,iy+1,mom2)];
       fim1=(w[fencode(p,ix,iy-1,mom3)]/w[fencode(p,ix,iy-1,rho)])*w[fencode(p,ix,iy-1,mom2)];
     //  fip2=(w[fencode(p,ix,iy+2,mom3)]/w[fencode(p,ix,iy+2,rho)])*w[fencode(p,ix,iy+2,mom2)];
     //  fim2=(w[fencode(p,ix,iy-2,mom3)]/w[fencode(p,ix,iy-2,rho)])*w[fencode(p,ix,iy-2,mom2)];
       ddcy=evalgrad(fi,fim1,0,0,p,1);
    break;
  }
  
  ddc1=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);
//fip2=0, fim2=0;
  //calculate bfield current
  switch(direction)
  {
    case 0:
       fi=w[fencode(p,ix+1,iy,b1)]*w[fencode(p,ix+1,iy,b1)];
       fim1=w[fencode(p,ix-1,iy,b1)]*w[fencode(p,ix-1,iy,b1)];
     // fip2=w[fencode(p,ix+2,iy,b1)]*w[fencode(p,ix+2,iy,b1)];
     //  fim2=w[fencode(p,ix-2,iy,b1)]*w[fencode(p,ix-2,iy,b1)];
       ddcx=evalgrad(fi,fim1,0,0,p,0);
       fi=w[fencode(p,ix,iy+1,b1)]*w[fencode(p,ix,iy+1,b2)];
       fim1=w[fencode(p,ix,iy-1,b1)]*w[fencode(p,ix,iy-1,b2)];
     //  fip2=w[fencode(p,ix,iy+2,b1)]*w[fencode(p,ix,iy+2,b2)];
      // fim2=w[fencode(p,ix,iy-2,b1)]*w[fencode(p,ix,iy-2,b2)];
       ddcy=evalgrad(fi,fim1,0,0,p,1);
    break;
    case 1:
       fi=w[fencode(p,ix+1,iy,b2)]*w[fencode(p,ix+1,iy,b1)];
       fim1=w[fencode(p,ix-1,iy,b2)]*w[fencode(p,ix-1,iy,b1)];
     //  fip2=w[fencode(p,ix+2,iy,b2)]*w[fencode(p,ix+2,iy,b1)];
      // fim2=w[fencode(p,ix-2,iy,b2)]*w[fencode(p,ix-2,iy,b1)];
       ddcx=evalgrad(fi,fim1,0,0,p,0);
       fi=w[fencode(p,ix,iy+1,b2)]*w[fencode(p,ix,iy+1,b2)];
       fim1=w[fencode(p,ix,iy-1,b2)]*w[fencode(p,ix,iy-1,b2)];
      // fip2=w[fencode(p,ix,iy+2,b2)]*w[fencode(p,ix,iy+2,b2)];
      // fim2=w[fencode(p,ix,iy-2,b2)]*w[fencode(p,ix,iy-2,b2)];
       ddcy=evalgrad(fi,fim1,0,0,p,1);
    break;
    case 2:
       fi=w[fencode(p,ix+1,iy,b3)]*w[fencode(p,ix+1,iy,b1)];
       fim1=w[fencode(p,ix-1,iy,b3)]*w[fencode(p,ix-1,iy,b1)];
      // fip2=w[fencode(p,ix+2,iy,b3)]*w[fencode(p,ix+2,iy,b1)];
      // fim2=w[fencode(p,ix-2,iy,b3)]*w[fencode(p,ix-2,iy,b1)];
       ddcx=evalgrad(fi,fim1,0,0,p,0);
       fi=w[fencode(p,ix,iy+1,b3)]*w[fencode(p,ix,iy+1,b2)];
       fim1=w[fencode(p,ix,iy-1,b3)]*w[fencode(p,ix,iy-1,b2)];
      // fip2=w[fencode(p,ix,iy+2,b3)]*w[fencode(p,ix,iy+2,b2)];
     //  fim2=w[fencode(p,ix,iy-2,b3)]*w[fencode(p,ix,iy-2,b2)];
       ddcy=evalgrad(fi,fim1,0,0,p,1);
    break;
  }
  //ddc2=ddcx+ddcy;
  ddc2=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);

  //ddc=ddc1-ddc2;

  return ( ddc1-ddc2);
}

__device__ __host__
float ddotcurrentb (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  //float ddc=0;

  float fi, fim1;// fip2=0, fim2=0;
  float ddc1,ddc2;
  float ddcx,ddcy;

  switch(direction)
  {
	case 0:
	       fi=w[fencode(p,ix+1,iy,mom1)]*w[fencode(p,ix+1,iy,b1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,mom1)]*w[fencode(p,ix-1,iy,b1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,mom1)]*w[fencode(p,ix+2,iy,b1)]/w[fencode(p,ix+2,iy,rho)];
	       //fim2=w[fencode(p,ix-2,iy,mom1)]*w[fencode(p,ix-2,iy,b1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,mom1)]*w[fencode(p,ix,iy+1,b2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,mom1)]*w[fencode(p,ix,iy-1,b2)]/w[fencode(p,ix,iy-1,rho)];
	       //fip2=w[fencode(p,ix,iy+2,mom1)]*w[fencode(p,ix,iy+2,b2)]/w[fencode(p,ix,iy+2,rho)];
	       //fim2=w[fencode(p,ix,iy-2,mom1)]*w[fencode(p,ix,iy-2,b2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);
        break;
	case 1:
	       fi=w[fencode(p,ix+1,iy,mom2)]*w[fencode(p,ix+1,iy,b1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,mom2)]*w[fencode(p,ix-1,iy,b1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,mom2)]*w[fencode(p,ix+2,iy,b1)]/w[fencode(p,ix+2,iy,rho)];
	       //fim2=w[fencode(p,ix-2,iy,mom2)]*w[fencode(p,ix-2,iy,b1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,mom2)]*w[fencode(p,ix,iy+1,b2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,mom2)]*w[fencode(p,ix,iy-1,b2)]/w[fencode(p,ix,iy-1,rho)];
	       //fip2=w[fencode(p,ix,iy+2,mom2)]*w[fencode(p,ix,iy+2,b2)]/w[fencode(p,ix,iy+2,rho)];
	       //fim2=w[fencode(p,ix,iy-2,mom2)]*w[fencode(p,ix,iy-2,b2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);
        break;
	case 2:
	       fi=w[fencode(p,ix+1,iy,mom3)]*w[fencode(p,ix+1,iy,b1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,mom3)]*w[fencode(p,ix-1,iy,b1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,mom3)]*w[fencode(p,ix+2,iy,b1)]/w[fencode(p,ix+2,iy,rho)];
	       //fim2=w[fencode(p,ix-2,iy,mom3)]*w[fencode(p,ix-2,iy,b1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,mom3)]*w[fencode(p,ix,iy+1,b2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,mom3)]*w[fencode(p,ix,iy-1,b2)]/w[fencode(p,ix,iy-1,rho)];
	       //fip2=w[fencode(p,ix,iy+2,mom3)]*w[fencode(p,ix,iy+2,b2)]/w[fencode(p,ix,iy+2,rho)];
	       //fim2=w[fencode(p,ix,iy-2,mom3)]*w[fencode(p,ix,iy-2,b2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);

        break;
  }
  ddc1=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);



  switch(direction)
  {
	case 0:
	       fi=w[fencode(p,ix+1,iy,b1)]*w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,b1)]*w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,b1)]*w[fencode(p,ix+2,iy,mom1)]/w[fencode(p,ix+2,iy,rho)];
	      // fim2=w[fencode(p,ix-2,iy,b1)]*w[fencode(p,ix-2,iy,mom1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,b1)]*w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,b1)]*w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)];
	       //fip2=w[fencode(p,ix,iy+2,b1)]*w[fencode(p,ix,iy+2,mom2)]/w[fencode(p,ix,iy+2,rho)];
	       //fim2=w[fencode(p,ix,iy-2,b1)]*w[fencode(p,ix,iy-2,mom2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);
        break;
	case 1:
	       fi=w[fencode(p,ix+1,iy,b2)]*w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,b2)]*w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,b2)]*w[fencode(p,ix+2,iy,mom1)]/w[fencode(p,ix+2,iy,rho)];
	      // fim2=w[fencode(p,ix-2,iy,b2)]*w[fencode(p,ix-2,iy,mom1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,b2)]*w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,b2)]*w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)];
	      // fip2=w[fencode(p,ix,iy+2,b2)]*w[fencode(p,ix,iy+2,mom2)]/w[fencode(p,ix,iy+2,rho)];
	      // fim2=w[fencode(p,ix,iy-2,b2)]*w[fencode(p,ix,iy-2,mom2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);
        break;
	case 2:
	       fi=w[fencode(p,ix+1,iy,b3)]*w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix+1,iy,rho)];
	       fim1=w[fencode(p,ix-1,iy,b3)]*w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)];
	       //fip2=w[fencode(p,ix+2,iy,b3)]*w[fencode(p,ix+2,iy,mom1)]/w[fencode(p,ix+2,iy,rho)];
	       //fim2=w[fencode(p,ix-2,iy,b3)]*w[fencode(p,ix-2,iy,mom1)]/w[fencode(p,ix-2,iy,rho)];
	       ddcx=evalgrad(fi,fim1,0,0,p,0);
	       fi=w[fencode(p,ix,iy+1,b3)]*w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)];
	       fim1=w[fencode(p,ix,iy-1,b3)]*w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)];
	       //fip2=w[fencode(p,ix,iy+2,b3)]*w[fencode(p,ix,iy+2,mom2)]/w[fencode(p,ix,iy+2,rho)];
	       //fim2=w[fencode(p,ix,iy-2,b3)]*w[fencode(p,ix,iy-2,mom2)]/w[fencode(p,ix,iy-2,rho)];
	       ddcy=evalgrad(fi,fim1,0,0,p,1);
        break;
  }
  ddc2=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);
  return(ddc1-ddc2);

}

__device__ __host__
float ddotcurrentenergy (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

 // float ddc=0;
  float dd1,dd2,dd3;
 
  float ddcx,ddcy;
  //float fi, fim1;//fip2=0, fim2=0;
  //float dpi, dpim1;//, dpip2=0, dpim2=0;


  //int field=energy;

  //fi=w[fencode(p,ix+1,iy,energy)]*w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix,iy,rho)];
  //fim1=w[fencode(p,ix-1,iy,energy)]*w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)];
  //fip2=w[fencode(p,ix+2,iy,energy)]*w[fencode(p,ix+2,iy,mom1)]/w[fencode(p,ix+2,iy,rho)];
 // fim2=w[fencode(p,ix-2,iy,energy)]*w[fencode(p,ix-2,iy,mom1)]/w[fencode(p,ix-2,iy,rho)];
 // ddcx=evalgrad(fi,fim1,0,0,p,0);
  ddcx=evalgrad(w[fencode(p,ix+1,iy,energy)]*w[fencode(p,ix+1,iy,mom1)]/w[fencode(p,ix,iy,rho)],w[fencode(p,ix-1,iy,energy)]*w[fencode(p,ix-1,iy,mom1)]/w[fencode(p,ix-1,iy,rho)],0,0,p,0);

 // fi=w[fencode(p,ix,iy+1,energy)]*w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)];
 // fim1=w[fencode(p,ix,iy-1,energy)]*w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)];
 // fip2=w[fencode(p,ix,iy+2,energy)]*w[fencode(p,ix,iy+2,mom2)]/w[fencode(p,ix,iy+2,rho)];
  //fim2=w[fencode(p,ix,iy-2,energy)]*w[fencode(p,ix,iy-2,mom2)]/w[fencode(p,ix,iy-2,rho)];
  //ddcy=evalgrad(fi,fim1,0,0,p,1);
  ddcy=evalgrad(w[fencode(p,ix,iy+1,energy)]*w[fencode(p,ix,iy+1,mom2)]/w[fencode(p,ix,iy+1,rho)],w[fencode(p,ix,iy-1,energy)]*w[fencode(p,ix,iy-1,mom2)]/w[fencode(p,ix,iy-1,rho)],0,0,p,1);

  dd1=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);


 // dpi=(w[fencode(p,ix+1,iy,b1)]*w[fencode(p,ix+1,iy,mom1)]+w[fencode(p,ix+1,iy,b2)]*w[fencode(p,ix+1,iy,mom2)]+w[fencode(p,ix+1,iy,b3)]*w[fencode(p,ix+1,iy,mom3)])/w[fencode(p,ix+1,iy,rho)];
 // dpim1=(w[fencode(p,ix-1,iy,b1)]*w[fencode(p,ix-1,iy,mom1)]+w[fencode(p,ix-1,iy,b2)]*w[fencode(p,ix-1,iy,mom2)]+w[fencode(p,ix-1,iy,b3)]*w[fencode(p,ix-1,iy,mom3)])/w[fencode(p,ix-1,iy,rho)];
  //dpip2=(w[fencode(p,ix+2,iy,b1)]*w[fencode(p,ix+2,iy,mom1)]+w[fencode(p,ix+2,iy,b2)]*w[fencode(p,ix+2,iy,mom2)]+w[fencode(p,ix+2,iy,b3)]*w[fencode(p,ix+2,iy,mom3)])/w[fencode(p,ix+2,iy,rho)];
 // dpim2=(w[fencode(p,ix-2,iy,b1)]*w[fencode(p,ix-2,iy,mom1)]+w[fencode(p,ix-2,iy,b2)]*w[fencode(p,ix-2,iy,mom2)]+w[fencode(p,ix-2,iy,b3)]*w[fencode(p,ix-2,iy,mom3)])/w[fencode(p,ix-2,iy,rho)];

 // fi=dpi*w[fencode(p,ix+1,iy,b1)];
 // fim1=dpim1*w[fencode(p,ix-1,iy,b1)];
  //fip2=dpip2*w[fencode(p,ix+2,iy,b1)];
 // fim2=dpim2*w[fencode(p,ix-2,iy,b1)];
 // ddcx=evalgrad(fi,fim1,0,0,p,0);
 //  ddcx=evalgrad(((w[fencode(p,ix+1,iy,b1)]*w[fencode(p,ix+1,iy,mom1)]+w[fencode(p,ix+1,iy,b2)]*w[fencode(p,ix+1,iy,mom2)]+w[fencode(p,ix+1,iy,b3)]*w[fencode(p,ix+1,iy,mom3)])/w[fencode(p,ix+1,iy,rho)])*w[fencode(p,ix+1,iy,b1)],((w[fencode(p,ix-1,iy,b1)]*w[fencode(p,ix-1,iy,mom1)]+w[fencode(p,ix-1,iy,b2)]*w[fencode(p,ix-1,iy,mom2)]+w[fencode(p,ix-1,iy,b3)]*w[fencode(p,ix-1,iy,mom3)])/w[fencode(p,ix-1,iy,rho)])*w[fencode(p,ix-1,iy,b1)],0,0,p,0);
  ddcx=evalgrad(wd[fencode(p,ix+1,iy,bdotv)]*w[fencode(p,ix+1,iy,b1)],wd[fencode(p,ix-1,iy,bdotv)]*w[fencode(p,ix-1,iy,b1)],0,0,p,1);

 // dpi=(w[fencode(p,ix,iy+1,b1)]*w[fencode(p,ix,iy+1,mom1)]+w[fencode(p,ix,iy+1,b2)]*w[fencode(p,ix,iy+1,mom2)]+w[fencode(p,ix,iy+1,b3)]*w[fencode(p,ix,iy+1,mom3)])/w[fencode(p,ix,iy+1,rho)];
 // dpim1=(w[fencode(p,ix,iy-1,b1)]*w[fencode(p,ix,iy-1,mom1)]+w[fencode(p,ix,iy-1,b2)]*w[fencode(p,ix,iy-1,mom2)]+w[fencode(p,ix,iy-1,b3)]*w[fencode(p,ix,iy-1,mom3)])/w[fencode(p,ix,iy-1,rho)];  
  //dpip2=(w[fencode(p,ix,iy+2,b1)]*w[fencode(p,ix,iy+2,mom1)]+w[fencode(p,ix,iy+2,b2)]*w[fencode(p,ix,iy+2,mom2)]+w[fencode(p,ix,iy+2,b3)]*w[fencode(p,ix,iy+2,mom3)])/w[fencode(p,ix,iy+2,rho)];
  //dpim2=(w[fencode(p,ix,iy-2,b1)]*w[fencode(p,ix,iy-2,mom1)]+w[fencode(p,ix,iy-2,b2)]*w[fencode(p,ix,iy-2,mom2)]+w[fencode(p,ix,iy-2,b3)]*w[fencode(p,ix,iy-2,mom3)])/w[fencode(p,ix,iy-2,rho)];

 // fi=dpi*w[fencode(p,ix,iy+1,b2)];
 // fim1=dpim1*w[fencode(p,ix,iy-1,b2)];
  //fip2=dpip2*w[fencode(p,ix,iy+2,b2)];
  //fim2=dpim2*w[fencode(p,ix,iy-2,b2)];

//fi=w[fencode(p,ix,iy+1,b2)];
//  fim1=w[fencode(p,ix,iy-1,b2)];
  ddcy=evalgrad(wd[fencode(p,ix,iy+1,bdotv)]*w[fencode(p,ix,iy+1,b2)],wd[fencode(p,ix,iy-1,bdotv)]*w[fencode(p,ix,iy-1,b2)],0,0,p,1);
//ddcx=0;
//ddcy=evalgrad(((w[fencode(p,ix,iy+1,b1)]*w[fencode(p,ix,iy+1,mom1)]+w[fencode(p,ix,iy+1,b2)]*w[fencode(p,ix,iy+1,mom2)]+w[fencode(p,ix,iy+1,b3)]*w[fencode(p,ix,iy+1,mom3)])/w[fencode(p,ix,iy+1,rho)])*w[fencode(p,ix,iy+1,b2)],((w[fencode(p,ix,iy-1,b1)]*w[fencode(p,ix,iy-1,mom1)]+w[fencode(p,ix,iy-1,b2)]*w[fencode(p,ix,iy-1,mom2)]+w[fencode(p,ix,iy-1,b3)]*w[fencode(p,ix,iy-1,mom3)])/w[fencode(p,ix,iy-1,rho)])*w[fencode(p,ix,iy-1,b2)],0,0,p,1);

  dd2=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);



  ddcx=wd[fencode(p,ix,iy,pressuret)]*grad(w,p,ix,iy,mom1,0)/w[fencode(p,ix,iy,rho)];
  ddcy=wd[fencode(p,ix,iy,pressuret)]*grad(w,p,ix,iy,mom2,1)/w[fencode(p,ix,iy,rho)];


  dd3=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);

  
  return(dd1+dd2+dd3);
 //return dd1;
 // return ( ddc);
}

__device__ __host__
int derivrho (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=rho;
        dw[fencode(p,ix,iy,field)]=sourcerho(dw,wd,w,p,ix,iy)-ddotcurrentrho(dw,wd,w,p,ix,iy);
     	//dw[fencode(p,ix,iy,field)]=w[fencode(p,ix,iy,field)]+10;
  return ( status);
}

__device__ __host__
int derivmom (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  int status=0;
     	//dw[fencode(p,ix,iy,field)]=w[fencode(p,ix,iy,field)]+20+5*(2*direction+1);
        dw[fencode(p,ix,iy,field)]=sourcemom(dw,wd,w,p,ix,iy,field,direction)-ddotcurrentmom(dw,wd,w,p,ix,iy,field,direction);
        //dw[fencode(p,ix,iy,field)]=-ddotcurrentmom(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivb (float *dw, float *wd, float *w, struct params *p,int ix, int iy, int field, int direction) {

  int status=0;
        dw[fencode(p,ix,iy,field)]=sourceb(dw,wd,w,p,ix,iy,field,direction)-ddotcurrentb(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivenergy (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=energy;
        dw[fencode(p,ix,iy,field)]=sourceenergy(dw,wd,w,p,ix,iy)-ddotcurrentenergy(dw,wd,w,p,ix,iy);

  return ( status);
}

//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void deriv (float *dw, float *wd, float *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case rho:
      derivrho(dw,wd,w,p,ix,iy);
     break;
     case mom1:
      derivmom(dw,wd,w,p,ix,iy,field,0);
     break;
     case mom2:
      derivmom(dw,wd,w,p,ix,iy,field,1);
     break;
     case mom3:
      derivmom(dw,wd,w,p,ix,iy,field,2);
     break;
     case energy:
       derivenergy(dw,wd,w,p,ix,iy);
     break;
     case b1:
      derivb(dw,wd,w,p,ix,iy,field,0);
     break;
     case b2:
      derivb(dw,wd,w,p,ix,iy,field,1);
     break;
     case b3:
      derivb(dw,wd,w,p,ix,iy,field,2);
     break;
  }
  //return ( status);
}



__global__ void prop_parallel(struct params *p, float *b, float *w, float *wnew, float *wmod, 
    float *dwn1, float *wd)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i>1 && j >1 && i<((p->ni)-2) && j<((p->nj)-2))
	{		               
               for(int f=rho; f<=b3; f++)               
                  wmod[fencode(p,i,j,f)]=w[fencode(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               computebdotv(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++)
               {              
                  deriv(dwn1,wd,wmod,p,i,j,f);
                  //dwn1[fencode(p,i,j,f)]=1.0;
                  __syncthreads();
               }
               
               /*for(int f=rho; f<=b3; f++) 
                  wmod[fencode(p,i,j,f)]=w[fencode(p,i,j,f)]+0.5*dt*dwn1[fencode(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn2,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode(p,i,j,f)]=w[fencode(p,i,j,f)]+0.5*dt*dwn2[fencode(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn3,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode(p,i,j,f)]=w[fencode(p,i,j,f)]+dt*dwn3[fencode(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn4,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  {
                  wnew[fencode(p,i,j,f)]=w[fencode(p,i,j,f)]+(dt/6.0)*(
                     dwn1[fencode(p,i,j,f)]+2.0*dwn2[fencode(p,i,j,f)]
                         +2.0*dwn3[fencode(p,i,j,f)]+dwn4[fencode(p,i,j,f)]);
               }*/
                __syncthreads();
               for(int f=rho; f<=b3; f++)
                   wnew[fencode(p,i,j,f)]=w[fencode(p,i,j,f)]+dt*dwn1[fencode(p,i,j,f)];
               computej(wnew,wd,p,i,j);
               computepk(wnew,wd,p,i,j);
               computept(wnew,wd,p,i,j);


	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors(char *label)
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




int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;

//__global__ void prop_parallel(struct params *p, float *b, float *w, float *wnew, float *wmod, 
  //  float *dwn1, float *dwn2, float *dwn3, float *dwn4, float *wd)
     //init_parallel(struct params *p, float *b, float *u, float *v, float *h)
     prop_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(float), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}



