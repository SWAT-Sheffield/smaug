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

 // float ddc=0;
//  int field=rho;

     // ddc= grad(w,p,ix,iy,mom1,0)+grad(w,p,ix,iy,mom2,1);
  return ( grad(w,p,ix,iy,mom1,0)+grad(w,p,ix,iy,mom2,1));
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



__global__ void init_parallel(struct params *p, float *w, float *wnew, float *b, float *wmod, 
    float *dwn1, float *wd)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int index,k;
int ni=p->ni;
  int nj=p->nj;

// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  float *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

 int nli = 0.45*(p->ni-1)+1;
  int nui = 0.55*(p->ni-1)+1;
  int nlj = 0.45*(p->nj-1)+1;
  int nuj = 0.55*(p->nj-1)+1; 
  int i,j;
   
   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->ni && j<p->nj)
	{
		b[i+j*(p->ni)]=0;

                 //Define b	
		if((i*(p->dx)) >20001)
		      b[j*(p->ni)+i]=0;
		else if((i*(p->dx)) <20000)
			//b[j*(p->ni)+i]=(5000/20000)*(20000-(i*(p->dx)));
                        b[j*(p->ni)+i]=0;
                        // b[j*(p->ni)+i]=5000*(1.0-(((float)i)/30.0));		



		//initialise the arrays here
               for(k=0;k<1;++k)
      		{
                    index=j*(p->ni)+i+k*(p->ni)*(p->nj);
                    //index=i+j*(p->ni)+(k*(p->nj)*(p->ni));
		    u[index]=0;
		    v[index]=0;
		    h[index]=5;
                    w[index+mom3*(p->ni)*(p->nj)]=0;
                    w[index+energy*(p->ni)*(p->nj)]=0;
                    w[index+b1*(p->ni)*(p->nj)]=0;
                    w[index+b2*(p->ni)*(p->nj)]=0;
                    w[index+b3*(p->ni)*(p->nj)]=0;

//float *wmod, 
//    float *dwn1, float *dwn2, float *dwn3, float *dwn4, float *wd)


      		}
		//h[iindex]=5000;
	
        __syncthreads();
        if(i>=nli && i<=nui && j>=nlj && j<=nuj)
	{
	   //j*(p->ni)+i;
           h[j*(p->ni)+i]=5.030;	
	}

       for(int f=0; f<=5; f++)
        { 
                  wd[fencode(p,i,j,f)]=0;
        }

        for(int f=rho; f<=b3; f++)
        {               
                  wnew[fencode(p,i,j,f)]=w[fencode(p,i,j,f)];
                  dwn1[fencode(p,i,j,f)]=0;
                  //dwn2[fencode(p,i,j,f)]=0;
                 // dwn3[fencode(p,i,j,f)]=0;
                  //dwn4[fencode(p,i,j,f)]=0;
                 
        }

	 __syncthreads();

			}	
	 __syncthreads();
  
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

__global__ void boundary_parallel(struct params *p, float *b, float *w, float *wnew)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;

  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;

  float *u,  *v,  *h;
  float *un,  *vn,  *hn;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

  hn=wnew+(p->ni)*(p->nj)*rho;
  un=wnew+(p->ni)*(p->nj)*mom1;
  vn=wnew+(p->ni)*(p->nj)*mom2;

    j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->ni && j<p->nj)
	{

		if(i==0 )
		{
			un[j*ni] = 2.5*un[1+j*ni] - 2*un[2+j*ni] + 0.5*un[3+j*ni];
			un[ni+j*ni] = 2.5*un[ni-1+j*ni] - 2*un[ni-2+ni*j] + 0.5*un[ni-3+j*ni];
			vn[j*ni] = 2.5*vn[1+j*ni] - 2*vn[2+j*ni] + 0.5*vn[3+j*ni];
		 	vn[ni+j*ni] = 2.5*vn[ni-1+j*ni] - 2*vn[ni-2+ni*j] + 0.5*vn[ni-3+j*ni];
		 	hn[j*ni] = 2.5*hn[1+j*ni] - 2*hn[2+j*ni] + 0.5*hn[3+j*ni];
			hn[ni+j*ni] = 2.5*hn[ni-1+j*ni] - 2*hn[ni-2+ni*j] + 0.5*hn[ni-3+j*ni];
		}

		if(j==0)
		{
			un[i+ni] = 2.5*un[i+1*ni] - 2*un[i+2*ni] + 0.5*un[i+3*ni];
			un[i+(nj )*ni] = 2.5*un[i+(nj-1)*ni] - 2*un[i+(nj-2)*ni] + 0.5*un[i+(nj-3)*ni];
			vn[i+ni] = 2.5*vn[i+1*ni] - 2*vn[i+2*ni] + 0.5*vn[i+3*ni];
			vn[i+(nj)*ni] = 2.5*vn[i+(nj-1)*ni] - 2*vn[i+(nj-2)*ni] + 0.5*vn[i+(nj-3)*ni];
			hn[i+ni] = 2.5*hn[i+1*ni] - 2*hn[i+2*ni] + 0.5*hn[i+3*ni];
			hn[i+(nj)*ni] = 2.5*hn[i+(nj-1)*ni] - 2*hn[i+(nj-2)*ni] + 0.5*hn[i+(nj-3)*ni];
		}
	}
 __syncthreads();
  
}

__global__ void update_parallel(struct params *p, float *b, float *w, float *wnew)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;


  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;
  float *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

  float *un,  *vn,  *hn;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  hn=wnew+(p->ni)*(p->nj)*rho;
  un=wnew+(p->ni)*(p->nj)*mom1;
  vn=wnew+(p->ni)*(p->nj)*mom2;
     j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  //if(i>2 && j >2 && i<((p->ni)-3) && j<((p->nj)-3))

  if(i<p->ni && j<p->nj)
	{
             for(int f=rho; f<=b3; f++)               
                  w[fencode(p,i,j,f)]=wnew[fencode(p,i,j,f)];
            // u[i+j*ni]=un[i+j*ni];
           // v[i+j*ni]=vn[i+j*ni];
	   // h[i+j*ni]=hn[i+j*ni];
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



int cuinit(struct params **p, float **w, float **wnew,  float **b, struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
{



/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
   
 // if (deviceCount == 0)
 // {
 //   fprintf(stderr, "Sorry, no CUDA device fount");
 //   return 1;
//  }
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }
  //cudaSetDevice(selectedDevice);
  printf("device count %d selected %d\n", deviceCount,selectedDevice);
  checkErrors("initialisations");
  
	// Build empty u, v, b matrices

  printf("in cuinit\n");
  float *adb;
  float *adw, *adwnew;
  struct params *adp;

  cudaMalloc((void**)d_wmod, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)d_dwn1, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)d_wd, 6*((*p)->ni)* ((*p)->nj)*sizeof(float));

  cudaMalloc((void**)&adw, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adwnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adb, 1*(((*p)->ni)* ((*p)->nj))*sizeof(float));
  cudaMalloc((void**)&adp, sizeof(struct params));
  checkErrors("memory allocation");

printf("ni is %d\n",(*p)->nj);

    *d_b=adb;
    *d_p=adp;
    *d_w=adw;
    *d_wnew=adwnew;


    cudaMemcpy(*d_w, *w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_wnew, *wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_b, *b, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
    
    dim3 dimBlock(16, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;
   

    printf("calling initialiser\n");
     //init_parallel(struct params *p, float *b, float *u, float *v, float *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b, *d_wmod, *d_dwn1,  *d_wd);
     cudaThreadSynchronize();
	    printf("called initialiser\n");
	cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(float), cudaMemcpyDeviceToHost);




  return 0;



}


int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
 dim3 dimBlock(16, 1);
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
     update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
 cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(float), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}

int cufinish(struct params **p, float **w, float **wnew, float **b, struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
{
  

 cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(float), cudaMemcpyDeviceToHost);

  checkErrors("copy data from device");


  cudaFree(*d_p);

  cudaFree(*d_w);
  cudaFree(*d_wnew);
  cudaFree(*d_b);

  cudaFree(*d_wmod);
  cudaFree(*d_dwn1);
  cudaFree(*d_wd);



}
