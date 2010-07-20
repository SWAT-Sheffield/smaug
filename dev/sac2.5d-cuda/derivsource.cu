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
int encode_ds (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_ds (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}

__device__ __host__
float evalgrad_ds(float fi, float fim1, float fip2, float fim2,struct params *p,int dir)
{
 //float valgrad_ds;

 if(dir == 0)
 {
     //valgrad_ds=(2.0/(3.0*(p->dx)))*(fi-fim1)-(1.0/(12.0*(p->dx)))*(fip2-fim2);
   return((1.0/(1.0*(p->dx)))*(fi-fim1));
 }
 else if(dir == 1)
 {
    // valgrad_ds=(2.0/(3.0*(p->dy)))*(fi-fim1)-(1.0/(12.0*(p->dy)))*(fip2-fim2);
      return((1.0/(1.0*(p->dy)))*(fi-fim1));
 }

 return -1;
}


__device__ __host__
float grad_ds(float *wmod,struct params *p,int i,int j,int field,int dir)
{
 //float valgrad_ds;

 if(dir == 0)
 {
    // valgrad_ds=(2.0/(3.0*(p->dx)))*(wmod[fencode_ds(p,i,j,field)]-wmod[fencode_ds(p,i-1,j,field)])-(1.0/(12.0*(p->dx)))*(wmod[fencode_ds(p,i+2,j,field)]-wmod[fencode_ds(p,i-2,j,field)]);
return((1.0/(1.0*(p->dx)))*(wmod[fencode_ds(p,i+1,j,field)]-wmod[fencode_ds(p,i-1,j,field)]));
 }
 else if(dir == 1)
 {
    // valgrad_ds=(2.0/(3.0*(p->dy)))*(wmod[fencode_ds(p,i,j,field)]-wmod[fencode_ds(p,i,j-1,field)])-(1.0/(12.0*(p->dy)))*(wmod[fencode_ds(p,i,j+2,field)]-wmod[fencode_ds(p,i,j-2,field)]);
 return((1.0/(1.0*(p->dy)))*(wmod[fencode_ds(p,i,j+1,field)]-wmod[fencode_ds(p,i,j-1,field)]));

 }

 return -1;
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
         return(w[fencode_ds(p,ix,iy,rho)]*(p->g1))-grad_ds(wd,p,ix,iy,pressuret,0);
	break;
	case 1:
         return(w[fencode_ds(p,ix,iy,rho)]*(p->g2))-grad_ds(wd,p,ix,iy,pressuret,1);
	break;
	case 2:
         return(w[fencode_ds(p,ix,iy,rho)]*(p->g3))-grad_ds(wd,p,ix,iy,pressuret,2);
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
         return(p->eta)*grad_ds(wd,p,ix,iy,current3,1);
	break;
	case 1:
         return -(p->eta)*grad_ds(wd,p,ix,iy,current3,0);
	break;
	case 2:
         return (p->eta)*(grad_ds(wd,p,ix,iy,current2,0)-grad_ds(wd,p,ix,iy,current1,1));
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
      srcg=(p->g1)*w[fencode_ds(p,ix,iy,mom1)]+(p->g2)*w[fencode_ds(p,ix,iy,mom2)]+(p->g3)*w[fencode_ds(p,ix,iy,mom3)];

       fi=(w[fencode_ds(p,ix+1,iy,b2)]*wd[fencode_ds(p,ix+1,iy,current3)]-w[fencode_ds(p,ix+1,iy,b3)]*wd[fencode_ds(p,ix+1,iy,current2)]);
       fim1=(w[fencode_ds(p,ix-1,iy,b2)]*wd[fencode_ds(p,ix-1,iy,current3)]-w[fencode_ds(p,ix-1,iy,b3)]*wd[fencode_ds(p,ix-1,iy,current2)]);
      // fip2=(w[fencode_ds(p,ix+2,iy,b2)]*wd[fencode_ds(p,ix+2,iy,current3)]-w[fencode_ds(p,ix+2,iy,b3)]*wd[fencode_ds(p,ix+2,iy,current2)]);
     //  fim2=(w[fencode_ds(p,ix-2,iy,b2)]*wd[fencode_ds(p,ix-2,iy,current3)]-w[fencode_ds(p,ix-2,iy,b3)]*wd[fencode_ds(p,ix-2,iy,current2)]);
      // ddcx=evalgrad_ds(fi,fim1,fip2,fim2,p,0);
      ddcx=evalgrad_ds(fi,fim1,0,0,p,0);

       fi=(w[fencode_ds(p,ix+1,iy,b3)]*wd[fencode_ds(p,ix+1,iy,current1)]-w[fencode_ds(p,ix+1,iy,b1)]*wd[fencode_ds(p,ix+1,iy,current3)]);
       fim1=(w[fencode_ds(p,ix,iy-1,b3)]*wd[fencode_ds(p,ix,iy-1,current1)]-w[fencode_ds(p,ix,iy-1,b1)]*wd[fencode_ds(p,ix,iy-1,current3)]);
     //  fip2=(w[fencode_ds(p,ix,iy+2,b3)]*wd[fencode_ds(p,ix,iy+2,current1)]-w[fencode_ds(p,ix,iy+2,b1)]*wd[fencode_ds(p,ix,iy+2,current3)]);
     //  fim2=(w[fencode_ds(p,ix,iy-2,b3)]*wd[fencode_ds(p,ix,iy-2,current1)]-w[fencode_ds(p,ix,iy-2,b1)]*wd[fencode_ds(p,ix,iy-2,current3)]);
      // ddcx=evalgrad_ds(fi,fim1,fip2,fim2,p,0);
      ddcy=evalgrad_ds(fi,fim1,0,0,p,1);

      srcb=(isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy);

 // src=srcg+srcb;
  return ( srcg+srcb);
}


__device__ __host__
int derivsourcerho (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=rho;
        dw[fencode_ds(p,ix,iy,field)]=dw[fencode_ds(p,ix,iy,field)]+sourcerho(dw,wd,w,p,ix,iy);
     	//dw[fencode_ds(p,ix,iy,field)]=w[fencode_ds(p,ix,iy,field)]+10;
  return ( status);
}

__device__ __host__
int derivsourcemom (float *dw, float *wd, float *w, struct params *p,int ix, int iy,int field, int direction) {

  int status=0;
     	//dw[fencode_ds(p,ix,iy,field)]=w[fencode_ds(p,ix,iy,field)]+20+5*(2*direction+1);
        dw[fencode_ds(p,ix,iy,field)]=dw[fencode_ds(p,ix,iy,field)]+sourcemom(dw,wd,w,p,ix,iy,field,direction);
        //dw[fencode_ds(p,ix,iy,field)]=-ddotcurrentmom(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivsourceb (float *dw, float *wd, float *w, struct params *p,int ix, int iy, int field, int direction) {

  int status=0;
        dw[fencode_ds(p,ix,iy,field)]=dw[fencode_ds(p,ix,iy,field)]+sourceb(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivsourceenergy (float *dw, float *wd, float *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=energy;
        dw[fencode_ds(p,ix,iy,field)]=dw[fencode_ds(p,ix,iy,field)]+sourceenergy(dw,wd,w,p,ix,iy);

  return ( status);
}

//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void derivsource (float *dw, float *wd, float *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case rho:
      derivsourcerho(dw,wd,w,p,ix,iy);
     break;
     case mom1:
      derivsourcemom(dw,wd,w,p,ix,iy,field,0);
     break;
     case mom2:
      derivsourcemom(dw,wd,w,p,ix,iy,field,1);
     break;
     case mom3:
      derivsourcemom(dw,wd,w,p,ix,iy,field,2);
     break;
     case energy:
       derivsourceenergy(dw,wd,w,p,ix,iy);
     break;
     case b1:
      derivsourceb(dw,wd,w,p,ix,iy,field,0);
     break;
     case b2:
      derivsourceb(dw,wd,w,p,ix,iy,field,1);
     break;
     case b3:
      derivsourceb(dw,wd,w,p,ix,iy,field,2);
     break;
  }
  //return ( status);
}


__global__ void derivsource_parallel(struct params *p, float *b, float *w, float *wnew, float *wmod, 
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
               /*for(int f=rho; f<=b3; f++)               
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               computebdotv(wmod,wd,p,i,j);*/
               for(int f=rho; f<=b3; f++)
               {              
                  derivsource(dwn1,wd,wmod,p,i,j,f);
                  //dwn1[fencode_ds(p,i,j,f)]=1.0;
                  __syncthreads();
               }
               
               /*for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+0.5*dt*dwn1[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn2,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+0.5*dt*dwn2[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn3,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+dt*dwn3[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn4,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  {
                  wnew[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+(dt/6.0)*(
                     dwn1[fencode_ds(p,i,j,f)]+2.0*dwn2[fencode_ds(p,i,j,f)]
                         +2.0*dwn3[fencode_ds(p,i,j,f)]+dwn4[fencode_ds(p,i,j,f)]);
               }*/
                __syncthreads();
              /* for(int f=rho; f<=b3; f++)
                   wnew[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+dt*dwn1[fencode_ds(p,i,j,f)];
               computej(wnew,wd,p,i,j);
               computepk(wnew,wd,p,i,j);
               computept(wnew,wd,p,i,j);*/ 


	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_ds(char *label)
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





int cuderivsource(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
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
     derivsource_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
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







