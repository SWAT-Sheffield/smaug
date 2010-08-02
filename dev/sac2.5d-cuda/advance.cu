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
int encode_adv (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_adv (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}

__device__ __host__
float evalgrad_adv(float fi, float fim1, float fip2, float fim2,struct params *p,int dir)
{
 //float valgrad;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx)))*(fi-fim1)-(1.0/(12.0*(p->dx)))*(fip2-fim2);
   return((1.0/(1.0*(p->dx)))*(fi-fim1));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(fi-fim1)-(1.0/(12.0*(p->dy)))*(fip2-fim2);
      return((1.0/(1.0*(p->dy)))*(fi-fim1));
 }

 return -1;
}


__device__ __host__
float grad_adv(float *wmod,struct params *p,int i,int j,int field,int dir)
{
 //float valgrad;

 if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx)))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
return((1.0/(1.0*(p->dx)))*(wmod[fencode_adv(p,i+1,j,field)]-wmod[fencode_adv(p,i-1,j,field)]));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dy)))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
 return((1.0/(1.0*(p->dy)))*(wmod[fencode_adv(p,i,j+1,field)]-wmod[fencode_adv(p,i,j-1,field)]));

 }

 return -1;
}

__device__ __host__
void computej_adv(float *wmod,float *wd,struct params *p,int i,int j)
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

  wd[fencode_adv(p,i,j,0)]=(grad_adv(wmod,p,i,j,b3,1))/(p->mu);
  wd[fencode_adv(p,i,j,1)]=(grad_adv(wmod,p,i,j,b3,0))/(p->mu);
  wd[fencode_adv(p,i,j,2)]=(grad_adv(wmod,p,i,j,b2,0)-grad_adv(wmod,p,i,j,b1,1))/(p->mu);
 
  //return ( status);
}

__device__ __host__
void computebdotv_adv(float *wmod,float *wd,struct params *p,int i,int j)
{
 // int status=0;
 //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
//  wd[fencode(p,i,j,4)]=  wd[fencode(p,i,j,3)]+0.5*(wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)]);

wd[fencode_adv(p,i,j,bdotv)]=(wmod[fencode_adv(p,i,j,b1)]*wmod[fencode_adv(p,i,j,mom1)]+wmod[fencode_adv(p,i,j,b2)]*wmod[fencode_adv(p,i,j,mom2)]+wmod[fencode_adv(p,i,j,b3)]*wmod[fencode_adv(p,i,j,mom3)])/wmod[fencode_adv(p,i,j,rho)];
 // return ( status);
}


__device__ __host__
void computepk_adv(float *wmod,float *wd,struct params *p,int i,int j)
{
 // int status=0;
 //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
  wd[fencode_adv(p,i,j,4)]=  wd[fencode_adv(p,i,j,3)]+0.5*(wmod[fencode_adv(p,i,j,b1)]*wmod[fencode_adv(p,i,j,b1)]+wmod[fencode_adv(p,i,j,b2)]*wmod[fencode_adv(p,i,j,b2)]+wmod[fencode_adv(p,i,j,b3)]*wmod[fencode_adv(p,i,j,b3)]);
 // return ( status);
}
__device__ __host__
void computept_adv(float *wmod,float *wd,struct params *p,int i,int j)
{
  //int status=0;
  //float momsq=wmod[fencode(p,i,j,mom1)]*wmod[fencode(p,i,j,mom1)]+wmod[fencode(p,i,j,mom2)]*wmod[fencode(p,i,j,mom2)]+wmod[fencode(p,i,j,mom3)]*wmod[fencode(p,i,j,mom3)];
  //float bsq=wmod[fencode(p,i,j,b1)]*wmod[fencode(p,i,j,b1)]+wmod[fencode(p,i,j,b2)]*wmod[fencode(p,i,j,b2)]+wmod[fencode(p,i,j,b3)]*wmod[fencode(p,i,j,b3)];
  wd[fencode_adv(p,i,j,3)]=((p->gamma)-1)*(wmod[fencode_adv(p,i,j,energy)]- 0.5*(wmod[fencode_adv(p,i,j,mom1)]*wmod[fencode_adv(p,i,j,mom1)]+wmod[fencode_adv(p,i,j,mom2)]*wmod[fencode_adv(p,i,j,mom2)]+wmod[fencode_adv(p,i,j,mom3)]*wmod[fencode_adv(p,i,j,mom3)])/wmod[fencode_adv(p,i,j,rho)]-0.5*(wmod[fencode_adv(p,i,j,b1)]*wmod[fencode_adv(p,i,j,b1)]+wmod[fencode_adv(p,i,j,b2)]*wmod[fencode_adv(p,i,j,b2)]+wmod[fencode_adv(p,i,j,b3)]*wmod[fencode_adv(p,i,j,b3)]) );
  //return ( status);
}

__global__ void advance_parallel(struct params *p, float *w, float *wnew, float *wmod, 
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
              /* for(int f=rho; f<=b3; f++)               
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
               }*/
               
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
                   wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)]+dt*dwn1[fencode_adv(p,i,j,f)];
               computej_adv(wnew,wd,p,i,j);
               computepk_adv(wnew,wd,p,i,j);
               computept_adv(wnew,wd,p,i,j);


	}
 __syncthreads();
  
}
/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_adv(char *label)
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






int cuadvance(struct params **p, float **w, float **wnew,struct params **d_p, float **d_w, float **d_wnew, float **d_wmod, float **d_dwn1, float **d_wd)
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
     advance_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
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



