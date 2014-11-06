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
#include "gradops_dc2.cuh"

__device__ __host__
real transportflux_dc2 (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;

   //transport flux
    switch(direction)
  {
     case 0:
      flux= wd[fencode_dc2(p,ix,iy,vel1)]*w[fencode_dc2(p,ix,iy,field)];
     break;
     case 1:
      flux= wd[fencode_dc2(p,ix,iy,vel2)]*w[fencode_dc2(p,ix,iy,field)];
     break;
     case 2:
      flux= wd[fencode_dc2(p,ix,iy,vel3)]*w[fencode_dc2(p,ix,iy,field)];
     break;
   }
  return flux;


  //return ( ddc1-ddc2);
}




__device__ __host__
real fluxb1(real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;


        #ifdef USE_SAC
      		flux= -(w[fencode_dc2(p,ix,iy,field)]+w[fencode_dc2(p,ix,iy,field+5)])*w[fencode_dc2(p,ix,iy,mom1+direction)]/(w[fencode_dc2(p,ix,iy,rho)]+w[fencode_dc2(p,ix,iy,rhob)]);

         #else
      		flux= -w[fencode_dc2(p,ix,iy,field)]*w[fencode_dc2(p,ix,iy,mom1+direction)]/w[fencode_dc2(p,ix,iy,rho)];
         #endif


  return flux;
}



__device__ __host__
real fluxe1(real *dw, real *wd, real *w, struct params *p,int ix, int iy, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;


        #ifdef USE_SAC
      		flux= -w[fencode_dc2(p,ix,iy,b1+direction)]*wd[fencode_dc2(p,ix,iy,bdotv)]+(w[fencode_dc2(p,ix,iy,mom1+direction)]*wd[fencode_dc2(p,ix,iy,pressuret)]/w[fencode_dc2(p,ix,iy,rho)]);
         #else
      		flux= -w[fencode_dc2(p,ix,iy,b1+direction)]*wd[fencode_dc2(p,ix,iy,bdotv)]+(w[fencode_dc2(p,ix,iy,mom1+direction)]*wd[fencode_dc2(p,ix,iy,pressuret)]/w[fencode_dc2(p,ix,iy,rho)]);
         #endif

  return flux;


  //return ( ddc1-ddc2);
}








__device__ __host__
int computefluxe(real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  int field, direction;
  int status=0;
  for(direction=0;direction<3;direction++)
         #ifdef USE_SAC
	     wd[fencode_dc2(p,ix,iy,f1+direction)]= transportflux_dc2(dw,wd,w,p,ix,iy,energy,direction)+transportflux_dc2(dw,wd,w,p,ix,iy,energyb,direction)+fluxe1(dw,wd,w,p,ix,iy,direction);
         #else
             wd[fencode_dc2(p,ix,iy,f1+direction)]= transportflux_dc2(dw,wd,w,p,ix,iy,energy,direction)+fluxe1(dw,wd,w,p,ix,iy,direction);
         #endif
        
  return ( status);
}

__device__ __host__
int computefluxb (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  int direction;
  int status=0;
  for(direction=0;direction<3;direction++)
         #ifdef USE_SAC
     wd[fencode_dc2(p,ix,iy,f1+direction)]= transportflux_dc2(dw,wd,w,p,ix,iy,field,direction)+transportflux_dc2(dw,wd,w,p,ix,iy,field+5,direction)+fluxb1(dw,wd,w,p,ix,iy,field,direction);
         #else
     wd[fencode_dc2(p,ix,iy,f1+direction)]= transportflux_dc2(dw,wd,w,p,ix,iy,field,direction)+fluxb1(dw,wd,w,p,ix,iy,field,direction);
         #endif
        
  return ( status);
}

__device__ __host__
int divflux_dc2(real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field) {

  int direction;
  int status=0;
  real divflux=0;
  dw[fencode_dc2(p,ix,iy,field)]= grad_dc2(wd,p,ix,iy,f1,0)+grad_dc2(wd,p,ix,iy,f2,1);      

  return ( status);
}





//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux_dc2 (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case energy:
      computefluxe(dw,wd,w,p,ix,iy);
      // add the following terms for SAC
      // del((b bb+ bb b).v)+ptb del v - bb bb del v
     break;
     case b1:
      computefluxb(dw,wd,w,p,ix,iy,field);
     break;
     case b2:
       computefluxb(dw,wd,w,p,ix,iy,field);
     break;
     case b3:
      computefluxb(dw,wd,w,p,ix,iy,field);
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
  int ni=p->n[0];
  int nj=p->n[1];
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


               for(int f=energy; f<NVAR; f++)
               {

			if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                  		computeflux_dc2(dwn1+(NVAR*(p->n[0])*(p->n[1])*order),wd,wmod,p,i,j,f); 

                        //might need to set boundaries correctly 
                        __syncthreads();
			if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                                divflux_dc2(dwn1+(NVAR*(p->n[0])*(p->n[1])*order),wd,wmod,p,i,j,f); 
                }


  
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


