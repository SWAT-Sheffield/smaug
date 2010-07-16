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
int encode_u (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_u (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
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
                  w[fencode_u(p,i,j,f)]=wnew[fencode_u(p,i,j,f)];
            // u[i+j*ni]=un[i+j*ni];
           // v[i+j*ni]=vn[i+j*ni];
	   // h[i+j*ni]=hn[i+j*ni];
	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_u(char *label)
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


int cuupdate(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd)
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
    // prop_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     //cudaThreadSynchronize();
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

  checkErrors_u("copy data from device");


  cudaFree(*d_p);

  cudaFree(*d_w);
  cudaFree(*d_wnew);
  cudaFree(*d_b);

  cudaFree(*d_wmod);
  cudaFree(*d_dwn1);
  cudaFree(*d_wd);



}
