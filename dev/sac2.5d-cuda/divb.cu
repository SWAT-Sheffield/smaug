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
int encode_db (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_db (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}



__global__ void divb_parallel(struct params *p, struct state *s, real *w, real *wnew)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  __shared__ int ntot;

  int ni=p->ni;
  int nj=p->nj;
  real dt=p->dt;
  real dy=p->dy;
  real dx=p->dx;
  real g=p->g;
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

  real *un,  *vn,  *hn;
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
             {               
                  w[fencode_db(p,i,j,f)]=wnew[fencode_db(p,i,j,f)];

              }
            // u[i+j*ni]=un[i+j*ni];
           // v[i+j*ni]=vn[i+j*ni];
	   // h[i+j*ni]=hn[i+j*ni];
	}
 __syncthreads();



  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_db(char *label)
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

int cudivb(struct params **p, real **w, real **wnew,  struct state **state,struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;

//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
    // prop_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     //cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
    divb_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wnew);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
    //cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
   // cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}



