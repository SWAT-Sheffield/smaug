//#define MODID pre


#include "cudapars.h"
#include "paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "step.h"
#include "gradops_cdf.cuh"
#include "dervfields_cdf.cuh"
/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////


__global__ void computedervfields_parallel(struct params *p,  real *w,  real *wmod, 
    real *dwn1, real *wd, int order,int ordero)
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
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
if(order == 0)
if(i<((p->n[0])) && j<((p->n[1])))
	{		
 
               for(int f=rho; f<=b2; f++)
                  wmod[fencode_cdf(p,i,j,f)+((p->n[0]))*((p->n[1]))*NVAR]=wmod[fencode_cdf(p,i,j,f)]; 
        }
               __syncthreads();


if(i<((p->n[0])) && j<((p->n[1])))
	{		

               for(int f=vel1; f<NDERV; f++)
                  wd[fencode_cdf(p,i,j,f)]=0; 
               for(int f=rho; f<NVAR; f++)
                  dwn1[fencode_cdf(p,i,j,f)]=0; 
        }
               __syncthreads();

//if(i>20 && j >20 && i<90 && j<90)
//	{
//               computepk_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
//              computept_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
//}
//              __syncthreads();
#ifdef USE_VAC
 if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
                    computej_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
#endif

#ifdef USE_SAC
 if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
                    computej_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
#endif

__syncthreads();


  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
  if(i<((p->n[0])) && j<((p->n[1])))
	{		               
             #ifdef ADIABHYDRO
               computepk_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
               computept_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
             #else
               //computej_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
               computepk_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
               computept_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);

               computebdotv_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
               computedivb_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);

             #endif

         }
              __syncthreads();
  if(i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
        }
              __syncthreads();

   if( i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{ 
               computecmax_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);




	}
 __syncthreads();

 /*if(i<(p->n[0]) && j<(p->n[1]))
	{ 
              // for(int f=vel1; f<NDERV; f++)
              for(int f=current1; f<=current2; f++)
                  bc_cont_cdf(wd,p,i,j,f);

	}
 __syncthreads();*/
  
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cdf(char *label)
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




int cucomputedervfields(struct params **p, real **w,  struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;

//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
     computedervfields_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

     //following used for testing to check current soundspeeds etc
     //cudaMemcpy(*w, *d_wd, 7*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}






