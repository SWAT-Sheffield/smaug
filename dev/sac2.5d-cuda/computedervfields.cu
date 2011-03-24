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


__global__ void computedervfields_parallel(struct params *p,   real *wmod, real *wd, int order)
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


  

   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


if(order == 0)
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;




if(i<((p->n[0])) && j<((p->n[1])))
	{		
 
               for(int f=rho; f<=b2; f++)
                  wmod[fencode_cdf(p,i,j,f)+((p->n[0]))*((p->n[1]))*NVAR]=wmod[fencode_cdf(p,i,j,f)]; 
        }

}
               __syncthreads();

/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
if(i<((p->n[0])) && j<((p->n[1])))
	{		

               for(int f=vel1; f<NDERV; f++)
                 ;// wd[fencode_cdf(p,i,j,f)]=0; 
               for(int f=rho; f<NVAR; f++)
                 ;// dwn1[fencode_cdf(p,i,j,f)]=0; 
        }

}
               __syncthreads();*/

//if(i>20 && j >20 && i<90 && j<90)
//	{
//               computepk_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
//              computept_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
//}
//              __syncthreads();


/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
#ifdef USE_VAC
 if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
                    computej_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
#endif

#ifdef USE_SAC
 if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
                    computej_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
#endif

}
__syncthreads(); */


  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
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
               //computedivb_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);

             #endif

         }

}
              __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  if(i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);
               p->cmax=0.0;
        }

}
              __syncthreads();

if(iindex==0)
{
 //  for(ipg=0;ipg<(p->npgp[0]);ipg++)
 //  for(jpg=0;jpg<(p->npgp[1]);jpg++)
  // {

  //   i=ip*(p->npgp[0])+ipg;
 //    j=jp*(p->npgp[1])+jpg;
   //if( i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
    p->cmax=0.0;
    for(i>1;i<((p->n[0])-2);i++)
      for(j>1;j<((p->n[1])-2);j++)
	{ 
               computecmax_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);




	}

 //  }
}
 __syncthreads(); 




  /* for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
     i=2*i;
     j=2*j;

  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
               computecmax_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);

	

   }

 __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
     i=2*i+1;
     j=2*j+1;

  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
               computecmax_cdf(wmod+(order*((p->n[0]))*((p->n[1]))*NVAR),wd,p,i,j);

	

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




int cucomputedervfields(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order)
{


 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;


     computedervfields_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}






