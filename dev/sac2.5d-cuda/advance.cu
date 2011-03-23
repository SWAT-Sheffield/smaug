#include "cudapars.h"
#include "iotypes.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "step.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "gradops_adv.cuh"
#include "dervfields_adv.cuh"


__global__ void advance_parallel(struct params *p, real *wmod, real *w,  int order)
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
  real dy=p->dx[0];
  real dx=p->dx[1];
  
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
if( i<((p->n[0])) && j<((p->n[1])))
	{		               
 
               float big=9999.0;
               for(int f=rho; f<NVAR; f++)
               {
                  
                   
                  if((p->rkon)==1)
                  {
                     //wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)];
                    switch(order)
                     {
                        case 0:
                       wmod[fencode_adv(p,i,j,f)+(2*(p->n[0])*(p->n[1])*NVAR)]=wmod[fencode_adv(p,i,j,f)];
                       // wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)];

                         break;
                        case 1:
                       wmod[fencode_adv(p,i,j,f)+(3*(p->n[0])*(p->n[1])*NVAR)]=wmod[fencode_adv(p,i,j,f)];
                      // wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)+(2*(p->n[0])*(p->n[1])*NVAR)];

                         break;
                        case 2:
                       wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)]=(wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)]+2.0*wmod[fencode_adv(p,i,j,f)+(2*(p->n[0])*(p->n[1])*NVAR)]+wmod[fencode_adv(p,i,j,f)+(3*(p->n[0])*(p->n[1])*NVAR)]-4.0*wmod[fencode_adv(p,i,j,f)])/3;


                         break;
                        case 3:
                      // wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)]+wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)];
                        wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)]+wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)];

                         break;

                     }
                   }
                  else
                  {
                  //if((dwn1[fencode_adv(p,i,j,f)]<(big/100)) && ( dwn1[fencode_adv(p,i,j,f)]>(-big/100)) )
                  //  if( j!=2)
                       //wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)+(order*(p->n[0])*(p->n[1])*NVAR)];
                      wmod[fencode_adv(p,i,j,f)]=wmod[fencode_adv(p,i,j,f)+((p->n[0])*(p->n[1])*NVAR)];
                   //lax-friedrichs
                  //wmod[fencode_adv(p,i,j,f)]=((w[fencode_adv(p,i+1,j,f)]+w[fencode_adv(p,i-1,j,f)]+w[fencode_adv(p,i,j+1,f)]+w[fencode_adv(p,i,j-1,f)])/4.0)+(dt)*(dwn1[fencode_adv(p,i,j,f)]);
                   }
                  
                   if(isnan(wmod[fencode_adv(p,i,j,f)])) wmod[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];
                   if(wmod[fencode_adv(p,i,j,f)]>big)
                           wmod[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];
                   if(wmod[fencode_adv(p,i,j,f)]<-big)
                           wmod[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];

                     if(f==rho)
                            if(wmod[fencode_adv(p,i,j,f)]<0)
                               wmod[fencode_adv(p,i,j,f)]=1.00;
               }
               //computej_adv(wmod,wd,p,i,j);
               //computepk_adv(wmod,wd,p,i,j);
               //computept_adv(wmod,wd,p,i,j);


	}
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






int cuadvance(struct params **p, struct params **d_p,  real **d_wmod, real **d_w,  int order)
{

 dim3 dimBlock(dimblock, 1);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;

     advance_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wmod, *d_w, order);
     cudaThreadSynchronize();
}



