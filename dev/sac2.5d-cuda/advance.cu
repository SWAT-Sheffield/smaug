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


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;

  int index,i,j,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];

    int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  


int shift=order*NVAR*dimp;


   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     iia[0]=ip*(p->npgp[0])+ipg;
     iia[1]=jp*(p->npgp[1])+jpg;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp*(p->npgp[2])+kpg;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
	{		               
 
               float big=9999.0;
               for(int f=rho; f<NVAR; f++)
               {
                  
                   
                  if((p->rkon)==1)
                  {
                    switch(order)
                     {
                        case 0:
                       wmod[fencode3_adv(p,iia,f)+(2*dimp*NVAR)]=wmod[fencode3_adv(p,iia,f)];

                         break;
                        case 1:
                       wmod[fencode3_adv(p,iia,f)+(3*dimp*NVAR)]=wmod[fencode3_adv(p,iia,f)];
 
                         break;
                        case 2:
                       wmod[fencode3_adv(p,iia,f)+(dimp*NVAR)]=(wmod[fencode3_adv(p,iia,f)+(dimp*NVAR)]+2.0*wmod[fencode3_adv(p,iia,f)+(2*dimp*NVAR)]+wmod[fencode3_adv(p,iia,f)+(3*dimp*NVAR)]-4.0*wmod[fencode3_adv(p,iia,f)])/3;


                         break;
                        case 3:

                        wmod[fencode3_adv(p,iia,f)]=wmod[fencode3_adv(p,iia,f)]+wmod[fencode3_adv(p,iia,f)+(dimp*NVAR)];

                         break;

                     }
                   }
                  else
                  {
                  //if((dwn1[fencode3_adv(p,iia,f)]<(big/100)) && ( dwn1[fencode3_adv(p,iia,f)]>(-big/100)) )
                  //  if( j!=2)
                       //wmod[fencode3_adv(p,iia,f)]=wmod[fencode3_adv(p,iia,f)+(order*(p->n[0])*(p->n[1])*NVAR)];
                      wmod[fencode3_adv(p,iia,f)]=wmod[fencode3_adv(p,iia,f)+(dimp*NVAR)];
                   //lax-friedrichs
                  //wmod[fencode3_adv(p,iia,f)]=((w[fencode3_adv(p,i+1,j,f)]+w[fencode3_adv(p,i-1,j,f)]+w[fencode3_adv(p,iia+1,f)]+w[fencode3_adv(p,iia-1,f)])/4.0)+(dt)*(dwn1[fencode3_adv(p,iia,f)]);
                   }
                  
                   if(isnan(wmod[fencode3_adv(p,iia,f)])) wmod[fencode3_adv(p,iia,f)]=w[fencode3_adv(p,iia,f)];
                   if(wmod[fencode3_adv(p,iia,f)]>big)
                           wmod[fencode3_adv(p,iia,f)]=w[fencode3_adv(p,iia,f)];
                   if(wmod[fencode3_adv(p,iia,f)]<-big)
                           wmod[fencode3_adv(p,iia,f)]=w[fencode3_adv(p,iia,f)];

                     if(f==rho)
                            if(wmod[fencode3_adv(p,iia,f)]<0)
                               wmod[fencode3_adv(p,iia,f)]=1.00;
               }



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



