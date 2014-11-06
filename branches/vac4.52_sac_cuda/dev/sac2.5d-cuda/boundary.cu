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
#include "gradops_b.cuh"

__global__ void boundary_parallel(struct params *p, real *w, real *wnew, real *wd, real *wmod, int order)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


/*#ifdef ADIABHYDRO
;
#else
  //This second call makes sure corners are set correctly
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {
             for( f=rho; f<=b2; f++)
{
     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  if(i<p->n[0] && j<p->n[1])

                  bc_periodic1a_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);

 __syncthreads();
   } 
}

#endif*/




   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

               for( f=rho; f<=b2; f++)
               {
  if(i<p->n[0] && j<p->n[1])
	{

               //default continuous BC for all
               //gradient kept zero by copying variable values from edge of mesh to ghost cells
                //  bc_cont_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,rho);
               
                 // bc_fixed_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,rho,1.0);
               //   bc_fixed(wnew,p,i,j,rho,1.0);
               //   bc_periodic_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,rho);
               

               

#ifdef ADIABHYDRO
                  bc_cont_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);
#else
                 bc_periodic1_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);
#endif                

                //  bc_fixed_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f,0.0);
                 // bc_fixed(wnew,p,i,j,f,val);

                 // bc_periodic1_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);
                




               /*for(int f=vel1; f<NDERV; f++)
               {
                  bc_cont_b(wd,p,i,j,f);

                 //bc_fixed_b(wd,p,i,j,f,0.0);
                 //   bc_periodic(wd,p,i,j,f);

                  
               }*/

	}
 __syncthreads();
               }
}


#ifdef ADIABHYDRO
;
#else
  //This second call makes sure corners are set correctly
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {
             for( f=rho; f<=b2; f++)
{
     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  if(i<p->n[0] && j<p->n[1])

                  bc_periodic2_b(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);

 __syncthreads();
   } 
}

#endif



  
}

int cuboundary(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, int order)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   //int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;
//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
 	    //printf("called prop\n"); 
    // cudaThreadSynchronize();
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wd, *d_wmod, order);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
	    //printf("called update\n"); 
    cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}

