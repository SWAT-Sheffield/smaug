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
#include "gradops_hdv1a.cuh"

__device__ __host__
void bc_periodic1_temp2_hdv1a(real *wt, struct params *p,int i, int j, int f) {

                if(i==1 )                
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,6,j,f)];
                else if((i==((p->n[0]))) )                
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i-4,j,f)];
                else if(j==1  )                
                  wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,6,f)];
                else if((j==((p->n[1]))) )                
                  wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,j-4,f)];
}

__device__ __host__
void bc_periodic2_temp2_hdv1a(real *wt, struct params *p,int i, int j, int f) {


               if(i<1 && j<1)
                {
                  if(i==j)
                    //wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,6,f)];
                  else                  
                    //wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,6,j,f)];                                    
                }
                else if(i<1 && j>((p->n[1])-1))
                {
                  if(i==(j-(p->n[1])-1))                  
                    //wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,6,j,f)];                                     
                  else                  
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,j-6,f)];                                     
                }
                else if(i>((p->n[0])-1) && j<1)
                {
                  if((i-(p->n[0])+1)==j)                  
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i-5,j,f)];                                    
                  else                  
                   wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,4,f)];                                    
                }
                else if(i>((p->n[0])-1) && j>((p->n[1])-1))
                {
                  if(i==j)                  
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i,j-5,f)];                                    
                  else                  
                    wt[fencode_hdv1a(p,i,j,f)]=wt[fencode_hdv1a(p,i-5,j,f)];                                    
                }                       
                 
                




}



__global__ void hyperdifvisc1a_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  int i,j,iv;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
//int numBlocks = (ni*nj+numThreadsPerBlock-1) / numThreadsPerBlock;
  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*(p->n[0])*(p->n[1]);
  __shared__ real wts[512];
  __shared__ real wms[512];



   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if( i<((p->n[0])) && j<((p->n[1])))
   {
     //tmp6 is tmp_nuI
     wtemp2[fencode_hdv1a(p,i+1,j+1,tmpnui)]=wtemp[fencode_hdv1a(p,i,j,tmp6)];

   }

   }
   __syncthreads();








 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv1a(char *label)
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





int cuhyperdifvisc1a(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand)
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
     hyperdifvisc1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_hdv1a,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv1a,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv1a,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_hdv1a, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}







