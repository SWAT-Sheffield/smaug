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
#include "gradops_hdv2.cuh"




__global__ void hyperdifvisc2_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
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

  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*(p->n[0])*(p->n[1]);


   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l




   //tmp1  tmp_nuI
 
//compute d3r and d1r
   //tmp2  d3r
    //tmp3 d1r

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
 
   if(i>1 && j>1 && i<((p->n[0])) && j<((p->n[1])))
   //if(i<((p->n[0])-1) && j<((p->n[1])-1))            
   { 
     if(hand==1)
     {
           wtemp1[fencode_hdv2(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv2(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i,j,tmpnui)] ) - (wtemp2[fencode_hdv2(p,i+2*(dim==0),j+2*(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i-(dim==0),j-(dim==1),tmpnui)]    ));
     }
     else
     {
          // wtemp1[fencode_hdv2(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv2(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i,j,tmpnui)]) - (wtemp2[fencode_hdv2(p,i+2*(dim==0),j+2*(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i-(dim==0),j-(dim==1),tmpnui)]    ));
           wtemp1[fencode_hdv2(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv2(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i-(dim==0),j-(dim==1),tmpnui)]) - (wtemp2[fencode_hdv2(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i-2*(dim==0),j-2*(dim==1),tmpnui)]    ));
     }
   }
}
   __syncthreads();
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

   if(i>1 && j>1 && i<((p->n[0])+1) && j<((p->n[1])+1))            
   { 
     if(hand==1)
     {

           wtemp1[fencode_hdv2(p,i,j,d1)]=fabs((wtemp2[fencode_hdv2(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv2(p,i,j,tmpnui)] ));
     }
     else
     {
           wtemp1[fencode_hdv2(p,i,j,d1)]=fabs((wtemp2[fencode_hdv2(p,i,j,tmpnui)] - wtemp2[fencode_hdv2(p,i-(dim==0),j-(dim==1),tmpnui)] ));
     }
   }
}
   __syncthreads();



/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic1_hdv2(wtemp,p,i,j,tmp2);
                  bc_periodic1_hdv2(wtemp,p,i,j,tmp3);
             }

}
                __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic2_hdv2(wtemp,p,i,j,tmp2);
                  bc_periodic2_hdv2(wtemp,p,i,j,tmp3);
             }
}
                __syncthreads();*/








 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv2(char *label)
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





int cuhyperdifvisc2(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand)
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
     hyperdifvisc2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_hdv2,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv2,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv2,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_hdv2, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}







