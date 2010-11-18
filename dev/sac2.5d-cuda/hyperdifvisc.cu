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
#include "gradops_hdv.cuh"

__global__ void hyperdifvisc_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order, real *wtemp, int field, int dim)
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
  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);


    //set viscosities
   if(i<((p->n[0])) && j<((p->n[1])))
   {
        for(int f=tmp1; f<=tmp9; f++)
                 wtemp[fencode_hdv(p,i,j,f)]=0;


        //temp value for viscosity

#ifdef USE_SAC
        wtemp[fencode_hdv(p,i,j,tmp1)]=wmod[fencode_hdv(p,i,j,field)+order*NVAR*(p->n[0])*(p->n[1])]/((field==rho || field>mom3)+(field>rho && field<energy)*(wmod[fencode_hdv(p,i,j,rho)+order*NVAR*(p->n[0])*(p->n[1])]+wmod[fencode_hdv(p,i,j,rhob)]+order*NVAR*(p->n[0])*(p->n[1])));
        if(field=rho)
           wtemp[fencode_hdv(p,i,j,tmp1)]+=wmod[fencode_hdv(p,i,j,rhob)+order*NVAR*(p->n[0])*(p->n[1])];

       if(field=b1 || field==b2)
           wtemp[fencode_hdv(p,i,j,tmp1)]+=wmod[fencode_hdv(p,i,j,field+5)+order*NVAR*(p->n[0])*(p->n[1])];
#else
        wtemp[fencode_hdv(p,i,j,tmp1)]=wmod[fencode_hdv(p,i,j,field)+order*NVAR*(p->n[0])*(p->n[1])]/( (field==rho || field>mom3)+(field>rho && field<energy)*wmod[fencode_hdv(p,i,j,rho)+order*NVAR*(p->n[0])*(p->n[1])] );
#endif
        wd[fencode_hdv(p,i,j,hdnur)]=0;
        wd[fencode_hdv(p,i,j,hdnul)]=0;
   }

   __syncthreads();


   //boundaries
     if(i<((p->n[0])) && j<((p->n[1])))            
               {                                                      
                if(dim==0)
                {
		        if(i==0 )
		          wtemp[fencode_hdv(p,i,j,tmp1)]=wtemp[fencode_hdv(p,4,j,tmp1)];
		        if( (i==((p->n[0])-1)) )
		          wtemp[fencode_hdv(p,i,j,tmp1)]=wtemp[fencode_hdv(p,((p->n[0])-4),j,tmp1)];
                }

                if(dim==1)
                {
		        if(j==0 )
		          wtemp[fencode_hdv(p,i,j,tmp1)]=wtemp[fencode_hdv(p,i,4,tmp1)];
		        if( (j==((p->n[1])-1)) )
		          wtemp[fencode_hdv(p,i,j,tmp1)]=wtemp[fencode_hdv(p,i,((p->n[1])-4),tmp1)];                                  
                }
               }
               
   __syncthreads();

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
 
   if(i>1 && j>1 && i<((p->n[0])-1) && j<((p->n[1])-1))            
   { 
           wtemp[fencode_hdv(p,i,j,tmp2)]=fabs(3.0*(wtemp[fencode_hdv(p,i+(dim==0),j+(dim==1),tmp1)] - wtemp[fencode_hdv(p,i,j,tmp1)] ) - (wtemp[fencode_hdv(p,i+2*(dim==0),j+2*(dim==1),tmp1)] - wtemp[fencode_hdv(p,i-(dim==0),j-(dim==1),tmp1)]    ));



           wtemp[fencode_hdv(p,i,j,tmp3)]=fabs((wtemp[fencode_hdv(p,i+(dim==0),j+(dim==1),tmp1)] - wtemp[fencode_hdv(p,i,j,tmp1)] ));
   }
   __syncthreads();



  //compute md3r and md1r
//tmp4    md3r
//tmp5    md1r
   if(i>1 && j>1 && i<((p->n[0])-1) && j<((p->n[1])-1))            
   {
         maxt=0;
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                {
                   if(wtemp[fencode_hdv(p,i+is,j+js,tmp2)]>maxt)
                         wtemp[fencode_hdv(p,i+is,j+js,tmp2)]=maxt;

                }
          wtemp[fencode_hdv(p,i,j,tmp4)]=maxt;

         maxt=0;
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                {
                   if(wtemp[fencode_hdv(p,i+is,j+js,tmp3)]>maxt)
                        maxt=wtemp[fencode_hdv(p,i+is,j+js,tmp3)];

                }
          wtemp[fencode_hdv(p,i,j,tmp5)]=maxt;
   }



  //compute d3l and d1l
//tmp6    d3l
//tmp7    d1l






   if(i>1 && j>1 && i<((p->n[0])-1) && j<((p->n[1])-1))            
   { 
           wtemp[fencode_hdv(p,i,j,tmp6)]=fabs(3.0*(wtemp[fencode_hdv(p,i,j,tmp1)] - wtemp[fencode_hdv(p,i-(dim==0),j-(dim==0),tmp1)] - wtemp[fencode_hdv(p,i+(dim==0),j+(dim==1),tmp1)] - wtemp[fencode_hdv(p,i-2*(dim==0),j-2*(dim==1),tmp1)]    ));
           wtemp[fencode_hdv(p,i,j,tmp7)]=fabs((wtemp[fencode_hdv(p,i,j,tmp1)] - wtemp[fencode_hdv(p,i-(dim==0),j-(dim==1),tmp1)] ));
   }
   __syncthreads();



  //compute md3l and md1l
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l

   if(i>1 && j>1 && i<((p->n[0])-1) && j<((p->n[1])-1))            
   {
         maxt=0;
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                {
                   if(wtemp[fencode_hdv(p,i+is,j+js,tmp6)]>maxt)
                         maxt=wtemp[fencode_hdv(p,i+is,j+js,tmp6)];

                }
          wtemp[fencode_hdv(p,i,j,tmp8)]=maxt;

         maxt=0;
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                {
                   if(wtemp[fencode_hdv(p,i+is,j+js,tmp7)]>maxt)
                        maxt=wtemp[fencode_hdv(p,i+is,j+js,tmp7)];

                }
          wtemp[fencode_hdv(p,i,j,tmp9)]=maxt;
   }
 __syncthreads();

   p->maxviscoef=0;
    //finally update nur and nul
//tmp4    md3r
//tmp5    md1r

   if(i<((p->n[0])) && j<((p->n[1])))
   {
     if(wtemp[fencode_hdv(p,i,j,tmp5)]>0)
	wd[fencode_hdv(p,i,j,hdnur)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1]))*(p->cmax)*(p->chyp)*wtemp[fencode_hdv(p,i,j,tmp4)]/wtemp[fencode_hdv(p,i,j,tmp5)];
     else
        wd[fencode_hdv(p,i,j,hdnur)]=0;
   }
 __syncthreads();
   if(i<((p->n[0])) && j<((p->n[1])))
   {
       if(wd[fencode_hdv(p,i,j,hdnur)]>(p->maxviscoef))
          p->maxviscoef=wd[fencode_hdv(p,i,j,hdnur)];
   }

//tmp8    md3l
//tmp9    md1l
   if(i<((p->n[0])) && j<((p->n[1])))
   {
     if(wtemp[fencode_hdv(p,i,j,tmp9)]>0)
	wd[fencode_hdv(p,i,j,hdnul)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1]))*(p->cmax)*(p->chyp)*wtemp[fencode_hdv(p,i,j,tmp8)]/wtemp[fencode_hdv(p,i,j,tmp9)];
     else
        wd[fencode_hdv(p,i,j,hdnul)]=0;
   }
 __syncthreads();
   if(i<((p->n[0])) && j<((p->n[1])))
   {
       if(wd[fencode_hdv(p,i,j,hdnul)]>(p->maxviscoef))
          p->maxviscoef=wd[fencode_hdv(p,i,j,hdnul)];
   }
  __syncthreads();
 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv(char *label)
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





int cuhyperdifvisc(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, int field, int dim)
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
     hyperdifvisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order, *d_wtemp, field, dim);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_hdv,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_hdv, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}







