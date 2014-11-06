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
#include "gradops_hdm3.cuh"



__global__ void hyperdifmomsource3_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real rdx;
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1));

   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
//if(i<((p->n[0])) && j<((p->n[1])))
  if(i<((p->n[0])) && j<((p->n[1])))
	{		               

dwn1[fencode_hdm3(p,i,j,energy)]=wtemp[fencode_hdm3(p,i,j,tmp6)]*wd[fencode_hdm3(p,i,j,hdnur)]*wtemp[fencode_hdm3(p,i,j,tmp8)]-wtemp[fencode_hdm3(p,i,j,tmp5)]*wd[fencode_hdm3(p,i,j,hdnul)]*wtemp[fencode_hdm3(p,i,j,tmp7)]/(rdx)/2;

dwn1[fencode_hdm3(p,i,j,mom1+ii0)]=(wtemp[fencode_hdm3(p,i,j,tmp3)]*wd[fencode_hdm3(p,i,j,hdnur)]*wtemp[fencode_hdm3(p,i,j,tmp8)]-wtemp[fencode_hdm3(p,i,j,tmp2)]*wd[fencode_hdm3(p,i,j,hdnul)]*wtemp[fencode_hdm3(p,i,j,tmp7)])/(rdx)/2;

/*dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,energy)]=(

(  wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+field)]+wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i+(dim==0),j+(dim==1),mom1+field)])*wd[fencode_hdm3(p,i,j,hdnur)]*grad1r_hdm3(wtemp,p,i,j,tmp1,dim)

-(wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+field)]+wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i-(dim==0),j+(dim==1),mom1+field)])*wd[fencode_hdm3(p,i,j,hdnul)]*grad1l_hdm3(wtemp,p,i,j,tmp1,dim))/(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1))/2;

dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+ii0)]=(wtemp[fencode_hdm3(p,i,j,tmp2)]*wd[fencode_hdm3(p,i,j,hdnur)]*grad1r_hdm3(wtemp,p,i,j,tmp1,dim)-wtemp[fencode_hdm3(p,i,j,tmp3)]*wd[fencode_hdm3(p,i,j,hdnul)]*grad1l_hdm3(wtemp,p,i,j,tmp1,dim))/(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1))/2;*/


/*dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,energy)]=(

(  wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+field)]+wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i+(dim==0),j+(dim==1),mom1+field)])*wd[fencode_hdm3(p,i,j,hdnur)]*grad1r_hdm3(wtemp,p,i,j,tmp1,dim)

-(wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+field)]+wmod[(order*NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i-(dim==0),j+(dim==1),mom1+field)])*wd[fencode_hdm3(p,i,j,hdnul)]*grad1l_hdm3(wtemp,p,i,j,tmp1,dim))/2;

dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_hdm3(p,i,j,mom1+ii0)]=(wtemp[fencode_hdm3(p,i,j,tmp2)]*wd[fencode_hdm3(p,i,j,hdnur)]*grad1r_hdm3(wtemp,p,i,j,tmp1,dim)-wtemp[fencode_hdm3(p,i,j,tmp3)]*wd[fencode_hdm3(p,i,j,hdnul)]*grad1l_hdm3(wtemp,p,i,j,tmp1,dim))/2;*/


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

                  bc_periodic1_hdm3(dwn1,p,i,j,mom1+ii0);
                  bc_periodic1_hdm3(dwn1,p,i,j,energy);
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

                  bc_periodic2_hdm3(dwn1,p,i,j,mom1+ii0);
                  bc_periodic2_hdm3(dwn1,p,i,j,energy);
             }
}
                __syncthreads();*/



   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
			// if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                        if(i<((p->n[0])) && j<((p->n[1])))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode_hdm3(p,i,j,mom1+ii0)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hdm3(p,i,j,mom1+ii0)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hdm3(p,i,j,mom1+ii0)]; 
                             wmod[fencode_hdm3(p,i,j,energy)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hdm3(p,i,j,energy)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hdm3(p,i,j,energy)]; 

                         }
              //  }	
}
  __syncthreads();


  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdm3(char *label)
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





int cuhyperdifmomsource3(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt)
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
     hyperdifmomsource3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}







