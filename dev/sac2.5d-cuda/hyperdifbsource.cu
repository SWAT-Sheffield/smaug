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
#include "gradops_hdb.cuh"














__global__ void hyperdifbsource_parallel(struct params *p, real *w, real *wnew, real *wmod, 
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
  int m,ii,jj,ii1,ii0;
  real fip,fim1,tmp2,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real sb;
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);

  //init rhol and rhor
  if(i<((p->n[0])) && j<((p->n[1])))
    for(int f=tmp1; f<=tmprhor; f++)	
        wtemp[fencode_hdb(p,i,j,f)]=0.0;

 __syncthreads();

  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
  {
       wtemp[fencode_hdb(p,i,j,tmprhor)]=(wmod[fencode_hdb(p,i,j,rho)]+wmod[fencode_hdb(p,i+(field==0),j+(field==1),rho)])/2;
       wtemp[fencode_hdb(p,i,j,tmprhol)]=(wmod[fencode_hdb(p,i,j,rho)]+wmod[fencode_hdb(p,i-(field==0),j+(field==1),rho)])/2;

   }
__syncthreads();


  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
  {
     wtemp[fencode_hdb(p,i,j,tmp1)]=wmod[fencode_hdb(p,i,j,mom1+field)]/wmod[fencode_hdb(p,i,j,rho)];
     wtemp[fencode_hdb(p,i,j,tmp2)]=grad1_hdb(wtemp,p,i,j,tmp1,0);
     wtemp[fencode_hdb(p,i,j,tmp3)]=grad1_hdb(wtemp,p,i,j,tmp1,1);
  }

__syncthreads();


 



  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
             //ii1=0
             //case i=k, ii0=l
                      //   ii0=k
           //   m=l
           //   sB=-1.d0
           //   j=k

             ii1=0;
             ii0=dim;
             m=field;
             jj=dim;
             sb=-1.0;

             

                     fip=wmod[fencode_hdb(p,i+(field==0),j+(field==1),rho)]*((field==0)*wtemp[fencode_hdb(p,i+(field==0),j+(field==1),tmp2)] + (field==1)*wtemp[fencode_hdb(p,i+(field==0),j+(field==1),tmp3)])*(wtemp[fencode_hdb(p,i+(field==0),j+(field==1),hdnur)]+wtemp[fencode_hdb(p,i+(field==0),j+(field==1),hdnul)])/4.0;




                     fim1=wmod[fencode_hdb(p,i-(field==0),j-(field==1),rho)]*((field==0)*wtemp[fencode_hdb(p,i-(field==0),j-(field==1),tmp2)] + (field==1)*wtemp[fencode_hdb(p,i-(field==0),j-(field==1),tmp3)])*(wtemp[fencode_hdb(p,i-(field==0),j-(field==1),hdnur)]+wtemp[fencode_hdb(p,i-(field==0),j-(field==1),hdnul)]);
                     
		     //dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]+(evalgrad1_hdb(fip, fim1, p,field))/(((p->dx[0])*(field==0))+(p->dx[1])*(field==1));
                      //dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+field)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]-(p->chyp)*grad2_hdb(wmod,p,i,i,mom1+field,dim);

//dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1)]-(p->chyp)*grad2_hdb(wmod,p,i,i,mom1+field,dim);

             ii1=1;
             ii0=field;
             m=dim;
             jj=field;
             sb=1.0;

		     //dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]+( wtemp[fencode_hdb(p,i,j,tmprhor)]*(wd[fencode_hdb(p,i,j,hdnur)]*grad1r_hdb(wtemp,p,i,j,tmp1,field))-wtemp[fencode_hdb(p,i,j,tmprhol)]*(wd[fencode_hdb(p,i,j,hdnul)]*grad1l_hdb(wtemp,p,i,j,tmp1,field)) )/(((p->dx[0])*(field==0))+(p->dx[1])*(field==1));

dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,field)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,field)]+(evalgrad1_hdb(fip, fim1, p,field))/(((p->dx[0])*(field==0))+(p->dx[1])*(field==1));



         /*    for(ii1=0; ii1<2; i++)
             {
		     if(ii1==0)
                     {
                        ii0=field;
                        ii=dim;
                     }
                     else
                     {
			ii=field;
                        ii0=dim;
                     }

		     if(ii==field)
		     {
		     ;//dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]+( wtemp[fencode_hdb(p,i,j,tmprhor)]*(wd[fencode_hdb(p,i,j,hdnur)]*grad1r_hdb(wtemp,p,i,j,tmp1,field))-wtemp[fencode_hdb(p,i,j,tmprhol)]*(wd[fencode_hdb(p,i,j,hdnul)]*grad1l_hdb(wtemp,p,i,j,tmp1,field)) )/(((p->dx[0])*(field==0))+(p->dx[1])*(field==1));
		     }
		     else
		     {

                     fip=wmod[fencode_hdb(p,i+(field==0),j+(field==1),rho)]*((field==0)*wtemp[fencode_hdb(p,i+(field==0),j+(field==1),tmp2)] + (field==1)*wtemp[fencode_hdb(p,i+(field==0),j+(field==1),tmp3)])*(wtemp[fencode_hdb(p,i+(field==0),j+(field==1),hdnur)]+wtemp[fencode_hdb(p,i+(field==0),j+(field==1),hdnul)])/4.0;




                     fim1=wmod[fencode_hdb(p,i-(field==0),j-(field==1),rho)]*((field==0)*wtemp[fencode_hdb(p,i-(field==0),j-(field==1),tmp2)] + (field==1)*wtemp[fencode_hdb(p,i-(field==0),j-(field==1),tmp3)])*(wtemp[fencode_hdb(p,i-(field==0),j-(field==1),hdnur)]+wtemp[fencode_hdb(p,i-(field==0),j-(field==1),hdnul)]);
                     
		     ;//dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]=dwn1[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_hdb(p,i,j,mom1+ii0)]+(evalgrad1_hdb(fip, fim1, p,field))/(((p->dx[0])*(field==0))+(p->dx[1])*(field==1));
		     
}
             }*/


 
               
 


	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdb(char *label)
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





int cuhyperdifbsource(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, int field, int dim)
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
     hyperdifbsource_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order,*d_wtemp, field, dim);
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







