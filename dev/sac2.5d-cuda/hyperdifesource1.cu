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
#include "gradops_hde1.cuh"



__global__ void hyperdifesource3_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, real dt)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii,ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
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

if(i<((p->n[0])) && j<((p->n[1])))
//  if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {


 
//dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wtemp[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             )/rdx;
//dwn1[fencode_hde1(p,i,j,field)]=( wd[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wd[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             );

//wtemp[fencode_hde1(p,i,j,tmp2)]= grad1r_hde1(wtemp,p,i,j,tmp1,dim) ;
//wtemp[fencode_hde1(p,i,j,tmp3)]= grad1l_hde1(wtemp,p,i,j,tmp1,dim) ;
dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] *wtemp[fencode_hde1(p,i,j,tmp3)] - wtemp[fencode_hde1(p,i,j,hdnul)] *wtemp[fencode_hde1(p,i,j,tmp2)])/rdx;

    wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hde1(p,i,j,field)]; 

  }
}
__syncthreads();



   
/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


			 //if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         if(i<((p->n[0])) && j<((p->n[1])))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hde1(p,i,j,field)]; 
//wmod[fencode_hde1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]=dwn1[fencode_hde1(p,i,j,f2)];
                              //dwn1[fencode_hde1(p,i,j,f)]=0;
                         }
              //  }	
}
  __syncthreads();*/



 
}


__global__ void hyperdifesource2a_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii,ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
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

if(i<((p->n[0])) && j<((p->n[1])))
//  if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {


 
//dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wtemp[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             )/rdx;
//dwn1[fencode_hde1(p,i,j,field)]=( wd[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wd[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             );

//wtemp[fencode_hde1(p,i,j,tmp2)]= grad1r_hde1(wtemp,p,i,j,tmp1,dim) ;
//wtemp[fencode_hde1(p,i,j,tmp3)]= grad1l_hde1(wtemp,p,i,j,tmp1,dim) ;
dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] *wtemp[fencode_hde1(p,i,j,tmp3)] - wtemp[fencode_hde1(p,i,j,hdnul)] *wtemp[fencode_hde1(p,i,j,tmp2)])/rdx;


  }
}
__syncthreads();



   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


			 //if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         if(i<((p->n[0])) && j<((p->n[1])))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hde1(p,i,j,field)]; 
//wmod[fencode_hde1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]=dwn1[fencode_hde1(p,i,j,f2)];
                              //dwn1[fencode_hde1(p,i,j,f)]=0;
                         }
              //  }	
}
  __syncthreads();



 
}

__global__ void hyperdifesource2_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii,ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
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

//if(i<((p->n[0])) && j<((p->n[1])))
  if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {


 
//dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wtemp[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             )/rdx;
//dwn1[fencode_hde1(p,i,j,field)]=( wd[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wd[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             );

wtemp[fencode_hde1(p,i,j,tmp2)]= grad1l_hde1(wtemp,p,i,j,tmp1,dim) ;
wtemp[fencode_hde1(p,i,j,tmp3)]= grad1r_hde1(wtemp,p,i,j,tmp1,dim) ;


  }
}
__syncthreads();



/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

if(i<((p->n[0])) && j<((p->n[1])))
//  if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {


 
//dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wtemp[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             )/rdx;
//dwn1[fencode_hde1(p,i,j,field)]=( wd[fencode_hde1(p,i,j,hdnur)] * grad1r_hde1(wtemp,p,i,j,tmp1,dim) - wd[fencode_hde1(p,i,j,hdnul)] *grad1l_hde1(wtemp,p,i,j,tmp1,dim)             );

//wtemp[fencode_hde1(p,i,j,tmp2)]= grad1r_hde1(wtemp,p,i,j,tmp1,dim) ;
//wtemp[fencode_hde1(p,i,j,tmp3)]= grad1l_hde1(wtemp,p,i,j,tmp1,dim) ;
dwn1[fencode_hde1(p,i,j,field)]=( wtemp[fencode_hde1(p,i,j,hdnur)] *wtemp[fencode_hde1(p,i,j,tmp3)] - wtemp[fencode_hde1(p,i,j,hdnul)] *wtemp[fencode_hde1(p,i,j,tmp2)])/rdx;


  }
}
__syncthreads();*/



   
 /*  for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


			 //if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         if(i<((p->n[0])) && j<((p->n[1])))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_hde1(p,i,j,field)+(ordero*NVAR*(p->n[0])*(p->n[1]))]+dt*dwn1[fencode_hde1(p,i,j,field)]; 
//wmod[fencode_hde1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]=dwn1[fencode_hde1(p,i,j,f2)];
                              //dwn1[fencode_hde1(p,i,j,f)]=0;
                         }
              //  }	
}
  __syncthreads();*/



 
}



__global__ void hyperdifesource1_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii,ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
int shift=order*NVAR*(p->n[0])*(p->n[1]);
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

  //init rhol and rhor
  if(i<((p->n[0])) && j<((p->n[1])))
  {
    for(int f=tmp1; f<=tmp8; f++)	
        wtemp[fencode_hde1(p,i,j,f)]=0.0;
    dwn1[fencode_hde1(p,i,j,field)]=0.0;
   }
}
 __syncthreads();

  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1));

   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

if(i<((p->n[0])) && j<((p->n[1])))
 // if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
  {

#ifdef USE_SAC
     wtemp[fencode_hde1(p,i,j,tmp1)]=wmod[shift+fencode_hde1(p,i,j,energy)]-0.5*(


(wmod[shift+fencode_hde1(p,i,j,b1)]*wmod[shift+fencode_hde1(p,i,j,b1)]+wmod[shift+fencode_hde1(p,i,j,b2)]*wmod[shift+fencode_hde1(p,i,j,b2)])

+((wmod[shift+fencode_hde1(p,i,j,mom1)]*wmod[shift+fencode_hde1(p,i,j,mom1)]+wmod[shift+fencode_hde1(p,i,j,mom2)]*wmod[shift+fencode_hde1(p,i,j,mom2)])/(wmod[shift+fencode_hde1(p,i,j,rho)]+wmod[shift+fencode_hde1(p,i,j,rhob)])));
#else

     wtemp[fencode_hde1(p,i,j,tmp1)]=wmod[shift+fencode_hde1(p,i,j,energy)]-0.5*((wmod[shift+fencode_hde1(p,i,j,b1)]*wmod[shift+fencode_hde1(p,i,j,b1)]+wmod[shift+fencode_hde1(p,i,j,b2)]*wmod[shift+fencode_hde1(p,i,j,b2)])

+((wmod[shift+fencode_hde1(p,i,j,mom1)]*wmod[shift+fencode_hde1(p,i,j,mom1)]+wmod[shift+fencode_hde1(p,i,j,mom2)]*wmod[shift+fencode_hde1(p,i,j,mom2)])/(wmod[shift+fencode_hde1(p,i,j,rho)]))
);

#endif
 


  }
}
__syncthreads();




 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hde1(char *label)
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





int cuhyperdifesource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim,real dt)
{


 dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifesource1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();

     hyperdifesource2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();

     hyperdifesource2a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();

     hyperdifesource3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,dt);
      cudaThreadSynchronize();


}







