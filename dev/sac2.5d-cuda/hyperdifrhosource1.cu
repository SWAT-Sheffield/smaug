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
#include "gradops_hdr1.cuh"
__global__ void hyperdifrhosource2_parallel(struct params *p,  real *wmod, 
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
  int ii1,ii0;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  real rdx;

   int ip,jp,ipg,jpg;
  int ii[NDIM];
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
   
     #ifdef USE_SAC_3D
	  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1)+(p->dx[2])*(dim==2));
	#else
	  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1));
	#endif
  
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  {
     



dwn1[fencode3_hdr1(p,ii,field)]=( wd[fencode3_hdr1(p,ii,hdnur)] * wtemp[fencode3_hdr1(p,ii,tmp1)] - wd[fencode3_hdr1(p,ii,hdnul)] *wtemp[fencode3_hdr1(p,ii,tmp2)]            )/rdx;

                              wmod[fencode3_hdr1(p,ii,field)+(ordero*NVAR*dimp)]=wmod[fencode3_hdr1(p,ii,field)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdr1(p,ii,field)]; 
  }
}
__syncthreads();




 
}



__global__ void hyperdifrhosource1_parallel(struct params *p,  real *wmod, 
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
  int ii1,ii0;
  real fip,fim1,tmpc;
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
  real rdx;

   int ip,jp,ipg,jpg;
  int ii[NDIM];
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

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  //init rhol and rhor
  //if(i<((p->n[0])) && j<((p->n[1])))
  {
    //for(int f=tmp1; f<=tmprhor; f++)	
    //    wtemp[fencode_hdr1(p,i,j,f)]=0.0;
    dwn1[fencode3_hdr1(p,ii,field)]=0.0;
    wtemp[fencode3_hdr1(p,ii,tmp1)]=0.0;
    wtemp[fencode3_hdr1(p,ii,tmp2)]=0.0;
    //wtemp[fencode_hdr1(p,i,j,tmp3)]=0.0;
   }
}
 __syncthreads();

     #ifdef USE_SAC_3D
	  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1)+(p->dx[2])*(dim==2));
	#else
	  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1));
	#endif

 
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1) && k<((p->n[2])-1))
     #else
       if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif

  
  //if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {
     

    wtemp[fencode3_hdr1(p,ii,tmp1)]=grad1r3_hdr1(wmod+shift,p,ii,rho,dim);
    wtemp[fencode3_hdr1(p,ii,tmp2)]=grad1l3_hdr1(wmod+shift,p,ii,rho,dim);
  }
}
__syncthreads();




 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdr1(char *label)
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





int cuhyperdifrhosource1(struct params **p, struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero,real **d_wtemp, int field, int dim, real dt)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifrhosource1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
     cudaThreadSynchronize();
    hyperdifrhosource2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,dt);
     cudaThreadSynchronize();


}







