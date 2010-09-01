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

__device__ __host__
int dimproduct_adv (struct params *dp) {

  int tot=1;
  for(int i=0;i<NDIM;i++)
    tot*=dp->n[i];
  return tot; 
}


__device__ __host__
int encode_adv (struct params *dp,int ix, int iy) {

  return (iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int encode3_adv (struct params *dp,int ix, int iy, int iz) {

  return (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int fencode_adv (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return(( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1]))));
}


__device__ __host__
int fencode3_adv (struct params *dp,int ix, int iy, int iz, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return(  iz*((dp)->n[0])*((dp)->n[1])+ (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2]))  );
}
__device__ __host__
real evalgrad_adv(real fi, real fim1, real fip2, real fim2,struct params *p,int dir)
{
 //real valgrad;

   return(p->sodifon?((1.0/(2.0*(p->dx[dir])))*(fi-fim1)):((1.0/(12.0*(p->dx[dir])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 
}


__device__ __host__
real grad_adv(real *wmod,struct params *p,int *ix,int field,int dir)
{
 //real valgrad;

 if(dir == 0)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode_adv(p,ix[0]+1,ix[1],field)]-NVAR*wmod[fencode_adv(p,ix[0]-1,ix[1],field)]+wmod[fencode_adv(p,ix[0]-2,ix[1],field)]-wmod[fencode_adv(p,ix[0]+2,ix[1],field)])/6.0):wmod[fencode_adv(p,ix[0]+1,ix[1],field)]-wmod[fencode_adv(p,ix[0]-1,ix[1],field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode_adv(p,ix[0],ix[1]+1,field)]-NVAR*wmod[fencode_adv(p,ix[0],ix[1]-1,field)]+wmod[fencode_adv(p,ix[0],ix[1]-2,field)]-wmod[fencode_adv(p,ix[0],ix[1]+2,field)])/6.0):wmod[fencode_adv(p,ix[0],ix[1]+1,field)]-wmod[fencode_adv(p,ix[0],ix[1]-1,field)])/(2.0*(p->dx[1]))    );

 }

 return -1;
}

__device__ __host__
real grad3_adv(real *wmod,struct params *p,int *ix,int field,int dir)
{
 //real valgrad;

 if(dir == 0)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_adv(p,ix[0]+1,ix[1],ix[2],field)]-NVAR*wmod[fencode3_adv(p,ix[0]-1,ix[1],ix[2],field)]+wmod[fencode3_adv(p,ix[0]-2,ix[1],ix[2],field)]-wmod[fencode3_adv(p,ix[0]+2,ix[1],ix[2],field)])/6.0):wmod[fencode3_adv(p,ix[0]+1,ix[1],ix[2],field)]-wmod[fencode3_adv(p,ix[0]-1,ix[1],ix[2],field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_adv(p,ix[0],ix[1]+1,ix[2],field)]-NVAR*wmod[fencode3_adv(p,ix[0],ix[1]-1,ix[2],field)]+wmod[fencode3_adv(p,ix[0],ix[1]-2,ix[2],field)]-wmod[fencode3_adv(p,ix[0],ix[1]+2,ix[2],field)])/6.0):wmod[fencode3_adv(p,ix[0],ix[1]+1,ix[2],field)]-wmod[fencode3_adv(p,ix[0],ix[1]-1,ix[2],field)])/(2.0*(p->dx[1]))    );

 }
else if(dir == 2)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_adv(p,ix[0],ix[1],ix[2]+1,field)]-NVAR*wmod[fencode3_adv(p,ix[0],ix[1],ix[2]-1,field)]+wmod[fencode3_adv(p,ix[0],ix[1],ix[2]-2,field)]-wmod[fencode3_adv(p,ix[0],ix[1],ix[2]+2,field)])/6.0):wmod[fencode3_adv(p,ix[0],ix[1],ix[2]+1,field)]-wmod[fencode3_adv(p,ix[0],ix[1],ix[2]-1,field)])/(2.0*(p->dx[2]))    );

 }
 return -1;
}

__device__ __host__
void computej_adv(real *wmod,real *wd,struct params *p,int *ix)
{
  wd[fencode_adv(p,ix[0],ix[1],0)]=(grad_adv(wmod,p,ix,b3,1))/(p->mu);
  wd[fencode_adv(p,ix[0],ix[1],1)]=(grad_adv(wmod,p,ix,b3,0))/(p->mu);
  wd[fencode_adv(p,ix[0],ix[1],2)]=(grad_adv(wmod,p,ix,b2,0)-grad_adv(wmod,p,ix,b1,1))/(p->mu);
}

__device__ __host__
void computebdotv_adv(real *wmod,real *wd,struct params *p,int *ix)
{
wd[fencode_adv(p,ix[0],ix[1],bdotv)]=(wmod[fencode_adv(p,ix[0],ix[1],b1)]*wmod[fencode_adv(p,ix[0],ix[1],mom1)]+wmod[fencode_adv(p,ix[0],ix[1],b2)]*wmod[fencode_adv(p,ix[0],ix[1],mom2)]+wmod[fencode_adv(p,ix[0],ix[1],b3)]*wmod[fencode_adv(p,ix[0],ix[1],mom3)])/wmod[fencode_adv(p,ix[0],ix[1],rho)];
}

__device__ __host__
void computedivb_adv(real *wmod,real *wd,struct params *p,int *ix)
{
wd[fencode_adv(p,ix[0],ix[1],divb)]=grad_adv(wmod,p,ix,b1,0)+grad_adv(wmod,p,ix,b2,1);
}


__device__ __host__
void computepk_adv(real *wmod,real *wd,struct params *p,int *ix)
{
 // int status=0;

         #ifdef ADIABHYDRO
/*below used for adiabatic hydrodynamics*/
wd[fencode_adv(p,ix[0],ix[1],pressurek)]=(p->adiab)*pow(wmod[fencode_adv(p,ix[0],ix[1],rho)],p->gamma);


#else


 wd[fencode_adv(p,ix[0],ix[1],pressurek)]=((p->gamma)-1)*(wmod[fencode_adv(p,ix[0],ix[1],energy)]- 0.5*(wmod[fencode_adv(p,ix[0],ix[1],mom1)]*wmod[fencode_adv(p,ix[0],ix[1],mom1)]+wmod[fencode_adv(p,ix[0],ix[1],mom2)]*wmod[fencode_adv(p,ix[0],ix[1],mom2)]+wmod[fencode_adv(p,ix[0],ix[1],mom3)]*wmod[fencode_adv(p,ix[0],ix[1],mom3)])/wmod[fencode_adv(p,ix[0],ix[1],rho)]-0.5*(wmod[fencode_adv(p,ix[0],ix[1],b1)]*wmod[fencode_adv(p,ix[0],ix[1],b1)]+wmod[fencode_adv(p,ix[0],ix[1],b2)]*wmod[fencode_adv(p,ix[0],ix[1],b2)]+wmod[fencode_adv(p,ix[0],ix[1],b3)]*wmod[fencode_adv(p,ix[0],ix[1],b3)]) );

#endif


  if(wd[fencode_adv(p,ix[0],ix[1],pressurek)]<0)
              wd[fencode_adv(p,ix[0],ix[1],pressurek)]=0.001;


}
__device__ __host__
void computept_adv(real *wmod,real *wd,struct params *p,int *ix)
{
         #ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
wd[fencode_adv(p,ix[0],ix[1],pressuret)]=(p->adiab)*pow(wmod[fencode_adv(p,ix[0],ix[1],rho)],p->gamma);

#else


   wd[fencode_adv(p,ix[0],ix[1],pressuret)]=  wd[fencode_adv(p,ix[0],ix[1],pressurek)]+0.5*(wmod[fencode_adv(p,ix[0],ix[1],b1)]*wmod[fencode_adv(p,ix[0],ix[1],b1)]+wmod[fencode_adv(p,ix[0],ix[1],b2)]*wmod[fencode_adv(p,ix[0],ix[1],b2)]+wmod[fencode_adv(p,ix[0],ix[1],b3)]*wmod[fencode_adv(p,ix[0],ix[1],b3)]);

#endif

  if(wd[fencode_adv(p,ix[0],ix[1],pressuret)]<0)
              wd[fencode_adv(p,ix[0],ix[1],pressuret)]=0.001;

  //return ( status);
}

__global__ void advance_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd)
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


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
 
               float big=9999.0;
               for(int f=rho; f<=b3; f++)
               {
                   
                   
                  if((p->rkon)==1)
                  {
                  wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)]+(dt/6.0)*(dwn1[fencode_adv(p,i,j,f)]+2*dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]+2*dwn1[(2*NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]+dwn1[(3*NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]);
                 // wnew[fencode_adv(p,i,j,f)]=((w[fencode_adv(p,i+1,j,f)]+w[fencode_adv(p,i-1,j,f)]+w[fencode_adv(p,i,j+1,f)]+w[fencode_adv(p,i,j-1,f)])/4.0)+(dt/6.0)*(dwn1[fencode_adv(p,i,j,f)]+2*dwn1[(NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]+2*dwn1[(2*NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]+dwn1[(3*NVAR*(p->n[0])*(p->n[1]))+fencode_adv(p,i,j,f)]);

                   }
                  else
                  {
                  //if((dwn1[fencode_adv(p,i,j,f)]<(big/100)) && ( dwn1[fencode_adv(p,i,j,f)]>(-big/100)) )
                  //  if( j!=2)
                       wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)]+dt*dwn1[fencode_adv(p,i,j,f)];

                   //lax-friedrichs
                  //wnew[fencode_adv(p,i,j,f)]=((w[fencode_adv(p,i+1,j,f)]+w[fencode_adv(p,i-1,j,f)]+w[fencode_adv(p,i,j+1,f)]+w[fencode_adv(p,i,j-1,f)])/4.0)+(dt)*(dwn1[fencode_adv(p,i,j,f)]);
                   }
                  
                   if(isnan(wnew[fencode_adv(p,i,j,f)])) wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];
                   if(wnew[fencode_adv(p,i,j,f)]>big)
                           wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];
                   if(wnew[fencode_adv(p,i,j,f)]<-big)
                           wnew[fencode_adv(p,i,j,f)]=w[fencode_adv(p,i,j,f)];

                     if(f==rho)
                            if(wnew[fencode_adv(p,i,j,f)]<0)
                               wnew[fencode_adv(p,i,j,f)]=1.001;
               }
               //computej_adv(wnew,wd,p,i,j);
               //computepk_adv(wnew,wd,p,i,j);
               //computept_adv(wnew,wd,p,i,j);


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






int cuadvance(struct params **p, real **w, real **wnew,struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd)
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
     advance_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
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



