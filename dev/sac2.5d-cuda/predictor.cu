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



__device__ __host__
int encode_pre (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int fencode_pre (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
}

__device__ __host__
real evalgrad_pre(real fi, real fim1, real fip2, real fim2,struct params *p,int dir)
{
 //real valgrad_pre;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx[0])))*(fi-fim1)-(1.0/(12.0*(p->dx[0])))*(fip2-fim2);
   //return((1.0/(2.0*(p->dx[0])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[0])))*(fi-fim1)):((1.0/(12.0*(p->dx[0])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(fi-fim1)-(1.0/(12.0*(p->dx[1])))*(fip2-fim2);
     // return((2.0/(1.0*(p->dx[1])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[1])))*(fi-fim1)):((1.0/(12.0*(p->dx[1])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }

 return -1;
}


__device__ __host__
real grad_pre(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_pre;

 if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx[0])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx[0])))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
//return((1.0/(2.0*(p->dx[0])))*(wmod[fencode_pre(p,i+1,j,field)]-wmod[fencode_pre(p,i-1,j,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_pre(p,i+1,j,field)]-NVAR*wmod[fencode_pre(p,i-1,j,field)]+wmod[fencode_pre(p,i-2,j,field)]-wmod[fencode_pre(p,i+2,j,field)])/6.0):wmod[fencode_pre(p,i+1,j,field)]-wmod[fencode_pre(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dx[1])))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
// return((1.0/(2.0*(p->dx[1])))*(wmod[fencode_pre(p,i,j+1,field)]-wmod[fencode_pre(p,i,j-1,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_pre(p,i,j+1,field)]-NVAR*wmod[fencode_pre(p,i,j-1,field)]+wmod[fencode_pre(p,i,j-2,field)]-wmod[fencode_pre(p,i,j+2,field)])/6.0):wmod[fencode_pre(p,i,j+1,field)]-wmod[fencode_pre(p,i,j-1,field)])/(2.0*(p->dx[1]))    );  
}


 return 0;
}

__device__ __host__
void computej(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;

 // real dbzdy, dbydz;
 // real dbzdx, dbxdz;
 // real dbydx, dbxdy;

 // dbzdy=grad_pre(wmod,p,i,j,b3,1);
 // dbydz=0.0;
 // dbzdx=grad_pre(wmod,p,i,j,b3,0);
//  dbxdz=0.0;
 // dbydx=grad_pre(wmod,p,i,j,b2,0);
 // dbxdy=grad_pre(wmod,p,i,j,b1,1);

  wd[fencode_pre(p,i,j,0)]=(grad_pre(wmod,p,i,j,b3,1))/(p->mu);
  wd[fencode_pre(p,i,j,1)]=(grad_pre(wmod,p,i,j,b3,0))/(p->mu);
  wd[fencode_pre(p,i,j,2)]=(grad_pre(wmod,p,i,j,b2,0)-grad_pre(wmod,p,i,j,b1,1))/(p->mu);
  
  

 
  //return ( status);
}

__device__ __host__
void computebdotv(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;
 //real bsq=wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)];
//  wd[fencode_pre(p,i,j,4)]=  wd[fencode_pre(p,i,j,3)]+0.5*(wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)]);

wd[fencode_pre(p,i,j,bdotv)]=(wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,mom1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,mom2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,mom3)])/wmod[fencode_pre(p,i,j,rho)];
 // return ( status);
}

__device__ __host__
void computedivb(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;
 //real bsq=wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)];
//  wd[fencode_pre(p,i,j,4)]=  wd[fencode_pre(p,i,j,3)]+0.5*(wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)]);

wd[fencode_pre(p,i,j,divb)]=grad_pre(wmod,p,i,j,b1,0)+grad_pre(wmod,p,i,j,b2,1);
 // return ( status);
}


__device__ __host__
void computept(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
 wd[fencode_pre(p,i,j,pressuret)]=(p->adiab)*pow(wmod[fencode_pre(p,i,j,rho)],p->gamma);


#else

 //real bsq=wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)];
  wd[fencode_pre(p,i,j,pressuret)]=  wd[fencode_pre(p,i,j,pressurek)]+0.5*(wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)]);

#endif



  if(wd[fencode_pre(p,i,j,pressuret)]<0)
              wd[fencode_pre(p,i,j,pressuret)]=0.001;


 // return ( status);
}
__device__ __host__
void computepk(real *wmod,real *wd,struct params *p,int i,int j)
{
  //int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
wd[fencode_pre(p,i,j,pressurek)]=(p->adiab)*pow(wmod[fencode_pre(p,i,j,rho)],p->gamma);

#else

  //real momsq=wmod[fencode_pre(p,i,j,mom1)]*wmod[fencode_pre(p,i,j,mom1)]+wmod[fencode_pre(p,i,j,mom2)]*wmod[fencode_pre(p,i,j,mom2)]+wmod[fencode_pre(p,i,j,mom3)]*wmod[fencode_pre(p,i,j,mom3)];
  //real bsq=wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)];
  wd[fencode_pre(p,i,j,pressurek)]=((p->gamma)-1)*(wmod[fencode_pre(p,i,j,energy)]- 0.5*(wmod[fencode_pre(p,i,j,mom1)]*wmod[fencode_pre(p,i,j,mom1)]+wmod[fencode_pre(p,i,j,mom2)]*wmod[fencode_pre(p,i,j,mom2)]+wmod[fencode_pre(p,i,j,mom3)]*wmod[fencode_pre(p,i,j,mom3)])/wmod[fencode_pre(p,i,j,rho)]-0.5*(wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)]) );


#endif






  if(wd[fencode_pre(p,i,j,pressurek)]<0)
              wd[fencode_pre(p,i,j,pressurek)]=0.001;
  //return ( status);
}

__device__ __host__
void computec(real *wmod,real *wd,struct params *p,int i,int j)
{

  
#ifdef ADIABHYDRO
/*below used for adiabatic hydrodynamics*/
  wd[fencode_pre(p,i,j,soundspeed)]=sqrt((p->adiab)/wmod[fencode_pre(p,i,j,rho)]);

#else
wd[fencode_pre(p,i,j,soundspeed)]=sqrt(((p->gamma))*wd[fencode_pre(p,i,j,pressuret)]/wmod[fencode_pre(p,i,j,rho)]);
wd[fencode_pre(p,i,j,cfast)]=sqrt(((wmod[fencode_pre(p,i,j,b1)]*wmod[fencode_pre(p,i,j,b1)]+wmod[fencode_pre(p,i,j,b2)]*wmod[fencode_pre(p,i,j,b2)]+wmod[fencode_pre(p,i,j,b3)]*wmod[fencode_pre(p,i,j,b3)])/wmod[fencode_pre(p,i,j,rho)])+(wd[fencode_pre(p,i,j,soundspeed)]*wd[fencode_pre(p,i,j,soundspeed)]));
#endif



  
}

__device__ __host__
void computecmax(real *wmod,real *wd,struct params *p,int i,int j)
{

       if(wd[fencode_pre(p,i,j,soundspeed)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode_pre(p,i,j,soundspeed)]));
                    p->cmax=(wd[fencode_pre(p,i,j,soundspeed)]);
       if(wd[fencode_pre(p,i,j,cfast)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode_pre(p,i,j,soundspeed)]));
                    p->cmax=(wd[fencode_pre(p,i,j,cfast)]);

}


__global__ void predictor_parallel(struct params *p,  real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order)
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
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
if(i<((p->n[0])) && j<((p->n[1])))
	{		
               for(int f=rho; f<=b3; f++)
               {               
                  wmod[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)];
                  wnew[fencode_pre(p,i,j,f)]=0.0;
               }
               for(int f=current1; f<=divb; f++)
                  wd[fencode_pre(p,i,j,f)]=0; 
        }
               __syncthreads();


  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);

               computebdotv(wmod,wd,p,i,j);
               computedivb(wmod,wd,p,i,j);
         }
              __syncthreads();
  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec(wmod,wd,p,i,j);
        }
              __syncthreads();

  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{ 
               computecmax(wmod,wd,p,i,j);

               /*for(int f=rho; f<=b3; f++)
               {              
                  deriv(dwn1,wd,wmod,p,i,j,f);
                  //dwn1[fencode_pre(p,i,j,f)]=1.0;
                  __syncthreads();
               }*/
               
               /*for(int f=rho; f<=b3; f++) 
                  wmod[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)]+0.5*dt*dwn1[fencode_pre(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn2,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)]+0.5*dt*dwn2[fencode_pre(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn3,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)]+dt*dwn3[fencode_pre(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn4,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  {
                  wnew[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)]+(dt/6.0)*(
                     dwn1[fencode_pre(p,i,j,f)]+2.0*dwn2[fencode_pre(p,i,j,f)]
                         +2.0*dwn3[fencode_pre(p,i,j,f)]+dwn4[fencode_pre(p,i,j,f)]);
               }*/
           //     __syncthreads();
              /* for(int f=rho; f<=b3; f++)
                   wnew[fencode_pre(p,i,j,f)]=w[fencode_pre(p,i,j,f)]+dt*dwn1[fencode_pre(p,i,j,f)];
               computej(wnew,wd,p,i,j);
               computepk(wnew,wd,p,i,j);
               computept(wnew,wd,p,i,j);*/ 


	}
 __syncthreads();
  
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_pre(char *label)
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




int cupredictor(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, int order)
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
     predictor_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

     //following used for testing to check current soundspeeds etc
     //cudaMemcpy(*w, *d_wd, 7*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}






