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
#include "gradops_db.cuh"
#include "dervfields_db.cuh"


__device__ __host__
real dbsourcerho (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  real src=0;

  
 
  return src;
}

__device__ __host__
real dbsourcemom (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real src=0;
  switch(direction)
  {
	case 0:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,b1)];
	break;
	case 1:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,b2)];
	break;
	/*case 2:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,b3)];
	break;*/
  }

  return(isnan(src)?0:src);


}

__device__ __host__
real dbsourceb (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real src=0;
  switch(direction)
  {
	case 0:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,mom1)]/w[fencode_db(p,ix,iy,rho)];
	break;
	case 1:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,mom2)]/w[fencode_db(p,ix,iy,rho)];
	break;
	/*case 2:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,mom3)]/w[fencode_db(p,ix,iy,rho)];
	break;*/
  }
   return(isnan(src)?0:src);
}

__device__ __host__
real dbsourceenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

 real src=0;
    src= -wd[fencode_db(p,ix,iy,divb)]*wd[fencode_db(p,ix,iy,bdotv)];
 
  return ( src);
}


__device__ __host__
int dbderivsourcerho (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=rho;
        dw[fencode_db(p,ix,iy,field)]=dw[fencode_db(p,ix,iy,field)]+dbsourcerho(dw,wd,w,p,ix,iy);
     	//dw[fencode_db(p,ix,iy,field)]=w[fencode_db(p,ix,iy,field)]+10;
  return ( status);
}

__device__ __host__
int dbderivsourcemom (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  int status=0;
     	//dw[fencode_db(p,ix,iy,field)]=w[fencode_db(p,ix,iy,field)]+20+5*(2*direction+1);
        dw[fencode_db(p,ix,iy,field)]=dw[fencode_db(p,ix,iy,field)]+dbsourcemom(dw,wd,w,p,ix,iy,field,direction);
        //dw[fencode_db(p,ix,iy,field)]=-ddotcurrentmom(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int dbderivsourceb (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field, int direction) {

  int status=0;
        dw[fencode_db(p,ix,iy,field)]=dw[fencode_db(p,ix,iy,field)]+dbsourceb(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int dbderivsourceenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  int status=0;
  int field=energy;
        dw[fencode_db(p,ix,iy,field)]=dw[fencode_db(p,ix,iy,field)]+dbsourceenergy(dw,wd,w,p,ix,iy);

  return ( status);
}

//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void dbderivsource (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case rho:
      dbderivsourcerho(dw,wd,w,p,ix,iy);
     break;
     case mom1:
      dbderivsourcemom(dw,wd,w,p,ix,iy,field,0);
     break;
     case mom2:
      dbderivsourcemom(dw,wd,w,p,ix,iy,field,1);
     break;
     /*case mom3:
      dbderivsourcemom(dw,wd,w,p,ix,iy,field,2);
     break;*/
     case energy:
       dbderivsourceenergy(dw,wd,w,p,ix,iy);
     break;
     case b1:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,0);
     break;
     case b2:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,1);
     break;
    /* case b3:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,2);
     break;*/
  }
  //return ( status);
}


__global__ void divb_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real dt)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;

  int ni=p->n[0];
  int nj=p->n[1];

  j=iindex/ni;
  i=iindex-(j*ni);

  if(i<(ni) && j<(nj))
     for(int f=rho; f<=b2; f++)
                dwn1[fencode_db(p,i,j,f)]=0;
 __syncthreads();

  if(i>2 && j>2 && i<(ni-2) && j<(nj-2))
	{
           if(p->divbfix)
           {   

               wd[fencode_db(p,i,j,divb)]=grad_db(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,b1,0)+grad_db(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,b2,1);
               #ifdef USE_SAC
		wd[fencode_db(p,i,j,divb)]+=grad_db(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,b1b,0)+grad_db(wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,b2b,1);
                #endif

               for(int f=rho; f<=b2; f++) 
               //for(int f=rho; f<=b3; f++)
               {              
                  //dbderivsource(dwn1+(NVAR*(p->n[0])*(p->n[1])*order),wd,wmod,p,i,j,f);
                  dbderivsource(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);
 
               }
            }

	}
 __syncthreads();

    if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         {
                         if(p->divbfix)
                          { 
                             for(int f=rho; f<=b2; f++) 
                             //                                                  - sign here same as vac maybe a +
                              wmod[fencode_db(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_db(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]-dt*dwn1[fencode_db(p,i,j,f)]; 
                          }

                         }
              //  }	

  __syncthreads();



  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_db(char *label)
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

int cudivb(struct params **p, real **w,  struct state **state,struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, int order,int ordero, real dt)
{
int status=0;

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
    // prop_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     //cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
    divb_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
    //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
   // cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 return status;


}



