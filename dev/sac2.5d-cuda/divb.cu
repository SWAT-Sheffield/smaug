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
int encode_db (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_db (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}

__device__ __host__
real evalgrad_db(real fi, real fim1, real fip2, real fim2,struct params *p,int dir)
{
 //real valgrad_db;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx)))*(fi-fim1)-(1.0/(12.0*(p->dx)))*(fip2-fim2);
   //return((1.0/(2.0*(p->dx)))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx)))*(fi-fim1)):((1.0/(12.0*(p->dx)))*((8*fi-8*fim1+fim2-fip2))));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(fi-fim1)-(1.0/(12.0*(p->dy)))*(fip2-fim2);
     // return((2.0/(1.0*(p->dy)))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dy)))*(fi-fim1)):((1.0/(12.0*(p->dy)))*((8*fi-8*fim1+fim2-fip2))));
 }

 return -1;
}


__device__ __host__
real grad_db(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_db;

 if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx)))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
//return((1.0/(2.0*(p->dx)))*(wmod[fencode_db(p,i+1,j,field)]-wmod[fencode_db(p,i-1,j,field)]));
 return(  ( (p->sodifon)?((8*wmod[fencode_db(p,i+1,j,field)]-8*wmod[fencode_db(p,i-1,j,field)]+wmod[fencode_db(p,i-1,j,field)]-wmod[fencode_db(p,i+1,j,field)])/6.0):wmod[fencode_db(p,i+1,j,field)]-wmod[fencode_db(p,i-1,j,field)])/(2.0*(p->dx))    );
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dy)))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dy)))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
// return((1.0/(2.0*(p->dy)))*(wmod[fencode_db(p,i,j+1,field)]-wmod[fencode_db(p,i,j-1,field)]));
 return(  ( (p->sodifon)?((8*wmod[fencode_db(p,i,j+1,field)]-8*wmod[fencode_db(p,i,j-1,field)]+wmod[fencode_db(p,i,j-1,field)]-wmod[fencode_db(p,i,j+1,field)])/6.0):wmod[fencode_db(p,i,j+1,field)]-wmod[fencode_db(p,i,j-1,field)])/(2.0*(p->dy))    );

 }

 return -1;
}

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
	case 2:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,b3)];
	break;
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
	case 2:
         src= -wd[fencode_db(p,ix,iy,divb)]*w[fencode_db(p,ix,iy,mom3)]/w[fencode_db(p,ix,iy,rho)];
	break;
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
     case mom3:
      dbderivsourcemom(dw,wd,w,p,ix,iy,field,2);
     break;
     case energy:
       dbderivsourceenergy(dw,wd,w,p,ix,iy);
     break;
     case b1:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,0);
     break;
     case b2:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,1);
     break;
     case b3:
      dbderivsourceb(dw,wd,w,p,ix,iy,field,2);
     break;
  }
  //return ( status);
}


__global__ void divb_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  __shared__ int ntot;

  int ni=p->ni;
  int nj=p->nj;
  real dt=p->dt;
  real dy=p->dy;
  real dx=p->dx;
  real g=p->g;
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

  real *un,  *vn,  *hn;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  hn=wnew+(p->ni)*(p->nj)*rho;
  un=wnew+(p->ni)*(p->nj)*mom1;
  vn=wnew+(p->ni)*(p->nj)*mom2;
     j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  //if(i>2 && j >2 && i<((p->ni)-3) && j<((p->nj)-3))



  if(i<p->ni && j<p->nj)
	{
           if(p->divbfix)
           {    
               for(int f=rho; f<=b3; f++)
               {              
                  dbderivsource(dwn1,wd,wmod,p,i,j,f);
                  //dwn1[fencode_ds(p,i,j,f)]=1.0;
                  __syncthreads();
               }
            }
            // u[i+j*ni]=un[i+j*ni];
           // v[i+j*ni]=vn[i+j*ni];
	   // h[i+j*ni]=hn[i+j*ni];
	}
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

int cudivb(struct params **p, real **w, real **wnew,  struct state **state,struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;

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
    divb_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
    //cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
   // cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}



