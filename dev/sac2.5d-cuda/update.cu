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
#include "gradops_u.cuh"


__device__ __host__
int updatestate (struct params *p, struct state *s, real *w ,int i, int j, int field) {

  int status=0;
                      // atomicExch(&(p->cmax),(wd[fencode_pre(p,i,j,soundspeed)]));
                    switch(field)
                    {
                      case rho:
                    	s->rho=s->rho+(w[fencode_u(p,i,j,field)]);
		      break;
                      case mom1:
                    	s->m1=s->m1+(w[fencode_u(p,i,j,field)]);
		      break;
                      case mom2:
                    	s->m2=s->m2+(w[fencode_u(p,i,j,field)]);
		      break;
                      /*case mom3:
                    	s->m3=s->m3+(w[fencode_u(p,i,j,field)]);
		      break;*/
                      case energy:
                    	s->e=s->e+(w[fencode_u(p,i,j,field)]);
		      break;
                      case b1:
                    	s->b1=s->b1+(w[fencode_u(p,i,j,field)]);
		      break;
                      case b2:
                    	s->b2=s->b2+(w[fencode_u(p,i,j,field)]);
		      break;
                      /*case b3:
                    	s->b3=s->b3+(w[fencode_u(p,i,j,field)]);
		      break;*/
                    };
  return status;
}



__global__ void update_parallel(struct params *p, struct state *s, real *w, real *wmod)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  __shared__ int ntot;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  //real g=p->g;
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->n[0])*(p->n[1])*rho;
  u=w+(p->n[0])*(p->n[1])*mom1;
  v=w+(p->n[0])*(p->n[1])*mom2;


     j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  //if(i>2 && j >2 && i<((p->n[0])-3) && j<((p->n[1])-3))






if (threadIdx.x == 0) 
{
 ntot=(p->n[0])*(p->n[1]);
 for(int f=rho; f<=NVAR; f++) 
 {
                    switch(f)
                    {
                      case rho:
                    	s->rho=0;
		      break;
                      case mom1:
                    	s->m1=0;
		      break;
                      case mom2:
                    	s->m2=0;
		      break;
                      /*case mom3:
                    	s->m3=0;
		      break;*/
                      case energy:
                    	s->e=0;
		      break;
                      case b1:
                    	s->b1=0;
		      break;
                      case b2:
                    	s->b2=0;
		      break;
                      /*case b3:
                    	s->b3=0;
		      break;*/
                    };

  }              
                 
}
__syncthreads();
 // if(i>1 && j>1 && i<((p->n[0])-2) && j<((p->n[1])-2))
 //if(i>0 && j>0 && i<((p->n[0])) && j<((p->n[1])))
if( i<((p->n[0])) && j<((p->n[1])))
	{
             for(int f=rho; f<NVAR; f++)
             {               
                  w[fencode_u(p,i,j,f)]=wmod[fencode_u(p,i,j,f)];
                  updatestate (p, s, w ,i, j, f);
              }
            // u[i+j*ni]=un[i+j*ni];
           // v[i+j*ni]=vn[i+j*ni];
	   // h[i+j*ni]=hn[i+j*ni];
	}
 __syncthreads();

if (threadIdx.x == 0) 
{
 for(int f=rho; f<NVAR; f++) 
 {
                    switch(f)
                    {
                      case rho:
                    	s->rho=(s->rho)/ntot;
		      break;
                      case mom1:
                    	s->m1=(s->m1)/ntot;
		      break;
                      case mom2:
                    	s->m2=(s->m2)/ntot;
		      break;
                      /*case mom3:
                    	s->m3=(s->m3)/ntot;
		      break;*/
                      case energy:
                    	s->e=(s->e)/ntot;
		      break;
                      case b1:
                    	s->b1=(s->b1)/ntot;
		      break;
                      case b2:
                    	s->b2=(s->b2)/ntot;
		      break;
                      /*case b3:
                    	s->b3=(s->b3)/ntot;
		      break;*/
                    };

  }              
                 
}
__syncthreads();






  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_u(char *label)
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


int cuupdate(struct params **p, real **w, real **wnew, struct state **state,struct params **d_p, real **d_w, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state)
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
    // prop_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_u,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_u,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     //cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_u,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wmod);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
    cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);

    //cudaMemcpy(*w, *d_wd, 6*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);

   cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}


int cufinish(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd)
{
  

 //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  checkErrors_u("copy data from device");


  cudaFree(*d_p);
//  cudaFree(*d_state);

  cudaFree(*d_w);
  cudaFree(*d_wnew);
 // cudaFree(*d_u);

  cudaFree(*d_wmod);
  cudaFree(*d_dwn1);
  cudaFree(*d_wd);



}
