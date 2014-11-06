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
#include "gradops_cd1a.cuh"
#include "dervfields_cd1a.cuh"



__device__ __host__
int divflux1(real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;
  //dw[fencode_cd1a(p,ix,iy,field)]= grad_cd1a(wd,p,ix,iy,f1,dir);//+grad_cd1a(wd,p,ix,iy,f2,1); 
dw[fencode_cd1a(p,ix,iy,field)]= grad_cd1a(wd,p,ix,iy,flux,dir);//+grad_cd1a(wd,p,ix,iy,f2,1); 
 /*switch(field)
  {
     case mom1:
       dw[fencode_cd1a(p,ix,iy,field)]+= grad_cd1a(wd,p,ix,iy,pressuret,0);
      break;

    case mom2:
      dw[fencode_cd1a(p,ix,iy,field)]+= grad_cd1a(wd,p,ix,iy,pressuret,1);
      break;


  } */   
 // dw[fencode_cd1a(p,ix,iy,field)]= gradd0_cd1a(wd,p,ix,iy,f1,0)+gradd1_cd1a(wd,p,ix,iy,f2,1);    
  return ( status);
}






__global__ void centdiff1a_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
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

   int ip,jp,ipg,jpg;
   fid=0;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
			 if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         //if( i<(ni) && j<(nj))
                               divflux1(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f,dir);  
 __syncthreads();
}
               // }
    

         // if( i<(ni) && j<(nj))
          //        bc_cont_cd1a(dwn1,p,i,j,f);
            //    __syncthreads();



             // for(int f=rho; f<=mom3; f++)
              // {

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


                       // switch(dir)
                       // {
                       //  case 0:
                         //if(i<(ni)  && j >1 &&  j<(nj-2))
                         //if( i >1 &&  i<(ni-2)  && j >1 &&  j<(nj-2))
                         if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         //if(i>2 && j >2 && i<(ni-3) && j<(nj-3))
                         //if(i<(ni)  && j >3 &&  j<(nj-4))
                         //if(i<(ni)  && j >1 &&  j<(nj-2))
                              wmod[fencode_cd1a(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_cd1a(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]-dt*dwn1[fencode_cd1a(p,i,j,f)]; 
                      //   break;
                      //   case 1:
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                         //if(i >1 &&  i<(ni-2)  && j >1 &&  j<(nj-2))
                         //if(i>3 && j >3 && i<(ni-4) && j<(nj-4))
                         //if(i>3 &&  i<(ni-4) && j<(nj))
                      //   if(i>1 &&  i<(ni-2) && j<(nj))
                      //        wmod[fencode_cd1a(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_cd1a(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]-dt*dwn1[fencode_cd1a(p,i,j,f)];
                       //  break;
                      //  }


              //  }
 __syncthreads();	
}
 


}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cd1a(char *label)
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




int cucentdiff1a(struct params **p, real **w, struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
 //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
 // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
//cudaThreadSynchronize();

     centdiff1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
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


