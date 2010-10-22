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
#include "gradops_cd1.cuh"
__device__ __host__
real transportflux (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;

   //transport flux
    switch(direction)
  {
     case 0:
     flux= wd[fencode_cd1(p,ix,iy,vel1)]*w[fencode_cd1(p,ix,iy,field)];
     break;
     case 1:
     flux= wd[fencode_cd1(p,ix,iy,vel2)]*w[fencode_cd1(p,ix,iy,field)];
     break;
     case 2:
     flux= wd[fencode_cd1(p,ix,iy,vel3)]*w[fencode_cd1(p,ix,iy,field)];
     break;
   }
  return flux;


  //return ( ddc1-ddc2);
}






__device__ __host__
real fluxmom1 (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;

   //transport flux
    switch(direction)
  {
     case 0:
        #ifdef USE_SAC
      		flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b1)]-w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b1b)]-w[fencode_cd1(p,ix,iy,field+9)]*w[fencode_cd1(p,ix,iy,b1)];
        #endif
        #ifdef USE_VAC
      		flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b1)];
         #endif
       // if((field==mom1 ) )
        //               flux+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     case 1:
        #ifdef USE_SAC
	      flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b2)]-w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b2b)]-w[fencode_cd1(p,ix,iy,field+9)]*w[fencode_cd1(p,ix,iy,b2)];
         #endif
        #ifdef USE_VAC
	      flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b2)];
         #endif
        // if((field==mom2 && direction==1)   )
        //             flux+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     case 2:
        #ifdef USE_SAC
     	 	flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b3)]-w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b3b)]-w[fencode_cd1(p,ix,iy,field+9)]*w[fencode_cd1(p,ix,iy,b3)];
         #endif
        #ifdef USE_VAC
      		flux= -w[fencode_cd1(p,ix,iy,field+4)]*w[fencode_cd1(p,ix,iy,b3)];	
         #endif
        // if((field==mom3 && direction==2))
         //           flux+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
   }




  return flux;


  //return ( ddc1-ddc2);
}










__device__ __host__
int computefluxrho (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

  int field, direction;
  int status=0;
  
   for(direction=0;direction<3;direction++)
         #ifdef USE_SAC
	      wd[fencode_cd1(p,ix,iy,f1+direction)]= transportflux(dw,wd,w,p,ix,iy,rho,direction)+transportflux(dw,wd,w,p,ix,iy,rhob,direction);
         #else
             wd[fencode_cd1(p,ix,iy,f1+direction)]= transportflux(dw,wd,w,p,ix,iy,rho,direction);
         #endif
  
  return ( status);
}

__device__ __host__
int computefluxmom (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  int direction;
  int status=0;
  for(direction=0;direction<3;direction++)
  {
   // switch(direction)
 // {
   //  case 0:
        #ifdef ADIABHYDRO
     		wd[fencode_cd1(p,ix,iy,f1+direction)]= transportflux(dw,wd,w,p,ix,iy,field,direction);
        #else
    		wd[fencode_cd1(p,ix,iy,f1+direction)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
        #endif
   //  break;
    // case 1:
    /*    #ifdef ADIABHYDRO
	     wd[fencode_cd1(p,ix,iy,f2)]= transportflux(dw,wd,w,p,ix,iy,field,direction);
        #else
	     wd[fencode_cd1(p,ix,iy,f2)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
        #endif*/
  //   break;
  //   case 2:
    /*    #ifdef ADIABHYDRO
	     wd[fencode_cd1(p,ix,iy,f3)]= transportflux(dw,wd,w,p,ix,iy,field,direction);
        #else
	     wd[fencode_cd1(p,ix,iy,f3)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
        #endif*/
    // break;
  // }
}
        
  return ( status);
}

__device__ __host__
int divflux1(real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field) {

  int direction;
  int status=0;
  real divflux=0;
  dw[fencode_cd1(p,ix,iy,field)]= grad_cd1(wd,p,ix,iy,f1,0)+grad_cd1(wd,p,ix,iy,f2,1);      
 // dw[fencode_cd1(p,ix,iy,field)]= gradd0_cd1(wd,p,ix,iy,f1,0)+gradd1_cd1(wd,p,ix,iy,f2,1);    
  return ( status);
}





//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field) {

  //int status=0;
  switch(field)
  {
     case rho:
      computefluxrho(dw,wd,w,p,ix,iy);
     break;
     case mom1:
      computefluxmom(dw,wd,w,p,ix,iy,field);
      wd[fencode_cd1(p,ix,iy,f1)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     case mom2:
       computefluxmom(dw,wd,w,p,ix,iy,field);
       wd[fencode_cd1(p,ix,iy,f2)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     case mom3:
      computefluxmom(dw,wd,w,p,ix,iy,field);
      wd[fencode_cd1(p,ix,iy,f3)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
  }
  //return ( status);
}



__global__ void centdiff1_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f)
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


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);

             //  for(int f=rho; f<=mom3; f++)
             //  {
			if(i<(ni) && j<(nj))
                        {
                  	    for(fid=0;fid<3;fid++)
                               dwn1[fencode_cd1(p,i,j,f1+fid)]=0.0;

                        }
                        __syncthreads();

			//if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                        if(i<(ni) && j<(nj))
                        {
                            computeflux(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f); 
                        }
              //  }
                        //might need to set boundaries correctly 
                        __syncthreads();

          if( i<(ni) && j<(nj))
             for(fid=0;fid<3;fid++)
                  bc_cont_cd1(dwn1,p,i,j,f1+fid);
                __syncthreads();
			//if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                        //        divflux1(dwn1+(NVAR*(p->n[0])*(p->n[1])*order),wd,wmod,p,i,j,f);
            //  for(int f=rho; f<=mom3; f++)
             //  {
			 if(i>1 && j >1 && i<(ni-2) && j<(nj-2))

                               divflux1(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);  

               // }
     __syncthreads();

         // if( i<(ni) && j<(nj))
          //        bc_cont_cd1(dwn1,p,i,j,f);
            //    __syncthreads();



             // for(int f=rho; f<=mom3; f++)
              // {
			 if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode_cd1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]=wmod[fencode_cd1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]-dt*dwn1[fencode_cd1(p,i,j,f)]; 
//wmod[fencode_cd1(p,i,j,f)+ordero*NVAR*(p->n[0])*(p->n[1])]=dwn1[fencode_cd1(p,i,j,f2)];
                              //dwn1[fencode_cd1(p,i,j,f)]=0;
                         }
              //  }	

  __syncthreads();


         //if( i<(ni) && j<(nj))
         //         bc_cont_cd1(wmod+ordero*NVAR*(p->n[0])*(p->n[1]),p,i,j,f);
         //       __syncthreads();


}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cd1(char *label)
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




int cucentdiff1(struct params **p, real **w, struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field)
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
     centdiff1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field);
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


