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
#include "dervfields_cd1.cuh"

__device__ __host__
int divflux1(real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;

dw[fencode_cd1(p,ix,iy,field)]= grad_cd1(wd,p,ix,iy,flux,dir); 
//dw[fencode_cd1(p,ix,iy,field)]=0.0;
 switch(field)
  {
     case mom1:
       dw[fencode_cd1(p,ix,iy,field)]+= (p->g[dir])*w[fencode_cd1(p,ix,iy,rho)];
      break;

    case mom2:
      dw[fencode_cd1(p,ix,iy,field)]+= (p->g[dir])*w[fencode_cd1(p,ix,iy,rho)];
      break;
    case rho:
     ;// dw[fencode_cd1(p,ix,iy,field)]+= ix/800;
      break;

  }    
 // dw[fencode_cd1(p,ix,iy,field)]= gradd0_cd1(wd,p,ix,iy,f1,0)+gradd1_cd1(wd,p,ix,iy,f2,1);    
  return ( status);
}




__device__ __host__
void bc_periodic1_cd1t1(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 || i==1 )                
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,(p->n[0])-4+i,j,f)];
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) )                
                   ;// wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,4-(p->n[0])+i,j,f)];
                else if(j==0 || j==1 )                
                 ;// wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,(p->n[1])-4+j,f)];
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                
                 ;// wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,4-(p->n[1])+j,f)];

 


}


__device__ __host__
void bc_periodic2_cd1t1(real *wt, struct params *p,int i, int j, int f) {


               if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,(p->n[1])-4+j,f)];
                  else                  
                    //wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,(p->n[0])-4+i,j,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,(p->n[0])-4+i,j,f)];                                     
                  else                  
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,4-(p->n[1])+j,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,4-(p->n[0])+i,j,f)];                                    
                  else                  
                   wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,(p->n[1])-4+j,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,i,4-(p->n[1])+j,f)];                                    
                  else                  
                    wt[fencode_cd1(p,i,j,f)]=wt[fencode_cd1(p,4-(p->n[0])+i,j,f)];                                    
                }                       
                 
                




}



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
     //flux= wd[fencode_cd1(p,ix,iy,vel1)]*w[fencode_cd1(p,ix,iy,field)];
        #ifdef USE_SAC
     flux= w[fencode_cd1(p,ix,iy,mom1)]*w[fencode_cd1(p,ix,iy,field)]/(w[fencode_cd1(p,ix,iy,rho)]+w[fencode_cd1(p,ix,iy,rhob)]);
    // flux= w[fencode_cd1(p,ix,iy,mom1)]*w[fencode_cd1(p,ix,iy,field)]/w[fencode_cd1(p,ix,iy,rho)];

        #else
     flux= w[fencode_cd1(p,ix,iy,mom1)]*w[fencode_cd1(p,ix,iy,field)]/w[fencode_cd1(p,ix,iy,rho)];

        #endif
     break;
     case 1:
        #ifdef USE_SAC
     flux= w[fencode_cd1(p,ix,iy,mom2)]*w[fencode_cd1(p,ix,iy,field)]/(w[fencode_cd1(p,ix,iy,rho)]+w[fencode_cd1(p,ix,iy,rhob)]);
     //flux= w[fencode_cd1(p,ix,iy,mom2)]*w[fencode_cd1(p,ix,iy,field)]/w[fencode_cd1(p,ix,iy,rho)];

        #else
     //flux= wd[fencode_cd1(p,ix,iy,vel2)]*w[fencode_cd1(p,ix,iy,field)];
     flux= w[fencode_cd1(p,ix,iy,mom2)]*w[fencode_cd1(p,ix,iy,field)]/w[fencode_cd1(p,ix,iy,rho)];

        #endif
     break;
    /* case 2:
     flux= wd[fencode_cd1(p,ix,iy,vel3)]*w[fencode_cd1(p,ix,iy,field)];
     break;*/
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
     		flux= -w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b1b)]-w[fencode_cd1(p,ix,iy,field+7)]*w[fencode_cd1(p,ix,iy,b1)]-w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b1)];
        #endif
        #ifdef USE_VAC
                flux= -w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b1)];
         #endif

     break;
     case 1:
        #ifdef USE_SAC
                flux= -w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b2b)]-w[fencode_cd1(p,ix,iy,field+7)]*w[fencode_cd1(p,ix,iy,b2)]-w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b2)];
         #endif
        #ifdef USE_VAC
              flux= -w[fencode_cd1(p,ix,iy,field+3)]*w[fencode_cd1(p,ix,iy,b1)];
         #endif
     break;

   }




  return flux;


  //return ( ddc1-ddc2);
}










__device__ __host__
int computefluxrho (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int direction) {

  int field;
  int status=0;

  // for(direction=0;direction<2;direction++)
         #ifdef USE_SAC
	      wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,rho,direction)+(w[fencode_cd1(p,ix,iy,rhob)]*w[fencode_cd1(p,ix,iy,mom1+direction)])/(w[fencode_cd1(p,ix,iy,rhob)]+w[fencode_cd1(p,ix,iy,rho)]);
         #else
             wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,rho,direction);
         #endif
  
  return ( status);
}

__device__ __host__
int computefluxmom (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field,int direction) {

 
  int status=0;
  //for(direction=0;direction<2;direction++)
  //{
    switch(field)
  {
     case mom1:
        #ifdef ADIABHYDRO
     		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction);
        #endif
        #ifdef USE_VAC
    		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
 
        #endif
        #ifdef USE_SAC
    		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
               //if(direction==0)
               //   wd[fencode_cd1(p,ix,iy,f1)]+=wd[fencode_cd1(p,ix,iy,ptb)];
 
        #endif
               if(direction==0)
               {

                // computept_cd1(w,wd,p,ix,iy);
                 //commented out to compare with vac 
                 wd[fencode_cd1(p,ix,iy,flux)]+=wd[fencode_cd1(p,ix,iy,pressuret)];

        #ifdef USE_SAC

                  wd[fencode_cd1(p,ix,iy,flux)]+=wd[fencode_cd1(p,ix,iy,ptb)];
       #endif
               }
 
     break;
     case mom2:
        #ifdef ADIABHYDRO
     		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction);
        #endif
        #ifdef USE_VAC
    		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
 
        #endif
        #ifdef USE_SAC
    		wd[fencode_cd1(p,ix,iy,flux)]= transportflux(dw,wd,w,p,ix,iy,field,direction)+fluxmom1(dw,wd,w,p,ix,iy,field,direction);
               //if(direction==1)
               //   wd[fencode_cd1(p,ix,iy,f1)]+=wd[fencode_cd1(p,ix,iy,ptb)];
 
        #endif
               if(direction==1)
               {
                //computept_cd1(w,wd,p,ix,iy);
                //commented out to compare with vac
                 wd[fencode_cd1(p,ix,iy,flux)]+=wd[fencode_cd1(p,ix,iy,pressuret)];

        #ifdef USE_SAC

                  wd[fencode_cd1(p,ix,iy,flux)]+=wd[fencode_cd1(p,ix,iy,ptb)];
        #endif

               }
 
     break;
 
  // }
}
        
  return ( status);
}







//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field,int dir) {

  //int status=0;
  switch(field)
  {
     case rho:
      computevel_cd1(w,wd,p,ix,iy);
      computept_cd1(w,wd,p,ix,iy);
      computefluxrho(dw,wd,w,p,ix,iy,dir);
     break;
     case mom1:
      computefluxmom(dw,wd,w,p,ix,iy,field,dir);
      //wd[fencode_cd1(p,ix,iy,f1)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     case mom2:
       computefluxmom(dw,wd,w,p,ix,iy,field,dir);
       //wd[fencode_cd1(p,ix,iy,f2)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;
     /*case mom3:
      computefluxmom(dw,wd,w,p,ix,iy,field);
      //wd[fencode_cd1(p,ix,iy,f3)]+=wd[fencode_cd1(p,ix,iy,pressuret)];
     break;*/
  }
  //return ( status);
}



__global__ void centdiff1_parallel(struct params *p, real *w, real *wmod, 
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


             //  for(int f=rho; f<=mom3; f++)
             //  {
			if(i<(ni) && j<(nj))
                        {
                            dwn1[fencode_cd1(p,i,j,f)]=0.0;
                  	    //for(fid=0;fid<2;fid++)
                               wd[fencode_cd1(p,i,j,flux)]=0.0;
                               //wmod[fencode_cd1(p,i,j,flux)+order*NVAR*(p->n[0])*(p->n[1])]=0.0;
                        }

   }
 __syncthreads();                       



   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

//1. 11/1/11 could swap cases below
                        switch(dir)
                        {
                         case 0:
                         if(i<(ni)  && j >1 &&  j<(nj-1))
                         
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f,dir); 
                         break;
                         case 1:
                         if(i>1 &&  i<(ni-1) && j<(nj))
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                        
                            computeflux(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f,dir); 
                         break;
                        }
              //  }
                        //might need to set boundaries correctly
 
}
__syncthreads();                        

/*  for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
             //for(fid=0;fid<2;fid++)
              #ifdef ADIABHYDRO
                  bc_cont_cd1(wd,p,i,j,flux);
              #else
                  bc_periodic1_cd1t1(wd,p,i,j,flux);


              #endif
                __syncthreads();
}



#ifndef ADIABHYDRO
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj) )
             //for(fid=0;fid<2;fid++)
                  //bc_cont_cd1(dwn1,p,i,j,f1+fid);
             
                 bc_periodic2_cd1t1(wd,p,i,j,flux);
                __syncthreads();
}

#endif*/


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
			 //if(i>2 && j >2 && i<(ni-3) && j<(nj-3))
                         //if( i<(ni) && j<(nj))
                               divflux1(dwn1,wd,wmod+order*NVAR*(p->n[0])*(p->n[1]),p,i,j,f,dir);  

}
 __syncthreads();
               // }
    

         // if( i<(ni) && j<(nj))
          //        bc_cont_cd1(dwn1,p,i,j,f);
            //    __syncthreads();



             // for(int f=rho; f<=mom3; f++)
              // {

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


                        switch(dir)
                        {
                         case 0:

                         if(i<(ni)  && j >1 &&  j<(nj-2))
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode_cd1(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_cd1(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]-dt*dwn1[fencode_cd1(p,i,j,f)]; 
                         break;
                         case 1:

                         //if(i<(ni)  && j >1 &&  j<(nj-2))
                         if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode_cd1(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]=wmod[fencode_cd1(p,i,j,f)+(ordero*NVAR*(p->n[0])*(p->n[1]))]-dt*dwn1[fencode_cd1(p,i,j,f)];
                         break;
                        }


              //  }
	
}
  __syncthreads();


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

int cucentdiff1(struct params **p, struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir)
{

 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
 //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
 // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

     centdiff1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     centdiff1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     cudaThreadSynchronize();
}


