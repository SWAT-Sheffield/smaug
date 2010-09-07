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
#include "gradops_ds.cuh"


__device__ __host__
real sourcerho (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

 // real src=0;
 // int field=rho;
  real src=0;
               src= -(p->chyp)*(grad2_ds(w,p,ix,iy,rho,0)+grad2_ds(w,p,ix,iy,rho,1));
        #ifdef USE_SAC
           src= -(p->chyp)*(grad2_ds(w,p,ix,iy,rhob,0)+grad2_ds(w,p,ix,iy,rhob,1));
         #endif
 
  return src;
}

__device__ __host__
real sourcemom (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real src=0;

  switch(direction)
  {
	case 0:
         src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[0]))-grad_ds(wd,p,ix,iy,pressuret,0)-(p->chyp)*(grad2_ds(w,p,ix,iy,mom1,0)+grad2_ds(w,p,ix,iy,mom1,1));

        // src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[0]));
	break;
	case 1:
         src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[1]))-grad_ds(wd,p,ix,iy,pressuret,1)-(p->chyp)*(grad2_ds(w,p,ix,iy,mom2,1)+grad2_ds(w,p,ix,iy,mom2,0));

         //src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[1]));
	break;
	case 2:
         //src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[2]))-grad_ds(wd,p,ix,iy,pressuret,2);
         src=(w[fencode_ds(p,ix,iy,rho)]*(p->g[2]));

	break;
  }

  return(isnan(src)?0:src);


}

__device__ __host__
real sourceb (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction) {

  real src=0;
  switch(direction)
  {
	case 0:
         src=(p->eta)*grad_ds(wd,p,ix,iy,current3,1);
	break;
	case 1:
         src= -(p->eta)*grad_ds(wd,p,ix,iy,current3,0);
	break;
	case 2:
         src= (p->eta)*(grad_ds(wd,p,ix,iy,current2,0)-grad_ds(wd,p,ix,iy,current1,1));
	break;
  }
   return(isnan(src)?0:src);
}

__device__ __host__
real sourceenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy) {

 // real src=0;
  real srcg,srcb,srcc;
  int field=energy;
  real ddcx,ddcy;
  real fi,fim1,fip2,fim2;
  fip2=0;
  fim2=0;
  srcc=0.0;

         #ifdef USE_SAC
  	    srcg=w[fencode_ds(p,ix,iy,rho)]*(((p->g[0])*w[fencode_ds(p,ix,iy,mom1)]+(p->g[1])*w[fencode_ds(p,ix,iy,mom2)]+(p->g[2])*w[fencode_ds(p,ix,iy,mom3)]))/(w[fencode_ds(p,ix,iy,rho)]+w[fencode_ds(p,ix,iy,rhob)]);
         #else
     		 srcg=(p->g[0])*w[fencode_ds(p,ix,iy,mom1)]+(p->g[1])*w[fencode_ds(p,ix,iy,mom2)]+(p->g[2])*w[fencode_ds(p,ix,iy,mom3)];
         #endif


         #ifdef USE_SAC
	       fi=((w[fencode_ds(p,ix+1,iy,b2)]+w[fencode_ds(p,ix+1,iy,b2b)])*wd[fencode_ds(p,ix+1,iy,current3)]-(w[fencode_ds(p,ix+1,iy,b3)]+w[fencode_ds(p,ix+1,iy,b3b)])*wd[fencode_ds(p,ix+1,iy,current2)]);
	       fim1=((w[fencode_ds(p,ix-1,iy,b2)]+w[fencode_ds(p,ix-1,iy,b2b)])*wd[fencode_ds(p,ix-1,iy,current3)]-(w[fencode_ds(p,ix-1,iy,b3)]+w[fencode_ds(p,ix-1,iy,b3b)])*wd[fencode_ds(p,ix-1,iy,current2)]);
	    if(p->sodifon)
	     {
	       fip2=((w[fencode_ds(p,ix+2,iy,b2)]+w[fencode_ds(p,ix+2,iy,b2b)])*wd[fencode_ds(p,ix+2,iy,current3)]-(w[fencode_ds(p,ix+2,iy,b3)]+w[fencode_ds(p,ix+2,iy,b3b)])*wd[fencode_ds(p,ix+2,iy,current2)]);
	       fim2=((w[fencode_ds(p,ix-2,iy,b2)]+w[fencode_ds(p,ix-2,iy,b2b)])*wd[fencode_ds(p,ix-2,iy,current3)]-(w[fencode_ds(p,ix-2,iy,b3)]+w[fencode_ds(p,ix-2,iy,b3b)])*wd[fencode_ds(p,ix-2,iy,current2)]);
	      }
         #else
	       fi=(w[fencode_ds(p,ix+1,iy,b2)]*wd[fencode_ds(p,ix+1,iy,current3)]-w[fencode_ds(p,ix+1,iy,b3)]*wd[fencode_ds(p,ix+1,iy,current2)]);
	       fim1=(w[fencode_ds(p,ix-1,iy,b2)]*wd[fencode_ds(p,ix-1,iy,current3)]-w[fencode_ds(p,ix-1,iy,b3)]*wd[fencode_ds(p,ix-1,iy,current2)]);
	    if(p->sodifon)
	     {
	       fip2=(w[fencode_ds(p,ix+2,iy,b2)]*wd[fencode_ds(p,ix+2,iy,current3)]-w[fencode_ds(p,ix+2,iy,b3)]*wd[fencode_ds(p,ix+2,iy,current2)]);
	       fim2=(w[fencode_ds(p,ix-2,iy,b2)]*wd[fencode_ds(p,ix-2,iy,current3)]-w[fencode_ds(p,ix-2,iy,b3)]*wd[fencode_ds(p,ix-2,iy,current2)]);
	      }

         #endif

       ddcx=evalgrad_ds(fi,fim1,fip2,fim2,p,0);
      //ddcx=evalgrad_ds(fi,fim1,0,0,p,0);

         #ifdef USE_SAC
	      fi=((w[fencode_ds(p,ix,iy+1,b3)]+w[fencode_ds(p,ix,iy+1,b3b)])*wd[fencode_ds(p,ix,iy+1,current1)]-(w[fencode_ds(p,ix,iy+1,b1)]+w[fencode_ds(p,ix,iy+1,b1b)])*wd[fencode_ds(p,ix,iy+1,current3)]);
	       fim1=((w[fencode_ds(p,ix,iy-1,b3)]+w[fencode_ds(p,ix,iy-1,b3b)])*wd[fencode_ds(p,ix,iy-1,current1)]-(w[fencode_ds(p,ix,iy-1,b1)]+w[fencode_ds(p,ix,iy-1,b1b)])*wd[fencode_ds(p,ix,iy-1,current3)]);
	    if(p->sodifon)
	     {
	       fip2=((w[fencode_ds(p,ix,iy+2,b3)]+w[fencode_ds(p,ix,iy+2,b3b)])*wd[fencode_ds(p,ix,iy+2,current1)]-(w[fencode_ds(p,ix,iy+2,b1)]+w[fencode_ds(p,ix,iy+2,b1b)])*wd[fencode_ds(p,ix,iy+2,current3)]);
	       fim2=((w[fencode_ds(p,ix,iy-2,b3)]+w[fencode_ds(p,ix,iy-2,b3b)])*wd[fencode_ds(p,ix,iy-2,current1)]-(w[fencode_ds(p,ix,iy-2,b1)]+w[fencode_ds(p,ix,iy-2,b1b)])*wd[fencode_ds(p,ix,iy-2,current3)]);
         #else
	      fi=(w[fencode_ds(p,ix,iy+1,b3)]*wd[fencode_ds(p,ix,iy+1,current1)]-w[fencode_ds(p,ix,iy+1,b1)]*wd[fencode_ds(p,ix,iy+1,current3)]);
	       fim1=(w[fencode_ds(p,ix,iy-1,b3)]*wd[fencode_ds(p,ix,iy-1,current1)]-w[fencode_ds(p,ix,iy-1,b1)]*wd[fencode_ds(p,ix,iy-1,current3)]);
	    if(p->sodifon)
	     {
	       fip2=(w[fencode_ds(p,ix,iy+2,b3)]*wd[fencode_ds(p,ix,iy+2,current1)]-w[fencode_ds(p,ix,iy+2,b1)]*wd[fencode_ds(p,ix,iy+2,current3)]);
	       fim2=(w[fencode_ds(p,ix,iy-2,b3)]*wd[fencode_ds(p,ix,iy-2,current1)]-w[fencode_ds(p,ix,iy-2,b1)]*wd[fencode_ds(p,ix,iy-2,current3)]);
         #endif
 
      }
      // ddcx=evalgrad_ds(fi,fim1,fip2,fim2,p,0);
      ddcy=evalgrad_ds(fi,fim1,fip2,fim2,p,1);

      srcb=(p->eta)*((isnan(ddcx)?0:ddcx)+(isnan(ddcy)?0:ddcy));

        #ifdef USE_SAC
         //calc gradv
	      fip2=0;
	      fim2=0;
	       fi=(w[fencode_ds(p,ix+1,iy,mom1)]/(w[fencode_ds(p,ix+1,iy,rho)]+w[fencode_ds(p,ix+1,iy,rhob)]));
	       fim1=(w[fencode_ds(p,ix-1,iy,mom1)]/(w[fencode_ds(p,ix-1,iy,rho)]+w[fencode_ds(p,ix-1,iy,rhob)]));
	    if(p->sodifon){
	       fip2=(w[fencode_ds(p,ix+2,iy,mom1)]/(w[fencode_ds(p,ix+2,iy,rho)]+w[fencode_ds(p,ix+2,iy,rhob)]));
	       fim2=(w[fencode_ds(p,ix-2,iy,mom1)]/(w[fencode_ds(p,ix-2,iy,rho)]+w[fencode_ds(p,ix-2,iy,rhob)]));
	      }
            ddcx=evalgrad_ds(fi,fim1,fip2,fim2,p,0);
	      fip2=0;
	      fim2=0;
	       fi=(w[fencode_ds(p,ix,iy+1,mom2)]/(w[fencode_ds(p,ix,iy+1,rho)]+w[fencode_ds(p,ix,iy+1,rhob)]));
	       fim1=(w[fencode_ds(p,ix,iy-1,mom2)]/(w[fencode_ds(p,ix,iy-1,rho)]+w[fencode_ds(p,ix,iy-1,rhob)]));
	    if(p->sodifon){
	       fip2=(w[fencode_ds(p,ix,iy+2,mom2)]/(w[fencode_ds(p,ix,iy+2,rho)]+w[fencode_ds(p,ix,iy+2,rhob)]));
	       fim2=(w[fencode_ds(p,ix,iy-2,mom2)]/(w[fencode_ds(p,ix,iy-2,rho)]+w[fencode_ds(p,ix,iy-2,rhob)]));
	      }
            ddcy=evalgrad_ds(fi,fim1,fip2,fim2,p,1);

                  srcc=-(ddcx+ddcy)*wd[fencode_ds(p,ix,iy,ptb)]+(ddcx*(w[fencode_ds(p,ix,iy,b1b)])+ddcy*(w[fencode_ds(p,ix,iy,b2b)]))*(w[fencode_ds(p,ix,iy,b1b)]+w[fencode_ds(p,ix,iy,b2b)]);

         #endif


 // src=srcg+srcb;
  return ( srcg+srcb+srcc);
}


__device__ __host__
int derivsourcerho (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int order) {

  int status=0;
  int field=rho;
        dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]=dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]+sourcerho(dw,wd,w,p,ix,iy);
     	//dw[fencode_ds(p,ix,iy,field)]=w[fencode_ds(p,ix,iy,field)]+10;
  return ( status);
}

__device__ __host__
int derivsourcemom (real *dw, real *wd, real *w, struct params *p,int ix, int iy,int field, int direction, int order) {

  int status=0;
     	//dw[fencode_ds(p,ix,iy,field)]=w[fencode_ds(p,ix,iy,field)]+20+5*(2*direction+1);
        dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]=dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]+sourcemom(dw,wd,w,p,ix,iy,field,direction);
        //dw[fencode_ds(p,ix,iy,field)]=-ddotcurrentmom(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivsourceb (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field, int direction, int order) {

  int status=0;
        dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]=dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]+sourceb(dw,wd,w,p,ix,iy,field,direction);

  return ( status);
}

__device__ __host__
int derivsourceenergy (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int order) {

  int status=0;
  int field=energy;
        dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]=dw[(NVAR*(p->n[0])*(p->n[1])*order)+fencode_ds(p,ix,iy,field)]+sourceenergy(dw,wd,w,p,ix,iy);

  return ( status);
}


//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void derivsource (real *dw, real *wd, real *w, struct params *p,int ix, int iy, int field, int order) {

  //int status=0;
  switch(field)
  {
     case rho:
      derivsourcerho(dw,wd,w,p,ix,iy,order);
     break;
     case mom1:
      derivsourcemom(dw,wd,w,p,ix,iy,field,0,order);
     break;
     case mom2:
      derivsourcemom(dw,wd,w,p,ix,iy,field,1,order);
     break;
     case mom3:
      derivsourcemom(dw,wd,w,p,ix,iy,field,2,order);
     break;
     case energy:
       derivsourceenergy(dw,wd,w,p,ix,iy,order);
     break;
     case b1:
    derivsourceb(dw,wd,w,p,ix,iy,field,0,order);
     break;
     case b2:
      derivsourceb(dw,wd,w,p,ix,iy,field,1,order);
     break;
     case b3:
      derivsourceb(dw,wd,w,p,ix,iy,field,2,order);
     break;
  }
  //return ( status);
}


__global__ void derivsource_parallel(struct params *p, real *w, real *wnew, real *wmod, 
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
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
               /*for(int f=rho; f<=b3; f++)               
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               computebdotv(wmod,wd,p,i,j);*/
               for(int f=rho; f<=b3; f++)
               {              
                  //if( (f==mom2) && (j==2))
                  //   ;//derivsource(dwn1,wd,wmod,p,i,j,f);
                  //else
                    derivsource(dwn1,wd,wmod,p,i,j,f,order);
                  //dwn1[fencode_ds(p,i,j,f)]=1.0;
                 // __syncthreads();
               }
               
               /*for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+0.5*dt*dwn1[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn2,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+0.5*dt*dwn2[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn3,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  wmod[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+dt*dwn3[fencode_ds(p,i,j,f)];
               computej(wmod,wd,p,i,j);
               computepk(wmod,wd,p,i,j);
               computept(wmod,wd,p,i,j);
               for(int f=rho; f<=b3; f++) 
                  deriv(dwn4,wd,wmod,p,i,j,f);
               
               for(int f=rho; f<=b3; f++) 
                  {
                  wnew[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+(dt/6.0)*(
                     dwn1[fencode_ds(p,i,j,f)]+2.0*dwn2[fencode_ds(p,i,j,f)]
                         +2.0*dwn3[fencode_ds(p,i,j,f)]+dwn4[fencode_ds(p,i,j,f)]);
               }*/
              //  __syncthreads();
              /* for(int f=rho; f<=b3; f++)
                   wnew[fencode_ds(p,i,j,f)]=w[fencode_ds(p,i,j,f)]+dt*dwn1[fencode_ds(p,i,j,f)];
               computej(wnew,wd,p,i,j);
               computepk(wnew,wd,p,i,j);
               computept(wnew,wd,p,i,j);*/ 


	}
 __syncthreads();
  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_ds(char *label)
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





int cuderivsource(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order)
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
     derivsource_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order);
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







