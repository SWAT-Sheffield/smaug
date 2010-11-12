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
#include "gradops_i.cuh"






__device__ __host__
void init_alftest (real *w, struct params *p,int i, int j) {
  int seg1,seg2,seg3,seg4;
  int width=10;
  real m2max=0.001;
  real start=((p->n[0])-width)/2;
  //seg1=((p->n[0])/3)-1;
  seg1=(p->n[0])/6;
  seg2=((p->n[0])/3);
  seg3=(2*(p->n[0])/3)-1;
  //seg4=(2*(p->n[0])/3);
  seg4=(p->n[0])-seg1;

	#ifdef USE_SAC

		    w[fencode_i(p,i,j,rhob)]=1.0;
		    w[fencode_i(p,i,j,b1b)]=1.0;
		    w[fencode_i(p,i,j,energy)]=0.01;

       #else
		    w[fencode_i(p,i,j,rho)]=1.0;
		    w[fencode_i(p,i,j,b1)]=1.0;
		    w[fencode_i(p,i,j,energy)]=0.01;

       #endif

		    //w[fencode_i(p,i,j,b1)]=15*j;
		    //w[fencode_i(p,i,j,b3)]=150*j;
		    
		   //if (i > seg2)
		    //if (i < seg3)
                   // if (i < seg1)
		   //   w[fencode_i(p,i,j,mom2)]=0.0;

		   if (i > seg1)
		    if (i < seg2)
		      w[fencode_i(p,i,j,mom2)]=m2max*(i-seg1)/(seg2-seg1);

		   if (i > seg2)
		    if (i < seg3)
		      //w[fencode_i(p,i,j,mom2)]=m2max*(i-seg2)/(seg3-seg2);
                      w[fencode_i(p,i,j,mom2)]=m2max;
		   if (i > seg3)
		    if (i < seg4)
		      w[fencode_i(p,i,j,mom2)]=m2max*(seg4-i)/(seg4-seg3);






}


__device__ __host__
void init_ozttest (real *w, struct params *p,int i, int j) {
                    
                    //real b0=1.0/sqrt(4.0*PI);
                    real b0=1.0;
                    //real ptot=5.0/(12.0*PI);
                    real ptot=5.0/3.0;
                    real rrho;
                    real rgamm1;


	#ifdef USE_SAC
                    b0=1.0/sqrt(4.0*PI);
                    ptot=5.0/(12.0*PI);
		    w[fencode_i(p,i,j,rhob)]=25.0/(36.0*PI);
                    //w[fencode_i(p,i,j,rhob)]=25.0/9.0;

                    rgamm1=1.0/((p->gamma)-1);
		    
		    w[fencode_i(p,i,j,b2b)]=-b0*sin((2.0*PI*p->dx[1])*j);
		    w[fencode_i(p,i,j,b1b)]=b0*sin(4.0*PI*(p->dx[0])*i);

                    //vx=-sin(2pi y)
                    //vy=sin(2pi x)
		    //w[fencode_i(p,i,j,mom1)]=-w[fencode_i(p,i,j,rhob)]*sin(2.0*PI*j*(p->dx[1]));
                    //w[fencode_i(p,i,j,mom2)]=w[fencode_i(p,i,j,rhob)]*sin(2.0*PI*j*(p->dx[0]));

		    w[fencode_i(p,i,j,mom1)]=-w[fencode_i(p,i,j,rhob)]*sin(2.0*PI*i*(p->dx[1]));
                    w[fencode_i(p,i,j,mom2)]=w[fencode_i(p,i,j,rhob)]*sin(2.0*PI*j*(p->dx[0]));


                    w[fencode_i(p,i,j,energyb)]=ptot*rgamm1;

       #else
		    //w[fencode_i(p,i,j,rho)]=25.0/(36.0*PI);
                    w[fencode_i(p,i,j,rho)]=25.0/9.0;
		    //w[fencode_i(p,i,j,b1)]=-b0*sin((p->dx[0])*i);
		    //w[fencode_i(p,i,j,b2)]=b0*sin(2.0*(p->dx[1])*j);
		    //w[fencode_i(p,i,j,b1)]=b0*sin((2.0*p->dx[1])*(p->n[0] -j));
		    //w[fencode_i(p,i,j,b2)]=-b0*sin(1.0*(p->dx[0])*(p->n[1] -i));
		    w[fencode_i(p,i,j,b1)]=b0*sin((2.0*p->dx[1])*j);
		    w[fencode_i(p,i,j,b2)]=-b0*sin(1.0*(p->dx[0])*i);

		    //w[fencode_i(p,i,j,b3)]=0.0;

                    //vx=-sin(2pi y)
                    //vy=sin(2pi x)
		    //w[fencode_i(p,i,j,mom1)]=-w[fencode_i(p,i,j,rho)]*sin(2.0*PI*j*(p->dx[1]));
                    //w[fencode_i(p,i,j,mom2)]=w[fencode_i(p,i,j,rho)]*sin(2.0*PI*j*(p->dx[0]));
		    //w[fencode_i(p,i,j,mom1)]=-w[fencode_i(p,i,j,rho)]*sin(i*(p->dx[0]));
                    //w[fencode_i(p,i,j,mom2)]=w[fencode_i(p,i,j,rho)]*sin(j*(p->dx[1]));
		    w[fencode_i(p,i,j,mom1)]=w[fencode_i(p,i,j,rho)]*sin(1.0*j*(p->dx[1]));
                    w[fencode_i(p,i,j,mom2)]=-w[fencode_i(p,i,j,rho)]*sin(1.0*i*(p->dx[0]));

		    //w[fencode_i(p,i,j,mom3)]=0;

                    //p=5/12pi  use this to determine the energy
                    //p=(gamma -1)*(e-0.5 rho v**2 - b**2/2)
                    rrho=1.0/w[fencode_i(p,i,j,rho)];
                    rgamm1=1.0/((p->gamma)-1);
		    //w[fencode_i(p,i,j,energy)]=(ptot/((p->gamma)-1))+0.5*rrho*(w[fencode_i(p,i,j,mom1)]*w[fencode_i(p,i,j,mom1)]+w[fencode_i(p,i,j,mom2)]*w[fencode_i(p,i,j,mom2)])+0.5*(w[fencode_i(p,i,j,b1)]*w[fencode_i(p,i,j,b1)]+w[fencode_i(p,i,j,b2)]*w[fencode_i(p,i,j,b2)]);
//w[fencode_i(p,i,j,energy)]=(ptot/((p->gamma)-1))+0.5*rrho;
		    //w[fencode_i(p,i,j,energy)]=(ptot-(p->gamma)*0.5*(w[fencode_i(p,i,j,b1)]*w[fencode_i(p,i,j,b1)]+w[fencode_i(p,i,j,b2)]*w[fencode_i(p,i,j,b2)]))*rgamm1+0.5*rrho*(w[fencode_i(p,i,j,mom1)]*w[fencode_i(p,i,j,mom1)]+w[fencode_i(p,i,j,mom2)]*w[fencode_i(p,i,j,mom2)]);
//w[fencode_i(p,i,j,energy)]=(ptot/((p->gamma)-1))+0.5*rrho*(w[fencode_i(p,i,j,mom1)]*w[fencode_i(p,i,j,mom1)]+w[fencode_i(p,i,j,mom2)]*w[fencode_i(p,i,j,mom2)])+0.5*(w[fencode_i(p,i,j,b1)]*w[fencode_i(p,i,j,b1)]+w[fencode_i(p,i,j,b2)]*w[fencode_i(p,i,j,b2)]);

//gives agreement with vac ozt
w[fencode_i(p,i,j,energy)]=-sin(1.0*i*(p->dx[0]));

       #endif





}


//*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd

__global__ void init_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];

// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;

   
   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->n[0] && j<p->n[1])
	{
		//b[i+j*(p->n[0])]=0;

                 //Define b	

 


	//apply this special condition
	//initiate alfven wave propagtion 
	//if no initial config read
	if(p->readini==0)
	{
	    for(int f=0; f<NVAR; f++)
            { 
		          w[fencode_i(p,i,j,f)]=0;
                          for(ord=0;ord<(2+3*(p->rkon==1));ord++)
                              wmod[fencode_i(p,i,j,f)+ord*NVAR*(p->n[0])*(p->n[1])]=0;
	    }
            w[fencode_i(p,i,j,rho)]=1.0;
            #ifdef ADIABHYDRO
		    if(i> (((p->n[0])/2)-2) && i<(((p->n[0])/2)+2) && j>(((p->n[1])/2)-2) && j<(((p->n[1])/2)+2) ) 
				w[fencode_i(p,i,j,rho)]=1.3;
            #else
                   // init_alftest (real *w, struct params *p,int i, int j)
                   // init_alftest(w,p,i,j);
                   // init_ozttest (real *w, struct params *p,int i, int j)
                    init_ozttest(w,p,i,j);
           #endif

	}


//	 __syncthreads();

			}	
	 __syncthreads();

  if(i<p->n[0] && j<p->n[1])
	{
        for(int f=rho; f<NVAR; f++)
        {               
                  wmod[fencode_i(p,i,j,f)]=w[fencode_i(p,i,j,f)];              
                  dwn1[fencode_i(p,i,j,f)]=0;
                            
        }

        for(int f=tmp1; f<NTEMP; f++)
                 wtemp[fencode_i(p,i,j,f)]=0;


}

 __syncthreads();
        if(i<p->n[0] && j<p->n[1])
               for(int f=vel1; f<NDERV; f++)
                    wd[fencode_i(p,i,j,f)]=0.0;

 __syncthreads(); 
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_i(char *label)
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



int cuinit(struct params **p, real **w, real **wnew, struct state **state, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp)
{



/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
   
 // if (deviceCount == 0)
 // {
 //   fprintf(stderr, "Sorry, no CUDA device fount");
 //   return 1;
//  }
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }
  //cudaSetDevice(selectedDevice);
  printf("device count %d selected %d\n", deviceCount,selectedDevice);
  checkErrors_i("initialisations");
  
	// Build empty u, v, b matrices

  printf("in cuinit\n");
 // real *adb;
  real *adw, *adwnew;
  struct params *adp;
  struct state *ads;

if(((*p)->rkon)==1)
  cudaMalloc((void**)d_wmod, 5*NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real));
else
  cudaMalloc((void**)d_wmod, 2*NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real));

  cudaMalloc((void**)d_dwn1, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real));
  cudaMalloc((void**)d_wd, NDERV*((*p)->n[0])* ((*p)->n[1])*sizeof(real));
  cudaMalloc((void**)d_wtemp, NTEMP*((*p)->n[0])* ((*p)->n[1])*sizeof(real));

  cudaMalloc((void**)&adw, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real));
  cudaMalloc((void**)&adwnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real));
  
  cudaMalloc((void**)&adp, sizeof(struct params));
  cudaMalloc((void**)&ads, sizeof(struct state));
  checkErrors_i("memory allocation");

printf("ni is %d\n",(*p)->n[1]);

   // *d_b=adb;
    *d_p=adp;
    *d_w=adw;
    *d_wnew=adwnew;
    *d_state=ads;

     
printf("allocating\n");
    cudaMemcpy(*d_w, *w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyHostToDevice);
   // cudaMemcpy(*d_wnew, *wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyHostToDevice);
    
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_state, *state, sizeof(struct state), cudaMemcpyHostToDevice);
    
    dim3 dimBlock(16, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
   

    printf("calling initialiser\n");
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp);
     cudaThreadSynchronize();
	    printf("called initialiser\n");
	cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);

	cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
        cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

        // printf("mod times step %f %f\n",(*p)->dt, ((*wnew)[10+16*((*p)->n[0])+((*p)->n[0])*((*p)->n[1])*b1]));



  return 0;



}


