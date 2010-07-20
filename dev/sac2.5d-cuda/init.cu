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
int encode_i (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_i (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}



__global__ void init_parallel(struct params *p, float *w, float *wnew, float *b, float *wmod, 
    float *dwn1, float *wd)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int index,k;
int ni=p->ni;
  int nj=p->nj;

// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  float *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

 int nli = 0.45*(p->ni-1)+1;
  int nui = 0.55*(p->ni-1)+1;
  int nlj = 0.45*(p->nj-1)+1;
  int nuj = 0.55*(p->nj-1)+1; 
  int i,j;
   
   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->ni && j<p->nj)
	{
		b[i+j*(p->ni)]=0;

                 //Define b	
       for(int f=0; f<=6; f++)
        { 
                  wd[fencode_i(p,i,j,f)]=0;
        }

        for(int f=rho; f<=b3; f++)
        {               
                  wnew[fencode_i(p,i,j,f)]=w[fencode_i(p,i,j,f)];
                  dwn1[fencode_i(p,i,j,f)]=0;
                  //dwn2[fencode(p,i,j,f)]=0;
                 // dwn3[fencode(p,i,j,f)]=0;
                  //dwn4[fencode(p,i,j,f)]=0;
                 
        }

	 __syncthreads();

			}	
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



int cuinit(struct params **p, float **w, float **wnew,  float **b, struct state **state, struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd, struct state **d_state)
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
  float *adb;
  float *adw, *adwnew;
  struct params *adp;
  struct state *ads;


  cudaMalloc((void**)d_wmod, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)d_dwn1, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)d_wd, 7*((*p)->ni)* ((*p)->nj)*sizeof(float));

  cudaMalloc((void**)&adw, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adwnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adb, 1*(((*p)->ni)* ((*p)->nj))*sizeof(float));
  cudaMalloc((void**)&adp, sizeof(struct params));
  cudaMalloc((void**)&ads, sizeof(struct state));
  checkErrors_i("memory allocation");

printf("ni is %d\n",(*p)->nj);

    *d_b=adb;
    *d_p=adp;
    *d_w=adw;
    *d_wnew=adwnew;
    *d_state=ads;


    cudaMemcpy(*d_w, *w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_wnew, *wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_b, *b, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_state, *state, sizeof(struct state), cudaMemcpyHostToDevice);
    
    dim3 dimBlock(16, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;
   

    printf("calling initialiser\n");
     //init_parallel(struct params *p, float *b, float *u, float *v, float *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b, *d_wmod, *d_dwn1,  *d_wd);
     cudaThreadSynchronize();
	    printf("called initialiser\n");
	cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);

	//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(float), cudaMemcpyDeviceToHost);




  return 0;



}


