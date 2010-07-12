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
__global__ void init_parallel(struct params *p, float *b, float *u, float *v, float *h)
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
    


 int nli = 0.45*(p->ni-1)+1;
  int nui = 0.55*(p->ni-1)+1;
  int nlj = 0.45*(p->nj-1)+1;
  int nuj = 0.55*(p->nj-1)+1; 
  int i,j;
   
   j=iindex/ni;
   i=iindex-j*ni;
  if(i<p->ni && j<p->nj)
	{
		b[j+i*(p->nj)]=0;
		//initialise the arrays here
               for(k=0;k<2;k++)
      		{
                    index=i+j*(p->ni)+(k*(p->nj)*(p->ni));
		    u[index]=0;
		    v[index]=0;
		    h[index]=0;
      		}
		h[i+j*(p->ni)]=5000;
	}
        __syncthreads();
        if(i>=nli && i<=nui && j>=nlj && j<=nuj)
	{
	   
           h[iindex]=5030;	
	}
	 __syncthreads();
//Define b	
  if(i<p->ni && j<p->nj)
	{
		
		if(i*p->dx >20001)
		      b[iindex]=0;
		else if(i*p->dx <20000)
			b[iindex]=(5000/20000)*(20000-i*p->dx);		
	}	
	 __syncthreads();
  
}



__global__ void prop_parallel(struct params *p, float *b, float *u, float *v, float *h)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;

  if(i>0 && j >0 && i<((p->ni)-1) && j<((p->nj)-1))
	{
		//update the arrays here
               u[i+j*ni+ni*nj] = ((u[i+1+j*ni] + u[i-1+j*ni] + u[i+(j+1)*ni] + u[i+(j-1)*ni])/4)- 0.5*(dt/dx)*(((u[(i+1)+ni*j])*u[(i+1)+ni*j])/2) - ((u[i-1+j*ni]*u[i-1+j*ni])/2) - 0.5*(dt/dy)*(v[i+j*ni])*(u[i+(j+1)*ni] - u[i+(j-1)*ni]) - 0.5*g*(dt/dx)*(h[i+1+j*ni]-h[i-1+j*ni]);

v[i+j*ni+ni*nj] = ((v[i+1+j*ni] + v[i-1+j*ni] + v[i+(j+1)*ni] + v[i+(j-1)*ni])/4)- 0.5*(dt/dy)*((v[i+ni*(j+1)]*v[(i)+ni*(j+1)])/2 - (v[i+(j-1)*ni]*v[i+(j-1)*ni])/2) - 0.5*(dt/dx)*(u[i+j*ni])*(v[i+1+j*ni] - v[i-1+j*ni]) - 0.5*g*(dt/dy)*(h[i+(j+1)*ni]-h[i+(j-1)*ni]);

h[i+j*ni+ni*nj] = ((h[i+1+j*ni] + h[i-1+j*ni] + h[i+(j+1)*ni] + h[i+(j-1)*ni])/4)- 0.5*(dt/dx)*(u[i+j*ni])*((h[i+1+j*ni]-b[i+1+j*ni]) - (h[i-1+j*ni]-b[i-1+j*ni])) - 0.5*(dt/dy)*(v[i+j*ni])*((h[i+(j+1)*ni]-b[i+(j+1)*ni]) - (h[i+(j-1)*ni]-b[i+(j-1)*ni])) - 0.5*(dt/dx)*(h[i+j*ni]-b[i+j*ni])*(u[i+1+j*ni]- u[i-1+j*ni])- 0.5*(dt/dy)*(h[i+j*ni]-b[i+j*ni])*(v[i+(j+1)*ni] - v[i+(j-1)*ni]);

	}
  
}

__global__ void boundary_parallel(struct params *p, float *b, float *u, float *v, float *h)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;

  if(i<p->ni && j<p->nj)
	{

		if(i=0)
		{
			u[1+j*ni+ni*nj] = 2.5*u[2+j*ni+ni*nj] - 2*u[3+j*ni+ni*nj] + 0.5*u[4+j*ni+ni*nj];
			u[ni-1+j*ni+ni*nj] = 2.5*u[ni-1+j*ni+ni*nj] - 2*u[ni-2+ni*j+ni*nj] + 0.5*u[ni-3+j*ni+ni*nj];
			v[1+j*ni+ni*nj] = 2.5*v[2+j*ni+ni*nj] - 2*v[3+j*ni+ni*nj] + 0.5*v[4+j*ni+ni*nj];
		 	v[ni-1+j*ni+ni*nj] = 2.5*v[ni-1+j*ni+ni*nj] - 2*v[ni-2+ni*j+ni*nj] + 0.5*v[ni-3+j*ni+ni*nj];
		 	h[1+j*ni+ni*nj] = 2.5*h[2+j*ni+ni*nj] - 2*h[3+j*ni+ni*nj] + 0.5*h[4+j*ni+ni*nj];
			h[ni-1+j*ni+ni*nj] = 2.5*h[ni-1+j*ni+ni*nj] - 2*h[ni-2+ni*j+ni*nj] + 0.5*h[ni-3+j*ni+ni*nj];
		}

		if(j=0)
		{
			u[i+ni+ni*nj] = 2.5*u[i+2*ni+ni*nj] - 2*u[i+3*ni+ni*nj] + 0.5*u[i+4*ni+ni*nj];
			u[i+(nj-1)*ni+ni*nj] = 2.5*u[i+(nj-2)*ni+ni*nj] - 2*u[i+(nj-3)*ni+ni*nj] + 0.5*u[i+(nj-4)*ni+ni*nj];
			v[i+ni+ni*nj] = 2.5*v[i+2*ni+ni*nj] - 2*v[i+3*ni+ni*nj] + 0.5*v[i+4*ni+ni*nj];
			v[i+(nj-1)*ni+ni*nj] = 2.5*v[i+(nj-2)*ni+ni*nj] - 2*v[i+(nj-3)*ni+ni*nj] + 0.5*v[i+(nj-4)*ni+ni*nj];
			h[i+ni+ni*nj] = 2.5*h[i+2*ni+ni*nj] - 2*h[i+3*ni+ni*nj] + 0.5*h[i+4*ni+ni*nj];
			h[i+(nj-1)*ni+ni*nj] = 2.5*h[i+(nj-2)*ni+ni*nj] - 2*h[i+(nj-3)*ni+ni*nj] + 0.5*h[i+(nj-4)*ni+ni*nj];
		}
	}
  
}

__global__ void update_parallel(struct params *p, float *b, float *u, float *v, float *h)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;


  int ni=p->ni;
  int nj=p->nj;
  float dt=p->dt;
  float dy=p->dy;
  float dx=p->dx;
  float g=p->g;



  if(i<p->ni && j<p->nj)
	{
            u[i+(nj-1)*ni]=u[i+(nj-1)*ni+ni*nj];
            v[i+(nj-1)*ni]=v[i+(nj-1)*ni+ni*nj];
	    h[i+(nj-1)*ni]=h[i+(nj-1)*ni+ni*nj];
	}
  
}
/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
__global__ void saxpy_parallel(int n, float alpha, float *x, float *y)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  // except for special cases, the total number of threads in all blocks
  // adds up to more than the vector length n, so this conditional is
  // EXTREMELY important to avoid writing past the allocated memory for
  // the vector y.
  if (i<n)
    y[i] = alpha*x[i] + y[i];
}

/////////////////////////////////////
// kernel function (CPU)
/////////////////////////////////////
void saxpy_serial(int n, float alpha, float *x, float *y)
{
  int i;
  for (i=0; i<n; i++)
    y[i] = alpha*x[i] + y[i];
}
/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors(char *label)
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

int stepfunc()
{
 /////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
  {
    fprintf(stderr, "Sorry, no CUDA device fount");
    return 1;
  }
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }
  cudaSetDevice(selectedDevice);
  checkErrors("initialisations");
  

  
  /////////////////////////////////////
  // (2) allocate memory on host (main CPU memory) and device,
  //     h_ denotes data residing on the host, d_ on device
  /////////////////////////////////////
  float *h_x = (float*)malloc(N*sizeof(float));
  float *h_y = (float*)malloc(N*sizeof(float));
  float *d_x;
  cudaMalloc((void**)&d_x, N*sizeof(float));
  float *d_y;
  cudaMalloc((void**)&d_y, N*sizeof(float));
  checkErrors("memory allocation");



  /////////////////////////////////////
  // (3) initialise data on the CPU
  /////////////////////////////////////
  int i;
  for (i=0; i<N; i++)
  {
    h_x[i] = 1.0f + i;
    h_y[i] = (float)(N-i+1);
  }



  /////////////////////////////////////
  // (4) copy data to device
  /////////////////////////////////////
  cudaMemcpy(d_x, h_x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y, N*sizeof(float), cudaMemcpyHostToDevice);
  checkErrors("copy data to device");



  /////////////////////////////////////
  // (5) perform computation on host (to enable result comparison later)
  /////////////////////////////////////
  saxpy_serial(N, 2.0f, h_x, h_y);



  /////////////////////////////////////
  // (6) perform computation on device
  //     - we use numThreadsPerBlock threads per block
  //     - the total number of blocks is obtained by rounding the
  //       vector length N up to the next multiple of numThreadsPerBlock
  /////////////////////////////////////
  int numBlocks = (N+numThreadsPerBlock-1) / numThreadsPerBlock;
  saxpy_parallel<<<numBlocks, numThreadsPerBlock>>>(N, 2.0, d_x, d_y);
  checkErrors("compute on device");



  /////////////////////////////////////
  // (7) read back result from device into temp vector
  /////////////////////////////////////
  float *h_z = (float*)malloc(N*sizeof(float));
  cudaMemcpy(h_z, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);
  checkErrors("copy data from device");

  
  /////////////////////////////////////
  // (8) perform result comparison
  /////////////////////////////////////
  int errorCount = 0;
  for (i=0; i<N; i++)
  {
    if (abs(h_y[i]-h_z[i]) > 1e-6)
      errorCount = errorCount + 1;
  }
  if (errorCount > 0)
    printf("Result comparison failed.\n");
  else
    printf("Result comparison passed.\n");



  /////////////////////////////////////
  // (9) clean up, free memory
  /////////////////////////////////////
  free(h_x);
  free(h_y);
  free(h_z);
  cudaFree(d_x);
  cudaFree(d_y);
  return 0;

}

int cuinit(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h)
{
/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
  {
    fprintf(stderr, "Sorry, no CUDA device fount");
    return 1;
  }
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }
  //cudaSetDevice(selectedDevice);
  checkErrors("initialisations");
  
	// Build empty u, v, b matrices

  printf("in cuinit\n");
  float *adu,*adv,*adh,*adb;
  struct params *adp;

  cudaMalloc((void**)&adu, ((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adv, ((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adh, ((*p)->ni)* ((*p)->nj)*sizeof(float));
  cudaMalloc((void**)&adb, ((*p)->ni)*sizeof(float));
  cudaMalloc((void**)&adp, sizeof(struct params));
  checkErrors("memory allocation");

printf("ni is %d\n",(*p)->nj);

    *d_u=adu;
    *d_v=adv;
    *d_h=adh;
    *d_b=adb;
    *d_p=adp;


    cudaMemcpy(*d_u, *u, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_v, *v, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_h, *h, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_b, *b, ((*p)->ni)*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
    
    dim3 dimBlock(16, 1);
    //dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);
   int numBlocks = (((*p)->ni)*((*p)->nj)+numThreadsPerBlock-1) / numThreadsPerBlock;
   

    printf("calling initialiser\n");
     //init_parallel(struct params *p, float *b, float *u, float *v, float *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
     cudaThreadSynchronize();
	    printf("called initialiser\n");

	cudaMemcpy(*u, *d_u, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*v, *d_v, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*b, *d_b, ((*p)->ni)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*h, *d_h, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);




  return 0;



}


int cuprop(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h)
{


printf("calling propagate solution\n");

    dim3 dimBlock(blocksize, blocksize);
    dim3 dimGrid(((*p)->ni)/dimBlock.x,((*p)->nj)/dimBlock.y);


     //init_parallel(struct params *p, float *b, float *u, float *v, float *h)
     prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    printf("called prop\n"); 
     cudaThreadSynchronize();
     boundary_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    printf("called boundary\n");  
     cudaThreadSynchronize();
     update_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    printf("called update\n"); 
    cudaThreadSynchronize();
 cudaMemcpy(*u, *d_u, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*v, *d_v, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*b, *d_b, ((*p)->ni)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*h, *d_h, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);

  checkErrors("copy data from device");


 


}

int cufinish(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h)
{
  

 cudaMemcpy(*u, *d_u, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*v, *d_v, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*b, *d_b, ((*p)->ni)*sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(*h, *d_h, ((*p)->ni)* ((*p)->nj)*sizeof(float), cudaMemcpyDeviceToHost);

  checkErrors("copy data from device");


  cudaFree(*d_p);

  cudaFree(*d_u);
  cudaFree(*d_v);
  cudaFree(*d_h);
  cudaFree(*d_b);


}
