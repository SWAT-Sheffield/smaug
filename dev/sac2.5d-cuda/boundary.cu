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
int encode_b (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( iy * ((dp)->ni) + ix);
}

__device__ __host__
int fencode_b (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return ( (iy * ((dp)->ni) + ix)+(field*((dp)->ni)*((dp)->nj)));
}


__global__ void boundary_parallel(struct params *p, real *w, real *wnew, real *wd, real *wmod)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;

  int ni=p->ni;
  int nj=p->nj;
  real dt=p->dt;
  real dy=p->dy;
  real dx=p->dx;
  real g=p->g;

  real *u,  *v,  *h;
  real *un,  *vn,  *hn;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(p->ni)*(p->nj)*rho;
  u=w+(p->ni)*(p->nj)*mom1;
  v=w+(p->ni)*(p->nj)*mom2;

  hn=wnew+(p->ni)*(p->nj)*rho;
  un=wnew+(p->ni)*(p->nj)*mom1;
  vn=wnew+(p->ni)*(p->nj)*mom2;

    j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->ni && j<p->nj)
	{

               //default continuous BC for all
               //gradient kept zero by copying variable values from edge of mesh to ghost cells
               
               for(int f=rho; f<=b3; f++)
               {
                   
                if(i==0 || i==1)
                {
                  wnew[fencode_b(p,i,j,f)]=wnew[fencode_b(p,2,j,f)];
                  wmod[fencode_b(p,i,j,f)]=wmod[fencode_b(p,2,j,f)];
                }
                if((i==((p->ni)-1)) || (i==((p->ni)-2)))
                {
                  wnew[fencode_b(p,i,j,f)]=wnew[fencode_b(p,((p->ni)-3),j,f)];
                  wmod[fencode_b(p,i,j,f)]=wmod[fencode_b(p,((p->ni)-3),j,f)];
                }
                if(j==0 || j==1)
                {
                  wnew[fencode_b(p,i,j,f)]=wnew[fencode_b(p,i,2,f)];
                  wmod[fencode_b(p,i,j,f)]=wmod[fencode_b(p,i,2,f)];
                }
                if((j==((p->nj)-1)) || (j==((p->nj)-2)))
                {
                  wnew[fencode_b(p,i,j,f)]=wnew[fencode_b(p,i,((p->nj)-3),f)];
                  wmod[fencode_b(p,i,j,f)]=wnew[fencode_b(p,i,((p->nj)-3),f)];
                }

                  
               }

               for(int f=current1; f<=divb; f++)
               {
                                      
                if(i==0 || i==1)
                  wd[fencode_b(p,i,j,f)]=wd[fencode_b(p,2,j,f)];
                if((i==((p->ni)-1)) || (i==((p->ni)-2)))
                  wd[fencode_b(p,i,j,f)]=wd[fencode_b(p,((p->ni)-3),j,f)];
                if(j==0 || j==1)
                  wd[fencode_b(p,i,j,f)]=wd[fencode_b(p,i,2,f)];
                if((j==((p->nj)-1)) || (j==((p->nj)-2)))
                  wd[fencode_b(p,i,j,f)]=wd[fencode_b(p,i,((p->nj)-3),f)];
                
                  
               }

		
               /*if(i==0 )
		{
			un[j*ni] = 2.5*un[1+j*ni] - 2*un[2+j*ni] + 0.5*un[3+j*ni];
			un[ni+j*ni] = 2.5*un[ni-1+j*ni] - 2*un[ni-2+ni*j] + 0.5*un[ni-3+j*ni];
			vn[j*ni] = 2.5*vn[1+j*ni] - 2*vn[2+j*ni] + 0.5*vn[3+j*ni];
		 	vn[ni+j*ni] = 2.5*vn[ni-1+j*ni] - 2*vn[ni-2+ni*j] + 0.5*vn[ni-3+j*ni];
		 	hn[j*ni] = 2.5*hn[1+j*ni] - 2*hn[2+j*ni] + 0.5*hn[3+j*ni];
			hn[ni+j*ni] = 2.5*hn[ni-1+j*ni] - 2*hn[ni-2+ni*j] + 0.5*hn[ni-3+j*ni];
		}

		if(j==0)
		{
			un[i+ni] = 2.5*un[i+1*ni] - 2*un[i+2*ni] + 0.5*un[i+3*ni];
			un[i+(nj )*ni] = 2.5*un[i+(nj-1)*ni] - 2*un[i+(nj-2)*ni] + 0.5*un[i+(nj-3)*ni];
			vn[i+ni] = 2.5*vn[i+1*ni] - 2*vn[i+2*ni] + 0.5*vn[i+3*ni];
			vn[i+(nj)*ni] = 2.5*vn[i+(nj-1)*ni] - 2*vn[i+(nj-2)*ni] + 0.5*vn[i+(nj-3)*ni];
			hn[i+ni] = 2.5*hn[i+1*ni] - 2*hn[i+2*ni] + 0.5*hn[i+3*ni];
			hn[i+(nj)*ni] = 2.5*hn[i+(nj-1)*ni] - 2*hn[i+(nj-2)*ni] + 0.5*hn[i+(nj-3)*ni];
		}*/
	}
 __syncthreads();
  
}

int cuboundary(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd)
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
 	    //printf("called prop\n"); 
    // cudaThreadSynchronize();
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wd, *d_wmod);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
	    //printf("called update\n"); 
    cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->ni)* ((*p)->nj)*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->ni)* ((*p)->nj))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}

