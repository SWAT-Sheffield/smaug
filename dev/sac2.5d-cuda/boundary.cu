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
#include "gradops_b.cuh"
__device__ __host__
void bc_cont(real *wt, struct params *p,int i, int j, int f) {

                if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i+2,j,f)];
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,j+2,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i+2,j,f)];                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,(j-3),f)];                  
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(i-3),j,f)];                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,j+2,f)];                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(i-3),j,f)];                   
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,(j-3),f)];                  
                }                       
                else if(i==0 || i==1)                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i+2,j,f)];                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(i-3),j,f)];                
                else if(j==0 || j==1)                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,j+2,f)];                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,(j-3),f)];
                




}

__device__ __host__
void bc_fixed(real *wt, struct params *p,int i, int j, int f, real val) {


                //(UPPER or LOWER)*NDIM*NVAR+dim*NVAR+varnum = picks out correct value for fixed BC
                //for array of values for fixed BC's

                if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_b(p,i,j,f)]=val;
                  else                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_b(p,i,j,f)]=val;                   
                  else                  
                    wt[fencode_b(p,i,j,f)]=val;                  
                }                       
                else if(i==0 || i==1)                
                  wt[fencode_b(p,i,j,f)]=val;                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[fencode_b(p,i,j,f)]=val;                
                else if(j==0 || j==1)                
                  wt[fencode_b(p,i,j,f)]=val;                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[fencode_b(p,i,j,f)]=val;
                




}

__device__ __host__
void bc_periodic(real *wt, struct params *p,int i, int j, int f) {

               if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(p->n[0])-3+i,(p->n[1])-3+j,f)];
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(p->n[0])-3+i,(p->n[1])-3+j,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(p->n[0])-3+i,2+((p->n[1])-j),f)];                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(p->n[0])-3+i,2+((p->n[1])-j),f)];                  
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,2+((p->n[0])-i),(p->n[1])-3+j,f)];                  
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,2+((p->n[0])-i),(p->n[1])-3+j,f)];                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,2+((p->n[0])-i),2+((p->n[1])-j),f)];                   
                  else                  
                    wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,2+((p->n[0])-i),2+((p->n[1])-j),f)];                  
                }                       
                else if(i==0 || i==1)                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,(p->n[0])-3+i,j,f)];                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,2+((p->n[0])-i),j,f)];                
                else if(j==0 || j==1)                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,(p->n[1])-3+j,f)];                
               else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[fencode_b(p,i,j,f)]=wt[fencode_b(p,i,2+((p->n[1])-j),f)];
                




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

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  



    j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
  if(i<p->n[0] && j<p->n[1])
	{

               //default continuous BC for all
               //gradient kept zero by copying variable values from edge of mesh to ghost cells
                  bc_cont(wmod,p,i,j,rho);
                  bc_cont(wnew,p,i,j,rho);
               //   bc_fixed(wmod,p,i,j,rho,1.0);
               //   bc_fixed(wnew,p,i,j,rho,1.0);
               //   bc_periodic(wmod,p,i,j,rho);
               //   bc_periodic(wnew,p,i,j,rho);

               
               for(int f=rho+1; f<NVAR; f++)
               {

                  bc_cont(wmod,p,i,j,f);
                  bc_cont(wnew,p,i,j,f);

                 // bc_fixed(wmod,p,i,j,f,val);
                 // bc_fixed(wnew,p,i,j,f,val);

                //  bc_periodic(wmod,p,i,j,f);
                //  bc_periodic(wnew,p,i,j,f);


               }

               for(int f=vel1; f<NDERV; f++)
               {
                  bc_cont(wd,p,i,j,f);

                // bc_fixed(wd,p,i,j,f,val);
                 //   bc_periodic(wd,p,i,j,f);

                  
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
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   //int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;
//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
 	    //printf("called prop\n"); 
    // cudaThreadSynchronize();
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wd, *d_wmod);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
	    //printf("called update\n"); 
    cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}

