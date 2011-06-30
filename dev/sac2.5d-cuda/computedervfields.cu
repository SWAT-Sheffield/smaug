//#define MODID pre


#include "cudapars.h"
#include "paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "step.h"
#include "gradops_cdf.cuh"
#include "dervfields_cdf.cuh"
/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
__global__ void computevels_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  








  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         #endif
                        }


         }

}
              __syncthreads();











  
}


__global__ void computepres_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  








  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

	     #ifdef ADIABHYDRO
	       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #else
	       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #endif         
              /* switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))
				     {

				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

				     }
                         break;
                         #endif
                        }*/


         }

}
              __syncthreads();











  
}


__global__ void computemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  



 /*  for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
        }

}*/
              __syncthreads();



if(iindex==0)
{
   
 //  for(ipg=0;ipg<(p->npgp[0]);ipg++)
 //  for(jpg=0;jpg<(p->npgp[1]);jpg++)
  // {

  //   i=ip*(p->npgp[0])+ipg;
 //    j=jp*(p->npgp[1])+jpg;
   //if( i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
    //p->cmax=0.0;
    for(ii[0]=2;ii[0]<((p->n[0])-2);ii[0]++)
      for(ii[1]=2;ii[1]<((p->n[1])-2);ii[1]++)
     #ifdef USE_SAC_3D
        for(ii[2]>1;ii[2]<((p->n[2])-2);ii[2]++)
     #endif
	{ 
               computecmax3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);




	}

 //  }
}
 __syncthreads(); 





  
}

__global__ void computec_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  






 p->cmax=0.0;
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
        }

}
              __syncthreads();












  
}


__global__ void computedervfields_parallel(struct params *p,   real *wmod, real *wd, int order)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  





if(order == 0)
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		
              #ifdef ADIABHYDRO
              for(int f=vel1; f<NDERV; f++)
             #else 
               for(int f=vel1; f<=pkb; f++)
             #endif
                        wd[fencode3_cdf(p,ii,f)]=0; 
		#ifdef USE_SAC_3D
		  for(int f=rho; f<NVAR; f++)
                  	wmod[fencode3_cdf(p,ii,f)+dimp*NVAR]=wmod[fencode3_cdf(p,ii,f)]; 

		#else
		  for(int f=rho; f<NVAR; f++)
                  	wmod[fencode3_cdf(p,ii,f)+dimp*NVAR]=wmod[fencode3_cdf(p,ii,f)]; 
		#endif               

        }

}
               __syncthreads();



  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if( ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               
             #ifdef ADIABHYDRO
               //computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
             #else
               //computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computej3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

               computebdotv3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computedivb3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

             #endif

         }

}
              __syncthreads();

  
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cdf(char *label)
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




int cucomputedervfields(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   ////cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  

 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computedervfields_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputevels(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   ////cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computevels_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputemaxc(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}


int cucomputec(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computec_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputepres(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{

 int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computepres_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}







