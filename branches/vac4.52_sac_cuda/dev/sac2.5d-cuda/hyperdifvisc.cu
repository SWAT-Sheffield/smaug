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
#include "gradops_hdv.cuh"

__device__ __host__
void bc_periodic1_temp2(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 )                
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,5,j,f)];
                else if((i==((p->n[0])+1)) )                
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i-6,j,f)];
                else if(j==0  )                
                  wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,5,f)];
                else if((j==((p->n[1])+1)) )                
                  wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,j-6,f)];

 


}

__device__ __host__
void bc_periodic2_temp2(real *wt, struct params *p,int i, int j, int f) {


               if(i<1 && j<1)
                {
                  if(i==j)
                    //wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,5,f)];
                  else                  
                    //wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,5,j,f)];                                    
                }
                else if(i<1 && j>((p->n[1])+1))
                {
                  if(i==(j-(p->n[1])+1))                  
                    //wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,5,j,f)];                                     
                  else                  
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,j-6,f)];                                     
                }
                else if(i>((p->n[0])+1) && j<1)
                {
                  if((i-(p->n[0])+1)==j)                  
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i-6,j,f)];                                    
                  else                  
                   wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,5,f)];                                    
                }
                else if(i>((p->n[0])+1) && j>((p->n[1])+1))
                {
                  if(i==j)                  
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i,j-6,f)];                                    
                  else                  
                    wt[fencode_hdv(p,i,j,f)]=wt[fencode_hdv(p,i-6,j,f)];                                    
                }                       
                 
                




}



__global__ void hyperdifvisc_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int is,js;
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

  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*(p->n[0])*(p->n[1]);

switch(field)
{
    case rho:
      bfac1=1.0;
      bfac2=1.0;
      bfac3=0.0;
    break;
    case mom1:
    case mom2:
      bfac1=1.0;
      bfac2=0.0;
      bfac3=1.0;
    break;
    case energy:
      bfac1=1.0;
      bfac2=1.0;
      bfac3=0.0;
    break;
    case b1:
    case b2:
      bfac1=1.0;
      bfac2=1.0;
      bfac3=0.0;
    break;
}

//init temp1 and temp2 to zero 
//the compute element initialising n[0] or n[1] element must do +1 and +2
//this is because we fit the problem geometrically to nixnj elements 
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if(i<((p->n[0])) && j<((p->n[1])))
   {

        for(int f=d1r; f<=d3l; f++)
                 wtemp1[fencode_hdv(p,i,j,f)]=0;
      wtemp2[fencode_hdv(p,i,j,tmpnui)]=0;
      if(i==((p->n[0])-1))
      {
        for(int f=d1r; f<=d3l; f++)
                 wtemp1[fencode_hdv(p,i+1,j,f)]=0;
        wtemp2[fencode_hdv(p,i+1,j,tmpnui)]=0;
        wtemp2[fencode_hdv(p,i+2,j,tmpnui)]=0;
      }
      if(j==((p->n[1])-1))
      {
          for(int f=d1r; f<=d3l; f++)
                 wtemp1[fencode_hdv(p,i,j+1,f)]=0;
      }
      wtemp2[fencode_hdv(p,i,j+1,tmpnui)]=0;
      wtemp2[fencode_hdv(p,i,j+2,tmpnui)]=0;


   }

  }
   __syncthreads();




   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if(i<((p->n[0])) && j<((p->n[1])))
   {
        for(int f=tmp1; f<=tmp6; f++)
                 wtemp[fencode_hdv(p,i,j,f)]=0;


        //temp value for viscosity

       //tmp6  tmpnu
#ifdef USE_SAC
     if(field !=energy)
        wtemp[fencode_hdv(p,i,j,tmp6)]=wmod[fencode_hdv(p,i,j,field)+shift]/(bfac2+bfac3*((wmod[fencode_hdv(p,i,j,rho)+shift] +wmod[fencode_hdv(p,i,j,rhob)+shift])));

     else
        wtemp[fencode_hdv(p,i,j,tmp6)]=wmod[fencode_hdv(p,i,j,energy)+shift]-0.5*(wmod[fencode_hdv(p,i,j,b1)+shift]*wmod[fencode_hdv(p,i,j,b1)+shift]+wmod[fencode_hdv(p,i,j,b2)+shift]*wmod[fencode_hdv(p,i,j,b2)+shift])+(wmod[fencode_hdv(p,i,j,mom1)+shift]*wmod[fencode_hdv(p,i,j,mom1)+shift]+wmod[fencode_hdv(p,i,j,mom2)+shift]*wmod[fencode_hdv(p,i,j,mom2)+shift])/(wmod[fencode_hdv(p,i,j,rho)+shift]+wmod[fencode_hdv(p,i,j,rhob)+shift] );

#else
     if(field !=energy)
        wtemp[fencode_hdv(p,i,j,tmp6)]=wmod[fencode_hdv(p,i,j,field)+shift]/(bfac2+bfac3*(wmod[fencode_hdv(p,i,j,rho)+shift] ));

     else
        wtemp[fencode_hdv(p,i,j,tmp6)]=wmod[fencode_hdv(p,i,j,energy)+shift]-0.5*(wmod[fencode_hdv(p,i,j,b1)+shift]*wmod[fencode_hdv(p,i,j,b1)+shift]+wmod[fencode_hdv(p,i,j,b2)+shift]*wmod[fencode_hdv(p,i,j,b2)+shift])+(wmod[fencode_hdv(p,i,j,mom1)+shift]*wmod[fencode_hdv(p,i,j,mom1)+shift]+wmod[fencode_hdv(p,i,j,mom2)+shift]*wmod[fencode_hdv(p,i,j,mom2)+shift])/(wmod[fencode_hdv(p,i,j,rho)+shift] );

#endif
        wd[fencode_hdv(p,i,j,hdnur)]=0;
        wd[fencode_hdv(p,i,j,hdnul)]=0;
   }

}
   __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if( i<((p->n[0])) && j<((p->n[1])))
   {
     //tmp6 is tmp_nuI
     wtemp2[fencode_hdv(p,i+1,j+1,tmpnui)]=wtemp[fencode_hdv(p,i,j,tmp6)];

   }

   }
   __syncthreads();
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
          {

                  bc_periodic1_temp2(wtemp2,p,i,j,tmpnui);
                  if(i==((p->n[0])-1))
                  {
                  bc_periodic1_temp2(wtemp2,p,i+1,j,tmpnui);
                  bc_periodic1_temp2(wtemp2,p,i+2,j,tmpnui);


                  }

                  if(j==((p->n[1])-1))
                  {
                  bc_periodic1_temp2(wtemp2,p,i,j+1,tmpnui);
                  bc_periodic1_temp2(wtemp2,p,i,j+2,tmpnui);
                   }

      
          }

}
                __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
          {
                  //bc_cont_cd1(dwn1,p,i,j,f1+fid);
                  bc_periodic2_temp2(wtemp2,p,i,j,tmpnui);
                  if(i==((p->n[0])-1))
                  {
                  bc_periodic2_temp2(wtemp2,p,i+1,j,tmpnui);
                  bc_periodic2_temp2(wtemp2,p,i+2,j,tmpnui);


                  }

                  if(j==((p->n[1])-1))
                  {
                  bc_periodic1_temp2(wtemp2,p,i,j+1,tmpnui);
                  bc_periodic1_temp2(wtemp2,p,i,j+2,tmpnui);
                   }
           }
}
                __syncthreads();

   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l




   //tmp1  tmp_nuI
 
//compute d3r and d1r
   //tmp2  d3r
    //tmp3 d1r

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
 
   //if(i>1 && j>1 && i<((p->n[0])-1) && j<((p->n[1])-1))
   if(i<((p->n[0])-1) && j<((p->n[1])-1))            
   { 
     if(hand==0)
     {
           wtemp1[fencode_hdv(p,i+2,j+2,d3r)]=fabs(3.0*(wtemp2[fencode_hdv(p,i+2+(dim==0),j+2+(dim==1),tmpnui)] - wtemp2[fencode_hdv(p,i+2,j+2,tmpnui)] ) - (wtemp2[fencode_hdv(p,i+2*(dim==0),j+2*(dim==1),tmpnui)] - wtemp2[fencode_hdv(p,i+2-(dim==0),j+2-(dim==1),tmpnui)]    ));
     }
     else
     {
           wtemp1[fencode_hdv(p,i+2,j+2,d3l)]=fabs(3.0*(wtemp2[fencode_hdv(p,i+2+(dim==0),j+2+(dim==1),tmpnui)] - wtemp2[fencode_hdv(p,i+2,j+2,tmpnui)]) - (wtemp2[fencode_hdv(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv(p,i+2-(dim==0),j-(dim==1),tmpnui)]    ));
     }
   }
}
   __syncthreads();
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

   if(i<((p->n[0])+(dim==0)-1) && j<((p->n[1])+(dim==1)-1))            
   { 
     if(hand==0)
     {

           wtemp1[fencode_hdv(p,i,j,d1r)]=fabs((wtemp2[fencode_hdv(p,i+1+(dim==0),j+1+(dim==1),tmpnui)] - wtemp2[fencode_hdv(p,i+1,j+1,tmpnui)] ));
     }
     else
     {
           wtemp1[fencode_hdv(p,i,j,d1l)]=fabs((wtemp2[fencode_hdv(p,i+1,j+1,tmpnui)] - wtemp2[fencode_hdv(p,i+1-(dim==0),j+1-(dim==1),tmpnui)] ));
     }
   }
}
   __syncthreads();



/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic1_hdv(wtemp,p,i,j,tmp2);
                  bc_periodic1_hdv(wtemp,p,i,j,tmp3);
             }

}
                __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic2_hdv(wtemp,p,i,j,tmp2);
                  bc_periodic2_hdv(wtemp,p,i,j,tmp3);
             }
}
                __syncthreads();*/



  //compute md3r and md1r
//tmp4    md3r
//tmp5    md1r
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  if( i<((p->n[0])) && j<((p->n[1])))            
   {
         maxt=0;
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                {
                   if(wtemp1[fencode_hdv(p,i+is,j+js,d3r+hand)]>maxt)
                         maxt=wtemp1[fencode_hdv(p,i+is,j+js,d3r+hand)];

                }
          wtemp[fencode_hdv(p,i,j,tmp4)]=maxt;

         maxt=0;
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                {
                   if(wtemp1[fencode_hdv(p,i+is,j+js,d1r+hand)]>maxt)
                        maxt=wtemp1[fencode_hdv(p,i+is,j+js,d1r+hand)];

                }
          wtemp[fencode_hdv(p,i,j,tmp5)]=maxt;
   }
}
   __syncthreads();





 /*  for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic1_hdv(wtemp,p,i,j,tmp4);
                  bc_periodic1_hdv(wtemp,p,i,j,tmp5);
             }

}
                __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic2_hdv(wtemp,p,i,j,tmp4);
                  bc_periodic2_hdv(wtemp,p,i,j,tmp5);
             }
}
                __syncthreads();*/



   p->maxviscoef=0;


    //finally update nur and nul
//tmp4    md3r
//tmp5    md1r
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;


   if(i<((p->n[0])) && j<((p->n[1])))
   {
     if(wtemp[fencode_hdv(p,i,j,tmp5)]>0)
	wd[fencode_hdv(p,i,j,hdnur+hand)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1]))*(p->cmax)*(p->chyp[field])*wtemp[fencode_hdv(p,i,j,tmp4)]/wtemp[fencode_hdv(p,i,j,tmp5)];


     else
        wd[fencode_hdv(p,i,j,hdnur+hand)]=0;


   }
}
 __syncthreads();



 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv(char *label)
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





int cuhyperdifvisc(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand)
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
     hyperdifvisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_hdv,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_hdv,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
// cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_hdv, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}







