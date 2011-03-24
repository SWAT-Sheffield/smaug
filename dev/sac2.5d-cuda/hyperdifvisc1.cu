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
#include "gradops_hdv1.cuh"

__device__ __host__
void bc_hyperdif(real *wt, struct params *p,int i, int j, int f,int dir) {




 
   if(  (dir == 0) && (i==(p->n[0])-1)   && j>0   && j<(p->n[1])           )
         wt[fencode_hdv1(p,i+2,j,f)]=wt[fencode_hdv1(p,(p->n[0])-5,j+1,f)];
   else if((dir == 1) && (j==(p->n[1])-1)    && i>0   && i<((p->n[0]))  )
       wt[fencode_hdv1(p,i,j+2,f)]=wt[fencode_hdv1(p,i+1,(p->n[1])-5,f)];
  else if((dir == 0) && (i==0)    && j>0   && j<((p->n[1]))   )
       wt[fencode_hdv1(p,0,j+1,f)]=wt[fencode_hdv1(p,6,j+1,f)];
   else if((dir == 1) && (j==0)    && i>0   && i<((p->n[0]))   )
       wt[fencode_hdv1(p,i+1,0,f)]=wt[fencode_hdv1(p,i+1,6,f)];
 
}


__device__ __host__
void bc_periodic1_temp2(real *wt, struct params *p,int i, int j, int f) {

                if(i==1 )                
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];
                else if((i==((p->n[0]))) )                
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-4,j,f)];
                else if(j==1  )                
                  wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,6,f)];
                else if((j==((p->n[1]))) )                
                  wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-4,f)];
}

__device__ __host__
void bc_periodic2_temp2(real *wt, struct params *p,int i, int j, int f) {


               if(i<1 && j<1)
                {
                  if(i==j)
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,6,f)];
                  else                  
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];                                    
                }
                else if(i<1 && j>((p->n[1])-1))
                {
                  if(i==(j-(p->n[1])-1))                  
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];                                     
                  else                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-6,f)];                                     
                }
                else if(i>((p->n[0])-1) && j<1)
                {
                  if((i-(p->n[0])+1)==j)                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-5,j,f)];                                    
                  else                  
                   wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,4,f)];                                    
                }
                else if(i>((p->n[0])-1) && j>((p->n[1])-1))
                {
                  if(i==j)                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-5,f)];                                    
                  else                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-5,j,f)];                                    
                }                       
                 
                




}


__global__ void hyperdifvisc4_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
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


int shift=order*NVAR*(p->n[0])*(p->n[1]);


   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l





 /*  for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
          if( i<(ni) && j<(nj))
            {
                  bc_periodic1_hdv1(wtemp,p,i,j,tmp4);
                  bc_periodic1_hdv1(wtemp,p,i,j,tmp5);
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
                  bc_periodic2_hdv1(wtemp,p,i,j,tmp4);
                  bc_periodic2_hdv1(wtemp,p,i,j,tmp5);
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


   if(i>1 && i<((p->n[0])-2) && j>1 && j<((p->n[1])-2))
   {
     //wd[fencode_hdv1(p,i,j,hdnur+hand)]=wtemp2[fencode_hdv1(p,i+1,j+1,tmpnui)];
     if(wtemp[fencode_hdv1(p,i,j,tmp5)]>0)
{
//p->cmax=1.0;
	wd[fencode_hdv1(p,i,j,hdnur+hand)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1]))*(p->cmax)*(p->chyp[field])*wtemp[fencode_hdv1(p,i,j,tmp4)]/wtemp[fencode_hdv1(p,i,j,tmp5)];

          //wd[fencode_hdv1(p,i,j,hdnur+hand)]=wtemp[fencode_hdv1(p,i,j,tmp4)];
	//wd[fencode_hdv1(p,i,j,hdnul+hand)]=0.01;
}
     else
        wd[fencode_hdv1(p,i,j,hdnur+hand)]=0;

    

     //temporary trap for debugging
     //if(wd[fencode_hdv1(p,i,j,hdnul+hand)]>0.02 ||wd[fencode_hdv1(p,i,j,hdnul+hand)] < -0.02)
     //                                                     wd[fencode_hdv1(p,i,j,hdnul+hand)]=0.0;

   }
}
 __syncthreads();



 
}





__global__ void hyperdifvisc3_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
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

  real maxt1=0,max3=0, maxt2=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));


int shift=order*NVAR*(p->n[0])*(p->n[1]);


   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l





  //compute md3r and md1r
//tmp4    md3r
//tmp5    md1r
  //js=0;
 // is=0;
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
  if( i>1 && j>1 && i<((p->n[0])-2) && j<((p->n[1])-2))            
   {
         maxt1=0;
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                {
                   if(wtemp1[fencode_hdv1(p,i+1+is,j+1+js,d3)]>maxt1)
                         maxt1=wtemp1[fencode_hdv1(p,i+1+is,j+1+js,d3)];

                }
          wtemp[fencode_hdv1(p,i,j,tmp4)]=maxt1;

         maxt2=0;
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                {
                   if(wtemp1[fencode_hdv1(p,i+1+is,j+1+js,d1)]>maxt2)
                        maxt2=wtemp1[fencode_hdv1(p,i+1+is,j+1+js,d1)];

                }
          wtemp[fencode_hdv1(p,i,j,tmp5)]=maxt2;
   }
}
   __syncthreads();







 
}




__global__ void hyperdifvisc2_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
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


int shift=order*NVAR*(p->n[0])*(p->n[1]);






   //tmp1  tmp_nuI
 
//compute d3r and d1r
   //tmp2  d3r
    //tmp3 d1r

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
 
   if(i>1 && j>1 && i<((p->n[0])) && j<((p->n[1])))
   //if(i<((p->n[0])-1) && j<((p->n[1])-1))            
   { 
     if(hand==0)
     {
           wtemp1[fencode_hdv1(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv1(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i,j,tmpnui)] ) - (wtemp2[fencode_hdv1(p,i+2*(dim==0),j+2*(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i-(dim==0),j-(dim==1),tmpnui)]    ));
     }
     else
     {
          // wtemp1[fencode_hdv1(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv1(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i,j,tmpnui)]) - (wtemp2[fencode_hdv1(p,i+2*(dim==0),j+2*(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i-(dim==0),j-(dim==1),tmpnui)]    ));
           wtemp1[fencode_hdv1(p,i,j,d3)]=fabs(3.0*(wtemp2[fencode_hdv1(p,i,j,tmpnui)] - wtemp2[fencode_hdv1(p,i-(dim==0),j-(dim==1),tmpnui)]) - (wtemp2[fencode_hdv1(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i-2*(dim==0),j-2*(dim==1),tmpnui)]    ));
     }
   }
}
   __syncthreads();
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;

   if(i>0 && j>0 && i<=((p->n[0])) && j<=((p->n[1])))            
   { 
     if(hand==0)
     {

           wtemp1[fencode_hdv1(p,i,j,d1)]=fabs((wtemp2[fencode_hdv1(p,i+(dim==0),j+(dim==1),tmpnui)] - wtemp2[fencode_hdv1(p,i,j,tmpnui)] ));
     }
     else
     {
           wtemp1[fencode_hdv1(p,i,j,d1)]=fabs((wtemp2[fencode_hdv1(p,i,j,tmpnui)] - wtemp2[fencode_hdv1(p,i-(dim==0),j-(dim==1),tmpnui)] ));
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
                  bc_periodic1_hdv1(wtemp,p,i,j,tmp2);
                  bc_periodic1_hdv1(wtemp,p,i,j,tmp3);
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
                  bc_periodic2_hdv1(wtemp,p,i,j,tmp2);
                  bc_periodic2_hdv1(wtemp,p,i,j,tmp3);
             }
}
                __syncthreads();*/








 
}



__global__ void hyperdifvisc1a_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  int i,j,iv;
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
//int numBlocks = (ni*nj+numThreadsPerBlock-1) / numThreadsPerBlock;
  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*(p->n[0])*(p->n[1]);
  __shared__ real wts[512];
  __shared__ real wms[512];



   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if( i<((p->n[0])) && j<((p->n[1])))
   {
     //tmp6 is tmp_nuI
     wtemp2[fencode_hdv1(p,i+1,j+1,tmpnui)]=wtemp[fencode_hdv1(p,i,j,tmp6)];

   }

   }
   __syncthreads();








 
}


__global__ void hyperdifvisc1_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  int i,j,iv;
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
//int numBlocks = (ni*nj+numThreadsPerBlock-1) / numThreadsPerBlock;
  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
   jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*(p->n[0])*(p->n[1]);
  __shared__ real wts[512];
  __shared__ real wms[512];




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


        for(int f=tmp1; f<=tmp8; f++)
                 wtemp[fencode_hdv1(p,i,j,f)]=0;

        for(int f=d1; f<=d3; f++)
                 wtemp1[fencode_hdv1(p,i,j,f)]=0;
      wtemp2[fencode_hdv1(p,i,j,tmpnui)]=0;
      if(i==((p->n[0])-1))
      {
        for(int f=d1; f<=d3; f++)
                 wtemp1[fencode_hdv1(p,i+1,j,f)]=0;
        wtemp2[fencode_hdv1(p,i+1,j,tmpnui)]=0;
        wtemp2[fencode_hdv1(p,i+2,j,tmpnui)]=0;
      }
      if(j==((p->n[1])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[fencode_hdv1(p,i,j+1,f)]=0;
          wtemp2[fencode_hdv1(p,i,j+1,tmpnui)]=0;
          wtemp2[fencode_hdv1(p,i,j+2,tmpnui)]=0;
      }
      if(j==((p->n[1])-1)  && i==j)
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[fencode_hdv1(p,i+1,j+1,f)]=0;
          for(int di=0; di<2; di++)
             for(int dj=0; dj<2; dj++)
                wtemp2[fencode_hdv1(p,i+1+di,j+1+dj,tmpnui)]=0;

      }



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

        //for(iv=0;iv<NVAR;iv++)
        //               wms[tid+iv*blockdim]=wmod[fencode_hdv1(p,i,j,iv)+shift];
        //wts[tid]=wtemp[fencode_hdv1(p,i,j,tmp6)];
        //temp value for viscosity

       //tmp6  tmpnu
#ifdef USE_SAC
	if((field ==mom1 || field == mom2))
		wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,field)+shift]/(((wmod[fencode_hdv1(p,i,j,rho)+shift] +wmod[fencode_hdv1(p,i,j,rhob)+shift])));
               //wts[tid]=wms[tid+field*blockdim]/(((wms[tid+rho*blockdim] +wms[tid+rhob*blockdim])));
     	else if(field !=energy)
        	wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,field)+shift];///(bfac2+bfac3*((wmod[fencode_hdv1(p,i,j,rho)+shift] +wmod[fencode_hdv1(p,i,j,rhob)+shift])));
               //wts[tid]=wms[tid+field*blockdim];///(bfac2+bfac3*((wmod[fencode_hdv1(p,i,j,rho)+shift] +wmod[fencode_hdv1(p,i,j,rhob)+shift])));
     	else
        wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,energy)+shift]-0.5*(wmod[fencode_hdv1(p,i,j,b1)+shift]*wmod[fencode_hdv1(p,i,j,b1)+shift]+wmod[fencode_hdv1(p,i,j,b2)+shift]*wmod[fencode_hdv1(p,i,j,b2)+shift])+(wmod[fencode_hdv1(p,i,j,mom1)+shift]*wmod[fencode_hdv1(p,i,j,mom1)+shift]+wmod[fencode_hdv1(p,i,j,mom2)+shift]*wmod[fencode_hdv1(p,i,j,mom2)+shift])/(wmod[fencode_hdv1(p,i,j,rho)+shift]+wmod[fencode_hdv1(p,i,j,rhob)+shift] );
//wts[tid]=wms[tid+energy*blockdim]-0.5*(wms[tid+b1*blockdim]*wms[tid+b1*blockdim]+wms[tid+b2*blockdim]*wms[tid+b2*blockdim])+(wms[tid+mom1*blockdim]*wms[tid+mom1*blockdim]+wms[tid+mom2*blockdim]*wms[tid+mom2*blockdim])/(wms[tid+rho*blockdim]+wms[tid+rhob*blockdim] );

#else
	if((field ==mom1 || field == mom2))
		wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,field)+shift]/(((wmod[fencode_hdv1(p,i,j,rho)+shift] )));
     else if(field !=energy)
        wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,field)+shift]/(bfac2+bfac3*(wmod[fencode_hdv1(p,i,j,rho)+shift] ));

     else
        wtemp[fencode_hdv1(p,i,j,tmp6)]=wmod[fencode_hdv1(p,i,j,energy)+shift]-0.5*(wmod[fencode_hdv1(p,i,j,b1)+shift]*wmod[fencode_hdv1(p,i,j,b1)+shift]+wmod[fencode_hdv1(p,i,j,b2)+shift]*wmod[fencode_hdv1(p,i,j,b2)+shift])+(wmod[fencode_hdv1(p,i,j,mom1)+shift]*wmod[fencode_hdv1(p,i,j,mom1)+shift]+wmod[fencode_hdv1(p,i,j,mom2)+shift]*wmod[fencode_hdv1(p,i,j,mom2)+shift])/(wmod[fencode_hdv1(p,i,j,rho)+shift] );

#endif

        //for(iv=0;iv<NVAR;iv++)
        //               wmod[fencode_hdv1(p,i,j,iv)+shift]=wms[tid+iv*blockdim];
        //          wtemp[fencode_hdv1(p,i,j,tmp6)]=wts[tid];

       // wd[fencode_hdv1(p,i,j,hdnur)]=0;
        wd[fencode_hdv1(p,i,j,hdnul+hand)]=0;
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
	
        bc_hyperdif(wtemp2, p,i, j, tmpnui,dim);

   }


    }
   __syncthreads();

/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   {

     i=ip*(p->npgp[0])+ipg;
     j=jp*(p->npgp[1])+jpg;
    //set viscosities
   if( i<((p->n[0])) && j<((p->n[1])))
   {
     //tmp6 is tmp_nuI
     wtemp2[fencode_hdv1(p,i+1,j+1,tmpnui)]=wtemp[fencode_hdv1(p,i,j,tmp6)];

   }

   }
   __syncthreads();
*/







 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv1(char *label)
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





int cuhyperdifvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand)
{



 dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifvisc1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();

}







