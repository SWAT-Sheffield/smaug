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
#include "usersource_cd1.cuh"

__device__ __host__
int divflux1(real *dw, real *wd, real *w, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;

dw[fencode3_cd1(p,ii,field)]= grad3d_cd1(wd,p,ii,flux,dir); 
//dw[fencode3_cd1(p,ii,field)]=0.0;
 switch(field)
  {

     case mom1:
       dw[fencode3_cd1(p,ii,field)]+= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
    case mom2:
      dw[fencode3_cd1(p,ii,field)]+= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
#if defined USE_SAC_3D
    case mom3:
      dw[fencode3_cd1(p,ii,field)]+= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
#endif
    case rho:
     ;// dw[fencode3_cd1(p,ii,field)]+= ix/800;
      break;

  }    
 // dw[fencode3_cd1(p,ii,field)]= gradd0_cd1(wd,p,ii,f1,0)+gradd1_cd1(wd,p,ii,f2,1);    
  return ( status);
}






__device__ __host__
real transportflux (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;

   //transport flux
   //this will work without the switch as follows
        #if defined USE_SAC || defined USE_SAC_3D
     //flux= w[fencode3_cd1(p,ii,mom1+direction)]*w[fencode3_cd1(p,ii,field)]/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
flux= wd[fencode3_cd1(p,ii,vel1+direction)]*w[fencode3_cd1(p,ii,field)];
        #else
     //flux= w[fencode3_cd1(p,ii,mom1+direction)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];
flux= wd[fencode3_cd1(p,ii,vel1+direction)]*w[fencode3_cd1(p,ii,field)];
        #endif


  /*  switch(direction)
  {
     case 0:
     //flux= wd[fencode3_cd1(p,ii,vel1)]*w[fencode3_cd1(p,ii,field)];
        #if defined USE_SAC || defined USE_SAC_3D
     flux= w[fencode3_cd1(p,ii,mom1)]*w[fencode3_cd1(p,ii,field)]/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
    // flux= w[fencode3_cd1(p,ii,mom1)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #else
     flux= w[fencode3_cd1(p,ii,mom1)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #endif
     break;
     case 1:
        #if defined USE_SAC || defined USE_SAC_3D
     flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
     //flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #else
     //flux= wd[fencode3_cd1(p,ii,vel2)]*w[fencode3_cd1(p,ii,field)];
     flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #endif
     break;
     case 2:
        #ifdef USE_SAC
     flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
     //flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #else
     //flux= wd[fencode3_cd1(p,ii,vel2)]*w[fencode3_cd1(p,ii,field)];
     flux= w[fencode3_cd1(p,ii,mom2)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)];

        #endif
     break;
    /* case 2:
     flux= wd[fencode3_cd1(p,ii,vel3)]*w[fencode3_cd1(p,ii,field)];
     break;*/
   //}*/
  return flux;


  //return ( ddc1-ddc2);
}






__device__ __host__
real fluxmom1 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;


         #if defined USE_SAC || defined USE_SAC_3D
     		flux= -(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]-w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)];
        #endif




   //transport flux
  /*  switch(direction)
  {
     case 0:
         #if defined USE_SAC || defined USE_SAC_3D
     		flux= -w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b)]-w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1)]-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1)];
        #endif
     break;
     case 1:
         #if defined USE_SAC || defined USE_SAC_3D
                flux= -w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b2b)]-w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b2)]-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b2)];
         #endif
     break;
#ifdef USE_SAC_3D
     case 2:
         
                flux= -w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b3b)]-w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b3)]-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b3)];

     break;
#endif

   }*/




  return flux;


  //return ( ddc1-ddc2);
}










__device__ __host__
int computefluxrho (real *dw, real *wd, real *w, struct params *p,int *ii,int direction) {

  int field;
  int status=0;

  // for(direction=0;direction<2;direction++)
         #if defined USE_SAC || defined USE_SAC_3D
	      wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,rho,direction)+(w[fencode3_cd1(p,ii,rhob)]*w[fencode3_cd1(p,ii,mom1+direction)])/(w[fencode3_cd1(p,ii,rhob)]+w[fencode3_cd1(p,ii,rho)]);
         #else
             wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,rho,direction);
         #endif
  
  return ( status);
}


__device__ __host__
int computefluxmom3 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

#ifdef USE_SAC_3D


    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom1(dw,wd,w,p,ii,field,direction);
               //if(direction==1)
               //   wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,ptb)];
 
 
               if(direction==2)
               {
                //computept_cd1(w,wd,p,ii);
                //commented out to compare with vac
                 wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,pressuret)];


                //  wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,ptb)];


               }
 

#endif

  return ( status);
}



__device__ __host__
int computefluxmom2 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

        #ifdef ADIABHYDRO
     		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction);
        #endif
        #ifdef USE_SAC
    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom1(dw,wd,w,p,ii,field,direction);
               //if(direction==1)
               //   wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,ptb)];
 
        #endif
        #ifdef USE_SAC_3D
    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom1(dw,wd,w,p,ii,field,direction);
               //if(direction==1)
               //   wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,ptb)];
 
        #endif
               if(direction==1)
               {
                //computept_cd1(w,wd,p,ii);
                //commented out to compare with vac
                 wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,pressuret)];

        #ifdef USE_SAC

               //   wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,ptb)];
        #endif
        #ifdef USE_SAC_3D

               //   wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,ptb)];
        #endif

               }


  return ( status);
}



__device__ __host__
int computefluxmom1 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

        #ifdef ADIABHYDRO
     		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction);
        #endif
        #ifdef USE_SAC
    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom1(dw,wd,w,p,ii,field,direction);
               //if(direction==0)
               //   wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,ptb)];
 
        #endif
        #ifdef USE_SAC_3D
    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom1(dw,wd,w,p,ii,field,direction);
               //if(direction==0)
               //   wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,ptb)];
 
        #endif
               if(direction==0)
               {

                // computept_cd1(w,wd,p,ii);
                 //commented out to compare with vac 
                 wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,pressuret)];

        #ifdef USE_SAC

                //  wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,ptb)];
       #endif
        #ifdef USE_SAC_3D

                //  wd[fencode3_cd1(p,ii,flux)]+=wd[fencode3_cd1(p,ii,ptb)];
       #endif
               }


        
  return ( status);
}







//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int dir) {

  //int status=0;
  switch(field)
  {
     case rho:
      //computevel3_cd1(w,wd,p,ii);
      //computept3_cd1(w,wd,p,ii);
      computefluxrho(dw,wd,w,p,ii,dir);
     break;
     case mom1:
      computefluxmom1(dw,wd,w,p,ii,field,dir);
      //wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,pressuret)];
     break;
     case mom2:
       computefluxmom2(dw,wd,w,p,ii,field,dir);
       //wd[fencode3_cd1(p,ii,f2)]+=wd[fencode3_cd1(p,ii,pressuret)];
     break;
     #ifdef USE_SAC_3D
       case mom3:
        computefluxmom3(dw,wd,w,p,ii,field,dir);
        //wd[fencode3_cd1(p,ii,f3)]+=wd[fencode3_cd1(p,ii,pressuret)];
       break;
     #endif
  }
  //return ( status);
}



__global__ void centdiff1_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dy=p->dx[1];
  real dx=p->dx[0];

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

   fid=0;
   
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
                            dwn1[fencode3_cd1(p,ii,f)]=0.0;
                  	    //for(fid=0;fid<2;fid++)
                               wd[fencode3_cd1(p,ii,flux)]=0.0;
                               //wmod[fencode_cd1(p,i,j,flux)+order*NVAR*(p->n[0])*(p->n[1])]=0.0;
                        }

   }
 __syncthreads();                       



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


//1. 11/1/11 could swap cases below
                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                         #endif
                        }
              //  }
                        //might need to set boundaries correctly
 
}
__syncthreads();                        



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

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
    int nk=p->n[2];
    real dz=p->dx[2];
#endif
 #ifdef USE_SAC_3D
   int kp,kpg;
   
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

   fid=0;



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

			// if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
			     #ifdef USE_SAC
				   if(ii[0]>1 && ii[1] >1 && ii[0]<(ni-2) && ii[1]<(nj-2))
			     #endif
			     #ifdef USE_SAC_3D
				  if(ii[0]>1 && ii[1] >1 && ii[2] >1 && ii[0]<(ni-2) && ii[1]<(nj-2) && ii[2]<(nk-2))
			     #endif                        
                               divflux1(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir);  

}
 __syncthreads();

#if(defined(USE_USERSOURCE))
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
#endif
   #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
#if(defined(USE_USERSOURCE))
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
#endif
     #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif


     #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
       if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
     #endif
     #if(defined(USE_SAC) && defined(USE_USERSOURCE))
       if(ii[0]<(p->n[0])-2 && ii[1]<(p->n[1])-2)
     #endif

                     #ifdef USE_USERSOURCE
                                addsourceterms1_cd2(dwn1,wd,wmod+ordero*NVAR*dimp,p,ii,f,dir); 


                      }
                    __syncthreads();
                     #endif



               // }
    

         // if( i<(ni) && j<(nj))
          //        bc_cont_cd1(dwn1,p,i,j,f);
            //    __syncthreads();



             // for(int f=rho; f<=mom3; f++)
              // {

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


                        switch(dir)
                        {
                         case 0:

                         //if(i<(ni)  && j >1 &&  j<(nj-2))
			     #ifdef USE_SAC
				   if(ii[1]>1 && ii[1] <(nj-2) && ii[0]<(ni) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[1]>1 && ii[1] <(nj-2) && ii[0]<(ni) &&  ii[2]>1 && ii[2] <(nk-2) )
			     #endif                          
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)]; 
                         break;
                         case 1:
			     #ifdef USE_SAC
				   if(ii[0]>1 && ii[1] <(nj) && ii[0]<(ni-2) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[0]>1 && ii[1] <(nj) && ii[0]<(ni-2) &&  ii[2]>1 && ii[2] <(nk-2) )
			     #endif 
                         
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)];
                         break;
                         case 2:

 			     #ifdef USE_SAC
				   if(ii[1]>1 && ii[0] <(ni) && ii[1]<(nj-2) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[0]>1 &&  ii[0]<(ni-2)  && ii[1]>1 &&  ii[1]<(nj-2) && ii[2] <(nk) )
			     #endif                         
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)];
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


