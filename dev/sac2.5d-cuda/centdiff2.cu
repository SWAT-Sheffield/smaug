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
#include "gradops_cd2.cuh"
#include "dervfields_cd2.cuh"

__device__ __host__
real fluxe2(real *dw, real *wd, real *w, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real flux=0;


        #if defined USE_SAC  || defined USE_SAC_3D
computept3_cd2(w,wd,p,ii);



      		flux= -wd[fencode3_cd2(p,ii,ptb)]*grad3d_cd2(wd,p,ii,vel1+dir,dir);
               flux += +w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1,0)+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1+1,1);
         #endif


        #if defined USE_SAC_3D
               flux += +w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel3,0);
        #endif

  return flux;


  //return ( ddc1-ddc2);
}



__device__ __host__
int divflux_cd2(real *dw, real *wd, real *w, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;
  dw[fencode3_cd2(p,ii,field)]= grad3d_cd2(wd,p,ii,flux,dir);//+grad_cd2(wd,p,ii,f2,1); 


 #ifdef USE_SAC

  //commented out to test against vac
 /* if(field==energy)
  {    
     dw[fencode3_cd2(p,ii,field)]+=fluxe2(dw, wd, w, p,ix, iy,dir)+w[fencode3_cd2(p,ii,rho)]*((p->g[dir])*w[fencode3_cd2(p,ii,mom1+dir)]    )/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}


__device__ __host__
int addenergyterms_cd2(real *dw, real *wd, real *w, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;
  //dw[fencode3_cd2(p,ii,field)]= grad_cd2(wd,p,ii,flux,dir);//+grad_cd2(wd,p,ii,f2,1); 


 #if defined USE_SAC  ||  defined USE_SAC_3D

  
  if(field==energy)
  {    
     computept3_cd2(w,wd,p,ii);
     dw[fencode3_cd2(p,ii,field)]=fluxe2(dw, wd, w, p,ii,dir)+w[fencode3_cd2(p,ii,rho)]*((p->g[dir])*w[fencode3_cd2(p,ii,mom1+dir)]    )/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
   }


 #endif
  return ( status);
}



__device__ __host__
real transportflux_cd2 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
  real ddcx=0,ddcy=0;

   real flux=0;

   //transport flux
   //use versions with velocity less ops may improve performance
    switch(direction)
  {
     case 0:
        #if defined USE_SAC  || defined USE_SAC_3D
     flux= w[fencode3_cd2(p,ii,mom1)]*w[fencode3_cd2(p,ii,field)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
     //flux= w[fencode3_cd2(p,ii,mom1)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)];

        #else
     flux= w[fencode3_cd2(p,ii,mom1)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)];

        #endif
     break;
     case 1:
        #if defined USE_SAC  || defined USE_SAC_3D
      flux= w[fencode3_cd2(p,ii,mom2)]*w[fencode3_cd2(p,ii,field)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
    // flux= w[fencode3_cd2(p,ii,mom2)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)];

        #else
     //flux= wd[fencode3_cd2(p,ii,vel2)]*w[fencode3_cd2(p,ii,field)];
     flux= w[fencode3_cd2(p,ii,mom2)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)];

        #endif
     break;
        #if  defined USE_SAC_3D
     case 2:
     flux= w[fencode3_cd2(p,ii,mom3)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)];
     //flux= wd[fencode3_cd2(p,ii,vel3)]*w[fencode3_cd2(p,ii,field)];
     break;
     #endif
   }
  return flux;


  //return ( ddc1-ddc2);
}




__device__ __host__
real fluxb1(real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
  real ddcx=0,ddcy=0;

   real flux=0;

    switch(field)
    {
      case b1:
      //if(direction !=0)
        #if defined USE_SAC  || defined USE_SAC_3D

  flux= -(w[fencode3_cd2(p,ii,field+direction)]+w[fencode3_cd2(p,ii,field+(NDIM+2)+direction)])*w[fencode3_cd2(p,ii,mom1)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

flux+= (w[fencode3_cd2(p,ii,field+(NDIM+2))])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
         #endif

       break;

      case b2:
      //if(direction !=1)
        #if defined USE_SAC  || defined USE_SAC_3D
		flux= -(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,b1+(NDIM+2)+direction)])*w[fencode3_cd2(p,ii,mom2)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

               flux+= (w[fencode3_cd2(p,ii,field+(NDIM+2))])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
         #endif 
       break;
        #if defined USE_SAC_3D
      case b3:
      //if(direction !=2)

		flux= -(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,b1+(NDIM+2)+direction)])*w[fencode3_cd2(p,ii,mom2)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

               flux+= (w[fencode3_cd2(p,ii,field+(NDIM+2))])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

       break;
         #endif
     }


  return flux;
}



__device__ __host__
real fluxe1(real *dw, real *wd, real *w, struct params *p,int *ii, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
  real ddcx=0,ddcy=0;

   real flux=0;

computept3_cd2(w,wd,p,ii);

        #if defined USE_SAC || defined USE_SAC_3D

flux= -w[fencode3_cd2(p,ii,b1+direction)]*wd[fencode3_cd2(p,ii,bdotv)]+(w[fencode3_cd2(p,ii,mom1+direction)]*wd[fencode3_cd2(p,ii,energyb)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]));

flux += w[fencode3_cd2(p,ii,mom1+direction)]*(wd[fencode3_cd2(p,ii,pressuret)]+wd[fencode3_cd2(p,ii,ptb)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

flux -= w[fencode3_cd2(p,ii,b1b+direction)]*(w[fencode3_cd2(p,ii,b1)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2)]*w[fencode3_cd2(p,ii,mom2)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)])
            - w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,mom2)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

         #endif

#ifdef USE_SAC_3D

flux -= w[fencode3_cd2(p,ii,b1b+direction)]*(w[fencode3_cd2(p,ii,b3)]*w[fencode3_cd2(p,ii,mom3)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)])
            - w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,mom3)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

#endif

  return flux;

}








__device__ __host__
int computefluxe(real *dw, real *wd, real *w, struct params *p,int *ii,int direction) {

  int field;//, direction;
  int status=0;

         #if defined USE_SAC  || defined USE_SAC_3D
	     wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,energy,direction)+fluxe1(dw,wd,w,p,ii,direction);
         #endif

        
  return ( status);
}

__device__ __host__
int computefluxb (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {


  int status=0;


     switch(field)
     {
       case b1 :
         #if defined USE_SAC  || defined USE_SAC_3D
      if(direction==0)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
      else
wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,field,direction)+fluxb1(dw,wd,w,p,ii,field,direction);
         #endif

       break;

       case b2 :
         #if defined USE_SAC  || defined USE_SAC_3D
      if(direction==1)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
else
wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,field,direction)+fluxb1(dw,wd,w,p,ii,field,direction);
         #endif

       break;

 #ifdef USE_SAC_3D
       case b3 :

      if(direction==2)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
else
wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,field,direction)+fluxb1(dw,wd,w,p,ii,field,direction);


       break;
  #endif

    }
   
    
  return ( status);
}






//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux_cd2 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int dir) {

  //int status=0;
  switch(field)
  {
     case energy:
      computefluxe(dw,wd,w,p,ii,dir);
      
      // add the following terms for SAC
      // del((b bb+ bb b).v)+ptb del v - bb bb del v
     break;
     case b1:
      computefluxb(dw,wd,w,p,ii,field,dir);
     break;
     case b2:
       computefluxb(dw,wd,w,p,ii,field,dir);
     break;
#ifdef USE_SAC_3D
     case b3:
      computefluxb(dw,wd,w,p,ii,field,dir);
     break;
#endif
  }
  //return ( status);
}



__global__ void centdiff2a_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 // int index;
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
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

			//if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                          #ifdef USE_SAC_3D
       				if(ii[0]<((p->n[0])-2) && ii[0]>1 && ii[1]>1 && ii[1]<((p->n[1])-2) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]<((p->n[0]))-2 && ii[0]>1  && ii[1]>1 && ii[1]<((p->n[1])-1))
     			  #endif
                                divflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 


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

                        switch(dir)
                        {
                         case 0:

                         //if(i<(ni)  && j >1 &&  j<(nj-2))
                         #ifdef USE_SAC_3D
       				if(ii[0]<((p->n[0]))  && ii[1]>1 && ii[1]<((p->n[1])-2) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]<((p->n[0]))   && ii[1]>1 && ii[1]<((p->n[1])-2))
     			  #endif
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)]; 
                         break;
                         case 1:
                         #ifdef USE_SAC_3D
       				if(ii[0]>1 && ii[0]<((p->n[0])-2)  &&  ii[1]<((p->n[1])) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]>1 && ii[0]<((p->n[0])-2)   && ii[1]<((p->n[1])) )
     			  #endif
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)];
                         break;
                         #ifdef USE_SAC_3D
                         case 2:

                         //if(i>1 &&  i<(ni-2) && j<(nj))
      			if(ii[0]>1 && ii[0]<((p->n[0])-2)  && ii[1]>1 && ii[1]<((p->n[1])-2)  && ii[2]<((p->n[2])))
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)];
                         break;
                         #endif
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


     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                            dwn1[fencode3_cd2(p,ii,f)]=0.0;


			//if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
     #ifdef USE_SAC_3D
       if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                                addenergyterms_cd2(dwn1,wd,wmod+ordero*NVAR*dimp,p,ii,f,dir); 

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

                        switch(dir)
                        {
                         case 0:

                         //if(i<(ni)  && j >1 &&  j<(nj-2))
			     #ifdef USE_SAC
				   if(ii[0]<ni && ii[1] >1 && ii[1]<(nj-2))
			     #endif
			     #ifdef USE_SAC_3D
				  if(ii[0]<ni && ii[1] >1 && ii[2] >1  && ii[1]<(nj-2) && ii[2]<(nk-2))
			     #endif   
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)]; 
                         break;
                         case 1:

                         //if(i>1 &&  i<(ni-2) && j<(nj))
			     #ifdef USE_SAC
				   if(ii[1]<nj && ii[0] >1 && ii[0]<(ni-2))
			     #endif
			     #ifdef USE_SAC_3D
				  if(ii[1]<nj && ii[0] >1 && ii[2] >1  && ii[0]<(ni-2) && ii[2]<(nk-2))
			     #endif
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)];
                         break;
                   #ifdef USE_SAC_3D
                         case 2:
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                           if(ii[2]<nk && ii[0] >1 && ii[1] >1  && ii[1]<(nj-2) && ii[0]<(ni-2))
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)];
                         break;
                   #endif
                        }


}
__syncthreads(); 

                         
}





__global__ void centdiff2_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 // int index;
  int ni=p->n[0];
  int nj=p->n[1];

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
                            dwn1[fencode3_cd2(p,ii,f)]=0.0;

                               wd[fencode3_cd2(p,ii,flux)]=0.0;

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
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-1) && ii[2]>1 && ii[2]<(p->n[2]-1))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-1))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-1) && ii[2]>1 && ii[2]<(p->n[2]-1))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-1))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-1) && ii[1]>1 && ii[1]<(p->n[1]-1))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 
                         break;
                         #endif
                        }
              //  }
                        //might need to set boundaries correctly
 
}
__syncthreads();                        






                         
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cd2(char *label)
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




int cucentdiff2(struct params **p, struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt, int field,int dir)
{

    dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;
   //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
   // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);


     centdiff2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();

     centdiff2a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();

     // cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

     //checkErrors("copy data from device");

}


