

__device__ __host__
void computej_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;

 // real dbzdy, dbydz;
 // real dbzdx, dbxdz;
 // real dbydx, dbxdy;

 // dbzdy=grad_MODID(wmod,p,i,j,b3,1);
 // dbydz=0.0;
 // dbzdx=grad_MODID(wmod,p,i,j,b3,0);
//  dbxdz=0.0;
 // dbydx=grad_MODID(wmod,p,i,j,b2,0);
 // dbxdy=grad_MODID(wmod,p,i,j,b1,1);

  wd[fencode_MODID(p,i,j,0)]=(grad_MODID(wmod,p,i,j,b3,1))/(p->mu);
  wd[fencode_MODID(p,i,j,1)]=(grad_MODID(wmod,p,i,j,b3,0))/(p->mu);
  wd[fencode_MODID(p,i,j,2)]=(grad_MODID(wmod,p,i,j,b2,0)-grad_MODID(wmod,p,i,j,b1,1))/(p->mu);
  
          #ifdef USE_SAC
	  wd[fencode_MODID(p,i,j,0)]+=(grad_MODID(wmod,p,i,j,b3b,1))/(p->mu);
	  wd[fencode_MODID(p,i,j,1)]+=(grad_MODID(wmod,p,i,j,b3b,0))/(p->mu);
	  wd[fencode_MODID(p,i,j,2)]+=(grad_MODID(wmod,p,i,j,b2b,0)-grad_MODID(wmod,p,i,j,b1b,1))/(p->mu);


         #endif

 
  //return ( status);
}

__device__ __host__
void computebdotv_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;
 //real bsq=wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)];
//  wd[fencode_MODID(p,i,j,4)]=  wd[fencode_MODID(p,i,j,3)]+0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]);
        #ifdef USE_SAC

wd[fencode_MODID(p,i,j,bdotv)]=((wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b1b)])*wmod[fencode_MODID(p,i,j,mom1)]+(wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b2b)])*wmod[fencode_MODID(p,i,j,mom2)]+(wmod[fencode_MODID(p,i,j,b3)]+wmod[fencode_MODID(p,i,j,b3b)])*wmod[fencode_MODID(p,i,j,mom3)])/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]);
         #else
wd[fencode_MODID(p,i,j,bdotv)]=(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,mom1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,mom2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,mom3)])/wmod[fencode_MODID(p,i,j,rho)];
         #endif
 // return ( status);
}

__device__ __host__
void computedivb_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;
 //real bsq=wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)];
//  wd[fencode_MODID(p,i,j,4)]=  wd[fencode_MODID(p,i,j,3)]+0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]);

wd[fencode_MODID(p,i,j,divb)]=grad_MODID(wmod,p,i,j,b1,0)+grad_MODID(wmod,p,i,j,b2,1);
        #ifdef USE_SAC
		wd[fencode_MODID(p,i,j,divb)]+=grad_MODID(wmod,p,i,j,b1b,0)+grad_MODID(wmod,p,i,j,b2b,1);
         #endif
 // return ( status);
}


__device__ __host__
void computept_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
 // int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
 wd[fencode_MODID(p,i,j,pressuret)]=(p->adiab)*pow(wmod[fencode_MODID(p,i,j,rho)],p->gamma);
#elif defined(USE_SAC)
  wd[fencode_MODID(p,i,j,pressuret)]=  wd[fencode_MODID(p,i,j,pressurek)]+0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)])+(wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3)]);

  wd[fencode_MODID(p,i,j,ptb)]=  wd[fencode_MODID(p,i,j,pkb)]+0.5*(wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1b)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2b)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3b)]);


#else

 //real bsq=wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)];
  wd[fencode_MODID(p,i,j,pressuret)]=  wd[fencode_MODID(p,i,j,pressurek)]+0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]);

#endif



  if(wd[fencode_MODID(p,i,j,pressuret)]<0)
              wd[fencode_MODID(p,i,j,pressuret)]=0.001;


 // return ( status);
}
__device__ __host__
void computepk_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
  //int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
wd[fencode_MODID(p,i,j,pressurek)]=(p->adiab)*pow(wmod[fencode_MODID(p,i,j,rho)],p->gamma);
wd[fencode_MODID(p,i,j,vel1)]=wmod[fencode_MODID(p,i,j,mom1)]/(wmod[fencode_MODID(p,i,j,rho)]);
wd[fencode_MODID(p,i,j,vel2)]=wmod[fencode_MODID(p,i,j,mom2)]/(wmod[fencode_MODID(p,i,j,rho)]);
wd[fencode_MODID(p,i,j,vel3)]=wmod[fencode_MODID(p,i,j,mom3)]/(wmod[fencode_MODID(p,i,j,rho)]);
#elif defined(USE_SAC)

wd[fencode_MODID(p,i,j,vel1)]=wmod[fencode_MODID(p,i,j,mom1)]/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]);
wd[fencode_MODID(p,i,j,vel2)]=wmod[fencode_MODID(p,i,j,mom2)]/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]);
wd[fencode_MODID(p,i,j,vel3)]=wmod[fencode_MODID(p,i,j,mom3)]/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]);

 wd[fencode_MODID(p,i,j,pressurek)]=((p->gamma)-1)*(wmod[fencode_MODID(p,i,j,energy)]- 0.5*(wmod[fencode_MODID(p,i,j,mom1)]*wmod[fencode_MODID(p,i,j,mom1)]+wmod[fencode_MODID(p,i,j,mom2)]*wmod[fencode_MODID(p,i,j,mom2)]+wmod[fencode_MODID(p,i,j,mom3)]*wmod[fencode_MODID(p,i,j,mom3)])-0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]) -(wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3)]) );


wd[fencode_MODID(p,i,j,pkb)]=((p->gamma)-1)*(wmod[fencode_MODID(p,i,j,energyb)]- 0.5*(wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1b)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2b)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3b)]) );

#else
wd[fencode_MODID(p,i,j,vel1)]=wmod[fencode_MODID(p,i,j,mom1)]/(wmod[fencode_MODID(p,i,j,rho)]);
wd[fencode_MODID(p,i,j,vel2)]=wmod[fencode_MODID(p,i,j,mom2)]/(wmod[fencode_MODID(p,i,j,rho)]);
wd[fencode_MODID(p,i,j,vel3)]=wmod[fencode_MODID(p,i,j,mom3)]/(wmod[fencode_MODID(p,i,j,rho)]);
  //real momsq=wmod[fencode_MODID(p,i,j,mom1)]*wmod[fencode_MODID(p,i,j,mom1)]+wmod[fencode_MODID(p,i,j,mom2)]*wmod[fencode_MODID(p,i,j,mom2)]+wmod[fencode_MODID(p,i,j,mom3)]*wmod[fencode_MODID(p,i,j,mom3)];
  //real bsq=wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)];
  wd[fencode_MODID(p,i,j,pressurek)]=((p->gamma)-1)*(wmod[fencode_MODID(p,i,j,energy)]- 0.5*(wmod[fencode_MODID(p,i,j,mom1)]*wmod[fencode_MODID(p,i,j,mom1)]+wmod[fencode_MODID(p,i,j,mom2)]*wmod[fencode_MODID(p,i,j,mom2)]+wmod[fencode_MODID(p,i,j,mom3)]*wmod[fencode_MODID(p,i,j,mom3)])/wmod[fencode_MODID(p,i,j,rho)]-0.5*(wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]) );


#endif






  if(wd[fencode_MODID(p,i,j,pressurek)]<0)
              wd[fencode_MODID(p,i,j,pressurek)]=0.001;
  //return ( status);
}

__device__ __host__
void computec_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{

  
#ifdef ADIABHYDRO
/*below used for adiabatic hydrodynamics*/
  wd[fencode_MODID(p,i,j,soundspeed)]=sqrt((p->adiab)/wmod[fencode_MODID(p,i,j,rho)]);
#elif defined(USE_SAC)
wd[fencode_MODID(p,i,j,soundspeed)]=sqrt((  (p->gamma))*(wd[fencode_MODID(p,i,j,pressuret)]+wd[fencode_MODID(p,i,j,ptb)])/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]   ));
wd[fencode_MODID(p,i,j,cfast)]=sqrt((   ( (wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)]) + (wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1b)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2b)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3b)]) +2.0*(wmod[fencode_MODID(p,i,j,b1b)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2b)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3b)]*wmod[fencode_MODID(p,i,j,b3)])    )/(wmod[fencode_MODID(p,i,j,rho)]+wmod[fencode_MODID(p,i,j,rhob)]))+(wd[fencode_MODID(p,i,j,soundspeed)]*wd[fencode_MODID(p,i,j,soundspeed)]));
#else
wd[fencode_MODID(p,i,j,soundspeed)]=sqrt(((p->gamma))*wd[fencode_MODID(p,i,j,pressuret)]/wmod[fencode_MODID(p,i,j,rho)]);
wd[fencode_MODID(p,i,j,cfast)]=sqrt(((wmod[fencode_MODID(p,i,j,b1)]*wmod[fencode_MODID(p,i,j,b1)]+wmod[fencode_MODID(p,i,j,b2)]*wmod[fencode_MODID(p,i,j,b2)]+wmod[fencode_MODID(p,i,j,b3)]*wmod[fencode_MODID(p,i,j,b3)])/wmod[fencode_MODID(p,i,j,rho)])+(wd[fencode_MODID(p,i,j,soundspeed)]*wd[fencode_MODID(p,i,j,soundspeed)]));
#endif



  
}

__device__ __host__
void computecmax_MODID(real *wmod,real *wd,struct params *p,int i,int j)
{
#ifdef ADIABHYDRO
       if(wd[fencode_MODID(p,i,j,soundspeed)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode_MODID(p,i,j,soundspeed)]));
                    p->cmax=(wd[fencode_MODID(p,i,j,soundspeed)]);
#else
       if(wd[fencode_MODID(p,i,j,soundspeed)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode_MODID(p,i,j,soundspeed)]));
                    p->cmax=(wd[fencode_MODID(p,i,j,soundspeed)]);
       if(wd[fencode_MODID(p,i,j,cfast)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode_MODID(p,i,j,soundspeed)]));
                    p->cmax=(wd[fencode_MODID(p,i,j,cfast)]);
#endif

}


