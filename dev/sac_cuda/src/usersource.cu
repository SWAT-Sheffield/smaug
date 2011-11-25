

__device__ __host__
int addsourceterms2_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;

   real xc1,xc2,xc3;
   real xxmax,yymax;
   real dx,dy,dz;
   real aa;
   real s_period;
   real tdep;

   real vx,vy;

   real exp_x,exp_y,exp_z,exp_xyz;


   real x,y,z;
   int i,j,k;
 	  
	  i=ii[0];
	  j=ii[1];
	  k=ii[2];

     xc1=0.1e6;
    xc2=1.0e6;
    xc3=1.0e6;

          x=(p->xmin[1])+(j*(p->dx[1]))-xc2;
          z=(p->xmin[0])+(i*(p->dx[0]))-xc1;
          y=(p->xmin[2])+(k*(p->dx[2]))-xc3;
     // xx=x(ix_1,ix_2,ix_3,2)-xc2
     // yy=x(ix_1,ix_2,ix_3,3)-xc3
     // zz=x(ix_1,ix_2,ix_3,1)-xc1  
  


    xxmax=2.0e6;
    yymax=2.0e6;

    dx=0.1e6;
    dy=0.1e6;
    dz=0.05e6;

    aa=10000.0;
    s_period=30.0;
    tdep=1.00;


        //exp_z=exp(-zz**2.d0/(delta_z**2.d0))
        //exp_x=exp(-xx**2.d0/(delta_x**2.d0))
        //exp_y=exp(-yy**2.d0/(delta_y**2.d0))       
        //exp_xyz=exp_x*exp_y*exp_z
        exp_z=exp(-z*z/(dz*dz));
        exp_x=exp(-x*x/(dx*dx));
        exp_y=exp(-y*y/(dy*dy));       
        exp_xyz=exp_x*exp_y*exp_z;

        //vvx(ix_1,ix_2,ix_3)=AA*yy/yymax*exp_xyz*tdep    
        //vvy(ix_1,ix_2,ix_3)=-AA*xx/xxmax*exp_xyz*tdep 
        vx=aa*y/yymax*exp_xyz*tdep;    
        vy=-aa*x/xxmax*exp_xyz*tdep; 
 
 switch(field)
  {

    case mom2:
                           dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-vx*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
    break;
    case mom3:     
                           dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-vy*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
    break;
    case energy:
                          dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-(vx*vx+vy*vy)*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;
    break;
   } 

  return ( status);
}

__device__ __host__
int addsourceterms1_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;

   real xc1,xc2,xc3;
   real xxmax,yymax;
   real delx,dely,delz;
   real aa;
   real s_period;
   real tdep;

   real vx,vy;

   real exp_x,exp_y,exp_z,exp_xyz;


   real x,y,z;
   int i,j,k;
 	  
	  i=ii[0];
	  j=ii[1];
	  k=ii[2];

     xc1=0.1e6;
    xc2=1.0e6;
    xc3=1.0e6;

          x=(p->xmin[1])+(j*(p->dx[1]))-xc2;
          z=(p->xmin[0])+(i*(p->dx[0]))-xc1;
          y=(p->xmin[2])+(k*(p->dx[2]))-xc3;
     // xx=x(ix_1,ix_2,ix_3,2)-xc2
     // yy=x(ix_1,ix_2,ix_3,3)-xc3
     // zz=x(ix_1,ix_2,ix_3,1)-xc1  
  


    xxmax=2.0e6;
    yymax=2.0e6;

    delx=0.1e6;
    dely=0.1e6;
    delz=0.05e6;

    aa=10000.0;
    s_period=30.0;
    tdep=1.00;


        //exp_z=exp(-zz**2.d0/(delta_z**2.d0))
        //exp_x=exp(-xx**2.d0/(delta_x**2.d0))
        //exp_y=exp(-yy**2.d0/(delta_y**2.d0))       
        //exp_xyz=exp_x*exp_y*exp_z
        exp_z=exp(-z*z/(delz*delz));
        exp_x=exp(-x*x/(delx*delx));
        exp_y=exp(-y*y/(dely*dely));       
        exp_xyz=exp_x*exp_y*exp_z;

        //vvx(ix_1,ix_2,ix_3)=AA*yy/yymax*exp_xyz*tdep    
        //vvy(ix_1,ix_2,ix_3)=-AA*xx/xxmax*exp_xyz*tdep 
        vx=aa*y/yymax*exp_xyz*tdep;    
        vy=-aa*x/xxmax*exp_xyz*tdep; 
 
 switch(field)
  {

    case mom2:
                           dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-vx*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
    break;
    case mom3:     
                           dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-vy*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
    break;
    case energy:
                          dw[fencode3_MODID(p,ii,field)]=dw[fencode3_MODID(p,ii,field)]-(vx*vx+vy*vy)*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;
    break;
   }
 
   


  return ( status);
}

