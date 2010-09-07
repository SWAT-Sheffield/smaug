__device__ __host__
int dimproduct_MODID (struct params *dp) {

  int tot=1;
  for(int i=0;i<NDIM;i++)
    tot*=dp->n[i];
  return tot; 
}


__device__ __host__
int encode_MODID (struct params *dp,int ix, int iy) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int encode3_MODID (struct params *dp,int ix, int iy, int iz) {

  return (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int fencode_MODID (struct params *dp,int ix, int iy, int field) {

  //int kSizeX=(dp)->n[0];
  //int kSizeY=(dp)->n[1];
  
  return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
}

 
__device__ __host__
int fencode3_MODID (struct params *dp,int ix, int iy, int iz, int field) {

  //int kSizeX=(dp)->ni;
  //int kSizeY=(dp)->nj;
  
  return(  iz*((dp)->n[0])*((dp)->n[1])+ (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2]))  );
}


__device__ __host__
real evalgrad_MODID(real fi, real fim1, real fip2, real fim2,struct params *p,int dir)
{
 //real valgrad_MODID;

 if(dir == 0)
 {
     //valgrad=(2.0/(3.0*(p->dx[0])))*(fi-fim1)-(1.0/(12.0*(p->dx[0])))*(fip2-fim2);
   //return((1.0/(2.0*(p->dx[0])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[0])))*(fi-fim1)):((1.0/(12.0*(p->dx[0])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(fi-fim1)-(1.0/(12.0*(p->dx[1])))*(fip2-fim2);
     // return((2.0/(1.0*(p->dx[1])))*(fi-fim1));
   return(p->sodifon?((1.0/(2.0*(p->dx[1])))*(fi-fim1)):((1.0/(12.0*(p->dx[1])))*((NVAR*fi-NVAR*fim1+fim2-fip2))));
 }

 return -1;
}


__device__ __host__
real grad_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_MODID;

 if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx[0])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx[0])))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
//return((1.0/(2.0*(p->dx[0])))*(wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_MODID(p,i+1,j,field)]-NVAR*wmod[fencode_MODID(p,i-1,j,field)]+wmod[fencode_MODID(p,i-2,j,field)]-wmod[fencode_MODID(p,i+2,j,field)])/6.0):wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dx[1])))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
// return((1.0/(2.0*(p->dx[1])))*(wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)]));
 return(  ( (p->sodifon)?((NVAR*wmod[fencode_MODID(p,i,j+1,field)]-NVAR*wmod[fencode_MODID(p,i,j-1,field)]+wmod[fencode_MODID(p,i,j-2,field)]-wmod[fencode_MODID(p,i,j+2,field)])/6.0):wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)])/(2.0*(p->dx[1]))    );  
}


 return 0;
}



__device__ __host__
real evalgrad1_MODID(real fi, real fim1, struct params *p,int dir)
{
 //real valgrad_MODID;

 if(dir == 0)
 {

   return(((1.0/(2*(p->dx[0])))*(fi-fim1)));
 }
 else if(dir == 1)
 {

   return(((1.0/(2*(p->dx[1])))*(fi-fim1)));
 }

 return -1;
}
__device__ __host__
real grad1l_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
  if(dir == 0)
 {
    return(  ( wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i-1,j,field)]) /((p->dx[0]))    );
 }
 else if(dir == 1)
 {
    return(  ( wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i,j-1,field)])/((p->dx[1]))    );
  }
 return 0;

}

__device__ __host__
real grad1r_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
  if(dir == 0)
 {
    return(  ( wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i,j,field)]) /((p->dx[0]))    );
 }
 else if(dir == 1)
 {
    return(  ( wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j,field)])/((p->dx[1]))    );
  }
 return 0;

}



__device__ __host__
real grad1_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_MODID;

  if(dir == 0)
 {

 return(  (wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {

 return(  (wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)])/(2.0*(p->dx[1]))    );
  }
 return 0;
}



__device__ __host__
real grad2_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_MODID;

  if(dir == 0)
 {
    // valgrad=(2.0/(3.0*(p->dx[0])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i-1,j,field)])-(1.0/(12.0*(p->dx[0])))*(wmod[fencode(p,i+2,j,field)]-wmod[fencode(p,i-2,j,field)]);
//return((1.0/(2.0*(p->dx[0])))*(wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)]));
 return(  ( (p->sodifon)?((16*wmod[fencode_MODID(p,i+1,j,field)]+16*wmod[fencode_MODID(p,i-1,j,field)]-wmod[fencode_MODID(p,i-2,j,field)]-wmod[fencode_MODID(p,i+2,j,field)]-30*wmod[fencode_MODID(p,i,j,field)])/6.0):2.0*(wmod[fencode_MODID(p,i+1,j,field)]-2*wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i-1,j,field)]))/(2.0*(p->dx[0])*(p->dx[0]))    );
 }
 else if(dir == 1)
 {
    // valgrad=(2.0/(3.0*(p->dx[1])))*(wmod[fencode(p,i,j,field)]-wmod[fencode(p,i,j-1,field)])-(1.0/(12.0*(p->dx[1])))*(wmod[fencode(p,i,j+2,field)]-wmod[fencode(p,i,j-2,field)]);
// return((1.0/(2.0*(p->dx[1])))*(wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)]));
 return(  ( (p->sodifon)?((16*wmod[fencode_MODID(p,i,j+1,field)]+16*wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i,j-2,field)]-wmod[fencode_MODID(p,i,j+2,field)]-30*wmod[fencode_MODID(p,i,j,field)])/6.0):2.0*(wmod[fencode_MODID(p,i,j+1,field)]-2.0*wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)]))/(2.0*(p->dx[1])*(p->dx[1]))    );
  }
 return 0;
}


__device__ __host__
real grad3_MODID(real *wmod,struct params *p,int *ix,int field,int dir)
{
 //real valgrad;

 if(dir == 0)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_MODID(p,ix[0]+1,ix[1],ix[2],field)]-NVAR*wmod[fencode3_MODID(p,ix[0]-1,ix[1],ix[2],field)]+wmod[fencode3_MODID(p,ix[0]-2,ix[1],ix[2],field)]-wmod[fencode3_MODID(p,ix[0]+2,ix[1],ix[2],field)])/6.0):wmod[fencode3_MODID(p,ix[0]+1,ix[1],ix[2],field)]-wmod[fencode3_MODID(p,ix[0]-1,ix[1],ix[2],field)])/(2.0*(p->dx[0]))    );
 }
 else if(dir == 1)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_MODID(p,ix[0],ix[1]+1,ix[2],field)]-NVAR*wmod[fencode3_MODID(p,ix[0],ix[1]-1,ix[2],field)]+wmod[fencode3_MODID(p,ix[0],ix[1]-2,ix[2],field)]-wmod[fencode3_MODID(p,ix[0],ix[1]+2,ix[2],field)])/6.0):wmod[fencode3_MODID(p,ix[0],ix[1]+1,ix[2],field)]-wmod[fencode3_MODID(p,ix[0],ix[1]-1,ix[2],field)])/(2.0*(p->dx[1]))    );

 }
else if(dir == 2)
 {

 return(  ( (p->sodifon)?((NVAR*wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]+1,field)]-NVAR*wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]-1,field)]+wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]-2,field)]-wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]+2,field)])/6.0):wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]+1,field)]-wmod[fencode3_MODID(p,ix[0],ix[1],ix[2]-1,field)])/(2.0*(p->dx[2]))    );

 }
 return -1;
}


