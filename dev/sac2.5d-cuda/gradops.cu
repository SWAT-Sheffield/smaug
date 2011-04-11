__device__ __host__
int dimproduct_MODID (struct params *dp) {

  int tot=1;
  for(int i=0;i<NDIM;i++)
    tot*=dp->n[i];
  return tot; 
}






__device__ __host__
int fencode_MODID (struct params *dp,int ix, int iy, int field) {


    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));

}




__device__ __host__
real grad_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_MODID;
 real grad=0;

 
 

 switch(dir)
 {
   case 0:
 
// return(  ( (p->sodifon)?((8*wmod[fencode_MODID(p,i+1,j,field)]-8*wmod[fencode_MODID(p,i-1,j,field)]+wmod[fencode_MODID(p,i-2,j,field)]-wmod[fencode_MODID(p,i+2,j,field)])/6.0):wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
 if(i>1 && i<((p->n[0])-2) )
 grad=(  ( ((8*wmod[fencode_MODID(p,i+1,j,field)]-8*wmod[fencode_MODID(p,i-1,j,field)]+wmod[fencode_MODID(p,i-2,j,field)]-wmod[fencode_MODID(p,i+2,j,field)])/6.0))/(2.0*(p->dx[0]))    );


   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
       grad=0;



   break;

   case 1:

// return(  ( (p->sodifon)?((8*wmod[fencode_MODID(p,i,j+1,field)]-8*wmod[fencode_MODID(p,i,j-1,field)]+wmod[fencode_MODID(p,i,j-2,field)]-wmod[fencode_MODID(p,i,j+2,field)])/6.0):wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)])/(2.0*(p->dx[1]))    ); 
 if( j >1 &&  j<((p->n[1])-2))
	grad=(  ( ((8*wmod[fencode_MODID(p,i,j+1,field)]-8*wmod[fencode_MODID(p,i,j-1,field)]+wmod[fencode_MODID(p,i,j-2,field)]-wmod[fencode_MODID(p,i,j+2,field)])/6.0))/(2.0*(p->dx[1]))    );

   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;

   break;
}



 return grad;
}

__device__ __host__
real gradd0_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{

 return(  ( ((8*wmod[fencode_MODID(p,i+1,j,field)]-8*wmod[fencode_MODID(p,i-1,j,field)]+wmod[fencode_MODID(p,i-2,j,field)]-wmod[fencode_MODID(p,i+2,j,field)])/6.0))/(2.0*(p->dx[0]))    );

}

__device__ __host__
real gradd1_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 
return(  ( ((8*wmod[fencode_MODID(p,i,j+1,field)]-8*wmod[fencode_MODID(p,i,j-1,field)]+wmod[fencode_MODID(p,i,j-2,field)]-wmod[fencode_MODID(p,i,j+2,field)])/6.0))/(2.0*(p->dx[1]))    ); 

}



__device__ __host__
real grad1l_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 real grad;
 if((dir == 0) && i>0 && i<(p->n[0]))
 {
    grad=(  ( wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i-1,j,field)]) /((p->dx[0]))    );

   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
       grad=0;
 }
 else if((dir == 1)    && j>0 && j<(p->n[1]))
 {
    grad=(  ( wmod[fencode_MODID(p,i,j,field)]-wmod[fencode_MODID(p,i,j-1,field)])/((p->dx[1]))    );

   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;


  }
 return grad;

}

__device__ __host__
real grad1r_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
  real grad;

  if((dir == 0) && i>=0 && i<((p->n[0])-1))
 {
    grad=(  ( wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i,j,field)]) /((p->dx[0]))    );
   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
       grad=0;
 }
 else if((dir == 1)    && j>=0 && j<((p->n[1])-1))
 {
    grad=(  ( wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j,field)])/((p->dx[1]))    );
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;
  }
 return grad;

}



__device__ __host__
real grad1_MODID(real *wmod,struct params *p,int i,int j,int field,int dir)
{
 //real valgrad_MODID;
  real grad;
  if((dir == 0) && i>0 && i<(p->n[0]))
 {
  
 grad=(  (wmod[fencode_MODID(p,i+1,j,field)]-wmod[fencode_MODID(p,i-1,j,field)])/(2.0*(p->dx[0]))    );
   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
       grad=0;
 }
 else if((dir == 1)    && j>0 && j<(p->n[1]))
 {

 grad=(  (wmod[fencode_MODID(p,i,j+1,field)]-wmod[fencode_MODID(p,i,j-1,field)])/(2.0*(p->dx[1]))    );
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;
  }
 return grad;
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
void bc_cont_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)]; 
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];                     
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,((p->n[1])-3),f)];                   
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,((p->n[0])-3),j,f)];                  
                  else                  
                   // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                        
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];                   
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];                  
                }                       
                else if(i==0 || i==1)                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];   
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];              
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];    
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];                            
                else if(j==0 || j==1)                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)]; 
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                    
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
                




}



__device__ __host__
void bc_cont_cd4_MODID(real *wt, struct params *p,int i, int j, int f) {

                
                if(i==0)              
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4,j,f)];
                else if(i==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,3,j,f)];
                else if( i==((p->n[0])-1))               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-5,j,f)];
                else if (i==((p->n[0])-2))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4,j,f)];
               


                if(j==0)               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4,f)];
                else if(j==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                else if (j== ((p->n[1])-1))               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-5,f)];
               else if (j== ((p->n[1])-2))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];

}



__device__ __host__
void bc_fixed_MODID(real *wt, struct params *p,int i, int j, int f, real val) {


                //(UPPER or LOWER)*NDIM*NVAR+dim*NVAR+varnum = picks out correct value for fixed BC
                //for array of values for fixed BC's

                if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_MODID(p,i,j,f)]=val;
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=val;                   
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=val;                  
                }                       
                else if(i==0 || i==1)                
                  wt[fencode_MODID(p,i,j,f)]=val;                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[fencode_MODID(p,i,j,f)]=val;                
                else if(j==0 || j==1)                
                  wt[fencode_MODID(p,i,j,f)]=val;                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[fencode_MODID(p,i,j,f)]=val;
                




}

__device__ __host__
void bc_periodic_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 || i==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];                
                else if(j==0 || j==1)                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];                                
               else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];
                 //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];


               if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];                                     
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                  else                  
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                    
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                }                       
                 
                




}

__device__ __host__
void bc_periodic1_test_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 || i==1 )                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];
                //else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) || (i==((p->n[0])-3)))
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i-2,j,f)];
                else if(j==0 || j==1 )                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                //else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) || (j==((p->n[1])-3)))
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j-2,f)];

 


}

__device__ __host__
void bc_periodic2_test_MODID(real *wt, struct params *p,int i, int j, int f) {


               if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];                                     
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j-2,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i-2,j,f)];                                    
                  else                  
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j-2,f)];                                    
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i-2,j,f)];                                    
                }                       
                 
                




}

//bc's are not applied to ghost cells?
__device__ __host__
void bc_periodic1a_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i==2 || i==3 )                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];
                //else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) || (i==((p->n[0])-3)))
                else if((i==((p->n[0])-3)) || (i==((p->n[0])-4)))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];
                else if(j==2 || j==3 )                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];
                //else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) || (j==((p->n[1])-3)))
                else if((j==((p->n[1])-3)) || (j==((p->n[1])-4)) )                 
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];



}

//periodic bc's labelled ori below
//are the original ones I used
__device__ __host__
void bc_periodic1_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 || i==1 )                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];
                //else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) || (i==((p->n[0])-3)))
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];


                if(j==0 || j==1 )                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];
                //else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) || (j==((p->n[1])-3)))
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];

 


}


__device__ __host__
void bc_periodic2_MODID(real *wt, struct params *p,int i, int j, int f) {


               if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];                                     
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                  else                  
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                    
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                }                       
                 
                




}

//periodic bc's labelled ori below
//are the original ones I used
__device__ __host__
void bc_periodic1_original_MODID(real *wt, struct params *p,int i, int j, int f) {

                if(i==0 || i==1 )                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];
                //else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) || (i==((p->n[0])-3)))
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];
                else if(j==0 || j==1 )                
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];
                //else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) || (j==((p->n[1])-3)))
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];

 


}

__device__ __host__
void bc_periodic2_original_MODID(real *wt, struct params *p,int i, int j, int f) {


               if(i<2 && j<2)
                {
                  if(i==j)
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];
                  else                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4+i,j,f)];                                     
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                  else                  
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4+j,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4-(p->n[1])+j,f)];                                    
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4-(p->n[0])+i,j,f)];                                    
                }                       
                 
                




}

//*********************************************************************************************************
//generalised 3d routines here
//**********************************************************************************************




__device__ __host__
int fencode3_MODID (struct params *dp,int *ii, int field) {

#ifdef USE_SAC
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif
#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#endif
  //return (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix);
}

__device__ __host__
int encode3_MODID (struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

__device__ __host__
real grad3d_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{

 //real valgrad_MODID;
 real grad=0;

 
 

 switch(dir)
 {
   case 0:
 
#ifdef USE_SAC_3D
 if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],ii[2],field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],ii[2],field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],ii[2],field)])/6.0))/(2.0*(p->dx[0]))    );

   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
#else
 if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],0,field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],0,field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],0,field)])/6.0))/(2.0*(p->dx[0]))    );

   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
#endif



   break;

   case 1:

#ifdef USE_SAC_3D
 if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,ii[2],field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,ii[2],field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,ii[2],field)])/6.0))/(2.0*(p->dx[1]))    );

   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
#else
 if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,0,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,0,field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,0,field)])/6.0))/(2.0*(p->dx[1]))    );

   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
#endif
   break;


   case 2:

#ifdef USE_SAC_3D
 if( ii[2] >1 &&  ii[2]<((p->n[2])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+1,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-1,field)]+wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+2,field)])/6.0))/(2.0*(p->dx[2]))    );

   if((ii[2]==(p->n[2])-3) || (ii[2]==(p->n[2])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[2]==2 || ii[2]==3  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
#endif
   break;

}



 return grad;


}

__device__ __host__
real grad1l3_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
 real grad;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && i>0 && i<(p->n[0]))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)]) /((p->dx[0]))    );

   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2 && k>1   && k<(p->n[2])-2 )
	       grad=0;
	   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2 && k>1   && k<(p->n[2])-2 )
	       grad=0;
   #else
	   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
	       grad=0;
	   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
	       grad=0;
   #endif
 }
 else if((dir == 1)    && j>0 && j<(p->n[1]))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/((p->dx[1]))    );
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2 && k>1   && k<(p->n[2])-2 )
	       grad=0;
	   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2 && k>1   && k<(p->n[2])-2 )
	       grad=0;
   #else
	   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
	       grad=0;
	   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
	       grad=0;
   #endif


  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<(p->n[2]))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/((p->dx[2]))    );

   if((k==(p->n[2])-3) || (k==(p->n[2])-4)  && i>1   && i<(p->n[0])-2  && j>1   && j<(p->n[1])-2 )
       grad=0;
   else if(k==2 || k==3  && i>1   && i<(p->n[0])-2  && j>1   && j<(p->n[1])-2 )
       grad=0;


  }
  #endif
 return grad;

}

__device__ __host__
real grad1r3_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
  real grad;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


  if((dir == 0) && i>=0 && i<((p->n[0])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i,j,k,field)]) /((p->dx[0]))    );
   #ifdef USE_SAC_3D
   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  && k>1   && k<(p->n[2])-2)
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2 && k>1   && k<(p->n[2])-2 )
       grad=0;
   #else
   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
       grad=0;
   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
       grad=0;
   #endif
 }
 else if((dir == 1)    && j>=0 && j<((p->n[1])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((p->dx[1]))    );
   #ifdef USE_SAC_3D
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  && k>1   && k<(p->n[2])-2 )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  && k>1   && k<(p->n[2])-2)
       grad=0;
   #else
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;
   #endif
  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>=0 && k<((p->n[2])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((p->dx[2]))    );

   if((k==(p->n[2])-3) || (k==(p->n[2])-4)  && i>1   && i<(p->n[0])-2  && j>1   && j<(p->n[1])-2)
       grad=0;
   else if(k==2 || k==3  && i>1   && i<(p->n[0])-2  && j>1   && j<(p->n[1])-2)
       grad=0;
  }
   #endif
 return grad;

}



__device__ __host__
real grad13_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
 //real valgrad_MODID;
  real grad;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


  if((dir == 0) && i>0 && i<(p->n[0]))
 {
  
 grad=(  (wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)])/(2.0*(p->dx[0]))    );
   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2   && k>1   && k<(p->n[2])-2)
	       grad=0;
	   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2 && k>1   && k<(p->n[2])-2  )
	       grad=0;
	#else
	   if((i==(p->n[0])-3) || (i==(p->n[0])-4)  && j>1   && j<(p->n[1])-2  )
	       grad=0;
	   else if(i==2 || i==3  && j>1   && j<(p->n[1])-2  )
	       grad=0;
	#endif
 }
 else if((dir == 1)    && j>0 && j<(p->n[1]))
 {

 grad=(  (wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/(2.0*(p->dx[1]))    );
   #ifdef USE_SAC_3D
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2  && k>1   && k<(p->n[2])-2  )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2 && k>1   && k<(p->n[2])-2  )
       grad=0;
   #else
   if((j==(p->n[1])-3) || (j==(p->n[1])-4)  && i>1   && i<(p->n[0])-2   )
       grad=0;
   else if(j==2 || j==3  && i>1   && i<(p->n[0])-2  )
       grad=0;
    #endif
  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<(p->n[2]))
 {

 grad=(  (wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/(2.0*(p->dx[2]))    );
   if((k==(p->n[2])-3) || (k==(p->n[2])-4)  && i>1   && i<(p->n[0])-2  && j>1   && j<(p->n[1])-2)
       grad=0;
   else if(k==2 || k==3  && i>1   && i<(p->n[0])-2 && j>1   && j<(p->n[1])-2  )
       grad=0;
  }

   #endif
 return grad;
}




__device__ __host__
void bc3_cont_MODID(real *wt, struct params *p,int *ii, int f) {

   
int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];

                if(i<2 && j<2  && k<2)
                {
                 if(i==j==k )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,2,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  &&  k<2)
                {
                  if(i==(j-(p->n[1]))  && k==(j-(p->n[1])))                  
                     wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,2,f)];                     
                }
                else if(i>((p->n[0])-3) && j<2 && k<2)
                {
                  if(j==(i-(p->n[0]))  && k==(i-(p->n[0])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),2,2,f)];                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3) && k<2)
                {
                  if(i==j  && k==(i-(p->n[0])) )                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),((p->n[1])-3),2,f)];                                                  
                }
                else if(i<2 && j<2  && k>((p->n[2])-3))
                {
                 if(i==j && k==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j<2  && k>((p->n[2])-3))
                {
                 if(i==k && j==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,2,(p->n[2])-3,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(j==k && i==(j-(p->n[1]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(i==j==k  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,(p->n[1])-3,(p->n[2])-3,f)];
                }                     
                else if(i==0 || i==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,j,k,f)];              
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,j,k,f)];                            
                else if(j==0 || j==1)                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-3,k,f)];
                else if(k==0 || k==1)                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-3,f)];

        #else
             if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];                     
                  else                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,((p->n[1])-3),f)];                   
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,((p->n[0])-3),j,f)];                  
                  else                  
                   // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                        
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];                   
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];                  
                }                       
                else if(i==0 || i==1)                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];   
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];              
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];    
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];                            
                else if(j==0 || j==1)                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)]; 
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                    
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
                
         #endif



}



__device__ __host__
void bc3_cont_cd4_MODID(real *wt, struct params *p,int *ii, int f) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
                if(i==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4,j,k,f)];
                else if(i==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,3,j,k,f)];
                else if( i==((p->n[0])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-5,j,k,f)];
                else if (i==((p->n[0])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4,j,k,f)];
               


                if(j==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4,k,f)];
                else if(j==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,3,k,f)];
                else if (j== ((p->n[1])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-5,k,f)];
               else if (j== ((p->n[1])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4,k,f)];




                if(k==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4,f)];
                else if(k==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,3,f)];
                else if (k== ((p->n[2])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-5,f)];
               else if (k== ((p->n[2])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4,f)];
        #else
                
                if(i==0)              
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,4,j,f)];
                else if(i==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,3,j,f)];
                else if( i==((p->n[0])-1))               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-5,j,f)];
                else if (i==((p->n[0])-2))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4,j,f)];
               


                if(j==0)               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4,f)];
                else if(j==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                else if (j== ((p->n[1])-1))               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-5,f)];
               else if (j== ((p->n[1])-2))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
         #endif

}



__device__ __host__
void bc3_fixed_MODID(real *wt, struct params *p,int *ii, int f, real val) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
        #endif

                   
                if(i==0 || i==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=val;                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=val;                
                else if(j==0 || j==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=val;                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=val;
           #ifdef USE_SAC_3D
                else if(k==0 || k==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=val;                
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2)))                
                  wt[encode3_MODID(p,i,j,k,f)]=val;
           #endif
               




}


__device__ __host__
void bc3_periodic1_MODID(real *wt, struct params *p,int *ii, int f) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
                if(i==0 || i==1 )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];
             
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];


                if(j==0 || j==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];

                if(k==0 || k==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4+j,f)];

                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4-(p->n[2])+k,f)];
       #else

                if(i==0 || i==1 )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];

                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];


                if(j==0 || j==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];

       #endif


}


__device__ __host__
void bc3_periodic2_MODID(real *wt, struct params *p,int *ii, int f) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];

                if(i<2 && j<2  && k<2)
                {
                 if(i==j==k )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                }
                else if(i<2 && j>((p->n[1])-3)  &&  k<2)
                {
                  if(i==(j-(p->n[1]))  && k==(j-(p->n[1])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                     
             
                }
                else if(i>((p->n[0])-3) && j<2 && k<2)
                {
                  if(j==(i-(p->n[0]))  && k==(i-(p->n[0])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
           
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3) && k<2)
                {
                  if(i==j  && k==(i-(p->n[0])))   
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                    
                                 
                                            
                }
                else if(i<2 && j<2  && k>((p->n[2])-3))
                {
                 if(i==j && i==(k-(p->n[2]))  )                 
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4+k,f)];                                     

                } 
                else if(i>((p->n[0])-3) && j<2  && k>((p->n[2])-3))
                {
                 if(i==k && j==(i-(p->n[0]))  )
                     wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,4-(p->n[2])+k,f)];   
                }
                else if(i<2 && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(j==k && i==(j-(p->n[1]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,4-(p->n[2])+k,f)]; 
                } 
                else if(i>((p->n[0])-3) && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(i==j==k  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,4-(p->n[1])+j,4-(p->n[2])+k,f)]; 
                }   

        #else

               if(i<2 && j<2)
                {
                  if(i==j)
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];
                  else                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                     
                  else                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
                  else                  
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                    
                  else                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
                }                       
                 
       #endif         




}



