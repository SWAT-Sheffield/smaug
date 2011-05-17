// iosac2.5d.cpp : Defines the entry point for the console application.
//


#include "iosac2.5d.h"
#include "step.h"

/*----------------------*/ 
real second()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (real) Time / (real) CLOCKS_PER_SEC;
   return secs;*/
   real retval;
	static long zsec=0;
	static long zusec=0;
	real esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}




int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int it=0; //test integer to be returned 

//getintparam_( int elist.id,char *sname,int *iv,  int elist.port, char *selist.server );
//int elist.id=0;
//int elist.port=8080;

int i1,i2,i3,j1;
int i,j,k,iv;


char *portfile=(char *)calloc(500,sizeof(char));
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));



#include "iosac2.5dparams.h"

FILE *portf;

#ifdef USE_IOME
if(argc>1)
{ 
   sprintf(portfile,"%s0_port.txt",argv[1]) ;  
   //strcpy(name,argv[1]);
   portf=fopen(portfile,"r");
   fscanf(portf,"%d %s",&elist.port,elist.server);
   fclose(portf);

   printf("read file junk is %d %s\n",elist.port,elist.server);
}
#endif





//int cuprop(struct params **p, real **w, real **wnew, real **b,struct params **d_p, real **d_w, real **d_wnew, real **d_b, real **d_wmod, real **d_dwn1, real **d_dwn2, real **d_dwn3, real **d_dwn4, real **d_wd)




printf("rho %d mom1 %d mom2 %d\n",rho,mom1,mom2);


printf("calling cuinit\n");



 //   getintparam_( &elist.id,"i1",&it,  &elist.port, "localhost" );	
//	printf("Get integer %d\n",it);
    //Set input filename as first arg
	//if NULL use defaults
	char *method=NULL;
	//CIoSimulation *TestSimulation;
	//this should be executed by the iome start up application
	//exec('ioshallowwater.sce');

	//this application is started using the io  start scilab application
	//exec('paramssteeringtest1.sce');
	//stacksize('max');
	//stacksize(268435454)
	//open the file generated
	//sprintf(elist.portfile,"%s0_elist.port.txt",meta.name);
	//FILE *fd=fopen(elist.portfile,"r");
	//int elist.portelist.id;
	//fscanf(fd,"%d",&elist.portelist.id);
	//fclose(fd);
	//elist.elist.port=elist.portelist.id;

    #ifdef USE_IOME
        if(argc>2)
        {
          //simfile already read by 
          readsim(p,&meta,argv[2],elist);
          //if((p->readini)!=0)
          //   readconfig(meta.ini_file,*p,meta,w);
        }
        else
	  createsim(*p,meta,simfile,elist);

	sprintf(simfile,"%s.xml",meta.name);
        sprintf(newsimfile,"%s_update.xml",meta.name);
     #endif
	//NewSimulation(metadata.name,'test1.xsl',elist);

// Build empty u, v, b matrices
// Define h
printf("allocating w and wnew\n");

  #ifdef USE_SAC_3D
 w=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));
wd=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
 wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));

 #else
 w=(real *)calloc(ni*nj*NVAR,sizeof(real ));
wd=(real *)calloc(ni*nj*NDERV,sizeof(real ));
 wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));
#endif

if((p->readini)==0)
 initconfig(p, &meta, w);
else
 readasciivacconfig(cfgfile,*p,meta,w,hlines);




//writeasciivacconfig(cfgout,*p, meta , w,hlines,*state);
//writevacconfig(cfgout,0,*p, meta , w,*state);
  /*   for( j1=2;j1<5;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
        //j1=2;
	printf("%d %d %f %f %f\n",i1,j1,w[j1*ni+i1+(ni*nj*rho)],w[j1*ni+i1+(ni*nj*mom1)],w[j1*ni+i1+(ni*nj*mom2)]);
        }     
       //

      }
 printf("\n");*/  
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(ni)*(nj)*rho;
  u=w+(ni)*(nj)*mom1;
  v=w+(ni)*(nj)*mom2;

cuinit(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
cuboundary(&p,&d_p,&d_state,&d_w, 0);
cuboundary(&p,&d_p,&d_state,&d_wmod, 0);

printf("after cuinit\n");
/*for( j1=0;j1<nj;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
              for(iv=0;iv<NVAR;iv++)               
                      printf("%f ", w[j1*ni+i1+(ni*nj*iv)]);
               printf("\n");

           }
         }*/


//For a steerable simulation generate and save a dxformfile that saves a single data step
//used for the steering dx module
//printf("here in runsim2a\n");

#ifdef USE_IOME
getmetadata_(elist.id,"directory",&sdir,elist.port,elist.server);
//sdir=metadata.directory

//name=metadata.name;

getmetadata_(elist.id,"name",&name,elist.port,elist.server);
//disp(sdir,name)
//printf("here in runsim3\n");
sprintf(outfile,"%s/%s.out",sdir,name);
#endif



//createlog(meta.log_file);



//while(finishsteering == 0)
//{
 
  //  if( steeringenabled==0)
  //    finishsteering=1;
 int n;  
// nt=24; 
real t1,t2,ttot;
int order=0;
int ordero=0;
int order1;
int orderb=0;
int ii,ii0,ii1;
ttot=0;
real time=0.0;
   state->it=0;
   state->t=0;
   state->dt=p->dt;

//printf("dx dy%f %f",dx,dy);
for( n=1;n<=nt;n++)
//for( n=0;n<1;n++)
{

    if((n%(p->cfgsavefrequency))==0)
    {
      //writeconfig(name,n,*p, meta , w);
      writevtkconfig(name,n,*p, meta , w);
      //writeasciivacconfig(cfgout,*p, meta , w,hlines,*state);

      writevacconfig(cfgout,n,*p, meta , w,*state);
    }
   order=0;
   t1=second();

   //printf("cmax is %f old dt %f new dt %f\n",p->cmax,p->dt,0.68*((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax)));
   /*if(p->moddton==1)
   {
      //if((p->cmax)>10)
      //         p->dt=0.68*(((p->dx)+(p->dy))/(2.0*(p->cmax))); 
      if((p->dt)>(((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax))))
               p->dt=(p->dt)/2;
       else if(2.0*(p->dt)<(((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax))))
               p->dt=2.0*p->dt; 
   }
   printf("new dt %f\n",p->dt);*/

if((p->rkon)==0)
{
  ordero=0;
 
  cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
  order=1;

 for(int dir=0;dir<NDIM; dir++)
 {

  cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
  cucomputepres(&p,&d_p,&d_wmod, &d_wd,order,dir);
  for(int f=rho; f<=(mom1+NDIM-1); f++)
  {
      cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
      //cucentdiff1a(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

  }

#ifndef ADIABHYDRO
   for(int f=energy; f<=(b1+NDIM-1); f++)
   {
     cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
     //cucentdiff2a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);

   }
#endif
  }
  // cuboundary(&p,&d_p,&d_wmod, ordero);
   if(p->divbon==1)
	       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
   if(p->hyperdifmom==1)
   {
    dt=(p->dt);
    for(int dim=0; dim<=(NDIM-1); dim++)
     {
       ;//cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);

       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);


       cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);
       //cuhyperdifrhosource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);

     }

     for(int dim=0; dim<=(NDIM-1); dim++)
     {
       ;//cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);

      cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
      //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);

       cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
       //cuhyperdifesource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       //cuhyperdifesource2a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       //cuhyperdifesource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);

     }

for(int dim=0; dim<=(NDIM-1); dim++)
       for(int f=0; f<=(NDIM-1); f++)
           	                 
	     {
       ;//cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);

               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           ii=dim;
                           ii0=f;
                          }
                          else
                           {
                           ii=f;
                           ii0=dim;
                            }

                  if(ii==dim)
                  {
                    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
                   }
                   else
                   {
                    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                   //cuhyperdifmomsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                   //cuhyperdifmomsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
		    }
                }
             }
            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=(NDIM-1); dim++)
	     for(int f=0; f<=(NDIM-1); f++)            
	     {
	       ;//cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
	       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
	       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);


               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);

                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           jj=dim;
                           mm=f;
                           sb=-1.0;
                           ii0=dim;
                          }
                          else
                           {
                           ii0=f;
                           mm=dim;
                           sb=1.0;
                           jj=f;
                           
                            }

                  if(mm==dim)
                  {
                     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     //cuhyperdifbsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     //cuhyperdifbsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     //cuhyperdifbsource4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   }
                   else
                   {
                     cuhyperdifbsourcene1(&p,&d_p,&d_w, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    // cuhyperdifbsourcene1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                    // cuhyperdifbsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                    // cuhyperdifbsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                   // cuhyperdifbsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   // cuhyperdifbsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   // cuhyperdifbsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    }

                }
              } 


   }
   //cuadvance(&p,&d_p,&d_wmod,&d_w, order);
   cuboundary(&p,&d_p,&d_state,&d_wmod, ordero);
}

   if((p->rkon)==1)
     for(order=0; order<4; order++) 
   {	   
           ordero=order+1;
           dt=(p->dt)/2.0;
           orderb=order+2;

           if(order==2)
           {
              dt=(p->dt);
              orderb=1;
            }


           if(order==3)
           {
              dt=(p->dt)/6.0;
              ordero=0;
              orderb=0;
           }


           cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
 for(int dir=0;dir<(NDIM-1); dir++)
 {
           cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);

           for(int f=rho; f<=mom1+(NDIM-1); f++)
           {
	       cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);
	       //cucentdiff1a(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);

            }

#ifndef ADIABHYDRO
           for(int f=energy; f<=b1+(NDIM-1); f++)
           {
	       cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
	       //cucentdiff2a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

           }
#endif

}
	   //cuderivsource(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
	   if(p->divbon==1)
	       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);
           if(p->hyperdifmom==1)
           {
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
	       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);

	       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
      //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
	       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
	       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
	       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);

	       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
      //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
	       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
	       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
	       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);

	      // cuhyperdifviscmax(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,rho,dim);

	       cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);
	       //cuhyperdifrhosource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);

	     }

     for(int dim=0; dim<=(NDIM-1); dim++)
     {
       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);

       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);

      cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
      //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);

       cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
       //cuhyperdifesource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       //cuhyperdifesource2a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       //cuhyperdifesource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);

     }
             
for(int dim=0; dim<=(NDIM-1); dim++)
       for(int f=0; f<=(NDIM-1); f++)
           	                 
	     {
	       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
	       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);

               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
 //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);

               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           ii=dim;
                           ii0=f;
                          }
                          else
                           {
                           ii=f;
                           ii0=dim;
                            }

                  if(ii==dim)
                  {
                    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
                   }
                   else
                   {
                    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    //cuhyperdifmomsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);


		    }
                }
             }

            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=(NDIM-1); dim++)
	     for(int f=0; f<=(NDIM-1); f++)            
	     {
	       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
	       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);

               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
              //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);


               cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               //cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);

                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           jj=dim;
                           mm=f;
                           sb=-1.0;
                           ii0=dim;
                          }
                          else
                           {
                           ii0=f;
                           mm=dim;
                           sb=1.0;
                           jj=f;
                           
                            }

                  if(mm==dim)
                  {
                     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     //cuhyperdifbsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     //cuhyperdifbsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     //cuhyperdifbsource4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   }
                   else
                   {
                     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     //cuhyperdifbsourcene1a(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     //cuhyperdifbsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     //cuhyperdifbsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                    //cuhyperdifbsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    //cuhyperdifbsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    //cuhyperdifbsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    }

                }
              } 

           }
           //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,ordero);
           cuadvance(&p,&d_p,&d_wmod,&d_w,order);
           cuboundary(&p,&d_p,&d_state,&d_wmod, orderb);
	   

   }

   cuupdate(&p,&w,&wd,&state,&d_p,&d_w,&d_wmod,  &d_state,n);
   //printf("nummaxthreads %d\n",p->mnthreads);

   t2=second()-t1;
   ttot+=t2;
   printf("step %d total time %f\n",n,ttot);

   state->it=n;
   state->t=time+(p->dt);
   time=state->t;
   state->dt=p->dt;

   //appendlog(meta.log_file,*p, *state);

 
    
    /*getintparam_(&elist.id,"steeringenabled",&steeringenabled,&elist.port,elist.server);
    if(steeringenabled==1)
    {
      //disp('getting updatea params');
      //for steering get the modified control params
      double dg;
      getintparam_(&elist.id,"finishsteering",&finishsteering,&elist.port,elist.server);//source y location  
        // Constants
      getdoubleparam_(elist.id,"g",&dg,elist.port,elist.server);

      g=dg;
     
    }*/
    
   /* for( j1=ngj;j1<nj-ngj;j1++)
        for( i1=ngi;i1<ni-ngi;i1++)
{

;//w[j1*ni+i1+(ni*nj*b1)]=wd[j1*ni+i1+(ni*nj*(hdnur))];
;//w[j1*ni+i1+(ni*nj*b2)]=wd[j1*ni+i1+(ni*nj*(hdnul))];

}*/
           
      //save file containing current data
     // sprintf(configfile,"tmp/%ss%d.out",name,n);-
     // printf("check dims %d %d \n",ni,nj);
     // FILE *fdt=fopen(configfile,"w");
     // fprintf(fdt,"%d\n",n);
     //for( j1=0;j1<nj;j1++)
     // {
      //  for( i1=0;i1<ni;i1++)
	//{
               // printf("%d %d ", i1,j1);
	//	fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,*(h+(j1*ni+i1)),(u[j1*ni+i1]),(v[j1*ni+i1]),w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
       // }     
        //printf("\n");   
        //fprintf(fdt,"\n");
     // }
     // fclose(fdt);
   
 //disp('writing data');
    //if finishsteering==1
    //  fprintf(fd,"%d\n",n);
    //  for( j1=0;j1<nj;j1++)
    //  {
     //   for( i1=0;i1<ni;i1++)
	//{
     //     fprintf(fd,"%f %f %f ",u[j1*ni+i1],v[j1*ni+i1],h[j1*ni+i1]);
	   	//fprintf(fd,"%f ",h[j1*ni+i1]);
	//}
       // fprintf(fd,"\n");
    //  }


    }//end of testep
  printf("params %f %f %f %f\n",g,dx,dy,dt);  
    //end
  
   
    //end
    
    //disp('written data');
//}//end of steering test finished

  // force the final ouput file to be over written
 // if(steeringenabled==1)
 // {
  //  if(finishsteering==0)
 //   {
  //    fclose(fd);
  //    fd=fopen(outfile,"w");
  //  }
 //  }

//}
//}//disp('while finsish steering');
//}//end //while finishsteering loop
//cufinish(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd);
cufinish(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
free(hlines);
free(p);
free(sdir);
free(name);
free(outfile);
free(formfile);

//for the completed simulation
//int nsteps=nt;
//fclose(fd);






	//[consts,domain,source]=loadsim('test1_16_02_09.xml',elist);
	//chdir(metadata.directory);
        //readsimulation_(elist.elist.id,simfile,elist.elist.port,elist.elist.server);
	//runsim(consts,dom,src,meta,simfile,elist);
#ifdef USE_IOME
	writesimulation_(elist.id,newsimfile,elist.port,elist.server);
#endif


	return 0;
}

