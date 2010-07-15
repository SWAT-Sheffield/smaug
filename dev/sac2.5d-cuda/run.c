#include "paramssteeringtest1.h"
#include "step.h"
#include <iome/genericsimulationlib/IoGenericSimulationLib.h>
void gendxgen(char *dir,char *jobname,int nsteps,int n1,int n2);


/*----------------------*/ 
double second()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (double) Time / (double) CLOCKS_PER_SEC;
   return secs;*/
   double retval;
	static long zsec=0;
	static long zusec=0;
	double esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}


void runsim(constants k, domain dom,source src, meta metadata,char *simname, iome el)
{
int i,j;
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));
char configfile[300];
//elist=list();  parameter used by iome to contain port and server address
//elist=list(); 

float sf=src.freq;//source frequency
float sa=src.amp;//source amplitude
float sx=src.xloc;//source x location
float sy=src.yloc;//source y location

// Constants
float g  = k.g;
float u0 = k.u0;                               
float v0 = k.v0;
float b0  = k.b0;                               
float h0 = k.h0; 


//Domain definition
// Define the x domain
//ni = 151; 
int ni=dom.ni;
float xmax = dom.xmax;                      
float dx = xmax/(ni-1);
float *x=(float *)calloc(ni,sizeof(float));
for(i=0;i<ni;i++)
		x[i]=i*dx;
int i1,i2,i3,j1;
// Define the y domain
//nj = 151;  
int nj=dom.nj;
float ymax = dom.ymax;                      
float dy = ymax/(nj-1);
float *y=(float *)calloc(nj,sizeof(float));
for(i=0;i<nj;i++)
		y[i]=i*dy;


float tmax = dom.tmax;
int steeringenabled=dom.steeringenabled;
int finishsteering=dom.finishsteering;


// Define the wavespeed
float wavespeed = u0 + sqrt(g*(h0 - b0));

// Define time-domain
float dt = 0.68*dx/wavespeed;
int nt=(int)((tmax-1)/dt);

float *t=(float *)calloc(nt,sizeof(float));
printf("runsim 1\n");
//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;
dom.nt=nt;

float courant = wavespeed*dt/dx;


//int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_dwn2, float **d_dwn3, float **d_dwn4, float **d_wd)


// Build empty u, v, b matrices
// Define h
float *w=(float *)calloc(ni*nj*8,sizeof(float ));
float *wnew=(float *)calloc(ni*nj*8,sizeof(float ));
float *b=(float *)calloc(ni*nj,sizeof(float ));
  float *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(ni)*(nj)*rho;
  u=w+(ni)*(nj)*mom1;
  v=w+(ni)*(nj)*mom2;

printf("rho %d mom1 %d mom2 %d\n",rho,mom1,mom2);

float *d_w;
float *d_wnew;
float *d_b;
float *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

struct params *d_p;
struct params *p=(struct params *)malloc(sizeof(struct params));

p->ni=ni;
p->nj=nj;
p->dt=dt;
p->dx=dx;
p->dy=dy;
p->g=g;
p->gamma=0.666666;
p->mu=1.0;
p->eta=0.0;
p->g1=0.0;
p->g2=0.0;
p->g3=0.0;
printf("calling cuinit\n");
cuinit(&p,&w,&wnew,&b,&d_p,&d_w,&d_wnew,&d_b,&d_wmod, &d_dwn1,  &d_wd);


printf("here in runsim\n");


printf("here in runsim1\n");







int nli = 0.45*(ni-1)+1;
int nui = 0.55*(ni-1)+1;
int nlj = 0.45*(nj-1)+1;
int nuj = 0.55*(nj-1)+1; 

printf("limits %d %d %d %d\n",nli,nui,nlj,nuj);                           
int in,ind;

printf("here in runsim2\n");

//For a steerable simulation generate and save a dxformfile that saves a single data step
//used for the steering dx module
//printf("here in runsim2a\n");
getmetadata_(el.id,"directory",&sdir,el.port,el.server);
//sdir=metadata.directory

//name=metadata.name;

getmetadata_(el.id,"name",&name,el.port,el.server);
//disp(sdir,name)
//printf("here in runsim3\n");
sprintf(outfile,"%s/%s.out",sdir,name);

FILE *fd=fopen(outfile,"w");
//if steeringenabled==1
 printf("\n %s %s here in runsim4 %s\n",sdir,name,outfile); 
  //mkdir('tmp');
  gendxgen(sdir,name,nt,ni,nj);
printf("here in runsim5\n");

sprintf(formfile,"%s/form%s.out",sdir,name);
FILE *fdform=fopen(formfile,"w");
  fprintf(fdform, "%d %d %d\n",nt-1, ni, nj);
fclose(fdform);


// Employ Lax
//disp(length[t));


//while(finishsteering == 0)
//{
 
  //  if( steeringenabled==0)
  //    finishsteering=1;
 int n;  
 nt=24; 
double t1,t2,ttot;
ttot=0;
for( n=0;n<nt;n++)
{
  
   t1=second();
   cuprop(&p,&w,&wnew,&b,&d_p,&d_w,&d_wnew,&d_b,&d_wmod, &d_dwn1, &d_wd);
   t2=second()-t1;
   ttot+=t2;
   printf("step %d total time %f\n",n,ttot);



 
    
    getintparam_(&el.id,"steeringenabled",&steeringenabled,&el.port,el.server);
    if(steeringenabled==1)
    {
      //disp('getting updatea params');
      //for steering get the modified control params
      double dsf,dsa,dsx,dsy,dg;
      getdoubleparam_(el.id,"frequency",&dsf,el.port,el.server);//source frequency
      getdoubleparam_(el.id,"amplitude",&dsa,el.port,el.server);//source amplitude
      getdoubleparam_(el.id,"sx",&dsx,el.port,el.server);//source x location
      getdoubleparam_(el.id,"sy",&dsy,el.port,el.server);//source y location
      getintparam_(&el.id,"finishsteering",&finishsteering,&el.port,el.server);//source y location  
        // Constants
      getdoubleparam_(el.id,"g",&dg,el.port,el.server);
      sf=dsf;
      sa=dsa;
      sx=dsx;
      sy=dsy;
      g=dg;
     
    }
    
    
      //save file containing current data
      sprintf(configfile,"tmp/%ss%d.out",name,n);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"w");
      fprintf(fdt,"%d\n",n);
     for( j1=0;j1<nj;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
               // printf("%d %d ", i1,j1);
		fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,*(h+(j1*ni+i1)),(u[j1*ni+i1]),(v[j1*ni+i1]),w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);
   
 //disp('writing data');
    //if finishsteering==1
      fprintf(fd,"%d\n",n);
      for( j1=0;j1<nj;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
          fprintf(fd,"%f %f %f ",u[j1*ni+i1],v[j1*ni+i1],h[j1*ni+i1]);
	   	//fprintf(fd,"%f ",h[j1*ni+i1]);
	}
        fprintf(fd,"\n");
      }


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
cufinish(&p,&w,&wnew,&b,&d_p,&d_w,&d_wnew,&d_b,&d_wmod, &d_dwn1,  &d_wd);
free(p);
free(sdir);
free(name);
free(outfile);
free(formfile);

//for the completed simulation
//int nsteps=nt;
fclose(fd);
}


