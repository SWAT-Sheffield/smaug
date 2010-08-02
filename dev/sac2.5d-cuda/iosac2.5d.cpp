// IoSimulationTest.cpp : Defines the entry point for the console application.
//

/*
IOME LICENSE
IOME Version 1.1.1

IOME Development  Tools
Copyright (C) 2001-2004, Michael Kenneth Griffiths, All Rights Reserved.

--------------------------------------------------------------------------------
IOME public license.

The contents of this file are subject to the IOME Public License Version 1.3
(the "License"); you may not use this file except in compliance with the
License. You may obtain a copy of the License at
http://81.174.178.112/iome/licensing/iomelicense.html
Software distributed under the License is distributed on an "AS IS" basis,
WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
for the specific language governing rights and limitations under the License.

The Initial Developer of the Original Code is Michael Kenneth Griffiths.
Copyright (C) 2000-2004 Michael Kenneth Griffiths. All Rights Reserved.
--------------------------------------------------------------------------------
GPL license.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA

Author contact information:
mikeg@photon0.freeserve.co.uk
--------------------------------------------------------------------------------
*/
#include "iosac2.5d.h"
#include "step.h"
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

int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int it=0; //test integer to be returned 

//getintparam_( int elist.id,char *sname,int *iv,  int elist.port, char *selist.server );
//int elist.id=0;
//int elist.port=8080;

int i1,i2,i3,j1;
int i,j;
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));


double g  = 9.81;
double u0 = 0;                               
double v0 = 0;
double b0  = 0;                               
double h0 = 5030; 


//Domain definition
// Define the x domain
int ni = 150; 
//ni=41;
float xmax = 3.0;                      
float dx = xmax/(ni-1);




// Define the y domain
int nj = 150;  
//nj=41;
float ymax = 3.0;                      
float dy = ymax/(nj-1);

float *x=(float *)calloc(ni,sizeof(float));
for(i=0;i<ni;i++)
		x[i]=i*dx;

float *y=(float *)calloc(nj,sizeof(float));
for(i=0;i<nj;i++)
		y[i]=i*dy;



int step=0;
double tmax = 200;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];

// Define the wavespeed
float wavespeed = u0 + sqrt(g*(h0 - b0));

// Define time-domain
float dt = 0.68*dx/wavespeed;
//dt=0.015985;
dt=0.001;
int nt=(int)((tmax-1)/dt);

float *t=(float *)calloc(nt,sizeof(float));
printf("runsim 1\n");
//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;

float courant = wavespeed*dt/dx;



iome elist;
meta meta;









meta.directory=(char *)calloc(500,sizeof(char));
meta.author=(char *)calloc(500,sizeof(char));
meta.sdate=(char *)calloc(500,sizeof(char));
meta.platform=(char *)calloc(500,sizeof(char));
meta.desc=(char *)calloc(500,sizeof(char));
meta.name=(char *)calloc(500,sizeof(char));
meta.ini_file=(char *)calloc(500,sizeof(char));
meta.log_file=(char *)calloc(500,sizeof(char));
meta.out_file=(char *)calloc(500,sizeof(char));

strcpy(meta.directory,"out");
strcpy(meta.author,"MikeG");
strcpy(meta.sdate,"Nov 2009");
strcpy(meta.platform,"swat");
strcpy(meta.desc,"A simple test of SAAS");
strcpy(meta.name,"test1");
strcpy(meta.ini_file,"test1.ini");
strcpy(meta.log_file,"test1.log");
strcpy(meta.out_file,"test1.out");
//meta.directory="out";
//meta.author="MikeG";
//meta.sdate="Nov 2009";
//meta.platform="felix";
//meta.desc="A simple test of SAAS";
//meta.name="tsteer1";



elist.server="localhost";
elist.port=8080;
elist.id=0;



//int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_dwn2, float **d_dwn3, float **d_dwn4, float **d_wd)




printf("rho %d mom1 %d mom2 %d\n",rho,mom1,mom2);

float *d_w;
float *d_wnew;

float *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

float *w,*wnew;
struct params *d_p;
struct params *p=(struct params *)malloc(sizeof(struct params));

struct state *d_state;
struct state *state=(struct state *)malloc(sizeof(struct state));

p->ni=ni;
p->nj=nj;
p->dt=dt;
p->dx=dx;
p->dy=dy;
p->g=g;
p->gamma=1.4;
p->mu=1.0;
p->eta=0.0;
p->g1=0.0;
p->g2=0.0;
p->g3=0.0;
p->cmax=1.0;

p->rkon=0.0;
p->sodifon=0.0;
p->moddton=0.0;
p->divbon=0.0;



p->xmax=xmax;

p->ymax=ymax;
p->nt=nt;
p->tmax=tmax;

p->steeringenabled=steeringenabled;
p->finishsteering=finishsteering;

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
	//stacksize(268435454);
	//open the file generated
	//sprintf(elist.portfile,"%s0_elist.port.txt",meta.name);
	//FILE *fd=fopen(elist.portfile,"r");
	//int elist.portelist.id;
	//fscanf(fd,"%d",&elist.portelist.id);
	//fclose(fd);
	//elist.elist.port=elist.portelist.id;
        if(argc>1)
        {
          readsim(p,&meta,argv[1],elist);
          if((p->readini)!=0)
             readconfig(meta.ini_file,*p,meta,w);
        }
        else
	  createsim(*p,meta,simfile,elist);

	sprintf(simfile,"%s.xml",meta.name);
        sprintf(newsimfile,"%s_update.xml",meta.name);
	//NewSimulation(metadata.name,'test1.xsl',elist);

// Build empty u, v, b matrices
// Define h
if((p->readini)==0)
{
 w=(float *)calloc(ni*nj*8,sizeof(float ));
 wnew=(float *)calloc(ni*nj*8,sizeof(float ));

 initconfig(p, &meta, w);
}
  float *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(ni)*(nj)*rho;
  u=w+(ni)*(nj)*mom1;
  v=w+(ni)*(nj)*mom2;



cuinit(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state);


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
getmetadata_(elist.id,"directory",&sdir,elist.port,elist.server);
//sdir=metadata.directory

//name=metadata.name;

getmetadata_(elist.id,"name",&name,elist.port,elist.server);
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

createlog(meta.log_file);


// Employ Lax
//disp(length[t));


//while(finishsteering == 0)
//{
 
  //  if( steeringenabled==0)
  //    finishsteering=1;
 int n;  
 nt=24; 
double t1,t2,ttot;
int order=0;
ttot=0;
float time=0.0;
for( n=0;n<nt;n++)
{
  
   t1=second();
   cupredictor(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order);

   if(order==0 && p->moddton==1)
       p->dt=((p->dx)+(p->dy))/(2.0*(p->cmax));
   cuderivcurrent(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order);
   cuderivsource(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd);

   cuadvance(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd);

   cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd);
   cuupdate(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd, &d_state);

   if(p->divbon==1)
       cudivb(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd, &d_state);

   t2=second()-t1;
   ttot+=t2;
   printf("step %d total time %f\n",n,ttot);

   state->it=n;
   state->t=time+(p->dt);
   time=state->t;
   state->dt=p->dt;

   appendlog(meta.log_file,*p, *state);

 
    
    getintparam_(&elist.id,"steeringenabled",&steeringenabled,&elist.port,elist.server);
    if(steeringenabled==1)
    {
      //disp('getting updatea params');
      //for steering get the modified control params
      double dg;
      getintparam_(&elist.id,"finishsteering",&finishsteering,&elist.port,elist.server);//source y location  
        // Constants
      getdoubleparam_(elist.id,"g",&dg,elist.port,elist.server);

      g=dg;
     
    }
    
    writeconfig(name,n,*p, meta , w);
    
      //save file containing current data
     // sprintf(configfile,"tmp/%ss%d.out",name,n);
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
cufinish(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd);
free(p);
free(sdir);
free(name);
free(outfile);
free(formfile);

//for the completed simulation
//int nsteps=nt;
fclose(fd);






	//[consts,domain,source]=loadsim('test1_16_02_09.xml',elist);
	//chdir(metadata.directory);
        //readsimulation_(elist.elist.id,simfile,elist.elist.port,elist.elist.server);
	//runsim(consts,dom,src,meta,simfile,elist);
	writesimulation_(elist.id,newsimfile,elist.port,elist.server);


	return 0;
}

