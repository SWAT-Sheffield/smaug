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
int i,j;


char *portfile=(char *)calloc(500,sizeof(char));
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));


real g  = 9.81;
real u0 = 0;                               
real v0 = 0;
real b0  = 0;                               
real h0 = 5030; 

int ngi=2;
int ngj=2;




//Domain definition
// Define the x domain
//adiab hydro
#ifdef ADIABHYDRO
int ni = 106;
ni=ni+2*ngi;
real xmax = 1.0;  
real dx = 0.55*xmax/(ni-4);
#endif

#ifndef ADIABHYDRO
//vac ozt
int ni;
ni=100;
ni=ni+2*ngi;
//ni=512;
//real xmax = 6.2831853;  
real xmax=1.0;
real dx = xmax/(ni-4);
#endif



// Define the y domain
//adiab hydro
#ifdef ADIABHYDRO
int nj = 106;
//nj=196;
nj=nj+2*ngj;
real ymax = 1.0;  
real dy = 0.55*ymax/(nj-4);
#endif

#ifndef ADIABHYDRO
//vac ozt
int nj = 196;
nj=100;
nj=nj+2*ngj;
//nj=512;
//real ymax = 6.2831853; 
real ymax = 1.0;   
real dy = ymax/(nj-4);    
//nj=41;
#endif
                    


real *x=(real *)calloc(ni,sizeof(real));
for(i=0;i<ni;i++)
		x[i]=i*dx;

real *y=(real *)calloc(nj,sizeof(real));
for(i=0;i<nj;i++)
		y[i]=i*dy;



int step=0;
//real tmax = 200;
real tmax = 0.2;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];


// Define time-domain
real dt;

#ifdef ADIABHYDRO
dt=0.0002985;  //ADIABHYDRO
#endif
//dt=0.15;

#ifndef ADIABHYDRO
dt=0.00065;
//dt=0.000139;
#endif
//dt=0.00009;
//dt=0.25;
//dt=0.00015125;
int nt=(int)((tmax)/dt);
//nt=3000;
nt=400;
nt=400;
//nt=2;
real *t=(real *)calloc(nt,sizeof(real));
printf("runsim 1%d \n",nt);
//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;

//real courant = wavespeed*dt/dx;



iome elist;
meta meta;






elist.server=(char *)calloc(500,sizeof(char));


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

	strcpy(elist.server,"localhost1");
	elist.port=80801;
	elist.id=0;

FILE *portf;

if(argc>1)
{ 
   sprintf(portfile,"%s0_port.txt",argv[1]) ;  
   //strcpy(name,argv[1]);
   portf=fopen(portfile,"r");
   fscanf(portf,"%d %s",&elist.port,elist.server);
   fclose(portf);

   printf("read file junk is %d %s\n",elist.port,elist.server);
}





//int cuprop(struct params **p, real **w, real **wnew, real **b,struct params **d_p, real **d_w, real **d_wnew, real **d_b, real **d_wmod, real **d_dwn1, real **d_dwn2, real **d_dwn3, real **d_dwn4, real **d_wd)




printf("rho %d mom1 %d mom2 %d\n",rho,mom1,mom2);

real *d_w;
real *d_wnew;

real *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

real *w,*wnew,*wd;
real *d_wtemp,*d_wtemp1,*d_wtemp2;
struct params *d_p;
struct params *p=(struct params *)malloc(sizeof(struct params));

struct state *d_state;
struct state *state=(struct state *)malloc(sizeof(struct state));

p->n[0]=ni;
p->n[1]=nj;
p->ng[0]=ngi;
p->ng[1]=ngj;
#ifdef ADIABHYDRO
p->npgp[0]=1;
p->npgp[1]=1;
#else
p->npgp[0]=1;
p->npgp[1]=1;
#endif

p->dt=dt;
p->dx[0]=dx;
p->dx[1]=dy;
//p->g=g;



/*constants used for adiabatic hydrodynamics*/
/*
p->gamma=2.0;
p->adiab=0.5;
*/
#ifdef ADIABHYDRO
p->gamma=2.0;
p->adiab=1.0;
#else

//ozt test
p->gamma=5.0/3.0;

//alfven test
//p->gamma=1.4;

#endif




p->mu=1.0;
p->eta=0.0;
p->g[0]=0.0;
p->g[1]=0.0;
p->g[2]=0.0;
//p->cmax=1.0;
p->cmax=0.02;

p->rkon=0.0;
p->sodifon=1.0;
p->moddton=0.0;
p->divbon=0.0;
p->divbfix=0.0;
p->hyperdifmom=0.0;
p->readini=0;
p->cfgsavefrequency=1;


p->xmax[0]=xmax;

p->xmax[1]=ymax;
p->nt=nt;
p->tmax=tmax;

p->steeringenabled=steeringenabled;
p->finishsteering=finishsteering;

p->maxviscoef=0;
//p->chyp=0.0;       
//p->chyp=0.00000;
p->chyp3=0.00000;
p->mnthreads=1;

for(i=0;i<NVAR;i++)
  p->chyp[i]=0.0;

p->chyp[rho]=0.02;
p->chyp[energy]=0.02;
p->chyp[b1]=0.02;
p->chyp[b2]=0.02;
p->chyp[mom1]=0.4;
p->chyp[mom2]=0.4;

p->chyp[rho]=0.02;
p->chyp[mom1]=0.4;
p->chyp[mom2]=0.4;
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
	//NewSimulation(metadata.name,'test1.xsl',elist);

// Build empty u, v, b matrices
// Define h
printf("allocating w and wnew\n");
 w=(real *)calloc(ni*nj*NVAR,sizeof(real ));
wd=(real *)calloc(ni*nj*NDERV,sizeof(real ));
 wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));
char *cfgfile;
if((p->readini)==0)
 initconfig(p, &meta, w);
else
 readconfig(cfgfile,*p,meta,w);
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(ni)*(nj)*rho;
  u=w+(ni)*(nj)*mom1;
  v=w+(ni)*(nj)*mom2;

printf("about to call cuinit\n");

cuinit(&p,&w,&wnew,&state,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);


printf("here in runsim\n");




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

//FILE *fd=fopen(outfile,"w");
//if steeringenabled==1
 printf("\n %s %s here in runsim4 %s\n",sdir,name,outfile); 
  //mkdir('tmp');
  gendxgen(sdir,name,nt,ni,nj);
printf("here in runsim5\n");

//sprintf(formfile,"%s/form%s.out",sdir,name);
//FILE *fdform=fopen(formfile,"w");
//  fprintf(fdform, "%d %d %d\n",nt-1, ni, nj);
//fclose(fdform);

//createlog(meta.log_file);


// Employ Lax
//disp(length[t));


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
for( n=0;n<nt;n++)
//for( n=0;n<1;n++)
{

    if((n%(p->cfgsavefrequency))==0)
    {
      //writeconfig(name,n,*p, meta , w);
      writevtkconfig(name,n,*p, meta , w);
    }
   order=0;
   t1=second();
   //cupredictor(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order);
   //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd);

   printf("cmax is %f old dt %f new dt %f\n",p->cmax,p->dt,0.68*((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax)));
   if(p->moddton==1)
   {
      //if((p->cmax)>10)
      //         p->dt=0.68*(((p->dx)+(p->dy))/(2.0*(p->cmax))); 
      if((p->dt)>(((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax))))
               p->dt=(p->dt)/2;
       else if(2.0*(p->dt)<(((p->dx[0])+(p->dx[1]))/(2.0*(p->cmax))))
               p->dt=2.0*p->dt; 
   }
   printf("new dt %f\n",p->dt);

if((p->rkon)==0)
{
  ordero=0;
 
 cucomputedervfields(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero);
 order=1; 
   //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order);
   //for(int f=rho; f<=mom3; f++)

 for(int dir=0;dir<2; dir++)
 {
 // int dir=0;
  for(int f=rho; f<=mom2; f++)
  {
      cucentdiff1(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
      cucentdiff1a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

  }

#ifndef ADIABHYDRO
   for(int f=energy; f<=b2; f++)
   {
     cucentdiff2(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
     cucentdiff2a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);

   }
#endif
  }

   //cuderivsource(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
   if(p->divbon==1)
	       cudivb(&p,&w,&state,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, &d_state,order,ordero,p->dt);
   if(p->hyperdifmom==1)
   {
    dt=(p->dt);
    for(int dim=0; dim<=1; dim++)
     {
       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);

       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);


       cuhyperdifrhosource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim);
       cuhyperdifrhosource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);

     }

     for(int dim=0; dim<=1; dim++)
     {
       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);

      cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);

       cuhyperdifesource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       cuhyperdifesource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       cuhyperdifesource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);

     }

for(int dim=0; dim<=1; dim++)
       for(int f=0; f<=1; f++)
           	                 
	     {
               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);

               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
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
                    cuhyperdifmomsource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
                   }
                   else
                   {
                    cuhyperdifmomsourcene1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                   cuhyperdifmomsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                   cuhyperdifmomsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
		    }
                }
             }
            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=1; dim++)
	     for(int f=0; f<=1; f++)            
	     {
               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);


               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);

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
                     cuhyperdifbsource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     cuhyperdifbsource4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   }
                   else
                   {
                     cuhyperdifbsourcene1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                    cuhyperdifbsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    cuhyperdifbsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    cuhyperdifbsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    }

                }
              } 


   }
   //cuadvance(&p,&w,&wnew,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order);
  // cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,ordero);
}

   if((p->rkon)==1)
     for(order=0; order<4; order++) 
   {	   
           ordero=order+1;
           dt=(p->dt)/2.0;
           //if(order==1)
           //{
           //   dt=(p->dt);
           //   //ordero=0;
           //}
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


           cucomputedervfields(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero);
	   //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order);
           //for(int f=rho; f<=mom3; f++)
 for(int dir=0;dir<2; dir++)
 {
           for(int f=rho; f<=mom2; f++)
           {
	       cucentdiff1(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);
	       cucentdiff1a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);

            }

#ifndef ADIABHYDRO
           for(int f=energy; f<=b2; f++)
           {
	       cucentdiff2(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
	       cucentdiff2a(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

           }
#endif

}
	   //cuderivsource(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
	   if(p->divbon==1)
	       cudivb(&p,&w,&state,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, &d_state,order,ordero,p->dt);
           if(p->hyperdifmom==1)
           {
	     for(int dim=0; dim<=1; dim++)
	     {
	       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
	       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
	       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);
	       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,0);

	       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
	       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
	       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);
	       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim,1);

	      // cuhyperdifviscmax(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,rho,dim);

	       cuhyperdifrhosource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim);
	       cuhyperdifrhosource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);

	     }

     for(int dim=0; dim<=1; dim++)
     {
       cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);

      cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);
       cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,1);

       cuhyperdifesource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       cuhyperdifesource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim);
       cuhyperdifesource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);

     }
             
for(int dim=0; dim<=1; dim++)
       for(int f=0; f<=1; f++)
           	                 
	     {
               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,0);

               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim,1);
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
                    cuhyperdifmomsource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
 
                   }
                   else
                   {
                    cuhyperdifmomsourcene1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                    cuhyperdifmomsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);


		    }
                }
             }

            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=1; dim++)
	     for(int f=0; f<=1; f++)            
	     {
               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,0);


               cuhyperdifvisc1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);
               cuhyperdifvisc4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim,1);

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
                     cuhyperdifbsource1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsource2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsource3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                     cuhyperdifbsource4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   }
                   else
                   {
                     cuhyperdifbsourcene1(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsourcene2(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                     cuhyperdifbsourcene3(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb);
                    cuhyperdifbsourcene4(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    cuhyperdifbsourcene5(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    cuhyperdifbsourcene6(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                    }

                }
              } 

           }
           //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,ordero);
           cuadvance(&p,&w,&wnew,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order);
           cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,orderb);
	   

   }

   cuupdate(&p,&w,&wd,&state,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, &d_state);
   printf("nummaxthreads %d\n",p->mnthreads);

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
    
    for( j1=ngj;j1<nj-ngj;j1++)
        for( i1=ngi;i1<ni-ngi;i1++)
{

;//w[j1*ni+i1+(ni*nj*b1)]=wd[j1*ni+i1+(ni*nj*(hdnur))];
;//w[j1*ni+i1+(ni*nj*b2)]=wd[j1*ni+i1+(ni*nj*(hdnul))];

}
           
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
cufinish(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd);
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
	writesimulation_(elist.id,newsimfile,elist.port,elist.server);


	return 0;
}

