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


int main(int argc, char* argv[])
{
int i;	
int itype=-1;
	int status=1;
	int it=0; //test integer to be returned 
    //getintparam_( int id,char *sname,int *iv,  int port, char *sserver );
    int id=0;
    int port=8080;



double g  = 9.81;
double u0 = 0;                               
double v0 = 0;
double b0  = 0;                               
double h0 = 5030; 


//Domain definition
// Define the x domain
int ni = 151; 
//ni=41;
double xmax = 100000;                      
double dx = xmax/(ni-1);
double *x;



// Define the y domain
int nj = 151;  
//nj=41;
double ymax = 100000;                      
double dy = ymax/(nj-1);
double *y;

int nt=10;
int step=0;
double tmax = 200;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];





double sf=10;//source frequency
double sa=5;//source amplitude
double sx=20;//source x location
double sy=30;//source y location


source src;
constants consts;
domain dom;
iome elist;
meta meta;

src.freq=sf;
src.amp=sa;
src.xloc=sx;
src.yloc=sy;




consts.g=g;
consts.u0=u0;
consts.v0=v0;
consts.b0=b0;
consts.h0=h0;

dom.ni=ni;
dom.xmax=xmax;
dom.nj=nj;
dom.ymax=ymax;
dom.nt=nt;
dom.tmax=tmax;
dom.step=step;
dom.steeringenabled=steeringenabled;
dom.finishsteering=finishsteering;


meta.directory="out";
meta.author="MikeG";
meta.sdate="Nov 2009";
meta.platform="felix";
meta.desc="A simple test of SAAS";
meta.name="tsteer1";



elist.server="localhost";
elist.port=8080;
elist.id=0;









 //   getintparam_( &id,"i1",&it,  &port, "localhost" );	
//	printf("Get integer %d\n",it);
    //Set input filename as first arg
	//if NULL use defaults
	char *method=NULL;
	//CIoSimulation *TestSimulation;
      x=(double *)calloc(ni,sizeof(double));
      for(i=0;i<ni;i++)
		x[i]=i*dx;


      y=(double *)calloc(nj,sizeof(double));
      for(i=0;i<nj;i++)
		y[i]=i*dy;

	//this should be executed by the iome start up application
	//exec('ioshallowwater.sce');

	//this application is started using the io  start scilab application
	//exec('paramssteeringtest1.sce');
	//stacksize('max');
	//stacksize(268435454);
	//open the file generated
	sprintf(portfile,"%s0_port.txt",meta.name);
	FILE *fd=fopen(portfile,"r");
	int portid;
	fscanf(fd,"%d",&portid);
	fclose(fd);

	elist.port=portid;

	sprintf(simfile,"%s.xml",meta.name);
	//NewSimulation(metadata.name,'test1.xsl',elist);
	createsim(consts,dom,src,meta,simfile,elist);

	//[consts,domain,source]=loadsim('test1_16_02_09.xml',elist);
	//chdir(metadata.directory);
        //readsimulation_(elist.id,simfile,elist.port,elist.server);
	runsim(consts,dom,src,meta,simfile,elist);
	writesimulation_(elist.id,simfile,elist.port,elist.server);


	return 0;
}

