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
#include "IoTest.h"
void gendxgen(char *dir,char *jobname,int nsteps,int n1,int n2)
{
//generate dx general file 
//describing magnetic fields
//particle locations
  //file=out/jobform.out
  //grid = 1
  //format = ascii
  //interleaving = record
  //majority = row
  //field = nsteps, nx, ny
  //structure = scalar, scalar, scalar
  //type = int, int, int
  //dependency = positions, positions, positions
  //positions = regular, 0, 1

  //end
  char dxformgenfile[300];
  char dxgenfile[300];
  char basedxgenfile[300];
  sprintf(dxformgenfile,"dx/%s_form.general",jobname);
  FILE *fdform=fopen(dxformgenfile,"w");
    fprintf(fdform, "file=%s/form%s.out",dir,jobname);
    fprintf(fdform,"grid=1\n");
    fprintf(fdform,"format = ascii \n interleaving = record \n majority = row \n");
    fprintf(fdform, "field = nsteps, nx, ny \n structure = scalar, scalar, scalar \n type = int, int, int  \n dependency = positions, positions,positions  \n positions = regular, 0, 1 \n end \n ");
  fclose(fdform);    
//generate dx general file for this data set
  //file=out/job.out
  //grid 51 x 51
  //format = ascii
  //interleaving = field
  //majority = row
  //header = lines 1

  //series =  24 , 1, 1, separator=lines 1
  //field = field0, field1
  //structure = 2-vector, scalar
  //type = float, float
  //dependency = positions, positions
  //positions = regular,regular, 0, 1,0,1

  //end
  sprintf(dxgenfile,"dx/%s.general",jobname);

  fdform=fopen(dxgenfile,"w");
    fprintf(fdform, "file=%s/%s.out",dir,jobname);
    fprintf(fdform,"grid %d X %d\n",n1,n2);
    fprintf(fdform,"format = ascii \n interleaving = field \n majority = row \n header = lines 1 \n");
    fprintf(fdform, "series =  %d  , 1, 1, separator=lines 1\n",nsteps-1);
    fprintf(fdform, "field = field0, field1 \n structure = 2-vector, scalar \n type = float, float  \n dependency = positions, positions  \n positions = regular,regular, 0, 1,0,1 \n end \n ");
  fclose(fdform);


 sprintf(basedxgenfile,"dx/base%s.general",jobname);

  FILE *fdbform=fopen(basedxgenfile,"w");
   // mfprintf(fdbform, 'file=%s\n', directory+'/'+jobname+'.out');
    fprintf(fdbform,"grid %d X %d\n",n1,n2);
    fprintf(fdbform,"format = ascii \n interleaving = field \n majority = row \n header = lines 1 \n");
    fprintf(fdbform, "series =  1  , 1, 1, separator=lines 1\n");
    fprintf(fdbform, "field = field0, field1 \n structure = 2-vector, scalar \n type = float, float  \n dependency = positions, positions  \n positions = regular,regular, 0, 1,0,1 \n end \n ");
  fclose(fdbform);



//endfunction
}

void runsim(constants k, domain dom,source src, meta metadata,char *simname, iome el)
{

char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));
//elist=list();  parameter used by iome to contain port and server address
//elist=list(); 

double sf=src.freq;//source frequency
double sa=src.amp;//source amplitude
double sx=src.xloc;//source x location
double sy=src.yloc;//source y location

// Constants
double g  = k.g;
double u0 = k.u0;                               
double v0 = k.v0;
double b0  = k.b0;                               
double h0 = k.h0; 

//Domain definition
// Define the x domain
//ni = 151; 
int ni=dom.ni;
double xmax = dom.xmax;                      
double dx = xmax/(ni-1);
double *x=(double *)calloc(ni,sizeof(double));
for(i=0;i<ni;i++)
		x[i]=i*dx;
int i1,j1;
// Define the y domain
//nj = 151;  
int nj=dom.nj;
double ymax = dom.ymax;                      
double dy = ymax/(nj-1);
double *y=(double *)calloc(nj,sizeof(double));
for(i=0;i<nj;i++)
		y[i]=i*dy;


double tmax = dom.tmax;
int steeringenabled=dom.steeringenabled;
int finishsteering=dom.finishsteering;


// Define the wavespeed
double wavespeed = u0 + sqrt(g*(h0 - b0));

// Define time-domain
double dt = 0.68*dx/wavespeed;
int nt=(int)((tmax-1)/dt);

double *t=(double *)calloc(nt,sizeof(double));

//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;
dom.nt=nt;

double courant = wavespeed*dt/dx;

// Build empty u, v, b matrices
double *u=(double *)calloc(2*ni*nj,sizeof(double));
double *v=(double *)calloc(2*ni*nj,sizeof(double));
double *b=(double *)calloc(ni*nj,sizeof(double));


// Define h
double *h=(double *)calloc(2*ni*nj,sizeof(double));


for(i=0; i<2*ni*nj; i++)
{
	u[i]=0;
	v[i]=0;
	h[i]=0;
}

for(i=0; i<ni*nj; i++)
		b[i]=0;


for(i=0;i<ni*nj;i++) 
		h[i] = 5000; 

int nli = (45000/100000*(ni-1));
int nui = floor(55000/100000*(ni-1));
int nlj = (45000/100000*(ni-1));
int nuj = floor(55000/100000*(ni-1));                            
//h[(45000/100000*(ni-1)+1):floor(55000/100000*(ni-1)+1),(45000/100000*(nj-1)+1):floor(55000/100000*(nj-1)+1),1) = 5030;
int in,ind;
for(i=nli;i<nui;i++)
	for(j=nlj;j<nuj;j++)
	{
	   in=i+j*nj;
           h[in]=5030;	
	}


//Define b
for(i=0;i<ni;i++)
{
    if(x[i] > 20001){

      for(j=0;j<ni;j++)
      {        
	ind=j+i*nj;
        b[ind]=0;
        //b(:,i) = 0;
       }
    }
    else if( x[i] < 20000)
    {
        //b(:,i) = 5000/20000*(20000-x(i));
        for(j=0;j<ni;j++)
        {        
	ind=j+i*nj;
        b[ind]=5000/20000*(20000-x[i]);
        //b(:,i) = 0;
         }
    }
}

//For a steerable simulation generate and save a dxformfile that saves a single data step
//used for the steering dx module

getmetadata_(el.id,"directory",sdir,el.port,el.server);
//sdir=metadata.directory

//name=metadata.name;

getmetadata_(el.id,"name",name,el.port,el.server);
//disp(sdir,name)

sprintf(outfile,"%s/%s.out",sdir,name);
FILE *fd=fopen(outfile,"w");
//if steeringenabled==1
  
  //mkdir('tmp');
  gendxgen("tmp",name,1,ni,nj);


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
for( n=1;n<nt;n++)
{
  //disp('step n');
  //disp(n);

   //disp('n');
   //disp(n);
   
    //for i=2:(ni-1)
       // for j=2:(nj-1)
       // ind=i+j*ni;
        for(ind=(2+2*ni);ind<((ni-1)+(nj-i)*ni);ind++)
        {   


            j=ind/ni;		            
            i=ind-ni*(int)(ind/ni);		
		//u[i,j,2] = ((u[i+1,j,1) + u[i-1,j,1) + u[i,j+1,1) + u[i,j-1,1))/4)...
            //    - 0.5*(dt/dx)*((u[i+1,j,1)^2)/2 - (u[i-1,j,1)^2)/2)...
            //    - 0.5*(dt/dy)*(v[i,j,1))*(u[i,j+1,1) - u[i,j-1,1)) - 0.5*g*(dt/dx)*(h[i+1,j,1)-h[i-1,j,1));
		u[i+j*ni+ni*nj] = ((u[i+1+j*ni] + u[i-1+j*ni] + u[i+(j+1)*ni] + u[i+(j-1)*ni])/4)- 0.5*(dt/dx)*(((u[(i+1)+ni*j])*u[(i+1)+ni*j])/2) - ((u[i-1+j*ni]*u[i-1+j*ni])/2) - 0.5*(dt/dy)*(v[i+j*ni])*(u[i+(j+1)*ni] - u[i+(j-1)*ni]) - 0.5*g*(dt/dx)*(h[i+1+j*ni]-h[i-1+j*ni]);

            
            //v[i,j,2) = ((v[i+1,j,1) + v[i-1,j,1) + v[i,j+1,1) + v[i,j-1,1))/4)- 0.5*(dt/dy)*((v[i,j+1,1)^2)/2 - (v[i,j+1,1)^2)/2)- 0.5*(dt/dx)*(u[i,j,1))*(v[i+1,j,1) - v[i-1,j,1)) - 0.5*g*(dt/dy)*(h[i,j+1,1)-h[i,j-1,1));

	v[i+j*ni+ni*nj] = ((v[i+1+j*ni] + v[i-1+j*ni] + v[i+(j+1)*ni] + v[i+(j-1)*ni])/4)- 0.5*(dt/dy)*((v[i+ni*(j+1)]*v[(i)+ni*(j+1)])/2 - (v[i+(j-1)*ni]*v[i+(j-1)*ni])/2) - 0.5*(dt/dx)*(u[i+j*ni])*(v[i+1+j*ni] - v[i-1+j*ni]) - 0.5*g*(dt/dy)*(h[i+(j+1)*ni]-h[i+(j-1)*ni]);		
                        
            
              //h[i,j,2) = sa*sin(n*sf)+((h[i+1,j,1) + h[i-1,j,1) + h[i,j+1,1) + h[i,j-1,1))/4)- 0.5*(dt/dx)*(u[i,j,1))*((h[i+1,j,1)-b(i+1,j)) - (h[i-1,j,1)-b(i-1,j))) - 0.5*(dt/dy)*(v[i,j,1))*((h[i,j+1,1)-b(i,j+1)) - (h[i,j-1,1)-b(i,j-1))) - 0.5*(dt/dx)*(h[i,j,1)-b(i,j))*(u[i+1,j,1)- u[i-1,j,1))- 0.5*(dt/dy)*(h[i,j,1)-b(i,j))*(v[i,j+1,1) - v[i,j-1,1));
        h[i+j*ni+ni*nj] = ((h[i+1+j*ni] + h[i-1+j*ni] + h[i+(j+1)*ni] + h[i+(j-1)*ni])/4)- 0.5*(dt/dx)*(u[i+j*ni])*((h[i+1+j*ni]-b[i+1+j*ni]) - (h[i-1+j*ni]-b[i-1+j*ni])) - 0.5*(dt/dy)*(v[i+j*ni])*((h[i+(j+1)*ni]-b[i+(j+1)*ni]) - (h[i+(j-1)*ni]-b[i+(j-1)*ni])) - 0.5*(dt/dx)*(h[i+j*ni]-b[i+j*ni])*(u[i+1+j*ni]- u[i-1+j*ni])- 0.5*(dt/dy)*(h[i+j*ni]-b[i+j*ni])*(v[i+(j+1)*ni] - v[i+(j-1)*ni]);     
           }     


 

    // Define Boundary Conditions
    for(j=0;j<nj;j++)
    {
	u[1+j*ni+ni*nj] = 2.5*u[2+j*ni+ni*nj] - 2*u[3+j*ni+ni*nj] + 0.5*u[4+j*ni+ni*nj];
    	//u[1,:,2) = 2.5*u[2,:,2) - 2*u[3,:,2) + 0.5*u[4,:,2);
        //u[length[x),:,2) = 2.5*u[ni-1,:,2) - 2*u[ni-2,:,2) + 0.5*u[ni-3,:,2);
       u[ni-1+j*ni+ni*nj] = 2.5*u[ni-1+j*ni+ni*nj] - 2*u[ni-2+ni*j+ni*nj] + 0.5*u[ni-3+j*ni+ni*nj];


	//v[1,:,2) = 2.5*v[2,:,2) - 2*v[3,:,2) + 0.5*v[4,:,2);
	v[1+j*ni+ni*nj] = 2.5*v[2+j*ni+ni*nj] - 2*v[3+j*ni+ni*nj] + 0.5*v[4+j*ni+ni*nj];

       //v[length[x),:,2) = 2.5*v[ni-1,:,2) - 2*v[ni-2,:,2) + 0.5*v[ni-3,:,2);
	v[ni-1+j*ni+ni*nj] = 2.5*v[ni-1+j*ni+ni*nj] - 2*v[ni-2+ni*j+ni*nj] + 0.5*v[ni-3+j*ni+ni*nj];


    	//h[1,:,2) = 2.5*h[2,:,2) - 2*h[3,:,2) + 0.5*h[4,:,2);
	h[1+j*ni+ni*nj] = 2.5*h[2+j*ni+ni*nj] - 2*h[3+j*ni+ni*nj] + 0.5*h[4+j*ni+ni*nj];

    	//h[length[x),:,2) = 2.5*h[ni-1,:,2) - 2*h[ni-2,:,2) + 0.5*h[ni-3,:,2);
        h[ni-1+j*ni+ni*nj] = 2.5*h[ni-1+j*ni+ni*nj] - 2*h[ni-2+ni*j+ni*nj] + 0.5*h[ni-3+j*ni+ni*nj];

    }


    for(i=0;i<ni;i++)
    {
	
    //u[:,1,2) = 2.5*u[:,2,2) - 2*u[:,3,2) + 0.5*u[:,4,2);
	u[i+ni+ni*nj] = 2.5*u[i+2*ni+ni*nj] - 2*u[i+3*ni+ni*nj] + 0.5*u[i+4*ni+ni*nj];

    //u[:,length[y),2) = 2.5*u[:,nj-1,2) - 2*u[:,nj-2,2) + 0.5*u[:,nj-3,2);
       u[i+(nj-1)*ni+ni*nj] = 2.5*u[i+(nj-2)*ni+ni*nj] - 2*u[i+(nj-3)*ni+ni*nj] + 0.5*u[i+(nj-4)*ni+ni*nj];

    //v[:,1,2) = 2.5*v[:,2,2) - 2*v[:,3,2) + 0.5*v[:,4,2);
     v[i+ni+ni*nj] = 2.5*v[i+2*ni+ni*nj] - 2*v[i+3*ni+ni*nj] + 0.5*v[i+4*ni+ni*nj];
    
     //v[:,length[y),2) = 2.5*v[:,nj-1,2) - 2*v[:,nj-2,2) + 0.5*v[:,nj-3,2);
       v[i+(nj-1)*ni+ni*nj] = 2.5*v[i+(nj-2)*ni+ni*nj] - 2*v[i+(nj-3)*ni+ni*nj] + 0.5*v[i+(nj-4)*ni+ni*nj];


    //h[:,1,2) = 2.5*h[:,2,2) - 2*h[:,3,2) + 0.5*h[:,4,2);
     h[i+ni+ni*nj] = 2.5*h[i+2*ni+ni*nj] - 2*h[i+3*ni+ni*nj] + 0.5*h[i+4*ni+ni*nj];

    //h[:,length[y),2) = 2.5*h[:,nj-1,2) - 2*h[:,nj-2,2) + 0.5*h[:,nj-3,2);
       h[i+(nj-1)*ni+ni*nj] = 2.5*h[i+(nj-2)*ni+ni*nj] - 2*h[i+(nj-3)*ni+ni*nj] + 0.5*h[i+(nj-4)*ni+ni*nj];


    }

        for(ind=0;ind<(ni*nj);ind++)
        {   
            j=ind/ni;		            
            i=ind-ni*(int)(ind/ni);

	    h[i+j*ni]=h[i+j*ni+ni*nj];
	    u[i+j*ni]=u[i+j*ni+ni*nj];
	    v[i+j*ni]=v[i+j*ni+ni*nj];
	}
    
    getintparam_(&el.id,"steeringenabled",&steeringenabled,&el.port,el.server);
    if(steeringenabled==1)
    {
      //disp('getting updatea params');
      //for steering get the modified control params
      getdoubleparam_(el.id,"frequency",&sf,el.port,el.server);//source frequency
      getdoubleparam_(el.id,"amplitude",&sa,el.port,el.server);//source amplitude
      getdoubleparam_(el.id,"sx",&sx,el.port,el.server);//source x location
      getdoubleparam_(el.id,"sy",&sy,el.port,el.server);//source y location
      getintparam_(&el.id,"finishsteering",&finishsteering,&el.port,el.server);//source y location  
        // Constants
      getdoubleparam_(el.id,"g",&g,el.port,el.server);
    }
    
    
      //save file containing current data
      sprintf(configfile,"tmp/%ss%d.out",name,n);
      FILE *fdt=fopen(configfile,"w");
      fprintf(fdt,"%d\n",n);
      for(i1=0;i1<nj;i1++)
      {
        for(j1=0;j1<ni;j1++)
	{
          fprintf(fdt,"%f %f %f",u[i1+j1*ni],v[i1+j1*ni],h[i1+j1*ni]);
        }        
        fprintf(fdt,"\n");
      }
      fclose(fdt);
    }//end of testep
    
    //end
  
    //disp('writing data');
    //if finishsteering==1
      fprintf(fd,"%d\n",n);
      for( i1=0;i1<ni;i1++)
      {
        for( j1=0;j1<nj;j1++)
	{
          fprintf(fd,"%f %f %f",u[i1+j1*ni],v[i1+j1*ni],h[i1+j1*ni]);
	}
        fprintf(fd,"\n");
      }
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
free(sdir);
free(name);
free(outfile);
free(formfile);

//for the completed simulation
//int nsteps=nt;
fclose(fd);
}


int main(int argc, char* argv[])
{
	int itype=-1;
	int status=1;
	int it=0; //test integer to be returned 
    //getintparam_( int id,char *sname,int *iv,  int port, char *sserver );
    int id=0;
    int port=8080;



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









    getintparam_( &id,"i1",&it,  &port, "localhost" );	
	printf("Get integer %d\n",it);
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
	//createsim(consts,domain,source,metadata,simfile,elist);

	//[consts,domain,source]=loadsim('test1_16_02_09.xml',elist);
	//chdir(metadata.directory);
        readsimulation_(elist.id,simfile,elist.port,elist.server);
	runsim(consts,dom,src,meta,simfile,elist);
	writesimulation_(elist.id,simfile,elist.port,elist.server);


	return 0;
}

