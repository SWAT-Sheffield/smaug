#include "readwrite.h"


int createlog(char *logfile)
{
	int status=0;

      FILE *fdt=0;

      fdt=fopen(logfile,"a+");
      fprintf(fdt,"it   t   dt    rho m1 m2 m3 e bx by bz\n");
      fclose(fdt);	

	return status;
}

int appendlog(char *logfile, params p, state s)
{
  int status=0;
      FILE *fdt=0;

      fdt=fopen(logfile,"a+");
      fprintf(fdt,"%d %f %f %f %f %f %f %f %f %f %f %f\n",s.it,
               s.t,s.dt,s.rho,s.m1,s.m2,s.m3,s.e,s.b1,s.b2,s.b3);
      fclose(fdt);
  return status;
}

int writeconfig(char *name,int n,params p, meta md, real *w)
{
  int status=0;
  int i1,j1;
  int ni,nj;
  char configfile[300];


  ni=p.n[0];
  nj=p.n[1];



      //save file containing current data
      sprintf(configfile,"tmp/%ss%d.out",name,n);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"w");
      fprintf(fdt,"%d\n",n);
     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
      {
        for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
	{
               // printf("%d %d ", i1,j1);
 #ifdef ADIABHYDRO
		fprintf(fdt,"%d %d %f %f %f %f \n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)]);
#else
		fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
#endif
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


      //save file containing current data
      sprintf(configfile,"out/%s.out",name);
      printf("write out check dims %s %d %d \n",configfile ,ni,nj);
      fdt=fopen(configfile,"a+");
      fprintf(fdt,"%d\n",n);
     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
      {
        for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
	{
               // printf("%d %d ", i1,j1);
 #ifdef ADIABHYDRO
 fprintf(fdt,"%d %d %f %f %f %f \n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)]);
#else
fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
        #endif
		
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


  return status;
}



int writevacconfig(char *name,int n,params p, meta md, real *w, state st)
{
  int status=0;
  int i1,j1;
  int ni,nj;
  char configfile[300];
  char buffer[800];
  double dbuffer[10];

  ni=p.n[0];
  nj=p.n[1];

      //save file containing current data
      sprintf(configfile,"out/v%s.out",name);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"a+");

      fwrite(md.name,sizeof(char)*strlen(md.name),1,fdt);
      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      sprintf(buffer,"%d %f 3 4 8\n",st.it,st.t);
      fwrite(buffer,sizeof(char)*strlen(buffer),1,fdt);

      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      sprintf(buffer,"%d %d\n",ni,nj);
      fwrite(buffer,sizeof(char)*strlen(buffer),1,fdt);

      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      sprintf(buffer,"%f %f %f %f\n",p.eta,p.g[0],p.g[1],p.g[2]);
      fwrite(buffer,sizeof(char)*strlen(buffer),1,fdt);

      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)
      sprintf(buffer,"x y rho mx my mz e bx by bz gamma eta g1 g2 g3\n");
      fwrite(buffer,sizeof(char)*strlen(buffer),1,fdt);

       for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
	{
         
     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
      {

                dbuffer[0]=i1*p.dx[0];
                dbuffer[1]=j1*p.dx[1];
                dbuffer[2]=w[(j1*ni+i1)+(ni*nj*rho)];
                dbuffer[3]=w[(j1*ni+i1)+(ni*nj*mom1)];
                dbuffer[4]=w[(j1*ni+i1)+(ni*nj*mom2)];
                dbuffer[5]=w[(j1*ni+i1)+(ni*nj*mom3)];
                dbuffer[6]=w[(j1*ni+i1)+(ni*nj*energy)];
                dbuffer[7]=w[(j1*ni+i1)+(ni*nj*b1)];
                dbuffer[8]=w[(j1*ni+i1)+(ni*nj*b2)];
                dbuffer[9]=w[(j1*ni+i1)+(ni*nj*b3)];

                fwrite(dbuffer,10*sizeof(double),1,fdt);		

        }     
      }
      fclose(fdt);

  return status;
}




int writevtkconfig(char *name,int n,params p, meta md, real *w)
{
  int status=0;
  int i1,j1;
  int ni,nj;
  char configfile[300];
  char labels[4][4]={"rho","e","mom","b"};
  int is;
  ni=p.n[0];
  nj=p.n[1];



      //save file containing current data

      //scalar fields
//n+=10;
        #ifdef ADIABHYDRO
      for(int i=0; i<=3; i+=4)
      //for(int i=0; i<=4; i+=3)
       #else
      for(int i=0,is=0; i<=4; i+=4,is+=3)
      //for(int i=0; i<=4; i+=3)
        #endif
      {
	      if(n<=9)
                 sprintf(configfile,"vtk/%s%ss00%d.vtk",labels[i/4],name,n);
              else if(n<=99)
                 sprintf(configfile,"vtk/%s%ss0%d.vtk",labels[i/4],name,n);
              else
                 sprintf(configfile,"vtk/%s%ss%d.vtk",labels[i/4],name,n);

	      printf("check dims %s %s %d %d \n",configfile,labels[i/4],ni,nj);
	      FILE *fdt=fopen(configfile,"w");


	      fprintf(fdt,"# vtk DataFile Version 2.0\n");
	      fprintf(fdt,"Structured Grid\n");
	      fprintf(fdt,"ASCII\n");
	      fprintf(fdt," \n");
	      fprintf(fdt,"DATASET RECTILINEAR_GRID\n");
	      fprintf(fdt,"DIMENSIONS %d %d 1\n",ni-2*(p.ng[0]),nj-2*(p.ng[1]));


	      fprintf(fdt,"X_COORDINATES %d double\n",ni-2*(p.ng[0]));
              for(i1=0;i1<ni-(2*(p.ng[0]));i1++)
	        fprintf(fdt,"%f\n",i1*p.dx[0]);

	      fprintf(fdt,"Y_COORDINATES %d double\n",nj-2*(p.ng[1]));
              for(i1=0;i1<nj-(2*(p.ng[1]));i1++)
	        fprintf(fdt,"%f\n",i1*p.dx[1]);

	      fprintf(fdt,"Z_COORDINATES 1 double\n");
	      fprintf(fdt,"0\n");

	      fprintf(fdt,"POINT_DATA  %d\n",(ni-2*(p.ng[0]))*(nj-2*(p.ng[1])));
	      fprintf(fdt,"SCALARS %s double 1\n",labels[i/4]);

             fprintf(fdt,"LOOKUP_TABLE TableName \n");

	     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
		for( i1=(p.ng[0]);i1<(ni-(p.ng[0]));i1++)
			fprintf(fdt,"%f\n",w[(j1*ni+i1)+(ni*nj*is)]);

	      fclose(fdt);
      }

      //vector fields
      int iv;
        #ifdef ADIABHYDRO
      for(int i=2; i<3; i++)
       #else
      for(int i=2; i<=3; i++)
        #endif

      {
	      if(i==2)
                iv=1;
              else
                //iv=5;
                iv=4;
              if(n<=9)
                 sprintf(configfile,"vtk/%s%ss00%d.vtk",labels[i],name,n);
              else if(n<=99)
                 sprintf(configfile,"vtk/%s%ss0%d.vtk",labels[i],name,n);
              else
                 sprintf(configfile,"vtk/%s%ss%d.vtk",labels[i],name,n);

	      printf("check dims %s %s %d %d \n",configfile,labels[i],ni,nj);
	      FILE *fdt=fopen(configfile,"w");


	      fprintf(fdt,"# vtk DataFile Version 2.0\n");
	      fprintf(fdt,"Structured Grid\n");
	      fprintf(fdt,"ASCII\n");
	      fprintf(fdt," \n");
	      fprintf(fdt,"DATASET RECTILINEAR_GRID\n");
	      fprintf(fdt,"DIMENSIONS %d %d 1\n",ni-2*(p.ng[0]),nj-2*(p.ng[1]));


	      fprintf(fdt,"X_COORDINATES %d double\n",ni-2*(p.ng[0]));
              for(i1=0;i1<ni-(2*(p.ng[0]));i1++)
	        fprintf(fdt,"%f\n",i1*p.dx[0]);

	      fprintf(fdt,"Y_COORDINATES %d double\n",nj-2*(p.ng[1]));
              for(i1=0;i1<nj-(2*(p.ng[1]));i1++)
	        fprintf(fdt,"%f\n",i1*p.dx[1]);

	      fprintf(fdt,"Z_COORDINATES 1 double\n");
	      fprintf(fdt,"0\n");

	      fprintf(fdt,"POINT_DATA  %d\n",(ni-2*(p.ng[0]))*(nj-2*(p.ng[1])));
	      fprintf(fdt,"VECTORS %s double \n",labels[i]);

		for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
	      		for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
			 //fprintf(fdt,"%f %f %f\n",w[(j1*ni+i1)+(ni*nj*iv)],w[(j1*ni+i1)+(ni*nj*(iv+1))],w[(j1*ni+i1)+(ni*nj*(iv+2))]);    
                         fprintf(fdt,"%f %f %f\n",w[(j1*ni+i1)+(ni*nj*iv)],w[(j1*ni+i1)+(ni*nj*(iv+1))]);

                       //printing mag fields including backround for SAC
                       //if(iv==4)
                       //  fprintf(fdt,"%f %f %f\n",w[(j1*ni+i1)+(ni*nj*iv)]+w[(j1*ni+i1)+(ni*nj*(iv+4))],w[(j1*ni+i1)+(ni*nj*(iv+1))]+w[(j1*ni+i1)+(ni*nj*(iv+1+4))]);


	      fclose(fdt);
      }




  return status;
}



int readconfig(char *cfgfile, params p, meta md, real *w)
{
  int status=0;

  return status;
}
