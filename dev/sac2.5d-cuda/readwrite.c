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

int writeconfig(char *name,int n,params p, meta md, float *w)
{
  int status=0;
  int i1,j1;
  int ni,nj;
  char configfile[300];


  ni=p.ni;
  nj=p.nj;



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
		fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


      //save file containing current data
      sprintf(configfile,"out/%s.out",name);
      printf("check dims %d %d \n",ni,nj);
      fdt=fopen(configfile,"a+");
      fprintf(fdt,"%d\n",n);
     for( j1=0;j1<nj;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
               // printf("%d %d ", i1,j1);
		fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


  return status;
}

int readconfig(char *cfgfile, params p, meta md, float *w)
{
  int status=0;

  return status;
}
