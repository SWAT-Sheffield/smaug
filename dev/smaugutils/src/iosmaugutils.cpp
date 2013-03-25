// iosmaug.cpp : Main routine for GPU enabled SAC
#include "../include/iosmaug.h"
#include "../include/iobparams.h"


int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int mode=run;//run a model 1=scatter 2=gather
int it=0; //test integer to be returned 
int n;
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

char configinfile[300];

 real tcom,tcom1, tcom2,tv,tcal,tc;

tcom=0.0;
tcal=0.0;

#include "../include/defs.h"
#include "../include/iosmaugparams.h"


struct bparams *d_bp;
struct bparams *bp=(struct bparams *)malloc(sizeof(struct bparams));


FILE *portf;


if(argc>3  && strcmp(argv[2],"gather")==0 && (atoi(argv[3])>=0) && (atoi(argv[3])<=nt)) 
{    
  mode=gather;
}

if(argc>2  && strcmp(argv[2],"scatter")==0)
{
  mode=scatter;
}

if(argc>2  && strcmp(argv[2],"init")==0)
{
  mode=init;
  printf("init mode=3\n");
}

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


       /*********************************************************************************************************/
       /* Start of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/


	char ext[3];
	char tcfg[300];
	char stemp[300];
	char *pch1,*pch2;
	strcpy(stemp,cfgfile);
	pch1 = strtok (stemp,".");
	sprintf(tcfg,"%s",pch1);
	pch2 = strtok (NULL,".");
	sprintf(ext,"%s",pch2);
	sprintf(configfile,"%s",cfgout);
	#ifdef USE_MULTIGPU
	#ifdef USE_MPI
	     MPI::Init(argc, argv);
	#endif
	mgpuinit(p);
	ipe2iped(p);     
	mgpuneighbours(0,p);
	mgpuneighbours(1,p);

	//compute the max and min domain dimensions for each processor
	p->xmax[0]=xmin+(1+(p->pipe[0]))*(xmax-xmin)/(p->pnpe[0]);
	p->xmax[1]=ymin+(1+(p->pipe[1]))*(ymax-ymin)/(p->pnpe[1]);
	p->xmin[0]=xmin+(p->pipe[0])*(xmax-xmin)/(p->pnpe[0]);
	p->xmin[1]=ymin+(p->pipe[1])*(ymax-ymin)/(p->pnpe[1]);

	//store global values for max and min domain dimensions
	p->gxmax[0]=xmax;
	p->gxmin[0]=xmin;
	p->gxmax[1]=ymax;
	p->gxmin[1]=ymin;

	#ifdef USE_SAC3D
	mgpuneighbours(2,p);
	p->xmax[2]=zmin+(1+(p->pipe[2]))*(zmax-zmin)/(p->pnpe[2]);
	p->xmin[2]=zmin+(p->pipe[2])*(zmax-zmin)/(p->pnpe[2]);
	p->gxmax[2]=zmax;
	p->gxmin[2]=zmin;

	#endif
	  sprintf(configinfile,"%s",cfgfile);

	//adopt the sac MPI naming convention append the file name npXXYY where XX and YY are the
	//number of processors in the x and y directions
	#ifdef USE_MPI
	     #ifdef USE_SAC3D
		      if(p->ipe>99)
			sprintf(configinfile,"%s_np%d%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np0%d0%d0%d_0%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np00%d00%d00%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);  	     
	     #else
		      if(p->ipe>99)
			sprintf(configinfile,"%s_np%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np%d%d_%d.0%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np0%d0%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);  	     	     
	     #endif
	#endif

	//if doing a scatter or gather set the domain size correctly
	//take a distribution and distribute domain to processors
	if(mode==scatter )
	{
	  printf("Scatter %s \n",cfgfile);
	  sprintf(configinfile,"%s",cfgfile);
	  p->n[0]=ni*(p->pnpe[0]);
	  p->n[1]=nj*(p->pnpe[1]);
	   #ifdef USE_SAC3D
		    p->n[2]=nk*(p->pnpe[2]);
	   #endif
	}

	if( mode==gather)
	{
	   ni=ni*(p->pnpe[0]);
	   nj=nj*(p->pnpe[1]);
	   #ifdef USE_SAC3D
		   nk=nk*(p->pnpe[2]);
	   #endif
	}


	if(mode==init)
	{
	    p->n[0]=ni;
	    p->n[1]=nj;
	    #ifdef USE_SAC3D
	      p->n[2]=nk;
	    #endif
	}
	printf("config files\n%s \n %s %d %d\n",configinfile,configfile,p->n[0],p->n[1]);


	#else
	     sprintf(configinfile,"%s",cfgfile);
	#endif   //#ifdef USE_MULTIGPU

       /*********************************************************************************************************/
       /* End of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/




char *method=NULL;


       /*********************************************************************************************************/
       /* Start of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/


        //printf("cfgfile %s\n",configfile);
        //   getintparam_( &elist.id,"i1",&it,  &elist.port, "localhost" );	
        //	printf("Get integer %d\n",it);
        //Set input filename as first arg
	//if NULL use defaults
	
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


       /*********************************************************************************************************/
       /* End of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/










       /*********************************************************************************************************/
       /* Start of section creating arrays on the host*/
       /*********************************************************************************************************/
  	printf("Creating arrays on the host\n");

       #ifdef USE_MULTIGPU
       if(mode==0)
       {
		if((p->pipe[0])==0) (p->n[0])+=ngi;
		if((p->pipe[0])==((p->pnpe[0])-1)) (p->n[0])+=ngi;
		if((p->pipe[1])==0) (p->n[1])+=ngj;
		if((p->pipe[1])==((p->pnpe[1])-1)) (p->n[1])+=ngj;

		#ifdef USE_SAC_3D
			if((p->pipe[2])==0) (p->n[2])+=ngk;
			if((p->pipe[2])==((p->pnpe[2])-1)) (p->n[2])+=ngk;
		#endif
	}
       #endif

	//allocate arrays to store fields, updated fields, dervived quantities and updated derived qunatities
	#ifdef USE_SAC_3D
		wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*((p)->n[2])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*((p)->n[2])*NVAR,sizeof(real ));
	#else
		wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*NVAR,sizeof(real ));
	#endif

        #ifdef USE_MULTIGPU
          //parameters used to set sizes of MPI communications buffers and
	  //data storage areas
          int szw,szw0,szw1,szw2,szvisc0,szvisc1,szvisc2;
	  #ifdef USE_SAC
		  szw=4*(  ((p)->n[1])  +  ((p)->n[0])   );
		  szw0=4*NDERV*(  ((p)->n[1])     );
		  szw1=4*NDERV*(  ((p)->n[0])     );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2 )   );
		  szvisc1=4*NVAR*(    (((p)->n[0]) +2 )  );
	  #endif
	  #ifdef USE_SAC_3D	  
		  szw=4*NDERV*(  ((p)->n[1])*((p)->n[2])  +  ((p)->n[0])*((p)->n[2])  +  ((p)->n[0])*((p)->n[1])  );
		  szw0=4*NDERV*(  ((p)->n[1])*((p)->n[2])    );
		  szw1=4*NDERV*(    ((p)->n[0])*((p)->n[2])   );
		  szw2=4*NDERV*(    ((p)->n[0])*((p)->n[1])  );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)  ); 
		  szvisc1=4*NVAR*(   (((p)->n[0])+2)*(((p)->n[2])+2)    );    
		  szvisc2=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)   );   
	  #endif

	  #ifdef USE_SAC
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2),sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2)* (((p)->n[2])+2),sizeof(real));
	  #endif

          //Data areas to store values communicated using MPI
	  gmpiwmod0=(real *)calloc(szw0,sizeof(real));
	  gmpiw0=(real *)calloc(szw0,sizeof(real));
	  gmpiwmod1=(real *)calloc(szw1,sizeof(real));
	  gmpiw1=(real *)calloc(szw1,sizeof(real));

          gmpivisc0=(real *)calloc(szvisc0,sizeof(real));
          gmpivisc1=(real *)calloc(szvisc1,sizeof(real));

	  #ifdef USE_SAC_3D
		  gmpiwmod2=(real *)calloc(szw2,sizeof(real));
		  gmpiw2=(real *)calloc(szw2,sizeof(real));
                  gmpivisc2=(real *)calloc(szvisc2,sizeof(real));
	  #endif
        #endif
       /*********************************************************************************************************/
       /* End of section creating arrays on the host*/
       /*********************************************************************************************************/

	//set initial time step to a large value
	if(p->moddton==1.0)
	{
		p->dt=1.0e-8;
	}
       int its=p->it;


       /*********************************************************************************************************/
       /* Start of section initialising the configuration 
          on the host and on GPU host memory*/
       /*********************************************************************************************************/

       if(mode !=init)
       {
               if((p->readini)==0)
               {
                 printf("init config\n");
		 initconfig(p, &meta, wmod,wd);
                }
		else
                {
	         printf("reading configuration from %s\n",configinfile);
		 readasciivacconfig(configinfile,*p,meta, state,wmod,wd,hlines,mode);
                }
       }



       /*********************************************************************************************************/
       /* Start of section to scatter data
        /*********************************************************************************************************/
        #ifdef USE_MULTIGPU
        //scatter/distribute configuration across each CPU
        if(mode==scatter)
        {
	       gpusync();
               if(p->ipe==0) //currently only processor zero
	       {

		  for(i=0; i<p->npe; i++)
		  {
		    p->ipe=i;
                    ipe2iped(p);
		    
		    //copy segment
		    printf("copy segment %d\n",i);                    
		    createconfigsegment(*p, wnew,wdnew,wmod,wd);  //in readwrite.c

		    //writeas
                    //set domain size to size for each processor		   
                    p->n[0]=ni;
                    p->n[1]=nj;
                    #ifdef USE_SAC3D
                      p->n[2]=nk;
                    #endif
		    writeasciivacconfig(configinfile, *p, meta,  wnew,wdnew, hlines, *state,mode);
                    //set domain size to the global domain size
                    //this will be used when we extract a segment
                    p->n[0]=ni*(p->pnpe[0]);
                    p->n[1]=nj*(p->pnpe[1]);
                    #ifdef USE_SAC3D
                      p->n[2]=nk*(p->pnpe[2]);
                    #endif
		  }
		}
                gpusync();
        }
       /*********************************************************************************************************/
       /* End of section to scatter data
        /*********************************************************************************************************/
 

	/*********************************************************************************************************/
	/* Start of section to gather data
	/*********************************************************************************************************/
	//gather configuration to single output file
	if(mode==gather)
	{
		n=atoi(argv[3]);
		if(p->ipe==0)
		{
			int myipe=p->ipe;
			for(i=0; i<p->npe; i++)
			{
				printf(" here nt=%d pid=%d i=%d\n",n,p->ipe,i);

				p->ipe=i;
				ipe2iped(p);
				strcpy(stemp,cfgout);
				pch1 = strtok (stemp,".");
				sprintf(tcfg,"%s",pch1);

				#ifdef USE_SAC3D
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np%d%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np0%d0%d0%d_0%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else
						sprintf(configinfile,"%s%d_np00%d00%d00%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);  	     
				#else
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else
						sprintf(configinfile,"%s%d_np0%d0%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);  	     	     
				#endif

				//copy segment
				printf("copy segment %d %s\n",i,configinfile);

				#ifdef USE_MULTIGPU
					readasciivacconfig(configinfile,*p, meta, state, wmod,wd, hlines,mode);
				#else
					readbinvacconfig(configinfile,*p, meta, wmod,wd, *state );
				#endif
				gathersegment(*p, wnew,wdnew,wmod,wd);
				printf(" here read and gath nt=%d pid=%d i=%d\n",n,p->ipe,i);
			}

			p->n[0]=ni;
			p->n[1]=nj;
			#ifdef USE_SAC3D
				p->n[2]=nk;
			#endif
			#ifdef USE_SAC3D
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     
			#else
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     	     
			#endif           

			state->it=n;
			//sprintf(configfile,"%s",cfggathout);
			writevacgatherconfig(configfile,n,*p, meta , wnew,wdnew,*state);
			printf(" here configfile %s nt=%d pid=%d \n",configfile,n,p->ipe);
			p->ipe=myipe;

		}//if p->ipe==0

		gpusync();
		//}//loop over nt steps
		//printf("proc %d here \n", p->ipe);

	}//if mode==gather
	#endif
	/*********************************************************************************************************/
	/* End of section to gather data
	/*********************************************************************************************************/



	/*********************************************************************************************************/
	/* Start of section to run special user initialisation
	/*********************************************************************************************************/
        //special user initialisation for the configuration 
        //this is a parallel routine
        if(mode==init)
        {
		p->mode=mode;

		#ifdef USE_MULTIGPU
			gpusync();
		#endif
                printf("init_config\n");
		initconfig(p, &meta, wmod,wd);
                printf("user initialisation\n");
		initialisation_user1(wmod,wd,p);

		// initialisation_user2(wmod,wd,p);
		//write the config file to ascii
                printf("writing ini file\n");
		writeasciivacconfig(configinfile,*p, meta , wmod,wd,hlines,*state,mode);
		#ifdef USE_MULTIGPU
			gpusync();
		#endif
        }   
	/*********************************************************************************************************/
	/* End of section to run special user initialisation
	/*********************************************************************************************************/



	//p->it=0;
	int order=0;



#ifdef USE_MPI
	     printf("at cumpifinish end here %d\n",p->ipe);
#endif
	free(hlines);
	free(p);
	free(bp);
	free(sdir);
	free(name);
	free(outfile);
	free(formfile);


	#ifdef USE_IOME
		writesimulation_(elist.id,newsimfile,elist.port,elist.server);
	#endif
	#ifdef USE_MPI
          ;// mgpufinalize(p);
        #endif
	#ifdef USE_MPI
	     printf("at cumpifinish end here %d\n",p->ipe);

	    //cufinishmgpu(&p,&w, &wmod, &temp2,&gmpivisc0,&gmpivisc1,&gmpivisc2,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,   &d_w, &d_wmod,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2);
            ;// mgpufinalize(p);


	#endif
		return 0;
	}

