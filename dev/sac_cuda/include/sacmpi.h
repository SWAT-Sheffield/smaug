#include "mpi.h"
#include "iotypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

MPI::Intracomm comm;
double gwall_time;

/*void mpiinit(params *p);
void mpifinalize(params *p);
void mpisetnpediped(params *p, char *string);
void ipe2iped(params *p);
void iped2ipe(params *p);*/

//!=============================================================================
//subroutine mpiinit
//
//! Initialize MPI variables
//

//!----------------------------------------------------------------------------
//call MPI_INIT(ierrmpi)
//call MPI_COMM_RANK (MPI_COMM_WORLD, ipe, ierrmpi)
//call MPI_COMM_SIZE (MPI_COMM_WORLD, npe, ierrmpi)

//! unset values for directional processor numbers
//npe1=-1;npe2=-1;
//! default value for test processor
//ipetest=0
void mpiinit(params *p)
{
     
     //MPI::Intracomm comm;
     //MPI_Init(&argc, &argv);
     gwall_time = MPI_Wtime();
     comm=MPI::COMM_WORLD;
     p->npe=comm.Get_size();
     p->ipe=comm.Get_rank();	
}



//!==============================================================================
//subroutine mpifinalize

//call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
//call MPI_FINALIZE(ierrmpi)

void mpifinalize(params *p)
{
     gwall_time = MPI_Wtime() - gwall_time;
     if ((p->ipe) == 0)
	  printf("\n Wall clock time = %f secs\n", gwall_time);
     MPI_Finalize();
}





//!==============================================================================
//subroutine mpisetnpeDipeD(name)

//! Set directional processor numbers and indexes based on a filename.
//! The filename contains _np followed by np^D written with 2 digit integers.
//! For example _np0203 means np1=2, np2=3 for 2D.

//! Extract and check the directional processor numbers and indexes
//! and concat the PE number to the input and output filenames

void mpisetnpediped(params *p, char *string)
{    
  
}


//!==============================================================================
//subroutine ipe2ipeD(qipe,qipe1,qipe2)
//
//! Convert serial processor index to directional processor indexes

//integer:: qipe1,qipe2, qipe
//!-----------------------------------------------------------------------------
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 

void ipe2iped(params *p)
{

#ifdef USE_SAC_3D
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 
//qipe3 = qipe/(npe1*npe2)
(p->pipe[0])=(p->ipe)-(p->pnpe[0])*(p->ipe)/(p->pnpe[0]);
(p->pipe[1])=((p->ipe)/(p->pnpe[0]))-(p->pnpe[1])*(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));
(p->pipe[2])=(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));   


//set upper boundary flags
//mpiupperB(1)=ipe1<npe1-1
//mpilowerB(1)=ipe1>0 

//mpiupperB(2)=ipe2<npe2-1
//mpilowerB(2)=ipe2>0 

(p->mpiupperb[0])=(p->pipe[0])<((p->npe[0]-1);
(p->mpiupperb[1])=(p->pipe[1])<((p->npe[1]-1);
(p->mpiupperb[2])=(p->pipe[2])<((p->npe[2]-1);

(p->mpilowerb[0])=(p->pipe[0])>0;
(p->mpilowerb[1])=(p->pipe[1])>0;
(p->mpilowerb[2])=(p->pipe[2])>0;
#else
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 
(p->pipe[0])=(p->ipe)-(p->pnpe[0])*(p->ipe)/(p->pnpe[0]);
(p->pipe[1])=((p->ipe)/(p->pnpe[0]))-(p->pnpe[1])*(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));

(p->mpiupperb[0])=(p->pipe[0])<((p->pnpe[0])-1);
(p->mpiupperb[1])=(p->pipe[1])<((p->pnpe[1])-1);

(p->mpilowerb[0])=(p->pipe[0])>0;
(p->mpilowerb[1])=(p->pipe[1])>0;

#endif

}

//!==============================================================================
//subroutine ipeD2ipe(qipe1,qipe2,qipe)
//
//! Convert directional processor indexes to serial processor index

//include 'vacdef.f'

//integer:: qipe1,qipe2, qipe
//!-----------------------------------------------------------------------------
//qipe = qipe1  + npe1*qipe2

void iped2ipe(int *tpipe,int *tpnp, int *oipe)
{
  #ifdef USE_SAC_3D
  //qipe = qipe1  + npe1*qipe2  + npe1*npe2*qipe3
  (*oipe)=tpipe[0]+(tpnp[0])*(tpipe[1])+(tpnp[0])*(tpnp[1])*(tpipe[2]);
  #else
  (*oipe)=tpipe[0]+(tpnp[0])*(tpipe[1]);
  #endif

}

//!==============================================================================
//subroutine mpineighbors(idir,hpe,jpe)

//! Find the hpe and jpe processors on the left and right side of this processor 
//! in direction idir. The processor cube is taken to be periodic in every
//! direction.

//!-----------------------------------------------------------------------------
//hpe1=ipe1-kr(1,idir);hpe2=ipe2-kr(2,idir);
//jpe1=ipe1+kr(1,idir);jpe2=ipe2+kr(2,idir);

//if(hpe1<0)hpe1=npe1-1
//if(jpe1>=npe1)jpe1=0

//if(hpe2<0)hpe2=npe2-1
//if(jpe2>=npe2)jpe2=0

//call ipeD2ipe(hpe1,hpe2,hpe)
//call ipeD2ipe(jpe1,jpe2,jpe)
void mpineighbours(int dir, params *p)
{
     int i;
     for(i=0; i<NDIM;i++)
     {
             (p->phpe[i])=(p->pipe[i])-(dir==i);
             (p->pjpe[i])=(p->pipe[i])+(dir==i);             
     }

     for(i=0; i<NDIM;i++)
     {
              if((p->phpe[i])<0) (p->phpe[i])=(p->pnpe[i])-1; 
              if((p->pjpe[i])<0) (p->pjpe[i])=0;                    
     }
     
     iped2ipe(p->phpe,p->pnpe,&(p->hpe));
     iped2ipe(p->pjpe,p->pnpe,&(p->jpe));
}





