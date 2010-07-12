#include "paramssteeringtest1.h"
#include <iome/genericsimulationlib/IoGenericSimulationLib.h>

void createsim(constants k, domain dom,source src, meta metadata,char *simname, iome el)
{
int i;
int ni,nj;
double xmax,ymax;
double dx,dy,dt,tmax,wavespeed;
double courant;
double *statsu, *statsv, *statsh;
//elist=list();  parameter used by iome to contain port and server address
//elist=list();

addmetadata_(el.id,"author",metadata.author,el.port,el.server);
addmetadata_(el.id,"directory",metadata.directory,el.port,el.server);
addmetadata_(el.id,"date",metadata.sdate,el.port,el.server);
addmetadata_(el.id,"platform",metadata.platform,el.port,el.server);
addmetadata_(el.id,"description",metadata.desc,el.port,el.server);
addmetadata_(el.id,"name",metadata.name,el.port,el.server);

adddoubleparam_(el.id,"frequency",src.freq,7,el.port,el.server);
adddoubleparam_(el.id,"amplitude",src.amp,7,el.port,el.server);
adddoubleparam_(el.id,"sx",src.xloc,7,el.port,el.server);
adddoubleparam_(el.id,"sy",src.yloc,7,el.port,el.server);

// Constants
adddoubleparam_(el.id,"g",k.g,7,el.port,el.server);
adddoubleparam_(el.id,"u0",k.u0,7,el.port,el.server);
adddoubleparam_(el.id,"v0",k.v0,7,el.port,el.server);
adddoubleparam_(el.id,"b",k.b0,7,el.port,el.server);
adddoubleparam_(el.id,"h0",k.h0,7,el.port,el.server);

//Domain definition
// Define the x domain
//ni = 151; 
ni=dom.ni;
xmax = dom.xmax;                      
dx = xmax/(ni-1);
//x  = [0:dx:xmax];

// Define the y domain
//nj = 151;  
nj=dom.nj;
ymax = dom.ymax;                      
dy = ymax/(nj-1);
//y  = [0:dy:ymax];

tmax = dom.tmax;

// Define the wavespeed
wavespeed = k.u0 + sqrt(k.g*(k.h0 - k.b0));

// Define time-domain
dt = 0.68*dx/wavespeed;

//t = [0:dt:tdomain];
//t=[1:dt:tmax];
dom.nt=(int)((tmax-1)/dt);
courant = wavespeed*dt/dx;

adddoubleparam_(el.id,"ni",dom.ni,7,el.port,el.server);
adddoubleparam_(el.id,"nj",dom.ni,7,el.port,el.server);
adddoubleparam_(el.id,"xmax",dom.xmax,7,el.port,el.server);
adddoubleparam_(el.id,"ymax",dom.ymax,7,el.port,el.server);
adddoubleparam_(el.id,"tmax",dom.tmax,7,el.port,el.server);
addintparam_(el.id,"nt",dom.nt,7,el.port,el.server);
addintparam_(el.id,"steeringenabled",dom.steeringenabled,7,el.port,el.server);
addintparam_(el.id,"finishsteering",dom.finishsteering,7,el.port,el.server);
addintparam_(el.id,"step",dom.step,7,el.port,el.server);




statsu=(double *)calloc(3*dom.nt,sizeof(double));
statsv=(double *)calloc(3*dom.nt,sizeof(double));
statsh=(double *)calloc(3*dom.nt,sizeof(double));
for(i=0;i<3*dom.nt;i++)
{
	statsu[i]=0;
	statsv[i]=0;
	statsh[i]=0;
}

addmatparam_(el.id,"statsu",statsu,dom.nt,3,7,el.port,el.server);
addmatparam_(el.id,"statsv",statsv,dom.nt,3,7,el.port,el.server);
addmatparam_(el.id,"statsh",statsh,dom.nt,3,7,el.port,el.server);

addstringparam_(el.id,"resultsfile","results.zip",7,el.port,el.server);
//simfile=sprintf('%s.xml',simname)



//endfunction
}

