#include "paramssteeringtest1.h"
#include <iome/genericsimulationlib/IoGenericSimulationLib.h>

void createsim(params k,  meta metadata,char *simname, iome el)
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


// Constants
//adddoubleparam_(el.id,"g",k.g,7,el.port,el.server);
//adddoubleparam_(el.id,"u0",k.u0,7,el.port,el.server);
//adddoubleparam_(el.id,"v0",k.v0,7,el.port,el.server);
//adddoubleparam_(el.id,"b",k.b0,7,el.port,el.server);
//adddoubleparam_(el.id,"h0",k.h0,7,el.port,el.server);

//Domain definition
// Define the x domain
//ni = 151; 
ni=k.ni;
xmax = k.xmax;                      
k.dx = xmax/(ni-1);
//x  = [0:dx:xmax];

// Define the y domain
//nj = 151;  
nj=k.nj;
ymax = k.ymax;                      
k.dy = ymax/(nj-1);
//y  = [0:dy:ymax];

tmax = k.tmax;

// Define the wavespeed
dt=k.dt;
//wavespeed = k.u0 + sqrt(k.g*(k.h0 - k.b0));

// Define time-domain
//dt = 0.68*dx/wavespeed;

//t = [0:dt:tdomain];
//t=[1:dt:tmax];
k.nt=(int)((tmax-1)/dt);
courant = wavespeed*dt/dx;

adddoubleparam_(el.id,"ni",k.ni,7,el.port,el.server);
adddoubleparam_(el.id,"nj",k.ni,7,el.port,el.server);
adddoubleparam_(el.id,"xmax",k.xmax,7,el.port,el.server);
adddoubleparam_(el.id,"ymax",k.ymax,7,el.port,el.server);
adddoubleparam_(el.id,"tmax",k.tmax,7,el.port,el.server);
addintparam_(el.id,"nt",k.nt,7,el.port,el.server);
addintparam_(el.id,"steeringenabled",k.steeringenabled,7,el.port,el.server);
addintparam_(el.id,"finishsteering",k.finishsteering,7,el.port,el.server);
addintparam_(el.id,"step",k.dt,7,el.port,el.server);






addstringparam_(el.id,"resultsfile","results.zip",7,el.port,el.server);
//simfile=sprintf('%s.xml',simname)



//endfunction
}

