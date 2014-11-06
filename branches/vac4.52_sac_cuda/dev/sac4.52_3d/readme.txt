The Versatile Advection Code is used for much of the work at Sheffield, a parallel version of the code is used, VAC is a collection of FORTRAN source code that features a number of modules that may be devloped and included by the user, for example, one of the models developed by the SPARC group investigates MHD for gravitationally stratified media.
The resulting computational tasks for 2D models require at least 10 processors they take approximately 100hrs of wall clock time (on an AMD 64 bit opteron processor) and generate around 50GB of data. It can be seen that the nature of these computations require the introduction of more sophisticated management techniques. The group will soon run 3D models making the provision of these management techniques a more pressing need.
Given the requirement for storage and processing cores it seems that this project may benefit from access to additional compute clusters and the processing power and storage made available. Utilisation of resources in this way increases the complexity of the research process. This necessitates the use of management tools to simplify  the research process. We will provide prototype tools from the web accessible portal described in an earlier posting. To simplify the development task the application will be tested on the Sheffield node of the White Rose grid.
After writing the user module, describing the physics for the problem of interest. The researcher undertakes the following steps, we assume here using 10 processors.
Change parameters and source expressions hard coded in the user FORTRAN module.
Before compiling the models set the parameters for building the model e.g. switch MPI on (or off) and set the mesh size. The setvac program enables these settings to be made.

setvac -d=22 -phi=0 -z=0 -g=1976,44 -p=mhd -u=sim1 \
       -on=cd,rk,mpi -off=mc,fct,tvdlf,tvd,impl,poisson,ct,gencoord,resist



 ./setvac -on=mpi   switching mpi on
./setvac -g=1976,44
In the above case the model is a 1976x400 model it will run on 10 processors each processor will have 1976x(40+2+2) cells
this includes 2 layers of ghost cells at each boundary
Edit the initialisation parameters (vacini.par) for example
the domain for this will be 1976,400
set the path to the .ini file
set the atmosphere file  this is density profile through solar interior
Edit vac.par
  set filelist the filenamein and filename  need to be set correctly
Make the initilisation and the simulation routine using the model must be compiled using the appropriate parallel compiler
make -vacini
make -vac
Generate the vacin file as follows
./vacini < vacini.par
Distribuute the ini file across the processors using the distribution program, eg.
./distribution -D -s =0 /data/username/VAC_NN/2_6Mnzx1976400.ini /data/username/VAC_NN/2_6mnzx1976400_np0110.ini
 Notation is as follows
Model depth is 2_6Mm
nzx  1976x400  (z direction is 1976)
np is number of processors in z and x directions
01 processors in z direction
10 processors in x direction
Submit job to batch queue using the correct parallel queue
Gather all the generated data into same data files
vac4.52/data/distribution 2_6Mnzx1023400_cont2_np0110.out test.out
the processor rank in filename(if any) is ignored
if successful remove the old files that are unmerged
convertdata to dx format
Run the data visualiser
In future posts we will describe how these requirements and user activities are converted into proposed user application modules. Preliminary suggestions are as follows.
Web application generates parameter and initilisation files and compiles executables the generated module is returned for execution/submission by the user.
Submit a job (as generated above) to the grid, ensure the generated data is moved to the correct location.
Post process the data and prepare for use by the visualisation tools
Tool to generate images and movies for a specified model
Metadata management and generation
