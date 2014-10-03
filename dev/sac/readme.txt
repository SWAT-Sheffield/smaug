

Compile, Run and visualise results for the Orszag-Tang Test

 cd src
 cp vacusr.t.sim1_OT vacusr.t.sim1
 

 ./setvac -d=22 -g=256,256 -u=sim1 -p=mhd -off=mpi
 make vacini
 make vac

 cd ..
 cp par/vac_OT.par vac.par

 Generate starting configuration for the Orszag-Tang test
 ./vacini < inipar/vacini_OT.par

 To generate a parallel configuration for a 2x2 array of processors
 ./distribution configs/zero1_ot_bin_256.ini configs/zero1_ot_bin_256_np0202.ini

 ./vac < vac.par

or qsub run.sh


 Visualise the results


 cd data
 idl
 
  From IDL (Using IDL procedure visex22D.pro in the data folder)

    .r procedures
    .r tvframe
    .r visex22D

******************************************************************************************

Compile, Run and visualise results for the Brio-Wu Test

 cd src
 cp vacusr.t.sim1_BW vacusr.t.sim1
 

 ./setvac -d=22 -g=800,10 -u=sim1 -p=mhd -off=mpi
 make vacini
 make vac

 cd ..
 cp par/vac_nm_BW.par vac.par

 Generate starting configuration for the Orszag-Tang test
 ./vacini < inipar/vacini_BW.par

 To generate a parallel configuration for a 4x1 array of processors
 ./distribution configs/zero1_bw_bin_256.ini configs/zero1_bw_bin_256_np0401.ini

 ./vac < vac.par

or qsub run.sh


 Visualise the results


 cd data
 idl
 
  From IDL (Using IDL procedure visex22D.pro in the data folder)

    .r procedures
    .r tvframe
    .r visex22D_BW

******************************************************************************************







Compile, Run and visualise results for the Flux tube Test

 cd src
 cp vacusr.t.sim1_TUBE vacusr.t.sim1
 

 ./setvac -d=33 -g=200,104,104 -u=sim1 -p=mhd -off=mpi
 make vacini
 make vac

 cd ..
 cp par/vac_TUBE.par vac.par

 Generate starting configuration for the Flux tube test
 ./vacini < inipar/vacini_TUBE.par

 Generate the magnetic flux tube
 cd mag_field
 idl

 From the idl session
 .r procedures
 .r tvframe
 .r B_field_vertical_tube


 To generate a parallel configuration for a 2x2 array of processors
 cd ..
 ./distribution configs/3D_tube_modif_200_100_100.ini configs/3D_tube_modif_200_100_100_np010202.ini.ini


cd src
make clean
module add mpi/pgi/openmpi/1.6.4 
./setvac -d=33 -g=200,54,54 -u=sim1 -p=mhd -on=mpi

make vacini -f Makefile_MPI
make vac -f Makefile_MPI

Submit the job
 qsub mpi_run_mhd.sh


 Visualise the results


 cd data
 idl
 
  From IDL (Using IDL procedure visex22D.pro in the data folder)

    .r procedures
    .r tvframe
    .r visex22D








To Do

1. Add example files - models
                     - can directly compare with smaug

2. Add basic readme

3. Guidelines for using hdf5


4. convert and check pbc implementation and switch between fixed


5. check that it compiles on N8 using intel compiler
                      
