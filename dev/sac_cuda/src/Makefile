#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_inputs
### End User configurable options ###







BASICEXOBJS	=\
	iosac2.5d.o\
	initialisation.o\
	dxroutines.o\
	hyperdifmomsource1.o\
	hyperdifmomsourcene1.o\
	hyperdifesource1.o\
	hyperdifbsource1.o\
	hyperdifbsourcene1.o\
	hyperdifrhosource1.o\
	hyperdifvisc1.o\
	centdiff1.o\
	centdiff2.o\
	computedervfields.o\
	boundary.o\
	init.o\
	update.o\
        advance.o\
	readwrite.o\
	divb.o

PROCEXOBJS	=\
	gradops_pre.o
	
#	../../IoGenericSimulationLib.o\
#	../../../IoSimulation/IoGenericSteerSimulation.o\
#	../../../IoSimulation/stdsoap2.o\
#	../../../IoSimulation/soapC.o\
#	../../../IoSimulation/soapClient.o
	
ALLOBJS	=	$(obj)/$(BASICEXOBJS)
ALLBIN	=	$(Bin)/iosac

#use -DUSE_IOME to enable IOME
#CUDACCFLAGS = --ptxas-options=-v -maxrregcount=32 --device-emulation   -DUSE_REAL -DUSE_SAC -DADIABHYDRO -DUSE_VAC -DD1D -DUSE_IOME -DUSE_USERSOURCE -DUSE_SAC_3D -DUSE_DORDER3
#compiling with  -arch sm_13 enables double precision to be used
CUDACCFLAGS = --ptxas-options=-v -arch sm_13  -maxrregcount=32 -DUSE_SAC
CCFLAGS = -DUSE_SAC -DD1D

#MYLIBS	=	libAgentSimulation2.a

#@# Targets follow ---------------------------------

##all:	IoSimTestSimple

bach3d:
	cp ../models/iosac2.5d_params.h.bach3d ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.bach3d ../src/init_user.cu
	cp ../models/boundary.cu.bach3d ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu


ah:
	cp ../models/iosac2.5d_params.h.ah ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.ah ../src/init_user.cu
	cp ../models/boundary.cu.ah ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu


bw:
	cp ../models/iosac2.5d_params.h.bw ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.bw ../src/init_user.cu
	cp ../models/boundary.cu.bw ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu

rt:
	cp ../models/iosac2.5d_params.h.rt ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.rt ../src/init_user.cu
	cp ../models/boundary.cu.rt ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu

kh:
	cp ../models/iosac2.5d_params.h.kh ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.kh ../src/init_user.cu
	cp ../models/boundary.cu.kh ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu

ot:
	cp ../models/iosac2.5d_params.h.ot ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.ot ../src/init_user.cu
	cp ../models/boundary.cu.ot ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu
kink1:
	cp ../models/iosac2.5d_params.h.kink1 ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.kink1 ../src/init_user.cu
	cp ../models/boundary.cu.kink1 ../src/boundary.cu
	cp ../models/usersource.cu.kink1 ../src/usersource.cu

kink2d:
	cp ../models/iosac2.5d_params.h.kink2d ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.kink2d ../src/init_user.cu
	cp ../models/boundary.cu.kink2d ../src/boundary.cu
	cp ../models/usersource.cu.kink2d ../src/usersource.cu

alf:
	cp ../models/iosac2.5d_params.h.alf ../include/iosac2.5dparams.h
	cp ../models/init_user.cu.alf ../src/init_user.cu
	cp ../models/boundary.cu.alf ../src/boundary.cu
	cp ../models/usersource.cu.default ../src/usersource.cu

cleanobjs:
	rm -f $(BASICEXOBJS) ../include/*.cuh

cleanbin:
	rm -f $(ALLBIN)


clean:	cleanobjs cleanbin

cleanall:	cleanobjs cleanbin

all: allsac
	$(NVCC) -g $(LIBS)  $(incDirs) -o $(Bin)/iosac $(ALLOBJS) $(MYLIBS)

allsac:$(ALLOBJS)
	$(NVCC) -g $(LIBS)  $(incDirs) -o $(Bin)/iosac $(ALLOBJS) $(MYLIBS)
#@# Dependency rules follow -----------------------------


#ioshalloww:$(BASICEXOBJS)
#	$(NVCC) --allow-multiple-defintion -g $(LIBS)  $(incDirs) -o $(Bin)/ioshalloww $(BASICEXOBJS) $(MYLIBS)
sac:$(BASICEXOBJS)
	rm iosac2.5d.o
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o ../obj/iosac2.5d.o -c iosac2.5d.cpp
	$(NVCC) -g $(LIBS)  $(incDirs) -o $(Bin)/iosac $(ALLOBJS) $(MYLIBS)




#IoTest.o: IoTest.f90
#	$(FF) -o $@ -c $<
#step.o: step.cu
#	$(NVCC) -o $@  --device-emulation -c $<
#step.o: step.cu
#	$(NVCC) -o $@ --ptxas-options=-v -c $<
hyperdifvisc1.o: hyperdifvisc1.cu
	sed 's/MODID/hdv1/g' < gradops.cu > ../include/gradops_hdv1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifviscmax.o: hyperdifviscmax.cu
	sed 's/MODID/hdvmax/g' < gradops.cu > ../include/gradops_hdvmax.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifrhosource1.o: hyperdifrhosource1.cu
	sed 's/MODID/hdr1/g' < gradops.cu > ../include/gradops_hdr1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifesource1.o: hyperdifesource1.cu
	sed 's/MODID/hde1/g' < gradops.cu > ../include/gradops_hde1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsource1.o: hyperdifbsource1.cu
	sed 's/MODID/hdb1/g' < gradops.cu > ../include/gradops_hdb1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsourcene1.o: hyperdifbsourcene1.cu
	sed 's/MODID/hdbne1/g' < gradops.cu > ../include/gradops_hdbne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsource1.o: hyperdifmomsource1.cu
	sed 's/MODID/hdm1/g' < gradops.cu > ../include/gradops_hdm1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsourcene1.o: hyperdifmomsourcene1.cu
	sed 's/MODID/hdmne1/g' < gradops.cu > ../include/gradops_hdmne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
centdiff1.o: centdiff1.cu
	sed 's/MODID/cd1/g' < gradops.cu > ../include/gradops_cd1.cuh
	sed 's/MODID/cd1/g' < dervfields.cu > ../include/dervfields_cd1.cuh
	sed 's/MODID/cd1/g' < usersource.cu > ../include/usersource_cd1.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
centdiff2.o: centdiff2.cu
	sed 's/MODID/cd2/g' < gradops.cu > ../include/gradops_cd2.cuh
	sed 's/MODID/cd2/g' < dervfields.cu > ../include/dervfields_cd2.cuh
	sed 's/MODID/cd2/g' < usersource.cu > ../include/usersource_cd2.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
computedervfields.o: computedervfields.cu
	sed 's/MODID/cdf/g' < gradops.cu > ../include/gradops_cdf.cuh
	sed 's/MODID/cdf/g' < dervfields.cu > ../include/dervfields_cdf.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS)   -c $<
boundary.o: boundary.cu
	sed 's/MODID/b/g' < gradops.cu > ../include/gradops_b.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
init.o: init.cu
	sed 's/MODID/i/g' < gradops.cu > ../include/gradops_i.cuh
	sed 's/MODID/i/g' < init_user.cu > ../include/init_user_i.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
update.o: update.cu
	sed 's/MODID/u/g' < gradops.cu > ../include/gradops_u.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
advance.o: advance.cu
	sed 's/MODID/adv/g' < gradops.cu > ../include/gradops_adv.cuh
	sed 's/MODID/adv/g' < dervfields.cu > ../include/dervfields_adv.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
divb.o: divb.cu
	sed 's/MODID/db/g' < gradops.cu > ../include/gradops_db.cuh
	sed 's/MODID/db/g' < dervfields.cu > ../include/dervfields_db.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
.cpp.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<





