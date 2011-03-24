#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_suilven
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
	derivsource.o\
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
	
ALLOBJS	=	$(BASICEXOBJS)
ALLBIN	=	$(Bin)/iosac

#use -DUSE_IOME to enable IOME
#CUDACCFLAGS = --ptxas-options=-v -maxrregcount=32 --device-emulation   -DUSE_REAL -DUSE_SAC -DADIABHYDRO -DUSE_VAC -DD1D -DUSE_IOME
#compiling with  -arch sm_13 enables double precision to be used
CUDACCFLAGS = --ptxas-options=-v -arch sm_13  -maxrregcount=32 -DUSE_SAC -DD1D
CCFLAGS = -DUSE_SAC -DD1D

#MYLIBS	=	libAgentSimulation2.a

#@# Targets follow ---------------------------------

##all:	IoSimTestSimple




cleanobjs:
	rm -f $(BASICEXOBJS) *.cuh

cleanbin:
	rm -f $(ALLBIN)


clean:	cleanobjs cleanbin

cleanall:	cleanobjs cleanbin

all: iosac


#@# Dependency rules follow -----------------------------


#ioshalloww:$(BASICEXOBJS)
#	$(NVCC) --allow-multiple-defintion -g $(LIBS)  $(incDirs) -o $(Bin)/ioshalloww $(BASICEXOBJS) $(MYLIBS)
iosac:$(BASICEXOBJS)
	$(NVCC) -g $(LIBS)  $(incDirs) -o $(Bin)/iosac $(BASICEXOBJS) $(MYLIBS)




#IoTest.o: IoTest.f90
#	$(FF) -o $@ -c $<
#step.o: step.cu
#	$(NVCC) -o $@  --device-emulation -c $<
#step.o: step.cu
#	$(NVCC) -o $@ --ptxas-options=-v -c $<
hyperdifvisc1.o: hyperdifvisc1.cu
	sed 's/MODID/hdv1/g' < gradops.cu > gradops_hdv1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifviscmax.o: hyperdifviscmax.cu
	sed 's/MODID/hdvmax/g' < gradops.cu > gradops_hdvmax.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifrhosource1.o: hyperdifrhosource1.cu
	sed 's/MODID/hdr1/g' < gradops.cu > gradops_hdr1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifesource1.o: hyperdifesource1.cu
	sed 's/MODID/hde1/g' < gradops.cu > gradops_hde1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsource1.o: hyperdifbsource1.cu
	sed 's/MODID/hdb1/g' < gradops.cu > gradops_hdb1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsourcene1.o: hyperdifbsourcene1.cu
	sed 's/MODID/hdbne1/g' < gradops.cu > gradops_hdbne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsource1.o: hyperdifmomsource1.cu
	sed 's/MODID/hdm1/g' < gradops.cu > gradops_hdm1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<


hyperdifmomsourcene1.o: hyperdifmomsourcene1.cu
	sed 's/MODID/hdmne1/g' < gradops.cu > gradops_hdmne1.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
derivsource.o: derivsource.cu
	sed 's/MODID/ds/g' < gradops.cu > gradops_ds.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
centdiff1.o: centdiff1.cu
	sed 's/MODID/cd1/g' < gradops.cu > gradops_cd1.cuh
	sed 's/MODID/cd1/g' < dervfields.cu > dervfields_cd1.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
centdiff2.o: centdiff2.cu
	sed 's/MODID/cd2/g' < gradops.cu > gradops_cd2.cuh
	sed 's/MODID/cd2/g' < dervfields.cu > dervfields_cd2.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
computedervfields.o: computedervfields.cu
	sed 's/MODID/cdf/g' < gradops.cu > gradops_cdf.cuh
	sed 's/MODID/cdf/g' < dervfields.cu > dervfields_cdf.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS)   -c $<
boundary.o: boundary.cu
	sed 's/MODID/b/g' < gradops.cu > gradops_b.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
init.o: init.cu
	sed 's/MODID/i/g' < gradops.cu > gradops_i.cuh
	sed 's/MODID/i/g' < init_user.cu > init_user_i.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
update.o: update.cu
	sed 's/MODID/u/g' < gradops.cu > gradops_u.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
advance.o: advance.cu
	sed 's/MODID/adv/g' < gradops.cu > gradops_adv.cuh
	sed 's/MODID/adv/g' < dervfields.cu > dervfields_adv.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
divb.o: divb.cu
	sed 's/MODID/db/g' < gradops.cu > gradops_db.cuh
	sed 's/MODID/db/g' < dervfields.cu > dervfields_db.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
.cpp.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<





