#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_suilven
### End User configurable options ###







BASICEXOBJS	=\
	iosac2.5d.o\
	initialisation.o\
	dxroutines.o\
	hyperdifmomsource.o\
	hyperdifesource.o\
	hyperdifbsource.o\
	hyperdifrhosource.o\
	hyperdifvisc.o\
	derivsource.o\
	derivcurrent1.o\
	derivcurrent2.o\
	predictor.o\
	corrector.o\
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
#CUDACCFLAGS = --ptxas-options=-v -maxrregcount=32 --device-emulation   -DUSE_REAL -DUSE_SAC -DADIABHYDRO -DUSE_VAC
#compiling with  -arch sm_13 enables double precision to be used
CUDACCFLAGS = --ptxas-options=-v -arch sm_13  -maxrregcount=32 -DADIABHYDRO
CCFLAGS = -DADIABHYDRO

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



hyperdifvisc.o: hyperdifvisc.cu
	sed 's/MODID/hdv/g' < gradops.cu > gradops_hdv.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifrhosource.o: hyperdifrhosource.cu
	sed 's/MODID/hdr/g' < gradops.cu > gradops_hdr.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifesource.o: hyperdifesource.cu
	sed 's/MODID/hde/g' < gradops.cu > gradops_hde.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsource.o: hyperdifbsource.cu
	sed 's/MODID/hdb/g' < gradops.cu > gradops_hdb.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsource.o: hyperdifmomsource.cu
	sed 's/MODID/hdm/g' < gradops.cu > gradops_hdm.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
derivsource.o: derivsource.cu
	sed 's/MODID/ds/g' < gradops.cu > gradops_ds.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
derivcurrent1.o: derivcurrent1.cu
	sed 's/MODID/dc1/g' < gradops.cu > gradops_dc1.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
derivcurrent2.o: derivcurrent2.cu
	sed 's/MODID/dc2/g' < gradops.cu > gradops_dc2.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
predictor.o: predictor.cu
	sed 's/MODID/pre/g' < gradops.cu > gradops_pre.cuh
	sed 's/MODID/pre/g' < dervfields.cu > dervfields_pre.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS)   -c $<
corrector.o: corrector.cu
	sed 's/MODID/cor/g' < gradops.cu > gradops_cor.cuh
	sed 's/MODID/cor/g' < dervfields.cu > dervfields_cor.cuh
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
boundary.o: boundary.cu
	sed 's/MODID/b/g' < gradops.cu > gradops_b.cuh
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
init.o: init.cu
	sed 's/MODID/i/g' < gradops.cu > gradops_i.cuh
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





