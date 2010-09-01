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

	
	
#	../../IoGenericSimulationLib.o\
#	../../../IoSimulation/IoGenericSteerSimulation.o\
#	../../../IoSimulation/stdsoap2.o\
#	../../../IoSimulation/soapC.o\
#	../../../IoSimulation/soapClient.o
	
ALLOBJS	=	$(BASICEXOBJS)
ALLBIN	=	$(Bin)/iosac
#CUDACCFLAGS = --ptxas-options=-v -maxrregcount=32 --device-emulation   -DUSE_REAL -DUSE_SAC -DADIABHYDRO
#compiling with  -arch sm_13 enables double precision to be used
CUDACCFLAGS = --ptxas-options=-v -arch sm_13  -maxrregcount=32 
CCFLAGS = 

#MYLIBS	=	libAgentSimulation2.a

#@# Targets follow ---------------------------------

##all:	IoSimTestSimple




cleanobjs:
	rm -f $(BASICEXOBJS)

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
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifrhosource.o: hyperdifrhosource.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifesource.o: hyperdifesource.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifbsource.o: hyperdifbsource.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
hyperdifmomsource.o: hyperdifmomsource.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
derivsource.o: derivsource.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
derivcurrent1.o: derivcurrent1.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
derivcurrent2.o: derivcurrent2.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
predictor.o: predictor.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
corrector.o: corrector.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
boundary.o: boundary.cu
	$(NVCC) -o $@    $(CUDACCFLAGS) -c $<
init.o: init.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
update.o: update.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
advance.o: advance.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
divb.o: divb.cu
	$(NVCC) -o $@     $(CUDACCFLAGS) -c $<
.cpp.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<





