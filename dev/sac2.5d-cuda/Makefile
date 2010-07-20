#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_suilven
### End User configurable options ###







BASICEXOBJS	=\
	iosac2.5d.o\
	initialisation.o\
	dxroutines.o\
	derivsource.o\
	derivcurrent.o\
	predictor.o\
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

derivsource.o: derivsource.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
derivcurrent.o: derivcurrent.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
predictor.o: predictor.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
boundary.o: boundary.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
init.o: init.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
update.o: update.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
advance.o: advance.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
divb.o: divb.cu
	$(NVCC) -o $@ --ptxas-options=-v -c $<
.cpp.o:
	$(CC) $(C_FLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<





