#Makefile for StarSimulation IoSimTestSimple

##### User configurable options #####
include make_inputs_mpi
### End User configurable options ###



	
ALLOBJS	=	$(obj)/$(BASICEXOBJS)
ALLBIN	=	$(Bin)/smaugutils





#@# Targets follow ---------------------------------

spicule4:
	cp ../solarmodels/iosmaugparams.h.spicule4_nob ../include/iosmaugparams.h
	cp ../solarmodels/initialisation_user.h.spicule4_b ../include/initialisation_user.h

fluxtube1:
	cp ../models/iosmaugparams.h.fluxtube1 ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

spiculejet1:
	cp ../models/iosmaugparams.h.spiculejet1 ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

spicule1:
	cp ../models/iosmaugparams.h.spicule1 ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.spicule1 ../include/initialisation_user.h
spicule1_mpi:
	cp ../models/iosmaugparams.h.spicule1_mpi ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.spicule1_mpi ../include/initialisation_user.h
sedov:
	cp ../models/iosmaugparams.h.sedov ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

bach:
	cp ../models/iosmaugparams.h.bach ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

ah:
	cp ../models/iosmaugparams.h.ah ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

bw:
	cp ../models/iosmaugparams.h.bw ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
rt:
	cp ../models/iosmaugparams.h.rt ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
kh:
	cp ../models/iosmaugparams.h.kh ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
ot:
	cp ../models/iosmaugparams.h.ot ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
ot_mpi:
	cp ../models/iosmaugparams.h.ot_mpi ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
kink1:
	cp ../models/iosmaugparams.h.kink1 ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
kink2d:
	cp ../models/iosmaugparams.h.kink2d ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
alf:
	cp ../models/iosmaugparams.h.alf ../include/iosmaugparams.h
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h
cleanobjs:
	rm -f $(BASICEXOBJS)

cleanbin:
	rm -f $(ALLBIN)


clean:	cleanobjs cleanbin

cleanall:	cleanobjs cleanbin

all: smaugutils
	$(CC) -g $(LIBS)  $(incDirs) -o $(Bin)/smaugutils $(ALLOBJS) $(MYLIBS)


#@# Dependency rules follow -----------------------------



smaugutils:$(BASICEXOBJS)
	rm iosmaugutils.o
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs)  -o ../obj/iosmaugutils.o -c iosmaugutils.cpp
	$(CC) $(LINKFLAGS)  $(LIBS)  $(incDirs) -o $(Bin)/smaugutils $(ALLOBJS) $(MYLIBS)




.cpp.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)  -o $@ -c $<

.c.o:
	$(CC) $(C_FLAGS) $(CCFLAGS) $(incDirs) $(LINKFLAGS)    -o $@ -c $<





