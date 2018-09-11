######################################################################
#
#  GDPME Makefile for Linux (Intel compiler), netcdf interface included
#
#  Use "make release" for optimized code
#      "make debug"   for debug code
#      "make clean"   for removing object files and executables
#
#######################################################################

####  NETCDF PATH  #####
       NETCDF ?= /usr/local 
NETCDF_INCDIR ?= $NETCDF/include   
NETCDF_LIBDIR ?= $NETCDF/lib

.SUFFIXES: .obj .o .f .h .mod
           FC ?= ifort -nologo -O3 # -C -gen-interface -warn interface -traceback

#### Flags with a "D" at the end are for the debug target  #####

#####  COMPILING FLAGS  #####
       FFLAGS ?= -c -fpp

#####  LINKING FLAGS  #####
LFLAGS = -o $(rca_target)

OBJECTS =  rca.o rca01.o rca02.o rca03.o rca04.o rca05.o rca06.o \
	rca07.o rca08.o rca09.o rca10.o rca11.o rca12.o \
	rca13.o rcabyss.o rcadwrit.o rcaexp.o rcaexp1.o \
	rcaexps.o rcalf1.o rcalfs.o rcamess.o rcamprof.o \
	rcaprnt.o rcasplt.o rcasplt1.o rcasplts.o pco2_water.o calcium_carbonate.o \
	$(tuner_oFILE) \
	sinter.o iunitcheck.o \
	fmter.o $(nc_oFILE)

###########################
########  TARGETS  ########
###########################

release: $(OBJECTS) RCACM $(nc_fFILE) 
	$(FC) $(LFLAGS) -I$(NETCDF_INCDIR) $(OBJECTS) -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf 
	@echo $(rca_target) is now up-to-date

clean :
	rm -rf *.o $(rca_target) RCAF* 

###########################
#####  DEPENDENCIES  ######
###########################

OBJECTS: RCACM $(nc_fFILE) netcdf.mod $(nc_hFILE)

OBJECTSD: RCACM $(nc_fFILE)

###########################
########  RULES  ##########
###########################

.f.obj:                   # for release
	$(FC) $(FFLAGS) -I$(NETCDF_INCDIR) -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf $<

.f.o:                     # for release
	$(FC) $(FFLAGS) -I$(NETCDF_INCDIR) -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf $<
