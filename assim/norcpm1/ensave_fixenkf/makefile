MPI = YES

VPATH = ../EnKF/:./:TMP:../shared

include make.inc.$(MACH)

PROGS = fixenkf_cice ensave ensave_ice

all: $(PROGS)

ENSAVE_OBJ = qmpi.o\
distribute.o\
nfw.o\
p_ensave.o

ENSAVE_ICE_OBJ = qmpi.o\
nfw.o\
p_ensave_ice.o

FIXCICE_OBJ = qmpi.o\
nfw.o\
eosdat.o\
mod_eosfun.o\
p_fixenkf_cice.o

ensave: $(ENSAVE_OBJ)
	@echo "->ensave"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../../../ensave $(ENSAVE_OBJ) $(LIBS)

ensave_ice: $(ENSAVE_ICE_OBJ)
	@echo "->ensave_ice"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../../../ensave_ice $(ENSAVE_ICE_OBJ) $(LIBS)

fixenkf_cice: $(FIXCICE_OBJ)
	@echo "->fixenkf_cice"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../../../fixenkf_cice $(FIXCICE_OBJ) $(LIBS)

clean:
	@rm -f *.o *.mod TMP/*.* $(PROGS)

%.o: %.F90
	@echo "  $*".F90
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90
	@cd ./TMP; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90

%.o: %.F
	@echo "  $*".F
	@rm -f ./TMP/$*.f
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f
	@cd ./TMP; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f 

%.o: %.c
	@echo "  $*".c
	@cd ./TMP ; $(CC) -c $(CFLAGS) -o $*.o -I.. -I../../EnKF ../../EnKF/$*.c
