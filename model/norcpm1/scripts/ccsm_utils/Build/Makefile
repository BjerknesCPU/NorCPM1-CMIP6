#===============================================================================
# SVN $Id: Makefile 28555 2011-05-20 16:56:22Z fischer $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/scripts/branch_tags/cesm1_0_4_rel_tags/cesm1_0_4_n13_scripts4_110531b/ccsm_utils/Build/Makefile $
#===============================================================================
# Common Makefile: a framework for building all CCSM components and more
#
# Command-line variables
#   MACFILE=<file> ~ the macros definition file to use/include
#   <macro defns>  ~ any macro definitions found in this file or the included 
#                    MACFILE will be over-riden by cmd-line macro definitions
#   MODEL=<model>  ~ a standard macro definition, often found in the included 
#                    MACFILE, used to trigger special compilation flags
#
# Usage examples:
#   % gmake MACFILE=Macros.AIX complib
#===============================================================================

# name of macros file - but default this is generic
ifneq ($(strip $(MACFILE)),)
   MACFILE  := Macros
endif

# Load dependency search path.
dirs := . $(shell cat Filepath)
cpp_dirs := $(dirs)
# Add INCROOT to path for Depends and Include
MINCROOT := 
ifneq ($(strip $(INCROOT)),)
  cpp_dirs += $(INCROOT)
  MINCROOT := $(INCROOT)
endif


# Expand any tildes in directory names. Change spaces to colons.
VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH := $(subst $(space),:,$(VPATH))               

RM    := rm
CP    := cp

#-------------------------------------------------------------------------------
# include the file that provides macro definitions required by build rules
# note: the MACFILE may not be needed for certain goals
#-------------------------------------------------------------------------------

exec_se: $(EXEC_SE) touch_filepath $(CURDIR)/Depends
complib: $(COMPLIB) touch_filepath $(CURDIR)/Depends

# Determine whether to compile threaded or not
ifeq ($(strip $(BUILD_THREADED)),TRUE)
   compile_threaded = true
endif
ifeq ($(strip $(SMP)),TRUE)
   compile_threaded = true
endif

# Machine specific macros file 
# This must be included before any settings are overwritten 
# But must be AFTER any definitions it uses are defined.
# So be careful if moving this either earlier or later in the makefile!!!
  -include $(MACFILE)

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)

ifeq ($(findstring -DLINUX,$(CPPDEFS)),-DLINUX)
  ifneq ($(findstring DBG,$(CPPDEFS)),DBG)
    ifneq ($(findstring Darwin,$(CPPDEFS)),Darwin)
      CFLAGS += -DHAVE_SLASHPROC
    endif
  endif
endif



# Set esmf.mk location with ESMF_LIBDIR having precedent over ESMFMKFILE
CCSM_ESMFMKFILE := undefined_CCSM_ESMFMKFILE
ifneq ($(strip $(ESMFMKFILE)),)
   CCSM_ESMFMKFILE := $(ESMFMKFILE)
endif
ifneq ($(strip $(ESMF_LIBDIR)),)
   CCSM_ESMFMKFILE := $(ESMF_LIBDIR)/esmf.mk
endif

LD := $(FC)

# For compiling and linking with external ESMF.
# If linking to external ESMF library then include esmf.mk to provide the macros:
# ESMF_F90COMPILEPATHS
# ESMF_F90LINKPATHS
# ESMF_F90LINKRPATHS
# ESMF_F90ESMFLINKLIBS
ifeq ($(USE_ESMF_LIB), TRUE)
  include $(CCSM_ESMFMKFILE)
  FFLAGS  += $(ESMF_F90COMPILEPATHS)
  LDFLAGS += $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# System libraries (netcdf, mpi, pnetcdf, esmf) 
ifeq ($(strip $(SLIBS)),)
   SLIBS := -L$(LIB_NETCDF) -lnetcdf 
else
   SLIBS += -L$(LIB_NETCDF) -lnetcdf 
endif
ifneq ($(strip $(LIB_PNETCDF)),)
   SLIBS += -L$(LIB_PNETCDF) -lpnetcdf
endif
ifneq ($(strip $(LAPACK_LIBDIR)),)
   SLIBS += -L$(LAPACK_LIBDIR) -llapack -lblas
endif
ifneq ($(strip $(LIB_MPI)),)
   ifeq ($(strip $(MPI_LIB_NAME)),)
      SLIBS += -L$(LIB_MPI) -lmpi
   else  
      SLIBS += -L$(LIB_MPI) -l$(MPI_LIB_NAME)
   endif
endif

# User libraries (mct, mpeu, pio)
ifeq ($(strip $(ULIBS)),)
   ifneq ($(strip $(LIBROOT)),)
     ifeq ($(USE_MPISERIAL),TRUE)
       ULIBS += -L$(LIBROOT) -lcsm_share -lmct -lmpeu -lpio -lmpi-serial 
     else
       ULIBS += -L$(LIBROOT) -lcsm_share -lmct -lmpeu -lpio 
     endif
   endif
endif

# Component libraries - hard-wired for current ccsm components
ifeq ($(strip $(CLIBS)),)
   ifneq ($(strip $(LIBROOT)),)
     CLIBS  += -L$(LIBROOT) -latm -llnd -lice -locn -lglc
   endif
endif

# Machine/Macros stuff to appear last on the link step
ifeq ($(strip $(MLIBS)),)
     MLIBS  :=
endif

#------------------------------------------------------------------------------
# Drive configure scripts for support libraries (mct and pio)
#------------------------------------------------------------------------------
FORCE:

configure: FORCE
	$(CONFIG_SHELL) ./configure $(CONFIG_ARGS) 


#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------

touch_filepath: 
	touch $(CURDIR)/Filepath

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends

SOURCES := $(shell cat Srcfiles)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))
INCS    := $(foreach dir,$(cpp_dirs),-I$(dir)) 

CURDIR := $(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Deppath
	$(CASETOOLS)/mkDepends Deppath Srcfiles > $@

$(CURDIR)/Deppath: $(CURDIR)/Filepath
	$(CP) -f $(CURDIR)/Filepath $@
	@echo "$(MINCROOT)" >> $@

$(CURDIR)/Srcfiles: $(CURDIR)/Filepath
	$(CASETOOLS)/mkSrcfiles > $@

$(CURDIR)/Filepath:
	@echo "$(VPATH)" > $@


#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

db_files:
	@echo " "
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* INCS    := $(INCS)"
	@echo "* OBJS    := $(OBJS)"
db_flags:
	@echo " "
	@echo "* cc      := $(CC)  $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(FC)  $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC)  $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR)"

#-------------------------------------------------------------------------------
# Rules used for the tests run by "configure -test"
#-------------------------------------------------------------------------------

test_fc: test_fc.o
	$(LD) -o $@ test_fc.o $(LDFLAGS)
test_nc: test_nc.o
	$(LD) -o $@ test_nc.o -L$(LIB_NETCDF) -lnetcdf $(LDFLAGS)
test_mpi: test_mpi.o
	$(LD) -o $@ test_mpi.o $(LDFLAGS)
test_esmf: test_esmf.o
	$(LD) -o $@ test_esmf.o $(LDFLAGS)

#-------------------------------------------------------------------------------
# build rules: MACFILE, cmd-line, or env vars must provide the needed macros
#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .F .f90 .c .o

ifeq ($(strip $(ULIBDEP)),)
   ifneq ($(strip $(LIBROOT)),)
     ULIBDEP := $(LIBROOT)/libmct.a
     ULIBDEP += $(LIBROOT)/libcsm_share.a
     ULIBDEP += $(LIBROOT)/libpio.a
     ULIBDEP += $(LIBROOT)/libatm.a
     ULIBDEP += $(LIBROOT)/libice.a
     ULIBDEP += $(LIBROOT)/liblnd.a
     ULIBDEP += $(LIBROOT)/libocn.a
     ULIBDEP += $(LIBROOT)/libglc.a
     ifeq ($(USE_MPISERIAL),TRUE)
       ULIBDEP += $(LIBROOT)/libmpi-serial.a
     endif
   endif
endif

$(EXEC_SE): $(OBJS) $(ULIBDEP)
	$(LD) -o $(EXEC_SE) $(OBJS) $(CLIBS) $(ULIBS) $(SLIBS) $(MLIBS) $(LDFLAGS)

$(COMPLIB): $(OBJS)
	$(AR) -r $(COMPLIB) $(OBJS)

.c.o:
	$(CC) -c $(INCLDIR) $(INCS) $(CFLAGS)  $<

.F.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $<

.f90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<

.F90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<

mostlyclean:
	$(RM) -f *.f *.f90 

clean:
	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS)

realclean:
	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS) $(EXEC_SE)

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), mostlyclean)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), realclean)
    -include $(CURDIR)/Depends
endif
endif
endif
endif
endif
