SHELL = /bin/bash # For our sanity!

#-----
# Compiler configuration

# pre-processing.
CPP = cpp -P -C -traditional
CPPFLAGS = -D__Linux -DLAPACK -DFFT_DEFAULT -DPOINTER8 -DMAXMEM='$(MAXMEM)' -D_VCS_VER='$(VCS_VERSION)' $(WORKING_DIR_CHANGES) 

# use compiler with perl scripts to avoid cascade compilation.
compiler = INTEL-ifort9-on-LINUX

# fortran compiler and flags.
FC = gfortran
FFLAGS = -m64 -fcray-pointer
FNEWFLAGS = 
F90FLAGS = -ffree-line-length-none

# tell fortran compiler to search for module files in the dest directory.
MODULEFLAG = -J $(DEST) -I $(SRC)

# c compiler and flags.
CC = gcc
CFLAGS = -m64 -O2

# linker, linker flags and libraries.
LD = gfortran
LDFLAGS = -m64
LFLAGS = -L~/local/lib/fftw-3.2/64/gfortran/lib -lfftw3 -framework Accelerate

# For building neci library.
AR = ar
ARFLAGS = -rcs

#-----
# Directory structure and setup

# Directories containing source files (space separated list).
SRC = src

# Directory in which compiled objects are placed.
DEST = dest

# Directory for compiled executable.
EXE = bin

# Directory for neci libraries.
LIB = lib

# Directory containing scripts to avoid cascade compilation
TOOLS = tools

# Directories that make searches for pre-requisites in addition to the local directory (./).
# Colon separated list.
# We include $(DEST) in this so we can do pre-processing and compilation in 2 steps.
empty :=
space := $(empty) $(empty) # stupid make language...
VPATH = $(subst $(space),:,$(SRC) $(DEST))

# Create output directories if they don't exist.
make_dest := $(shell test -e $(DEST) || mkdir -p $(DEST))
make_exe := $(shell test -e $(EXE) || mkdir -p $(EXE))
make_lib := $(shell test -e $(LIB) || mkdir -p $(LIB))

#-----
# VCS info.

# Get the version control id.  Git only.  Outputs a string.
VCS_VERSION := $(shell set -o pipefail && echo -n \" && ( git log --max-count=1 --pretty=format:%H || echo -n 'Not under version control.' ) 2> /dev/null | tr -d '\r\n'  && echo -n \")

# Test to see if the working directory contains changes.  Git only.  If the
# working directory contains changes (or is not under version control) then
# the _WORKING_DIR_CHANGES flag is set.
WORKING_DIR_CHANGES := $(shell git diff --quiet --cached && git diff --quiet 2> /dev/null || echo -n "-D_WORKING_DIR_CHANGES")

#-----
# Additional info.

MAXMEM = 2048 # RAM available, in MB.  Used by the memory logger.

#-----
# Find source files and resultant object files.

# Source extensions.
EXTS = .F90 .F .C .c

# Source filenames.
find_files = $(wildcard $(1)/*$(2))
FSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.F))
F90SRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.F90))
CSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.C))
cSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.c))
SRCFILES := $(FSRCFILES) $(F90SRCFILES) $(CSRCFILES) $(cSRCFILES)

# Objects (strip path and replace extension of source files with .o).
FOBJ := $(addsuffix .o,$(basename $(notdir $(FSRCFILES))))
F90OBJ := $(addsuffix .o,$(basename $(notdir $(F90SRCFILES))))
COBJ := $(addsuffix .o,$(basename $(notdir $(CSRCFILES))))
cOBJ := $(addsuffix .o,$(basename $(notdir $(cSRCFILES))))

# Full path to all objects.
FOBJ := $(addprefix $(DEST)/, $(FOBJ))
F90OBJ := $(addprefix $(DEST)/, $(F90OBJ))
COBJ := $(addprefix $(DEST)/, $(COBJ))
cOBJ := $(addprefix $(DEST)/, $(cOBJ))
OBJECTS := $(FOBJ) $(F90OBJ) $(COBJ) $(cOBJ)

# Objects for standalone NECI.
# We don't need libstub.
OBJECTS_NECI := $(filter-out %libstub.o,$(OBJECTS))

# Objects for CPMD library.
# We don't need necimain.o, cpmdstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_CPMD := $(filter-out %necimain.o %cpmdstub.o % %init_coul.o %init_could2D.o, $(OBJECTS)) 

# Objects for VASP library.
# We don't need necimain.o, vaspstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_VASP := $(filter-out %necimain.o %vaspstub.o % %init_coul.o %init_could2D.o, $(OBJECTS)) 

#-----
# Dependency files.

# Fortran dependencies.
# We need these before compiling any fortran files.
FDEPEND = $(DEST)/.depend

FDEPEND_EXISTS := $(wildcard $(FDEPEND))

# If the dependency file does not exist, then it is generated.  This 
# pass of make will not have the correct dependency knowledge though,
# so we force an exit.
ifneq ($(FDEPEND_EXISTS),$(FDEPEND))
	TEST_DEPEND = no_depend
	DEPEND_TARGET = $(FDEPEND)
else
	DEPEND_TARGET = depend_target
endif

# C dependencies.
# We don't need these when we first compile, only when we recompile.
# We achieve this (most of the time) by recompiling the C dependencies every
# time we compile.
CDEPEND_FILES = $(COBJ:.o=.d)
cDEPEND_FILES = $(cOBJ:.o=.d)

#-----
# Goals

.PHONY: clean $(DEPEND_TARGET) depend help neci.x

$(EXE)/neci.x : $(TEST_DEPEND) $(OBJECTS_NECI)
	$(MAKE) -B environment_report.o
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS_NECI) $(LFLAGS)

neci.x: 
	$(MAKE) $(EXE)/neci.x

#neci-cpmd.a : $(OBJECTSCPMDLIB)
#	 rm -f environment_report.F90
#	 $(CPP) $(CPPFLAGS) $(SRC)/environment_report.F90 $(DEST)/environment_report.f90
#	 $(FC) $(FFLAGS) $(FNEWFLAGS) $(DEST)/environment_report.f90
#	 $(AR) $(ARFLAGS) $(DEST)/neci-cpmd.a environment_report.o $(OBJECTSCPMDLIB)
#
#neci-vasp.a : $(OBJECTSVASPLIB)
#	 rm -f environment_report.F90
#	 $(CPP) $(CPPFLAGS) $(SRC)/environment_report.F90 $(DEST)/environment_report.f90
#	 $(FC) $(FFLAGS) $(FNEWFLAGS) $(DEST)/environment_report.f90
#	 $(AR) $(ARFLAGS) $(DEST)/neci-vasp.a environment_report.o $(OBJECTSVASPLIB)

clean : 
	  rm -f $(DEST)/*.{f,f90,mod,o,c,x,a,d} $(EXE)/neci.x

# Generate dependency file.
$(DEPEND_TARGET):
	$(TOOLS)/sfmakedepend --file - --silent $(SRCFILES) --objdir \$$\(DEST\) --moddir \$$\(DEST\) > $(FDEPEND)

depend: $(DEPEND_TARGET)

# Force exit if dependency file didn't exist as make didn't pickup the correct dependencies on this pass.
no_depend:
	@echo "The required dependency file did not exist but has now been generated."
	@echo "Please re-run make."
	exit 2

#-----
# Compilation macros (explicit rules)

.SUFFIXES:
.SUFFIXES: $(EXTS) .f .f90

# Compiling fortran source files...
$(DEST)/%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< $@

$(F90OBJ): %.o: %.f90
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(F90FLAGS) $(MODULEFLAG) -c $< -o $@" -provides "$@" -requires "$^"

$(DEST)/%.f: %.F
	@echo CPP: $< $@
	$(CPP) $(CPPFLAGS) $< $@

$(FOBJ): %.o: %.f
	@echo compiling: $< $@
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(FNEWFLAGS) $(MODULEFLAG) -c $< -o $@" -provides "$@" -requires "$^"

# Compiling C source files...
$(COBJ): $(DEST)/%.o: %.C
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ 

$(cOBJ): $(DEST)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ 

# Update C dependency files.
$(cDEP): $(DEST)/%.d: %.c
	$(CC) $(CFLAGS) -MM -MT \$$\(DEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@

$(CDEP): $(DEST)/%.d: %.C
	$(CC) $(CFLAGS) -MM -MT \$$\(DEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@

#-----
# Generic dependencies.
# These are for files which don't depend on anything other than themselves.
# Other files are overridden in the dependency files which are included below.

#-----
# Include dependency files

# $(FDEPEND) will be generated if it doesn't exist.
include $(FDEPEND)
include $(CDEPEND_FILES)
include $(cDEPEND_FILES)
