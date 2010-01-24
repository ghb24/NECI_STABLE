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

# flag to tell fortran compiler where to search for module files.
MODULEFLAG = -J

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

# Directories in which compiled objects are placed.
DEST = dest/opt
# REAL (molecular and gamma-point) objects
GDEST = $(DEST)/real
# COMPLEX (k-point) objects
KDEST = $(DEST)/cmplx

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
make_gdest := $(shell test -e $(GDEST) || mkdir -p $(GDEST))
make_kdest := $(shell test -e $(KDEST) || mkdir -p $(KDEST))
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
FOBJ := $(addprefix $(GDEST)/, $(FOBJ))
F90OBJ := $(addprefix $(GDEST)/, $(F90OBJ))
COBJ := $(addprefix $(GDEST)/, $(COBJ))
cOBJ := $(addprefix $(GDEST)/, $(cOBJ))
OBJECTS := $(FOBJ) $(F90OBJ) $(COBJ) $(cOBJ)

# Objects for standalone NECI.
# We don't need libstub.
OBJECTS_NECI := $(filter-out %libstub.o,$(OBJECTS))

# Objects for CPMD library.
# We don't need necimain.o, cpmdstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_RCPMD := $(filter-out %necimain.o %cpmdstub.o %init_coul.o %init_could2D.o,$(OBJECTS)) 

# Objects for VASP library.
# We don't need necimain.o, vaspstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_RVASP := $(filter-out %necimain.o %vaspstub.o % %init_coul.o %init_could2D.o, $(OBJECTS)) 

#-----
# Dependency files.

# Fortran dependencies.
# We need these before compiling any fortran files.
FRDEPEND = $(GDEST)/.depend
FCDEPEND = $(KDEST)/.depend
FDEPEND = $(FRDEPEND) $(FCDEPEND)

# C dependencies.
# We don't need these when we first compile, only when we recompile.
# We achieve this (most of the time) by recompiling the C dependencies every
# time we compile.
CDEPEND_FILES = $(COBJ:.o=.d)
cDEPEND_FILES = $(cOBJ:.o=.d)

#-----
# Goals

.PHONY: clean depend help neci.x

$(EXE)/neci.x : $(OBJECTS_NECI)
	rm $(GDEST)/environment_report.* && $(MAKE) $(GDEST)/environment_report.o
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS_NECI) $(LFLAGS)

neci.x: 
	$(MAKE) $(EXE)/neci.x

$(LIB)/gneci-cpmd.a : $(OBJECTS_RCPMD)
	rm $(GDEST)/environment_report.* && $(MAKE) $(GDEST)/environment_report.o
	$(AR) $(ARFLAGS) $@ $^

$(LIB)/gneci-vasp.a : $(OBJECTS_RVASP)
	rm $(GDEST)/environment_report.* && $(MAKE) $(GDEST)/environment_report.o
	$(AR) $(ARFLAGS) $@ $^

clean : 
	  rm -f $(GDEST)/*.{f,f90,mod,o,c,x,a,d} $(EXE)/neci.x $(LIB)/*.a

# Generate dependency files.
$(FRDEPEND):
	$(TOOLS)/sfmakedepend --file - --silent $(SRCFILES) --objdir \$$\(GDEST\) --moddir \$$\(GDEST\) > $(FRDEPEND)

$(FCDEPEND):
	$(TOOLS)/sfmakedepend --file - --silent $(SRCFILES) --objdir \$$\(KDEST\) --moddir \$$\(KDEST\) > $(FCDEPEND)

depend: $(FRDEPEND) $(FCDEPEND)

#-----
# Compilation macros (explicit rules)

.SUFFIXES:
.SUFFIXES: $(EXTS) .f .f90

# Compiling fortran source files...
$(GDEST)/%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< $@

$(F90OBJ): %.o: %.f90
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(F90FLAGS) $(MODULEFLAG) $(dir $@) -I $(SRC) -c $< -o $@" -provides "$@" -requires "$^"

$(GDEST)/%.f: %.F
	$(CPP) $(CPPFLAGS) $< $@

$(FOBJ): %.o: %.f
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(FNEWFLAGS) $(MODULEFLAG) $(dir $@) -I $(SRC) -c $< -o $@" -provides "$@" -requires "$^"

# Compiling C source files...
$(COBJ): $(GDEST)/%.o: %.C
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ 

$(cOBJ): $(GDEST)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ 

# Update C dependency files.
$(cDEPEND_FILES): $(GDEST)/%.d: %.c
	$(CC) $(CFLAGS) -MM -MT \$$\(DEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@

$(CDEPEND_FILES): $(GDEST)/%.d: %.C
	$(CC) $(CFLAGS) -MM -MT \$$\(DEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@

#-----
# Include dependency files

# Dependency files will be generated if they don't exist.
include $(FDEPEND)
include $(CDEPEND_FILES)
include $(cDEPEND_FILES)
