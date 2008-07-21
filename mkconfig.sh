#!/bin/sh
#      INTEGER*8 PLAN

#--------------------------------------------------------------------#
#Create Makefile for neci.x program.                                 #
#--------------------------------------------------------------------#
# Options for CPPFLAG                                                #
#  -DLAPACK  Use LAPACK routine, you must have this option and       #
#            the LAPACK library                                      #
#  -DFFT_DEFAULT Use default FFT routine in the code (from Goedecker)#
#  -DFFT_ESSL    Use FFT ESSL (need ESSL library)                    #
#  -DFFT_HP      Use FFT HP (Hewlet-Packard)                         #
#  -DFFT_YMP     Use CRAY YMP, C94, T90 FFT routines                 #
#  -DFFT_T3D     Use CRAY T3D or T3E FFT routines                    #
#  -DPARALLEL    For PARALLEL computers, you have to indicate        #
#                the parallel library with :                         #
#  -DMP_LIBRARY=__MPI   Message Passing Interface library            #
#            or __SHMEM SHMEM library (CRAY)                         #
#            or __MPL                                                #
#  -D__VECTOR    For vectorial computers                             #
#  -DPOINTER8    If the addresses are coded in INTEGER*8             #
#                (for LOC and MALLOC)                                #
#  -DMALLOC8     If the argument of MALLOC is INTEGER*8              #
#  -D__NOINT8    If the compiler does not use INTEGER*8 or           #
#                if the function BTEST and IBTEST do not like it     #
#  -DJOBLIMIT    For batch job, try to know the remaining time       #
#                (For CRAY only)                                     #
#--------------------------------------------------------------------#
#Create Makefile for neci.x program.

#Display all configurations
Message () {
cat << END >&2
Configure  [options] Name_of_Configuration
where the name of the configuration is:
END
\ls CONFIGURE >&2
cat << END >&2

Use Configure in the directory where the SOURCE FILES are.
Ex: Configure SUN > Makefile
See the description of options with Configure -help

END
}


#Help
Help () {
Message
cat << END >&2

Description of options:
   -debug    (-d) Give special options for debugging
   -help     (-h) Print this message
   -makefile (-m) Create the file Makefile in DEST directory
                  instead of using standard output
   -verbose  (-v) Display the name of files for the dependencies
   -VAR=value     Use the value of VAR in the makefile
            Ex: -DEST=destdir specifies where the compilation will be
   You can use:
     FFLAGS   Fortran flags
     LFLAGS   Link editor flags
     CFLAGS   C compiler flags
     CPP      Preprocessor program
     CPPFLAGS CPP flags
     FC       Fortran compiler
     LD       Link editor
     AR       Archive library maintainer
     RANLIB   Converts archives to random libraries
     SRC      Source files directory (default=.)
     DEST     Object files directory (default=.)
              Use to write irat.inc in the right directory
     BIN      Use to put the neci.x executable (default=.)
Note: if you want to compile neci.x in a different directory
      than the source files, use:
         Configure -DEST=destdir SUN > destdir/Makefile
      or Configure -makefile -DEST=destdir SUN
      and then cd destdir; make
      In this case, SRC is equal to source_dir_pathname
      except if you specify SRC.
END
}

#No argument: Error
if [ $# -eq 0 ]; then
  Message
 exit 2
fi

#Is it help option or debug option
opt=0
i=1
PreprocSep=1
while [ $i -le $# ];
do
  #Select the i-th argument
  eval option=\$$i
  case $option in
    -debug|-d)
      debug=1
      ;;
    -help|-h)
      Help
      exit 0
      ;;
    -makefile|-m)
      makefile=1
      ;;
    -verbose|-v)
      verbose=1
      ;;
    -DEST=*)
      opt=1
      DEST=`echo $option | cut -c2- | cut -d= -f 2`
      ;;
    -*=*)
      opt=1
      ;;
    -*)
      echo "Unknow option: $option" >&2
      Message
      exit 1
      ;;
    *)
      Configuration=$option
      ;;
  esac
  i=`expr $i + 1`
done

#No configuration given
if [ -z "$Configuration" ]; then
  echo "Configure: No configuration name" >&2
  Message
  exit 2
fi

#Check if the Configuration does exist.
if [ -f CONFIGURE/${Configuration} ]; then
   . CONFIGURE/${Configuration}
else
    echo "mkconfig.sh: Unknown configuration '${Configuration}'" >&2
    Message
    exit 2
fi

# File descriptor usage:
# 0 standard input
# 1 file creation (standard output or makefile)
# 2 errors and warnings
# 3 Makefile if $makefile or &1

echo "Default configuration for $Configuration." >&2
if [ $debug ]; then
   echo "Debug option is used." >&2
fi

#Default DEST is .
DEST=${DEST:-'.'}

#Check if DEST directory exists.
if [ ! -d ${DEST} ]; then
   echo "Warning: the directory ${DEST} does not exist" >&2
fi
if [ ! -w ${DEST} ]; then
  echo "Warning: the directory ${DEST} has no write permission" >&2
fi

#Create Makefile if output
if [ $makefile ]; then
  echo "Creation of Makefile: ${DEST}/Makefile" >&2
  exec 3>${DEST}/Makefile
else
  echo "Use standard output for the Makefile." >&2
  exec 3>&1
fi

if [ -z "$LDLIB" ] ; then
  LDLIB=${LD}
  LDLIBFLAGS=${LDFLAGS}
fi

cat << END >&3
#----------------------------------------------------------------------------
# Makefile for neci.x (plane wave electronic calculation)
# Configuration: ${Configuration}
# Creation of Makefile: `date '+%b %e %Y'`
# on `uname -a`
# Author: `who am i | cut -d' ' -f 1`
#----------------------------------------------------------------------------
#
SHELL = /bin/bash
#
#--------------- Default Configuration for $Configuration ---------------
SVNVER := \$(shell grep Revision <(svn info || echo 'Revision:"not under svn"') |sed "s/Revision://")
MAXMEM := 1024 # RAM available, in MB.
compiler = ${compiler}
SRC  = .
DEST = ${DEST}
BIN  = .
FFLAGS = ${FFLAGS}
FNEWFLAGS = ${FNEWFLAGS}
F90FLAGS = ${F90FLAGS}
LFLAGS = ${LFLAGS}
LDFLAGS = ${LDFLAGS}
CFLAGS = ${CFLAGS}
CPP = ${CPP}
CPPFLAGS = ${CPPFLAGS}
FC = ${FC}
LD = ${LD}
AR = ${AR}
FCD = ${FCD}
LDD = ${LDD}
ARD = ${ARD}
LDLIB= ${LDLIB}
LDLIBFLAGS= ${LDLIBFLAGS}
GLOBALS= ${GLOBALS}
#----------------------------------------------------------------------------
END
#There is Personal Configuration.
if [ $opt -eq 1 ]; then
  printf "Personal configurations..." >&2

  #Check the options and personal variables 
  cat << END >&3

#----------------------------------------------------------------------------
# Personal Configuration
#----------------------------------------------------------------------------
END

  i=1
  while [ $i -le $# ];
  do
    eval option=\$$i
    case $option in
    -DEST=*)
      #Special case for DEST (only used for makefile)
      ;;
    -*=*)
      var=`echo $option | cut -c2- | cut -d= -f 1`
      val=`echo $option | cut -c2- | cut -d= -f 2`
      echo "$var = $val" >&3
      eval $var='$val'
      ;;
    *)
      ;;
    esac
    i=`expr $i + 1`
  done

  #Default personal option
  if [ $DEST != "." ]; then
    if [ -z "$SRC" ]; then
      SRC=`pwd`
      echo "SRC = $SRC" >&3
    fi
  fi
  
  if [ -n "$SRC" ]; then
    echo "FC = $FC -I\$(SRC)" >&3
  fi
  cat << END >&3
#----------------------------------------------------------------------------
# End of Personal Configuration
#----------------------------------------------------------------------------
END
  echo "done." >&2
  if [ -n "$BIN" ]; then
    echo "neci.x executable will be in: $BIN" >&2
  fi
  
fi #End of Personal Configuration.

#Default SRC is .
SRC=${SRC:-'.'}

echo "The source directory is: $SRC" >&2
echo "The object directory is: $DEST" >&2
  

#Now we change irat if necessary
if [ -f ${DEST}/irat.inc ]; then
  iratfile="${DEST}/irat.inc"
else
  iratfile="irat.inc"
fi
irat=`awk -F= '/PARAMETER/ { print substr($2,1,1) }' $iratfile`
if [ ${irat} != ${IRAT} -o ! -f ${DEST}/irat.inc ]; then
  echo "We change the file $iratfile " \
       "[old irat=${irat}, new irat=${IRAT}]" >&2
  if [ -w /dev/null ]; then
    exec 4> /dev/null
  else
    exec 4>&2
  fi
  ed -s $iratfile << END >&4
/PARAMETER
d
i
      PARAMETER(IRAT=${IRAT})
.
w ${DEST}/irat.inc
q
END
  echo "The file ${DEST}/irat.inc is created." >&2
else
  echo "The file ${DEST}/irat.inc is correct." >&2
fi

#Now DEST is always ./
DEST='.'

#Default SRC
SRC=${SRC:-'.'}

#Print OBJECTS
if [ -f OBJECTS ]; then
  printf "Add OBJECTS (object files)..." >&2
  cat OBJECTS >&3
  echo "done."  >&2
else
  printf "\nThe file OBJECTS does not exist.\n" >&2
  echo "The file OBJECTS does not exist" >&2
  exit 1
fi

#Include files
printf "Add INCFILES (include files)..." >&2
IncludeFile=`ls *.inc`
printf "INCFILES = " >&3
i=0
for name in ${IncludeFile}
do
  if [ $i -eq 6 ]; then
    printf "\\\\\n           " >&3
    i=0
  fi
  i=`expr $i + 1`
  printf "%s " ${name} >&3
done
printf "\n\n" >&3
echo "done." >&2

printf "Add explicit rules..." >&2
if [ $PreprocSep ] ; then
cat << END >&3
#----------------------------------------------------------------------------
# Compile neci.x
#----------------------------------------------------------------------------
#MAIN_OBJS=necimain.o
clean : 
	  rm -f \$(DEST)/*.f90 \$(DEST)/*.f \$(DEST)/*.o \$(DEST)/*.mod

neci.x : necimain.o \$(SRC)/timetag.F \$(OBJECTS)
	 rm -f timetag.f
	 \$(CPP) \$(CPPFLAGS) \$(SRC)/timetag.F \$(DEST)/timetag.f
	 \$(FC) \$(FFLAGS) \$(FNEWFLAGS) \$(DEST)/timetag.f
	 rm -f neci.x
	 if [ \$(BIN) != '.' ]; then ln -s \$(BIN)/neci.x neci.x; fi
	 \$(LD) -o \$(BIN)/neci.x necimain.o timetag.o \$(OBJECTS) \$(LFLAGS)

neci-cpmd.a : \$(OBJECTSCPMDLIB)
	 rm -f timetag.f
	 \$(CPP) \$(CPPFLAGS) \$(SRC)/timetag.F \$(DEST)/timetag.f
	 \$(FC) \$(FFLAGS) \$(FNEWFLAGS) \$(DEST)/timetag.f
	 \$(LDLIB) \$(LDLIBFLAGS) -o \$(DEST)/neci2.a timetag.o \$(OBJECTSCPMDLIB)
	 objcopy --keep-global-symbols=\$(SRC)/\$(GLOBALS) \$(DEST)/neci2.a \$(DEST)/neci-cpmd.a

neci-vasp.a : \$(OBJECTSVASPLIB)
	 rm -f timetag.f
	 \$(CPP) \$(CPPFLAGS) \$(SRC)/timetag.F \$(DEST)/timetag.f
	 \$(FC) \$(FFLAGS) \$(FNEWFLAGS) \$(DEST)/timetag.f
	 \$(LDLIB) \$(LDLIBFLAGS) -o \$(DEST)/neci2.a timetag.o \$(OBJECTSVASPLIB)
	 objcopy --keep-global-symbols=\$(SRC)/\$(GLOBALS) \$(DEST)/neci2.a \$(DEST)/neci-vasp.a

#----------------------------------------------------------------------------
# Generate library libcpmd.a
#----------------------------------------------------------------------------
lib : \$(OBJ_LIB)
	 rm -f timetag.f
	 \$(CPP) \$(CPPFLAGS) \$(SRC)/timetag.F \$(DEST)/timetag.f
	 \$(FC) \$(FFLAGS) \$(DEST)/timetag.f
	 \$(AR) libcpmd.a timetag.o \$(OBJ_LIB)
	 \$(RANLIB) libcpmd.a

#----------------------------------------------------------------------------
# Explicit rules
#----------------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .f .F .f90 .F90 .mod .inc

\$(OBJECTSF:.o=.f) : 
	rm -f \$@
	\$(CPP) \$(CPPFLAGS) \$(SRC)/\$(@:.f=.F) \$(DEST)/\$@ 
\$(OBJECTSF90:.o=.f90) : 
	rm -f \$@
	\$(CPP) \$(CPPFLAGS) \$(SRC)/\$(@:.f90=.F90) \$(DEST)/\$@
.f.o :
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) \$(FFLAGS) \$(FNEWFLAGS) \$<" -provides "\$(*F).o " -requires "\$^"
.f90.o :
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) \$(FFLAGS) \$(F90FLAGS) \$<" -provides "\$(*F).o " -requires "\$^"

\$(OBJ_CC) :
	\$(CC) \$(CPPFLAGS) \$(CFLAGS) -c \$(SRC)/\$(@:.o=.c)

util-sp.o :
	\$(FCD) \$(FFLAGS) \$(FNEWFLAGS) \$(DEST)/util-sp.f

%.mod :
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) \$(FFLAGS) \$(F90FLAGS) \$<" -provides "\$(*F).mod" -requires "\$^"

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Dependencies
#----------------------------------------------------------------------------
END
else
cat << END >&3
#----------------------------------------------------------------------------
# Compile neci.x
#----------------------------------------------------------------------------
#MAIN_OBJS=necimain.o
DOBJECTS = \$(patsubst %,\$(DEST)/%,\$(OBJECTS))
DOBJECTSF90 = \$(patsubst %,\$(DEST)/%,\$(OBJECTSF90))
DOBJECTSF = \$(patsubst %,\$(DEST)/%,\$(OBJECTSF))
DOBJECTSNOCPMD = \$(patsubst %,\$(DEST)/%,\$(OBJECTSNOCPMD))
clean : 
	  rm -f \$(DEST)/*.f90 \$(DEST)/*.f \$(DEST)/*.o \$(DEST)/*.mod

neci.x : \$(DEST)/necimain.o \$(DOBJECTS) \$(OBJ_CC)
	 \$(FC) -o \$(DEST)/timetag.o \$(CPPFLAGS) \$(FFLAGS) timetag.F
	 rm -f \$(DEST)/neci.x
	 \$(LD) -o \$(DEST)/neci.x \$(DEST)/necimain.o \$(DEST)/timetag.o \$(DOBJECTS) \$(OBJ_CC) \$(LFLAGS)

neci-cpmd.a : \$(OBJECTSCPMDLIB)
	 \$(FC) \$(FFLAGS)\$(CPPFLAGS) -o \$(DEST)/timetag.f
	 \$(LDLIB) \$(LDLIBFLAGS) -o \$(DEST)/neci2.a \$(DEST)/timetag.o \$(OBJECTSCPMDLIB)
	 objcopy --keep-global-symbols=\$(SRC)/\$(GLOBALS) \$(DEST)/neci2.a \$(DEST)/neci-cpmd.a
     rm neci2.a

neci-vasp.a : \$(OBJECTSVASPLIB)
	 \$(FC) \$(FFLAGS)\$(CPPFLAGS) -o \$(DEST)/timetag.f
	 \$(LDLIB) \$(LDLIBFLAGS) -o \$(DEST)/neci2.a \$(DEST)/timetag.o \$(OBJECTSCPMDLIB)
	 objcopy --keep-global-symbols=\$(SRC)/\$(GLOBALS) \$(DEST)/neci2.a \$(DEST)/neci-vasp.a
     rm neci2.a

#----------------------------------------------------------------------------
# Explicit rules
#----------------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .F .F90 .mod .inc

\$(DOBJECTSF): 
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) -o \$(@) \$(FFLAGS) \$(CPPFLAGS) \$(shell basename \$@ .o).F" -provides "\$(*).o " -requires "\$^"
#	\$(FC) -o \$(@) \$(FFLAGS) \$(CPPFLAGS) \$(shell basename \$@ .o).F
\$(DOBJECTSF90):
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) -o \$(@) \$(FFLAGS) \$(F90FLAGS) \$(CPPFLAGS) \$(shell basename \$@ .o).F90" -provides "\$(*).o " -requires "\$^"
#	\$(FC) -o \$(@) \$(FFLAGS) \$(F90FLAGS) \$(CPPFLAGS) \$(shell basename \$@ .o).F90
\$(OBJ_CC) :
	\$(CC) \$(CPPFLAGS) \$(CFLAGS) -c \$(SRC)/\$(@:.o=.c)

\$(DEST)/util-sp.o :
	\$(FCD) -o \$(@) \$(FFLAGS) \$(CPPFLAGS) \$(SRC)/util-sp.F

%.mod :
#	rm -f \$(@)
	perl -w \$(SRC)/compile_mod.pl -cmp "perl -w \$(SRC)/compare_module_file.pl -compiler \$(compiler)" -fc "\$(FC) -o \$(DEST)/\$(<:.F90=.o) \$(FFLAGS) \$(CPPFLAGS) \$<" -provides "\$(*).mod" -requires "\$^"
#	\$(FC) -o \$(DEST)/\$(<:.F90=.o) \$(FFLAGS) \$(CPPFLAGS) \$<

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Dependencies
#----------------------------------------------------------------------------
END
fi
echo "done." >&2

printf "Create dependencies..." >&2
FortranFiles=`ls *.F *.c *.f90 *.F90`

for name in ${FortranFiles}
do
  if [ $verbose ]; then
    printf "[%s]" $name
  fi
  if [ $PreprocSep ] ; then
  awk 'NR==1 { Ninclude = 2;
               OldInclude["pimd.inc"] = 0;
               OldInclude["mpif.h"] = 0;
               MaxLength=60;
               linebuf="";
               ll = length(FILENAME);
               prefix = substr(FILENAME,1,ll-2);
               suffix = substr(FILENAME,ll-1,2);
               src    = "$(SRC)";
               dest = "$(DEST)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F") {
                 line = sprintf("%s.f:", prefix);
                 printf "%-16s %s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f";
               }
               if (suffix == "90") {
                 prefix = substr(FILENAME,1,ll-4);
                 suffix = substr(FILENAME,ll-3,4);
                 line = sprintf("%s.f90:", prefix);
                 printf "%-16s %s/%s%s\n", line,src,prefix,suffix;
                 suffix = ".f90";
               }
               line = sprintf("%s.o:", prefix);
               if (suffix == ".c") {
                 line = sprintf("%-16s %s/%s%s",line,src,prefix,suffix);
               } else {
                 line = sprintf("%-16s %s%s",line,prefix,suffix);
               }
             }
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               no = 0;
               if (Ninclude != 0) {
                 for ( Include in OldInclude ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (no != 1) {
                 Ninclude=Ninclude+1;
                 OldInclude[Name] = 0;
                 if (Name == "irat.inc") {
                   line=sprintf("%s %s/%s",line,dest,Name)
                 } else {
                   line=sprintf("%s %s/%s",line,src,Name)
                 }
               }
       }
       /^[ ]*USE / || /^[ ]*use / || /^[ ]*Use / {
               split($0,a,"[ ,]*");
               if(tolower(a[1])=="use") {
                  Name =  tolower(a[2]) ".mod";
               }
               else
               {
                  Name =  tolower(a[3]) ".mod";
               }
               no = 0;
               if (Ninclude != 0) {
                 for ( Include in OldInclude ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (nModules != 0) {
                 for ( Include in Modules ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (no != 1) {
                 Ninclude=Ninclude+1;
                 OldInclude[Name] = 0;
                 line=sprintf("%s %s",line,Name)
               }
       }
       (/^[ ]*module/ || /^[ ]*MODULE/){
            if(!/PROCEDURE/ && !/procedure/ && !/Procedure/) {
               split($0,a,"[ ,]*");
               if(tolower(a[1])=="module") {
                  Name =  tolower(a[2]) ".mod";
               }
               else
                {
                  Name =  tolower(a[3]) ".mod";
               }
               Modules[Name]=1;
               nModules=nModules+1
               line=sprintf("%s %s",Name,line);
            }
      }
       END   { printf "%s%s\n", linebuf,line }' ${name} >&3;
    else
  awk 'NR==1 { Ninclude = 2;
               OldInclude["pimd.inc"] = 0;
               OldInclude["mpif.h"] = 0;
               MaxLength=60;
               linebuf="";
               ll = length(FILENAME);
               prefix = substr(FILENAME,1,ll-2);
               suffix = substr(FILENAME,ll-1,2);
               src    = "$(SRC)";
               dest = "$(DEST)";
               #To have the quote.
               quote = sprintf("%c",39);
               if (suffix == ".F") {
                 line = sprintf("%s/%s.o:", dest,prefix);
                 line=sprintf("%-16s %s%s ", line,prefix,suffix);
               }
               if (suffix == "90") {
                 prefix = substr(FILENAME,1,ll-4);
                 suffix = substr(FILENAME,ll-3,4);
                 line = sprintf("%s/%s.o:",dest, prefix);
                 line=sprintf("%-16s %s%s ", line,prefix,suffix);
               }
             }
       /^[ ]*INCLUDE/ || /^[ ]*include/ {
               split($0,a,quote);
               Name = a[2];
               no = 0;
               if (Ninclude != 0) {
                 for ( Include in OldInclude ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (no != 1) {
                 Ninclude=Ninclude+1;
                 OldInclude[Name] = 0;
                 if (Name == "irat.inc") {
                   line=sprintf("%s %s/%s",line,dest,Name)
                 } else {
                   line=sprintf("%s %s",line,Name)
                 }
               }
       }
       /^[ ]*USE / || /^[ ]*use / || /^[ ]*Use / {
               split($0,a,"[ ,]*");
               if(tolower(a[1])=="use") {
                  Name =  tolower(a[2]) ".mod";
               }
               else
               {
                  Name =  tolower(a[3]) ".mod";
               }
               no = 0;
               if (Ninclude != 0) {
                 for ( Include in OldInclude ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (nModules != 0) {
                 for ( Include in Modules ) {
                   if (Include == Name ) { 
                     no = 1;
                     break;
                   }
                 }
               }
               if (no != 1) {
                 Ninclude=Ninclude+1;
                 OldInclude[Name] = 0;
                 line=sprintf("%s %s/%s",line,dest,Name)
               }
       }
       (/^[ ]*module/ || /^[ ]*MODULE/){
            if(!/PROCEDURE/ && !/procedure/ && !/Procedure/) {
               split($0,a,"[ ,]*");
               if(tolower(a[1])=="module") {
                  Name =  tolower(a[2]) ".mod";
               }
               else
                {
                  Name =  tolower(a[3]) ".mod";
               }
               Modules[Name]=1;
               nModules=nModules+1
               line=sprintf("%s/%s %s",dest,Name,line);
            }
      }
       END   { printf "%s%s\n", linebuf,line }' ${name} >&3 ; 
    fi
done

echo "done."  >&2
echo "O.K." >&2
exit 0

