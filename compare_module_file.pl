# Compare two .mod files which are given on the command line.
#
# Note that if a compiler (i.e. a compiler and operating system) is specified on
# the command line that the script does not know about, then it is not an error.
# In such a case the script is like "cmp -s".  This is to allow portability to
# untested systems without freeking out the end-user.
#
#*******************************************************************************
# (c) Daniel Grimwood, University of Western Australia, 2000-2003
#
# $Id: compare_module_file.pl,v 2.7 2003/12/18 07:55:52 reaper Exp $ *patched*
#
#*******************************************************************************

#*******************************************************************************
# Set up the hash array of compilers known, and the subroutine to call for each.
#
# In most cases the module file format will probably only depend on the
# compiler, but we allow it to depend on the OS just in case.
#
# When writing the routines, remember that "-" is reserved, so change it to "_".
# Don't forget to write the routine at the bottom of the script.
my (%compiler_array);
%compiler_array = (
  "ABSOFT-f95-on-DARWIN"   =>  ABSOFT_f95_on_DARWIN,
  "ABSOFT-f95-on-LINUX"    =>  ABSOFT_f95_on_LINUX,
  "COMPAQ-f90-on-OSF1"     =>  COMPAQ_f90_on_OSF1,
  "COMPAQ-f95-on-OSF1"     =>  COMPAQ_f95_on_OSF1,
  "COMPAQ-f95-on-LINUX"    =>  COMPAQ_f95_on_LINUX,
  "DEC-f90-on-OSF1"        =>  DEC_f90_on_OSF1,
  "DEC-f95-on-OSF1"        =>  DEC_f95_on_OSF1,
  "FUJITSU-f90-on-LINUX"   =>  FUJITSU_f90_on_LINUX,
  "FUJITSU-f95-on-LINUX"   =>  FUJITSU_f95_on_LINUX,
  "GCC-f95-on-LINUX"       =>  GCC_f95_on_LINUX,
  "IBM-xlf-on-AIX"         =>  IBM_xlf_on_AIX,
  "IBM-xlf90-on-AIX"       =>  IBM_xlf90_on_AIX,
  "INTEL-ifc-on-LINUX"     =>  INTEL_ifc_on_LINUX,
  "INTEL-ifort-on-LINUX"   =>  INTEL_ifort_on_LINUX,
  "INTEL-ifort9-on-LINUX"  =>  INTEL_ifort9_on_LINUX,
  "LAHEY-lf95-on-LINUX"    =>  LAHEY_lf95_on_LINUX,
  "MIPSPRO-f90-on-IRIX64"  =>  MIPSPRO_f90_on_IRIX64,
  "NEC-f90-on-SX8"         =>  NEC_f90_on_SX8,
  "PGI-pgf95-on-LINUX"     =>  PGI_pgf95_on_LINUX,
  "SUNSTUDIO-f95-on-LINUX" =>  SUNSTUDIO_f95_on_LINUX,
  "WORKSHOP-f90-on-SUNOS"  =>  WORKSHOP_f90_on_SUNOS,
  "WORKSHOP-f95-on-SUNOS"  =>  WORKSHOP_f95_on_SUNOS);

#*******************************************************************************
# Argument checking.
#
$argerr=0;
$n_arg=$#ARGV+1;
$n_arg >= 2 || do {print STDERR "Error : need at least two arguments\n"; $argerr=1; };
@ARGS = @ARGV[0 .. $n_arg-3];
$file1 = $ARGV[$n_arg-2];
$file2 = $ARGV[$n_arg-1];

while (@ARGS) {
  $arg=shift @ARGS;
  for ($arg) {
    /^-compiler/ && do {
      $fc=shift @ARGS;
      defined $fc || do {print STDERR "Error : no compiler specified\n"; $argerr=1};
      last;
    };
    print STDERR ("Error : unexpected argument $arg\n");
    $argerr=1;
  }
}

if ($argerr==1) {
  warn(
    "\nUsage :\n",
    "\t perl -w compare_module_file.pl [-compiler compiler_id] \\\n",
    "\t\t filename1 filename2 \n",
    "\n",
    "Where :\n",
    "\tfilename1 and filename2 are the two files to compare.\n",
    "\t\"-compiler compiler_id\" specifies the name of the compiler.\n");
   print STDERR "\nThe following list of compilers are recognised by -compiler :\n";
   foreach $word (keys %compiler_array) { print STDERR "\t$word\n"; }
   print STDERR "\n";
   exit 1;
}

#*******************************************************************************
# check whether files exist.
#
(-f $file1) or die "File $file1 does not exist\n";
(-f $file2) or die "File $file2 does not exist\n";

#*******************************************************************************
# quick filesize comparison.
#
$size1 = (stat($file1))[7];
$size2 = (stat($file2))[7];
($size1 == $size2) or exit 1;

#*******************************************************************************
# In this bit, we farm out to various routines to set an array of file offsets
# where differences between the files are to be ignored.  Then do the
# comparison.
@skip_array = ();
if (defined $fc) {
  $routine = $compiler_array{$fc};
  defined ($routine) && &$routine; # call the corresponding subroutine to $fc
}

open(FILE1,$file1) or die "Cannot open $file1\n";
open(FILE2,$file2) or die "Cannot open $file2\n";
binmode(FILE1);
binmode(FILE2);

$result = &do_compare();  # do the actual comparison
close(FILE2);
close(FILE1);
exit $result;

#*******************************************************************************
# Now for the subroutines.
#*******************************************************************************
sub do_compare {
  $i = 0;
  while ((! eof FILE1) && (! eof FILE2)) {
    $i++;
    if (getc(FILE1) ne getc(FILE2)) {
      # should we ignore this file offset?
      if (scalar(grep (/^$i$/,@skip_array)) eq 0) {
        return 1;
      }
    }
  }
  return 0;
}

#*******************************************************************************
# User defined subroutines for each compiler...

sub COMPAQ_f90_on_OSF1 {
  push @skip_array,(25,26,27,45,46,47);
# 45,46,47 in Compaq Fortran X5.4A
}

sub COMPAQ_f95_on_OSF1 {
  push @skip_array,(25,26,27,45,46,47);
# 45,46,47 in Compaq Fortran X5.4A
}

sub COMPAQ_f95_on_LINUX {
  push @skip_array,(48,49,50,51);
}

sub DEC_f90_on_OSF1 {
# don't skip anything.
}

sub DEC_f95_on_OSF1 {
# don't skip anything.
}

sub FUJITSU_f90_on_LINUX {
# don't skip anything.
}

sub FUJITSU_f95_on_LINUX {
# don't skip anything.
}

sub IBM_xlf_on_AIX {
# version 5.1, untested, copied from xlf90 below.
  @reverse = (29,30,31,32,33,34,35,36,37,38,39,40,41,42);
  foreach $i (@reverse) {
    push @skip_array,$size1-$i+1;
  }
}

sub IBM_xlf90_on_AIX {
# version 5.1.  Offsets start from end of file.
  @reverse = (29,30,31,32,33,34,35,36,37,38,39,40,41,42);
  foreach $i (@reverse) {
    push @skip_array,$size1-$i+1;
  }
}

sub GCC_f95_on_LINUX {
  # Gfortran and G95 have plain text module
  # First line has to be skipped. Put all characters
  # in this line into the skip_array
  my $i = 0;
  open(FILE1,$file1) or die "Cannot open $file1\n";
  binmode(FILE1);
  while (getc(FILE1) ne "\n") {
    $i++;
    push @skip_array,$i;
  }
  close(FILE1);
}

sub INTEL_ifc_on_LINUX {
  # version 7.1.  For version 08.00.00 if their module file format, ignore the
  # last record of the file, where the record separator is chr(0).
  my $lastreclength = 0;
  open(FILE1,$file1) or die "Cannot open $file1\n";
  binmode(FILE1);
  local $/ = chr(0);
  while (<FILE1>) {
    $lastreclength = length($_);
  }
  close(FILE1);
  for ($j=0; $j<$lastreclength; $j++) {
    push @skip_array,$size1-$j;
  }
  return 0;
}

sub INTEL_ifort_on_LINUX {
  # version 8.0.
  push @skip_array,(45,46,47,48);
}

sub INTEL_ifort9_on_LINUX {
  # version 9.x and 10.0
  push @skip_array,(49,50,51,52);
}

sub LAHEY_lf95_on_LINUX {
# don't skip anything.
}

sub MIPSPRO_f90_on_IRIX64 {
# don't skip anything.
}

sub NEC_f90_on_SX8 {
# don't skip anything.
}

sub PGI_pgf95_on_LINUX {
  # PGI has plain text module
  # Date and Time is found in third line
  my $i = 0;
  my $nl = 0;
  open(FILE1,$file1) or die "Cannot open $file1\n";
  binmode(FILE1);

  for ($nl=0; $nl<2; $nl++) {
    while (getc(FILE1) ne "\n") {
      $i++;
    }
    $i++;
  }
  while (getc(FILE1) ne "\n") {
    $i++;
    push @skip_array,$i;
  }

  close(FILE1);
}

sub SUNSTUDIO_f95_on_LINUX {
# don't skip anything.
}

sub WORKSHOP_f90_on_SUNOS {
# don't skip anything.
}

sub WORKSHOP_f95_on_SUNOS {
# don't skip anything.
}


sub ABSOFT_f95_on_LINUX {
# don't skip anything.
}

sub ABSOFT_f95_on_DARWIN {
# don't skip anything.
}

