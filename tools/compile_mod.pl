# compile_mod.pl
#
# This is a Perl script that will check to see if a module really needs
# recompiling, or if we can get away with pretending to compile it.  If the
# module needs recompiling, then the "-fc" argument is invoked.  If not, then
# "touch" is used on the files that would be outputted by the compiler, as
# specified by the "-provides" argument.
#
# The script records the last compilation date of a module, and stores this
# information as the timestamp of a file with extension .time .  If this file is
# not present, then it is assumed the module has not been recompiled before,
# which is the safest option in such a case.
#
#
# Usage :
# perl -w compile_mod.pl -fc command -provides filenames
#       -requires filenames -cmp command [-mod_ext ext]
#
# Where :
#       -fc command             : "command" is the complete compiler command.
#       -provides filenames     : "filenames" are the files produced by
#                                 running the compiler command.
#       -requires filenames     : "filenames" are all files required before
#                                 running the compiler command, including
#                                 module information files from "USE"
#                                 statements.
#       -cmp command            : "command" is the program used to compare two
#                                 binary files and return error code <> 0 if
#                                 they differ.  E.g. "cmp -s".
#       -mod_ext ext            : "ext is the module information file
#                                 extension.  Default is "mod".
#
# Notes :
#       Use quotes to enclose multiple words as a single argument.
#
# Example :
#       perl -w ./compile_mod.pl -fc "f90 -c a.f90" -provides "a.mod b.mod"
#               -requires "z.mod y.mod" -cmp "cmp -s"
#
#*******************************************************************************
# $Id: compile_mod.pl,v 2.3 2004/03/05 23:28:44 reaper Exp $
#
# Copyright 2002 Daniel Grimwood, <reaper@theochem.uwa.edu.au>.
#
# Freely use or hack this script, provided you don't sell it.
#*******************************************************************************

use File::Copy;  # use the native file copy/move mechanisms.
use File::Basename; # to extract the base part of a filename.

#*******************************************************************************
# Default value
$mod_ext = "mod";

#*******************************************************************************
# command line parsing
#
$argerr=0;
$argerr=1 if ($#ARGV+1==0);

while (@ARGV) {
  $arg=shift;
  for ($arg) {
    /^-cmp/ && do {$cmp=shift; last; };
    /^-fc/ && do {$fc=shift; last; };
    /^-provides/ && do {$provides=shift; last; };
    /^-requires/ && do {$requires=shift; last; };
    /^-mod_ext/ && do {$mod_ext=shift; last; };
    warn "Error : unexpected argument $arg\n";
    $argerr=1;
  }
}
defined $cmp || do {$argerr=1; warn "Error : -cmp flag not supplied.\n"};
defined $fc || do {$argerr=1; warn "Error : -fc flag not supplied.\n"};
defined $provides || do {$argerr=1; warn "Error : -provides flag not supplied.\n"};
defined $requires || do {$argerr=1; warn "Error : -requires flag not supplied.\n"};

if ($argerr==1) {
  die(
    "\nUsage :\n",
    "\t perl -w ./compile_mod.pl -fc command -provides filenames \\\n",
    "\t\t-requires filenames -cmp command [-mod_ext ext] \n",
    "\n",
    "Where :\n",
    "\t-fc command\t\t: \"command\" is the complete compiler command.\n",
    "\t-provides filenames\t: \"filenames\" are the files produced by\n",
    "\t\t\t\t  running the compiler command.\n",
    "\t-requires filenames\t: \"filenames\" are all files required before\n",
    "\t\t\t\t  running the compiler command, including\n",
    "\t\t\t\t  module information files from \"USE\"\n",
    "\t\t\t\t  statements.\n",
    "\t-cmp command\t\t: \"command\" is the program used to compare two\n",
    "\t\t\t\t  binary files and return error code <> 0 if\n",
    "\t\t\t\t  they differ.  E.g. \"cmp -s\".\n",
    "\t-mod_ext name\t\t: \"ext\" is the module information file\n",
    "\t\t\t\t  extension.  Default is \"mod\"\n",
    "\n",
    "Notes :\n",
    "\tUse quotes to enclose multiple words as a single argument.\n",
    "\n",
    "Example :\n",
    "\t perl -w ./compile_mod.pl -fc \"f90 -c a.f90\" -provides \"a.mod b.mod\"\n",
    "\t\t-requires \"z.mod y.mod\" -cmp \"cmp -s\"\n",
    "\n");
}

@provides = split / /,$provides;
@requires = split / /,$requires;

#*******************************************************************************
# The main part of the program.
#
$do_recompile = &TEST_RECOMPILE; # test whether to recompile or just touch.

if ($do_recompile==0) {
 &TOUCH_COMPILE;
} else {
 $exit_val = &REAL_COMPILE;
 exit $exit_val;
}


#*******************************************************************************
# In TOUCH_COMPILE, touch the resulting object and module file(s) to be up to
# date without recompiling them.
sub TOUCH_COMPILE {
  my $now;
  my $tmp;
  $now = time;
  utime $now, $now, @provides;
  $tmp = join(", ",@provides);
  print "No need to recompile $tmp\n";
}

#*******************************************************************************
# In REAL_COMPILE, compile the resulting object and module file(s).
# 
# If the .mod files already exist, then they are compared before and after the
# compile.  If they change upon recompilation, the corresponding .time files get
# touched.
#
sub REAL_COMPILE {
  my (@touch_array, @provides_mods, @touch_this);
  my %touch_mod;
  my ($now, $item, $exit_value, $timefile);

  # We're only interested in the .mod files, so make them a convenient array.
  # Set provides_mods=1 for the .mod if no need to do cmp test, otherwise
  # default to do it.
  foreach $item (@provides) {
    $_ = $item;
    if (m/\.$mod_ext/) {
      push @provides_mods, $item;
      $touch_mod{"$item"}=0;
    }
  }

  #check which .mod files don't exist.  No need to do cmp test on them.
  foreach $item (@provides_mods) {
    if ($touch_mod{$item}==0) {
      $touch_mod{$item}=1 if (! -e $item);
    };
  }

  # Any missing .time files?  Safest to assume their .mod files are new.
  foreach $item (@provides_mods) {
    if ($touch_mod{$item}==0) {
      ($timefile = $item) =~ s/\.$mod_ext/\.time/i;
      -e $timefile || do {$touch_mod{$item}=1};
    }
  }

  # backup any .mod files that already exist.
  foreach $item (@provides_mods) {
    $touch_mod{$item}==0 && do { copy($item,"$item~"); }
  }

  # now do the real compilation.
  $exit_value = &ECHO_AND_RUN($fc);
  $now = time;

  # perhaps the compiler put it into the current directory instead of the module
  # directory.
  foreach $item (@provides_mods) {
    $basemod = basename($item);
    -e $basemod && move($basemod,$item);
  }

  # if compilation not successful, then cleanup and exit.
  if ($exit_value != 0) {
    foreach $item (@provides_mods) {
      if ($touch_mod{$item} == 0) {
        unlink "$item\~";
        ($timefile = $item) =~ s/\.$mod_ext/\.time/i;
        unlink "$timefile";
      }
    }
    return $exit_value;
  }

  # For each .mod file that was backed up, do the comparison.
  # After each comparison, delete the corresponding backup file.
  foreach $item (@provides_mods) {
    if ($touch_mod{$item} == 0) {
      system("$cmp $item $item~");
      $exit_value = $?;
      $exit_value == 0 || do {$touch_mod{$item}=1};
      unlink "$item\~";
    }
  }

  # touch the .time files that need to be touched.
  foreach $item (@provides_mods) {
    if ($touch_mod{$item} == 1) {
      ($timefile = $item) =~ s/\.$mod_ext/\.time/i;
      if (! -e $timefile) {      # missing files get created.
        open(TOUCHFILE, "> $timefile");
        close(TOUCHFILE);
      }
      push @touch_this, $timefile;
    }
  }
  utime $now, $now, @touch_this;
  return $exit_value;
}

#*******************************************************************************
# In TEST_RECOMPILE, test whether resulting object and module file(s) are out of
# date.
#
# * If any of the resulting files or prerequisite files don't exist, then
# recompile, otherwise,
# * If any prerequisite modules were compiled more recently than the resulting
# object or module files then recompile, otherwise,
# * If no need to recompile, touch the resulting object and module files.
#
# Currently, if any file is missing we take the safe option and recompile.
#
sub TEST_RECOMPILE {
  my @provide_date;
  my ($require_date, $recompile, $n, $filename, $oldest_provided);

  # default to not be required to recompile.  Then test to see whether we do.
  $recompile = 0;

  # test - are there any resulting files?
  if ($recompile==0) {
    $recompile=1 if scalar(@provides)==0;
  }
  # test - do all resulting files exist?
  if ($recompile==0) {
    foreach $item (@provides) {
      $_ = $item;
      -e $item || do { $recompile = 1;
                       print STDERR "Target $item does not exist, compilation forced.\n";
                       last; };
      if (/\.$mod_ext/i) { # if .mod file, check whether .time file exists too.
        ($filename = $item) =~ s/\.$mod_ext/.time/is ;
        -e $filename || do { $recompile = 1;
                       print STDERR "Target $item ($filename) does not exist, compilation forced.\n";
                       last; };
      }
    }
  }

  # test - are there any prerequisites?
  if ($recompile==0) {
    $recompile=1 if scalar(@requires)==0;
  }

  # test - do all prerequisite files exist?
  if ($recompile==0) {
    foreach $item (@requires) {
      $_ = $item;
      -e $item || do { $recompile = 1;
                       print STDERR "Prerequisite $item does not exist, compilation forced.\n";
                       last; };
      if (/\.$mod_ext/i) { # if .mod file, check whether .time file exists too.
        ($filename = $item) =~ s/\.$mod_ext/.time/is ;
        -e $filename || do { $recompile = 1;
                       print STDERR "Prerequisite $item ($filename) does not exist, compilation forced.\n";
                       last; };
      }
    }
  }

  # test - compare file modification dates.  We can now safely assume all
  # files exist.
  if ($recompile==0) {
    # get modification dates of resulting files.
    foreach $item (@provides) {
      push @provide_date, (stat($item))[9];
    }
    $oldest_provided = &min(@provide_date);

    # get modification dates of prerequisites.
    foreach $item (@requires) {
      # for prerequisites that are modules, want .time not .mod.
      ($filename = $item) =~ s/\.$mod_ext/.time/is ;
      (stat($filename))[9] > $oldest_provided && do { $recompile = 1; last; };
    }
  }
  return $recompile;

}

#*******************************************************************************
# ECHO_AND_RUN prints out the command before executing it.  The return value of
# the routine is the exit value of the command.
#
sub ECHO_AND_RUN {
  my ($mycommand, $exit_value);
  $mycommand = join(' ',@_);
  $mycommand =~ s/{/\\{/g; # replace curlies to stop globbing
  print $mycommand,"\n";
  system($mycommand);
  if ($?==-1) {
    $exit_value = -1;
  } else {
    $exit_value = $?/256;
  }
  return $exit_value;
}

#*******************************************************************************
# routine to return maximum of an array.
sub max {
  my $max = shift(@_);
  foreach $item (@_) {
    $max = $item if $max < $item;
  }
  return $max;
}

#*******************************************************************************
# routine to return minimum of an array.
sub min {
  my $min = shift(@_);
  foreach $item (@_) {
    $min = $item if $min > $item;
  }
  return $min;
}

