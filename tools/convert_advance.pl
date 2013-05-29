#!/usr/bin/perl -pi 
# Usage:
# ./convert_advance.pl filenames
#
# Convert non-standard fortran write statements from (e.g.):
# write (unit,'(a3,f7.2$)')
# to:
# write (unit,'(a3,f7.2)',advance='no')
# If the file ends in .F, then it is assumed to be fixed format, and a line 
# which consequently exceeds 72 characters is split into two lines, with the
# appropriate continuation marks.
# Require "valid" format specifiers (format string enclosed in brackets).
# Does not work with something like write (6,'A,$').
# Works with/without a comma before the $.
# Format string can be quoted using single or double quotes.
# WARNING: fails if a format string is produced on the fly containing a $,
# where the $ really ought to be in the write statement using that format
# string.
if (/write/i && /\$/) {
    $_=~s/,?\$(\)("|'))/\1,advance='no'/;
    if (length($_)>72 && !/^ *!/ && !/^C/i && $ARGV=~/.*\.F$/) {
        # Assume fixed format.
        @w=split('(?<=advance=\'no\'\))');
        $_=sprintf "%s\n     &%66s",$w[0],$w[1];
    }
}
