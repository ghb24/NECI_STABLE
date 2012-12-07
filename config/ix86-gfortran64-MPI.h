[main]
fc = gfortran
cc = mpicxx
ld = mpicxx
ldflags = -m64
compiler = GCC-f95-on-LINUX
cpp = cpp -C -traditional
cppflags = -D__Linux -DPOINTER8 -DPARALLEL -D__INT64 -D__SHARED_MEM -DCBINDMPI -D__GFORTRAN__ 
libs = -lacml -lrt -lgfortran -lfftw3 -lm
module_flag = -J

[dbg]
fflags = -ffree-line-length-none -fcray-pointer -g -m64 -fbacktrace -fdefault-real-8 -Waggregate-return -Waliasing -Wampersand -Wcharacter-truncation -Wintrinsics-std -Wno-tabs -Wsurprising -Wunderflow -fdefault-integer-8
cflags = -g -m64 -Waddress -Wcast-align -Wchar-subscripts -Wcomment -Wformat -Wimplicit -Wimplicit-int -Wimplicit-function-declaration -Wmain -Wmissing-braces -Wmultichar -Wnested-externs -Wparentheses -Wpointer-arith -Wpointer-sign -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -DZLIB
f90flags = -fbounds-check

[opt]
fflags = -fcray-pointer -O4 -m64 -fdefault-real-8 -Waggregate-return -Waliasing -Wampersand -Wcharacter-truncation -Wintrinsics-std -Wno-tabs -Wsurprising -Wunderflow -fdefault-integer-8
cflags = -O -m64 -Waddress -Wcast-align -Wchar-subscripts -Wcomment -Wformat -Wimplicit -Wimplicit-int -Wimplicit-function-declaration -Wmain -Wmissing-braces -Wmultichar -Wnested-externs -Wparentheses -Wpointer-arith -Wpointer-sign -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -DZLIB
