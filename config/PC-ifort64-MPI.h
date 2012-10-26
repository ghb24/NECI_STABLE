[main]
fc = ifort
cc = mpicxx
ld = mpicxx
ldflags = -i-dynamic -L $(LD_LIBRARY_PATH)
compiler = INTEL-ifort9-on-LINUX
cpp = cpp -C -traditional
cppflags = -D__Linux -DPOINTER8 -DPARALLEL -DHAVE_SSE2 -D__INT64 -D__SHARED_MEM -DCBINDMPI
libs = -lfftw3 -lacml -lrt -lifcore -lifport -lm
module_flag = -module

[dbg]
fflags = -r8 -g -traceback -i8 -pc64 -auto -vec-report0 -warn nousage
cflags = -g -C -traceback -vec-report0 -DZLIB
f90flags = -check bounds -stand f03

[opt]
fflags = -r8 -O3 -i8 -pc64 -auto -vec-report0 -warn nousage
cflags = -O -vec-report0 -DZLIB
f90flags = -stand f03
