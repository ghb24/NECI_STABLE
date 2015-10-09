[main]
fc = ifort
cc = mpiicc
ld = mpiicc
ldflags = -L $(LD_LIBRARY_PATH) -rdynamic -xHost -i-dynamic
compiler = INTEL-ifort9-on-LINUX
cpp = cpp -traditional
cppflags = -D__Linux -DPOINTER8 -DPARALLEL -DHAVE_SSE2 -D__INT64 -D__SHARED_MEM -DDISABLE_FFTW -DCBINDMPI -D__IFORT
libs = -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lrt -lifcore -lifport
module_flag = -module

[dbg]
fflags = -r8 -g -traceback -i8 -pc64 -auto -warn nousage
cflags = -g -C -traceback -DZLIB
f90flags = -check bounds -stand f03

[opt]
# n.b. -xHost can be replaced with -xavx (or nothing) if processors are inhomogeneous in a cluster.
fflags = -r8 -O3 -i8 -pc64 -auto -warn nousage -xHost
cflags = -O -DZLIB -xHost
f90flags = -stand f03
