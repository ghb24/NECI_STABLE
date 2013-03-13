[main]
fc = ftn
cc = CC
ld = CC
ldflags = -m64 -rdynamic
libs = -lfftw3 -lacml -lrt -lm
cpp = cpp -C -traditional
cppflags =  -D__Linux -DPOINTER8 -DPARALLEL -D__PATHSCALE__ -D__INT64 -DCBINDMPI
compiler = pathf95-on-LINUX
module_flag = -module

[dbg]
fflags = -g -fno-second-underscore -m64
cflags = -g -m64
ldflags = -g -Bstatic
f90flags = -ffortran-bounds-check

[opt]
fflags = -O3 -OPT:Ofast -LNO:simd_verbose=ON -fno-second-underscore -m64
cflags = -O -m64
ldflags = -O3
