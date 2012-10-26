[main]
fc = gfortran
cc = g++
ld = gfortran
ldflags =
compiler = GCC-f95-on-LINUX
cpp = cpp -C -traditional
cppflags = -D__Linux -DPOINTER8 -DHAVE_SSE2 -D__INT64
libs = -lfftw3 -lacml -lstdc++ -lrt -lm
module_flag = -J

[dbg]
fflags = -fcray-pointer -ffree-line-length-none -g
cflags = -g -traceback

[opt]
fflags = -fcray-pointer -ffree-line-length-none -O3 
cflags = -O
