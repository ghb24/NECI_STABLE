[main]
fc = pgf90 
compiler = PGI-pgf95-on-LINUX
cc = g++
ccd = g++
cpp = cpp -C -traditional
cppflags =  -D__Linux -DPOINTER8 -DPARALLEL -DHAVE_SSE2 -D__INT64 -D__PGI -D__SHARED_MEM -DCBINDMPI
ld = mpic++
libs = -L ~/src/lib/fftw-3.1.2/lib -lfftw3 -lacml -pgf90libs -lrt -lm
module_flag = -module

[dbg]
fflags = -g -r8 pc=64 -Msignextend -Minform=warn
cflags = -g -C
f90flags = -Mfree -Mbounds

[opt]
fflags = -fastsse -tp k8-64 -r8 pc=64 -Msignextend -Minform=warn
cflags = -O
f90flags = -Mfree 
