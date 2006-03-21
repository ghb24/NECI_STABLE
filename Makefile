#----------------------------------------------------------------------------
# Makefile for sort.x (plane wave electronic calculation)
# Configuration: PC-PGI
# Creation of Makefile: Mar 21 2006
# on Linux slou 2.6.11.4-21.10-smp #1 SMP Tue Nov 29 14:32:49 UTC 2005 i686 i686 i386 GNU/Linux
# Author: ajwt3
#----------------------------------------------------------------------------
#
SHELL = /bin/sh
#
#--------------- Default Configuration for PC-PGI ---------------
SRC  = .
DEST = /home/ajwt3/NE-CI/SOURCE-new/dest
BIN  = .
FFLAGS = -I /usr/local/shared/fftw3/include -I /usr/local/fftw-3.0.1/include -Mr8 -Msignextend -Minform=warn
LFLAGS = -L /usr/local/shared/fftw3/lib -L /usr/local/fftw-3.0.1/lib -lfftw3 -llapack -lblas  $(QMMM_LIBS)
CFLAGS = -D__PGI
CPP = /lib/cpp -P -C -traditional
CPPFLAGS = -D__Linux -D__PGI -DLAPACK -DFFT_DEFAULT
FC = pgf90 -c -O4 -fast
LD = pgf90 -O4
AR = 
FCD = pgf90 -c -g
LDD = pgf90 -g
ARD = 
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Personal Configuration
#----------------------------------------------------------------------------
SRC = /home/ajwt3/NE-CI/SOURCE-new
FC = pgf90 -c -O4 -fast -I$(SRC)
#----------------------------------------------------------------------------
# End of Personal Configuration
#----------------------------------------------------------------------------
#
#  LIST OF FILES
#
OBJECTS = $(OBJ_AL) $(OBJ_CU) $(OBJ_MC)

OBJ_AL  = sort2.o fctrl.o sltcnd.o \
          scr.o scr1.o scr2.o matmul.o dotp.o \
          freem.o get_addr.o kb07ad.o memory.o prmem.o \
          readsr.o stopgm.o detham.o util.o frsblk.o \
          sort_1.o lineup.o gen_coul.o \
          gndts.o util-sp.o xcnrgy_dnsty.o xcener.o readinput.o

OBJ_CU  = init_coul.o fcoul.o rhoofr.o write_rho.o \
          xchole.o rlsxc_sc.o rlsxc1_sc.o rlsxc2_sc.o \
          xcholes.o fodmat.o blas_tuned_NECSX.o erf.o \
	  read_psi.o gen_coul_ueg.o init_coul2D.o \
          hfbasis.o scrtransf.o csf.o

OBJ_MC  = calcrho.o mcpaths.o montecarlo.o calcpath.o zsum.o \
          chebint.o sorti.o calcpathnci.o factorpoly.o rhodiag.o \
          hdiag.o gndts_blk.o mcpathsismc.o mcpathsmcmc.o \
          mcpathsch.o readint.o hubbard.o excit.o mcpathsold.o \
          rootfind.o

MAIN_OBJS = sort.o histogram.o

#----------------------------------------------------------------------------
# LIST OF INCLUDE FILES
#----------------------------------------------------------------------------
INCFILES = envj.inc geq0.inc irat.inc kb07ad.inc memc.inc system.h \
           time.inc quant.inc tabul.inc cons.inc fcoul.inc
INCFILES = calcp.inc cons.inc csf.inc envj.inc fcoul.inc geq0.inc \
           irat.inc memc.inc quant.inc time.inc uhfdet.inc vmc.inc 

#----------------------------------------------------------------------------
# Compile sort.x
#----------------------------------------------------------------------------
#MAIN_OBJS=sort.o histogram.o

sort.x : sort.o $(OBJECTS) $(OBJ_CC)
	 rm -f timetag.f
	 $(CPP) $(CPPFLAGS) $(SRC)/timetag.F $(DEST)/timetag.f
	 $(FC) $(FFLAGS) $(DEST)/timetag.f
	 rm -f sort.x
	 if [ $(BIN) != '.' ]; then ln -s $(BIN)/sort.x sort.x; fi
	 $(LD) -o $(BIN)/sort.x timetag.o sort.o $(OBJECTS) $(OBJ_CC) $(LFLAGS)

histogram.x : histogram.o $(OBJECTS) $(OBJ_CC)
	 rm -f timetag.f
	 $(CPP) $(CPPFLAGS) $(SRC)/timetag.F $(DEST)/timetag.f
	 $(FC) $(FFLAGS) $(DEST)/timetag.f
	 rm -f histogram.x
	 if [ $(BIN) != '.' ]; then ln -s $(BIN)/histogram.x histogram.x; fi
	 $(LD) -o $(BIN)/histogram.x histogram.o timetag.o $(OBJECTS) $(OBJ_CC) input.o $(LFLAGS)

#----------------------------------------------------------------------------
# Generate library libcpmd.a
#----------------------------------------------------------------------------
lib : $(OBJ_LIB)
	 rm -f timetag.f
	 $(CPP) $(CPPFLAGS) $(SRC)/timetag.F $(DEST)/timetag.f
	 $(FC) $(FFLAGS) $(DEST)/timetag.f
	 $(AR) libcpmd.a timetag.o $(OBJ_LIB)
	 $(RANLIB) libcpmd.a

#----------------------------------------------------------------------------
# Remove all *.o and *.f
#----------------------------------------------------------------------------
clean : 
	 rm -f $(MAIN_OBJS) $(OBJECTS) $(OBJ_CC) $(DEST)/$(OBJECTS:.o=.f)

#----------------------------------------------------------------------------
# Explicit rules
#----------------------------------------------------------------------------
ALLOBJECTS= $(MAIN_OBJS) $(OBJECTS) input.mod
.SUFFIXES:
.SUFFIXES: .o .f .F .mod .F90

$(ALLOBJECTS:.o=.f) : 
	rm -f $@
	$(CPP) $(CPPFLAGS) $(SRC)/$(@:.f=.F) $(DEST)/$@
.f.o :
	$(FC) $(FFLAGS) $<

$(OBJ_CC) :
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $(SRC)/$(@:.o=.c)

util-sp.o :
	$(FCD) $(FFLAGS) $(DEST)/util-sp.f

$(DEST)/input.mod :
	$(FCD) -c $(SRC)/input.F90

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Dependencies
#----------------------------------------------------------------------------
blas_tuned_NECSX.f:$(SRC)/blas_tuned_NECSX.F
blas_tuned_NECSX.o:blas_tuned_NECSX.f
calcpath.f:     $(SRC)/calcpath.F
calcpath.o:     calcpath.f $(SRC)/calcp.inc
calcpathnci.f:  $(SRC)/calcpathnci.F
calcpathnci.o:  calcpathnci.f
calcrho.f:      $(SRC)/calcrho.F
calcrho.o:      calcrho.f
chebint.f:      $(SRC)/chebint.F
chebint.o:      chebint.f
csf.f:          $(SRC)/csf.F
csf.o:          csf.f $(SRC)/csf.inc
detham.f:       $(SRC)/detham.F
detham.o:       detham.f
dotp.f:         $(SRC)/dotp.F
dotp.o:         dotp.f $(SRC)/geq0.inc
erf.f:          $(SRC)/erf.F
erf.o:          erf.f
excit.f:        $(SRC)/excit.F
excit.o:        excit.f
factorpoly.f:   $(SRC)/factorpoly.F
factorpoly.o:   factorpoly.f
fcoul.f:        $(SRC)/fcoul.F
fcoul.o:        fcoul.f
fctrl.f:        $(SRC)/fctrl.F
fctrl.o:        fctrl.f
fodmat.f:       $(SRC)/fodmat.F
fodmat.o:       fodmat.f
freem.f:        $(SRC)/freem.F
freem.o:        freem.f $(SRC)/system.h $(SRC)/memc.inc
frsblk.f:       $(SRC)/frsblk.F
frsblk.o:       frsblk.f $(SRC)/time.inc
gen_coul.f:     $(SRC)/gen_coul.F
gen_coul.o:     gen_coul.f
gen_coul_ueg.f: $(SRC)/gen_coul_ueg.F
gen_coul_ueg.o: gen_coul_ueg.f $(SRC)/cons.inc
get_addr.f:     $(SRC)/get_addr.F
get_addr.o:     get_addr.f
gndts_blk.f:    $(SRC)/gndts_blk.F
gndts_blk.o:    gndts_blk.f
gndts.f:        $(SRC)/gndts.F
gndts.o:        gndts.f
hdiag.f:        $(SRC)/hdiag.F
hdiag.o:        hdiag.f
hfbasis.f:      $(SRC)/hfbasis.F
hfbasis.o:      hfbasis.f
histogram.f:    $(SRC)/histogram.F
histogram.o:    histogram.f $(SRC)/vmc.inc $(SRC)/fcoul.inc \
                $(SRC)/cons.inc $(DEST)/irat.inc $(SRC)/csf.inc \
                $(SRC)/uhfdet.inc
hubbard.f:      $(SRC)/hubbard.F
hubbard.o:      hubbard.f $(SRC)/cons.inc
init_coul2D.f:  $(SRC)/init_coul2D.F
init_coul2D.o:  init_coul2D.f
init_coul.f:    $(SRC)/init_coul.F
init_coul.o:    init_coul.f
kb07ad.f:       $(SRC)/kb07ad.F
kb07ad.o:       kb07ad.f
lineup.f:       $(SRC)/lineup.F
lineup.o:       lineup.f
matmul.f:       $(SRC)/matmul.F
matmul.o:       matmul.f
mcpathsch.f:    $(SRC)/mcpathsch.F
mcpathsch.o:    mcpathsch.f
mcpaths.f:      $(SRC)/mcpaths.F
mcpaths.o:      mcpaths.f
mcpathsismc.f:  $(SRC)/mcpathsismc.F
mcpathsismc.o:  mcpathsismc.f $(SRC)/vmc.inc $(SRC)/uhfdet.inc
mcpathsmcmc.f:  $(SRC)/mcpathsmcmc.F
mcpathsmcmc.o:  mcpathsmcmc.f
mcpathsold.f:   $(SRC)/mcpathsold.F
mcpathsold.o:   mcpathsold.f
memory.f:       $(SRC)/memory.F
memory.o:       memory.f $(SRC)/system.h $(SRC)/memc.inc $(DEST)/irat.inc
montecarlo.f:   $(SRC)/montecarlo.F
montecarlo.o:   montecarlo.f $(DEST)/irat.inc
o.mcpaths.f:    $(SRC)/o.mcpaths.F
o.mcpaths.o:    o.mcpaths.f
prmem.f:        $(SRC)/prmem.F
prmem.o:        prmem.f $(SRC)/envj.inc
readinput.f:    $(SRC)/readinput.F
readinput.o:    readinput.f $(DEST)/input.mod $(SRC)/vmc.inc \
                $(SRC)/fcoul.inc
readint.f:      $(SRC)/readint.F
readint.o:      readint.f
read_psi.f:     $(SRC)/read_psi.F
read_psi.o:     read_psi.f
readsr.f:       $(SRC)/readsr.F
readsr.o:       readsr.f
rhodiag.f:      $(SRC)/rhodiag.F
rhodiag.o:      rhodiag.f $(DEST)/irat.inc
rhoofr.f:       $(SRC)/rhoofr.F
rhoofr.o:       rhoofr.f
rlsxc1_sc.f:    $(SRC)/rlsxc1_sc.F
rlsxc1_sc.o:    rlsxc1_sc.f
rlsxc2_sc.f:    $(SRC)/rlsxc2_sc.F
rlsxc2_sc.o:    rlsxc2_sc.f
rlsxc_sc.f:     $(SRC)/rlsxc_sc.F
rlsxc_sc.o:     rlsxc_sc.f
rootfind.f:     $(SRC)/rootfind.F
rootfind.o:     rootfind.f
scr1.f:         $(SRC)/scr1.F
scr1.o:         scr1.f
scr2.f:         $(SRC)/scr2.F
scr2.o:         scr2.f
scr.f:          $(SRC)/scr.F
scr.o:          scr.f
scrtransf.f:    $(SRC)/scrtransf.F
scrtransf.o:    scrtransf.f
sltcnd.f:       $(SRC)/sltcnd.F
sltcnd.o:       sltcnd.f $(SRC)/fcoul.inc
sort_1.f:       $(SRC)/sort_1.F
sort_1.o:       sort_1.f
sort2.f:        $(SRC)/sort2.F
sort2.o:        sort2.f
sort.f:         $(SRC)/sort.F
sort.o:         sort.f $(SRC)/cons.inc $(DEST)/irat.inc
sorti.f:        $(SRC)/sorti.F
sorti.o:        sorti.f
sortmods.f:     $(SRC)/sortmods.F
sortmods.o:     sortmods.f
stochdiag.f:    $(SRC)/stochdiag.F
stochdiag.o:    stochdiag.f $(SRC)/vmc.inc
stopgm.f:       $(SRC)/stopgm.F
stopgm.o:       stopgm.f $(SRC)/system.h
timec.f:        $(SRC)/timec.F
timec.o:        timec.f
timer.f:        $(SRC)/timer.F
timer.o:        timer.f $(SRC)/time.inc $(SRC)/envj.inc
timetag.f:      $(SRC)/timetag.F
timetag.o:      timetag.f
util.f:         $(SRC)/util.F
util.o:         util.f $(SRC)/geq0.inc $(SRC)/system.h
util-sp.f:      $(SRC)/util-sp.F
util-sp.o:      util-sp.f
write_rho.f:    $(SRC)/write_rho.F
write_rho.o:    write_rho.f
xcener.f:       $(SRC)/xcener.F
xcener.o:       xcener.f
xchole.f:       $(SRC)/xchole.F
xchole.o:       xchole.f
xcholes.f:      $(SRC)/xcholes.F
xcholes.o:      xcholes.f
xcnrgy_dnsty.f: $(SRC)/xcnrgy_dnsty.F
xcnrgy_dnsty.o: xcnrgy_dnsty.f
zsum.f:         $(SRC)/zsum.F
zsum.o:         zsum.f
input.mod:     $(SRC)/input.F90
