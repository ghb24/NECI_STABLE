! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#ifdef MOLPRO
      subroutine GetDipMomInts(pints,isize,znuc,zcor)
      implicit double precision(a-h,o-z)
      integer isize
      double precision pints(isize,3)
      double precision znuc(3),zcor(3)
      include "common/maxatm"
      include "common/tapes"
      include "common/dumpinfow"
      include "common/cstate"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      include "common/cbas"
      include "common/clseg"
      include "common/cref"
      include "common/ctran2"
      include "common/code"
      include "common/cmpp"
      include "common/d2gen_cvb"
      npr2 = 0
      do i = 1,3
          do j = 1,isize
              pints(j,i) = 0.0d0
          enddo
          znuc(i) = 0.0d0
          zcor(i) = 0.0d0
      enddo
      do i=1,nirr
          npr2 = npr2 + (nactt(i)*(nactt(i)+1))/2
      enddo
      if(npr2.ne.isize) then
          write(6,*) "Molpro active orbital space: ",nactt(1:nirr)
          write(6,*) "Size of molpro symmetry packed arrays: ",npr2
          write(6,*) "Size of neci symmetry packed arrays: ",isize
          stop 'Size of symmetry packed arrays wrong in GetDipMomInts'
      endif
      iprop = icorr(npr2)
      do ipr = 4,6
          call pget(q(iprop),ipr,zznuc,zzcor)
          do i = 1,isize
              pints(i,ipr-3)=q(iprop+i-1)
          enddo
          znuc(ipr-3) = zznuc
          zcor(ipr-3) = zzcor
      enddo

      return
      end subroutine GetDipMomInts

#else
      !In case it doesn't like empty files
      subroutine null_wrapper
          implicit none
          stop 'Should not be here'
      end subroutine null_wrapper

#endif


