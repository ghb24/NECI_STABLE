#ifdef MOLPRO
      subroutine GetDipMomInts(pints,isize,znuc,zcor,nt_frz,ntd_frz)
      implicit double precision(a-h,o-z)
      integer isize,nt_frz(8),ntd_frz(8)
      double precision pints(isize,3)
      double precision znuc(3),zcor(3)
      character(*), parameter :: this_routine = 'GetDipMomInts'
      include "common/maxatm"
      include "common/tapes"
!      include "common/dumpinfow"
#include "common/dumpinfow.fh"
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

      write(iout,*) "Extracting dipole moment integrals from molpro..."
          
      do i = 1,3
          do j = 1,isize
              pints(j,i) = 0.0d0
          enddo
          znuc(i) = 0.0d0
          zcor(i) = 0.0d0
      enddo

      npr2 = 0
      do i = 1,nirr
          npr2 = npr2 + nactt(i)
      enddo
      if(npr2.gt.0) then
          !Assume that we have done a prior mcscf/multi calculation for the dipole moment integrals.
          !See: multi/muprop.F 
          write(iout,*) "Active space found. "
          write(iout,*) "Assuming prior mcscf/multi calc for integrals"

          npr2 = 0
          do i=1,nirr
              npr2 = npr2 + (nactt(i)*(nactt(i)+1))/2
          enddo
          if(npr2.ne.isize) then
            write(iout,*) "Molpro active orbital space: ",nactt(1:nirr)
            write(iout,*) "Size of molpro symmetry packed arrays: ",npr2
            write(iout,*) "Size of neci symmetry packed arrays: ",isize
            call stop_all(this_routine,
     &            'Size of symmetry packed arrays wrong: GetDipMomInts')
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
          write(6,*) "z-component dipole moment integrals: "
          do i = 1,isize
              write(6,*) i,pints(i,3)
          enddo
          write(6,*) znuc(3)
          write(6,*) zcor(3)
      else
          !Assume that we have done a prior CISD calculation for the dipole moment integrals.
          !See: ccsd/ccmp2.F
          call GetDipMomInts_cisd(pints,isize,znuc,zcor,nt_frz,ntd_frz)

      endif

      return
      end subroutine GetDipMomInts
  
      !Assume that we have done a prior CISD calculation for the dipole moment integrals.
      !See: ccsd/ccmp2.F
      subroutine GetDipMomInts_cisd(pints,isize,znuc,zcor,nt_frz,       &
     &         ntd_frz)
      implicit double precision(a-h,o-z)
      integer isize,nt_frz(8),ntd_frz(8),nfrz_neci(8),tri_neci
      double precision pints(isize,3)
      double precision znuc(3),zcor(3)
      include "common/big"
      include "common/cbas"
      include "common/corb"
#include "common/cpair.fh"
      include "common/cvec"
      include "common/ceig"
      include "common/cpar"
#include "common/cpfil.fh"
      include "common/maxatm"
!      include "common/dumpinfow"
#include "common/dumpinfow.fh"
      include "common/tapes"
      include "common/clseg"
      logical exists
      double precision origin(3)
      character(*), parameter :: this_routine = 'GetDipMomInts_cisd'

      !First, get the MO coefficients, since we will need to transform the integrals ourselves
      ! get orbital coefficients
      ! see cidft/rpa_polar.F90? This is not great, because we want symmetry
!      allocate(coef_mo(nta,nta))
!      call get_matrix('ORB',coef_mo,nta,carec,cafil)
      origin(1) = 0.0d0
      origin(2) = 0.0d0
      origin(3) = 0.0d0

      !For some reason, the code below doesn't seem to be able to get the sym blocked MO coefficients correctly.
      !Why not?! For a start, ntqgcx is filled with rubbosh. I am just trying to reproduce the same 
      !call in multiple places in the code, ie. ccsd/ccmp2.F, line 1842
      inquire(file='MOInts',exist=exists)
      if(.not.exists) then
           call stop_all(this_routine,
     &             'Cannot find MO integrals-Run hacked CISD')
      endif
      iorbmo=icorr(ntqgcx(1))
      open(69,file='MOInts',status='old')
      read(69,*) ntgx
      do i = 1,8
          read(69,*) ntqgcx(i),ntqgc(i),nt(i),ntd(i)
      enddo
      do i = 1,ntqgcx(1)
          read(69,*) q(iorbmo+i-1)
      enddo
      close(69)

      !compress
      do i=1,ntqgcx(1)
          if(q(iorbmo+i-1).eq.0.0d0) then
              !Move all the next values down
              do j = i,ntqgcx(1)-1
                  q(iorbmo+j-1) = q(iorbmo+j)
              enddo
          endif
      enddo

!      write(6,*) "After compression: "
!      do i = 1,ntqgcx(1)
!          write(6,*) i,q(iorbmo+i-1)
!      enddo

      !Work out the number of frozen orbitals in each symmetry
      nfrz_neci(:) = 0
      do i = 1,8
          nfrz_neci(i) = nt(i)-nt_frz(i)
      enddo

!      call excom(2)
!      iorbmo=icorr(ntqgcx(1))
!      call cread(q(iorbmo),ntqgcx(1),kopfil,nwf,nwf)

      call druckm(q(iorbmo),1,0,'orbitals hack',0)   !This would print out the orbitals

      do i = 1,8
          if(nt(i).ne.nt_frz(i)) then
              write(6,'(a,i6,a,i6)') "number of orbs frozen in sym ",   &
     &          i," is: ",nt(i)-nt_frz(i)
          endif
      enddo

      iscr=icorr(2*ntqgcx(1)+2*ntgx)  !Memory for integral transforms

      iints=icorr(2*ntqgcx(1)+2*ntgx)   !Get memory pointer for integrals
      
      do i = 1,3
          !read the operators (in the AO basis)
          !1) Operator in AO basis
          !2) Nuclear value of operator
          !3) name of operator
          !4) Component of operator
          !5) Symmetry of operator
          !6) np = 0: return as square array
          !7) optional
          !8) Return sym of operator
          call read_op(q(iints),dipnu,'DM',i,1,1,origin,isyop)

!          write(6,*) "AO basis dip mom ints: ",i
!          do j = 1,2*ntqgcx(1)+2*ntgx
!              write(6,*) j,q(iints+j-1)
!          enddo

          !Transform the AO dipole moment integrals into the MO basis
          call tranop(q(iints),q(iorbmo),q(iscr),1,1,1,0)   !If the penultimate argument was zero, this would transform it from MO->AO

          !We now have the symmetry packed integrals in the MO basis
          !However, this is the *full* MO basis, not the frozen one.
          !Work out which integrals we want, and calculate the core contribution
          ind = 0
          do isym = 1,8
              do j = 1,nt_frz(isym) !isize
                  do k = 1,j
                     if(tri_neci(j,k,isym,nt_frz,ntd_frz).gt.isize) then
                          call stop_all(this_routine, 'indexing error')
                      endif
                      pints(tri_neci(j,k,isym,nt_frz,ntd_frz),i) =      &
     &                  q(iints+tri_neci(j+nfrz_neci(isym),             &
     &                      k+nfrz_neci(isym),isym,nt,ntd)-1)
                  enddo
              enddo
!              write(6,*) j, q(iints+j-1)
!              pints(j,i) = q(iints+j-1)
          enddo
          znuc(i) = dipnu
          !Now to calculate core contribution:
          zcor(i) = 0.0d0
          do isym=1,8
              do j=1,nfrz_neci(isym)
                  zcor(i) = zcor(i) + 2.0d0*q(iints+                    &
     &                  tri_neci(j,j,isym,nt,ntd)-1)
              enddo
          enddo
      enddo
!      write(6,*) "z-component dipole moment integrals: "
!      do i = 1,isize
!          write(6,*) i,pints(i,3)
!      enddo
!      write(6,*) znuc(3)
!      write(6,*) zcor(3)

      return

      end subroutine GetDipMomInts_cisd

      !Gives the index in a triangular packed array
      !i must be larger than or equal to j
      integer function tri_neci(i,j,sym,symlist,symlistoff)
          implicit none
          integer :: i,j,sym,symlist(8),symlistoff(8)

          tri_neci = symlistoff(sym)+(((i-1)*i)/2)+j

      end function


#else
      !In case it doesn't like empty files
      subroutine null_wrapper
          implicit none
          character(*), parameter :: this_routine = 'null_wrapper'
          call stop_all(this_routine, 'Should not be here')
      end subroutine null_wrapper

#endif


