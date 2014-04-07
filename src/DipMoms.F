#ifdef MOLPRO
      subroutine GetDipMomInts(pints,isize,znuc,zcor)
      implicit double precision(a-h,o-z)
      integer isize
      double precision pints(isize,3)
      double precision znuc(3),zcor(3)
      double precision origin(3)
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
            stop 'Size of symmetry packed arrays wrong: GetDipMomInts'
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
      else
          !Assume that we have done a prior CISD calculation for the dipole moment integrals.
          !See: ccsd/ccmp2.F

          !First, get the MO coefficients, since we will need to transform the integrals ourselves
          call excom(2)
          iorbmo=icorr(ntqgcx(1))
          call cread(q(iorbmo),ntqgcx(1),kopfil,nwf,nwf)
!          call druckm(q(iorbmo),1,0,'orbitals',0)   !This would print out the orbitals
          iscr=icorr(2*ntqgcx(1)+2*ntgx)  !Memory for integral transforms


          iints=icorr(2*ntqgcx(1)+2*ntgx)   !Get memory pointer for integrals
          
          do i = 1,3
              !read the operators (in the AO basis)
              call read_op(q(iints),dipnu,'DM',i,1,1,origin,isyop)

              call tranop(q(iints),q(iorbmo),q(iscr),1,1,1,0)   !If the penultimate argument was zero, this would transform it from MO->AO
              !We now have the symmetry packed integrals in the MO basis
              do j = 1,isize
                  pints(j,i) = q(iints+j-1)
              enddo
              znuc(i) = dipnu
              zcor(i) = 0.0d0
          enddo
          

      endif

      return
      end subroutine GetDipMomInts

#else
      !In case it doesn't like empty files
      subroutine null_wrapper
          implicit none
          stop 'Should not be here'
      end subroutine null_wrapper

#endif


