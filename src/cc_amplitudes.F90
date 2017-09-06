#include "macros.h" 

module cc_amplitudes 

    use SystemData, only: nbasis, nel, nOccAlpha, nOccBeta
    use detbitops, only: get_bit_excitmat
    use FciMCData, only: totwalkers, ilutref, currentdets, AllNoatHf
    use bit_reps, only: extract_sign
    use constants, only: dp, lenof_sign
    use cepa_shifts, only: calc_number_of_excitations

    implicit none 

    logical :: t_cc_amplitudes = .false. 

    ! we need the order of cluster operators..
    integer :: cc_order = 0 

    integer :: n_triples = 0, n_quads = 0

    ! maybe it would be nice to have a type which encodes this information 
    ! and which gives an easy and nice way to de/encode the indices 
    ! involved.. 
    type cc_amplitude 
        integer :: order = 0
        integer, allocatable :: operators(:,:,:)
        real(dp), allocatable :: amplitudes(:)

        integer :: n_ops = 0

    contains 

        procedure :: get_ex
        procedure :: get_ind

    end type cc_amplitude

    ! and make a global cc_ops 
    type(cc_amplitude), allocatable :: cc_ops(:)

contains 

    ! and now go to the routines to calculate the number of triples and 
    ! quadrupels 
    subroutine calc_n_triples() 
        ! this routine calculates the number of "important" triples 
        ! from the samples singles and doubles amplitudes: 
        ! essentiall calulating T1 * T2 
        character(*), parameter :: this_routine = "calc_n_triples"

        integer :: i, j, ia(2,1), jk_cd(2,2)

        do i = 1, cc_ops(1)%n_ops 
            ! for each t_i^a i have to check if a double excitation is 
            ! possible on top of it.. 
            ia = cc_ops(1)%get_ex(i)

            do j = 1, cc_ops(2)%n_ops

                jk_cd = cc_ops(2)%operators(j,:,:)

                ! check if the operators fit..
                if (any(ia(1,1) == jk_cd(1,:)) .or. any(ia(2,1) == jk_cd(2,:))) then
                    cycle
                end if
                ! if it fits increase the triples counter:
                n_triples = n_triples + 1

            end do
        end do

    end subroutine calc_n_triples

    subroutine calc_n_quads 
        ! this routines calculates the number of "important" quadrupels
        ! calculating T2*T2 essentially 
        character(*), parameter :: this_routine = "calc_n_quads" 

        integer :: i, j, ij_ab(2,2), kl_cd(2,2)

        do i = 1, cc_ops(2)%n_ops 
            ij_ab = cc_ops(2)%operators(i,:,:)

            do j = i + 1, cc_ops(2)%n_ops 

                kl_cd = cc_ops(2)%operators(j,:,:)

                ! can i just: any(a==b)
                if (.not. any(ij_ab == kl_cd)) then 

                    n_quads = n_quads + 1

                end if
            end do
        end do

    end subroutine calc_n_quads 

    ! do it other way.. only store the possible non-zero cluster operators! 
    subroutine fill_cc_amplitudes() 
        ! design decisions: since the singles are not so many in general 
        ! and because i could need them to correct the doubles amplitudes
        ! store all of the possible ones! and encode them specifically through 
        ! (i) and (a) 
        character(*), parameter :: this_routine = "fill_cc_amplitudes"
        
        integer :: idet, ic, ex(2,cc_order), n_singles, n_doubles, j, i
        integer :: ia, ib, ja, jb, ind
        integer, allocatable :: n_excits(:)

        HElement_t(dp) :: sign_tmp(lenof_sign) 

        ! for this it is helpful to have an upper limit of the number of 
        ! possible amplitudes, but just do it for the singles for now..
        allocate(n_excits(1)) 

        n_excits = calc_number_of_excitations(nOccAlpha, nOccBeta, 1, & 
            nbasis/2)

        allocate(cc_ops(cc_order))

        ! and do a nice initialization depending on the order 
        do i = 1, cc_order 

            cc_ops(i)%order = i 

        end do

        allocate(cc_ops(1)%amplitudes(n_excits(1)))

!         cc_ops(i)%operators = 0 
        cc_ops(1)%amplitudes = 0.0_dp
        cc_ops(1)%n_ops = n_excits(1)

        ! first figure out the number of double and fill the singles 
        ! amplitudes! 
        n_doubles = 0 
        n_singles = 0

        do idet = 1, int(totwalkers) 
            ic = cc_order 
            call get_bit_excitmat(iLutRef, CurrentDets(:,idet), ex, ic)
            select case (ic) 
            case (1) 
                n_singles = n_singles + 1 
                call extract_sign(CurrentDets(:,idet), sign_tmp)
                ind = cc_ops(1)%get_ind(ex(1,1),ex(2,1))

                ! for now only do it for single runs.. do i need the normalising?
                cc_ops(1)%amplitudes(ind) = sign_tmp(1) / AllNoatHf(1)

            case (2) 
                n_doubles = n_doubles + 1

            case (0)

            case default 
                call stop_all(this_routine, "cc_order > 2 not yet implemented!")

            end select
        end do

        if (allocated(cc_ops(2)%operators))     deallocate(cc_ops(2)%operators)
        if (allocated(cc_ops(2)%amplitudes))    deallocate(cc_ops(2)%amplitudes)

        allocate(cc_ops(2)%operators(n_doubles, 2, 2))
        allocate(cc_ops(2)%amplitudes(n_doubles))

        cc_ops(2)%operators = 0
        cc_ops(2)%amplitudes = 0.0_dp
        cc_ops(2)%n_ops = n_doubles

        j = 1 
        do idet = 1, int(totwalkers)

            ic = cc_order 
            call get_bit_excitmat(iLutRef, CurrentDets(:,idet), ex, ic)

            select case (ic) 
            case (2) 
                call extract_sign(CurrentDets(:,idet), sign_tmp)

                cc_ops(2)%operators(j,:,:) = ex 
                ia = cc_ops(1)%get_ind(ex(1,1),ex(2,1))
                ib = cc_ops(1)%get_ind(ex(1,1),ex(2,2))
                ja = cc_ops(1)%get_ind(ex(1,2),ex(2,1))
                jb = cc_ops(1)%get_ind(ex(1,2),ex(2,2))

                associate(cc_1 => cc_ops(1)%amplitudes)
                    cc_ops(2)%amplitudes(j) = sign_tmp(i)/AllNoatHf(1) + & 
                        cc_1(ja)*cc_1(ib) - cc_1(ia)*cc_1(jb)
                end associate

            end select 
        end do

    end subroutine fill_cc_amplitudes

    function get_ex(this, ind) result(ex) 
        ! this function gives the specific excitation, if present in the 
        ! cc_amplitudes, which are encoded in a linear fashion 
        ! with the convention that all the electron and orbital indices are 
        ! always provided in an ordered(from lowest to highest) fashion
        class(cc_amplitude), intent(in) :: this
        integer, intent(in) :: ind
        integer :: ex(2, this%order)
        character(*), parameter :: this_routine = "get_ex"

        integer :: ij, ab, cum, i, j, a, b, nij

        ASSERT(ind > 0) 
        select case (this%order) 
        case (1) 
            ASSERT(ind <= nel * (nbasis - nel))

            ! it is a single excition so this is easy to decompose 
            ! this decomposition depends on the encoding below!
            ex(2,1) = mod(ind-1,(nbasis - nel))
            ! lets hope the integer division is done correctly.. on all compilers
            ex(1,1) = int((ind - ex(2,1))/(nbasis - nel)) + 1
            ! and modify by nel the orbital index again to get the real 
            ! orbital index! 
            ex(2,1) = ex(2,1) + nel + 1

        case (2) 
            call stop_all(this_routine, "still buggy for double excitations!")
            ! first we have to get ij, ab back: 

            nij = nel * (nel - 1) / 2

            ab = mod(ind - 1, nij)
            ij = (ind - ab)/(nij) + 1 

            ab = ab + 1 

            ! then we have to get ij and ab from those indices in the same 
            ! way 
            ! but here it gets more tricky because we can just do the mod.. 
            ! maybe i have to store the indices in the end.. 
            i = 1
            cum = (nel - i)

            print *, ""
            print *, "nij: ", nij 
            print *, "ab: ", ab 
            print *, "ij: ", ij 
            print *, "cum: ", cum

            ! do a stupid sum..
            do while(cum < ij .and. i <= nel)
                i = i + 1
                cum = cum + (nel - i)
            end do

            j = ij + i - (cum - nel + i)

            ! and the same for ab 
            a = 1 
            cum = (nbasis - nel - a)
            do while (cum < ab .and. a <= (nbasis - nel)) 
                a = a + 1
                cum = cum + (nbasis - nel - a)
            end do

            b = ab + a - (cum - (nbasis - nel) + a)

            ex(1,:) = [i,j]
            ex(2,:) = [a,b] + nel

            print *, "(i,j): ", i,j 
            print *, "(a,b): ", a, b
        case default 
            call stop_all(this_routine, "higher than cc_order 2 not yet implemented!")

        end select

    end function get_ex

    function get_ind(this, elec_ind, orb_ind) result(ind)
        ! depending on the cc_order this encodes all the electron and 
        ! orbital indices of an excitation in a unique linear fashion 
        ! we can assume the electron and orbital indices to be ordered from 
        ! highest to lowest
        class(cc_amplitude), intent(in) :: this 
        integer, intent(in) :: elec_ind(this%order), orb_ind(this%order) 
        integer :: ind

        character(*), parameter :: this_routine = "get_ind"

        integer :: ij, ab, nij, orb(2)

        ind = -1 

        ! to i want to access this with invalid excitations?
        ASSERT(all(elec_ind > 0))
        ASSERT(all(elec_ind <= nel))
        ASSERT(all(orb_ind > nel))
        ASSERT(all(orb_ind <= nbasis))

        ! and assert ordered inputs.. 
        ! or do we want to order it here, if it is not ordered? 
        ! tbd!
        ASSERT(minloc(elec_ind,1) == 1) 
        ASSERT(maxloc(elec_ind,1) == this%order)

        ASSERT(minloc(orb_ind,1) == 1) 
        ASSERT(maxloc(orb_ind,1) == this%order)

        select case (this%order)
        case (1) 
            ! single excitation 
            ! the elec_ind goes from 1 to nel and the orb_ind goes from 
            ! nel + 1 to nbasis (has nbasis - nel) possible values 
            ! can we assume a closed shell ordered reference? 
            ! with the nel lowest orbitals occupied? 
            ! otherwise this is a bit more tricky.. 
            ind = (elec_ind(1) - 1) * (nbasis - nel) + (orb_ind(1) - nel)

        case (2) 
            call stop_all(this_routine, "still buggy for double excitations!")
            ! double excitation 
            ! first encode the ij electron index in a linear fashion
            ASSERT(elec_ind(1) < elec_ind(2))
            ASSERT(orb_ind(1) < orb_ind(2))

            orb = orb_ind - nel 

            ij = (elec_ind(1) - 1) * nel - elec_ind(1)*(elec_ind(1) - 1)/2  + (elec_ind(2) - elec_ind(1))
!             ij = (elec_ind(2) - 1)*elec_ind(2)/2 + elec_ind(1) - 1
            ! i shift the orb_indices by nel.. 
!             ab = (orb_ind(1) - nel - 1)*(nbasis - nel) + (orb_ind(2) - orb_ind(1))
            ab = (orb(1) - 1) *(nbasis - nel) - orb(1)*(orb(1) - 1)/2 + (orb(2) - orb(1))
!             ab = (orb_ind(2) - nel - 1)*(orb_ind(2)-nel)/2 + orb_ind(1) - nel -1

            nij = nel * (nel - 1) / 2 
            ind = (ij - 1) * nij + ab 

        case default 
            call stop_all(this_routine, "higher than cc_order 2 not yet implemented!")
        end select

    end function get_ind

!     subroutine calc_cc_amplitudes
!         character(*), parameter :: this_routine = "calc_cc_amplitudes" 
! 
!         integer :: idet, i
!         integer, allocatable :: n_excits(:)
!         type(cc_amplitude), allocatable :: cc_amp(:)
! 
!         ! i want to calculate the amplitudes up to a certain order given 
!         ! in the input 
! 
!         ! for this it is helpful to have an upper limit of the number of 
!         ! possible amplitudes 
!         allocate(n_excits(cc_order)) 
! 
!         n_excits = calc_number_of_excitations(nOccAlpha, nOccBeta, cc_order, & 
!             nbasis/2)
! 
!         allocate(cc_amp(cc_order))
! 
!     end subroutine calc_cc_amplitudes


    subroutine dongxia_amplitudes 
 
! dongxia test
        integer :: IC,idet,ierr,Nvirt,ind,idi,idj,idk,ida,idb,idc
        integer :: Ex(2,3)
        real(dp), allocatable :: C1(:),C2(:)
        real(dp), allocatable :: T1(:),T2(:)
        real(dp) :: sign_tmp(lenof_sign),C3,T3

          Ex = 0

! Dongxia 1.8.2017
! to obtain CC amplitudes from walker numbers
        Nvirt = Nbasis - Nel
        allocate(C1(Nvirt*nel),stat=ierr)
          c1 = 0.0_dp
        allocate(C2(Nvirt*(Nvirt+1)*nel*(nel+1)/4),stat=ierr)
          c2 = 0.0_dp
        allocate(T1(Nvirt*nel),stat=ierr)
          t1 = 0.0_dp
        allocate(T2(Nvirt*(Nvirt+1)*nel*(nel+1)/4),stat=ierr)
          t2 = 0.0_dp
        write(6,*) 'dongxia testing C1 and T1'
! look for C1 and obtain T1
        do idet=1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if (IC==1) then
              call extract_sign(CurrentDets(:,idet),sign_tmp)
              idi = Ex(1,1)
              ida = Ex(2,1) - Nel
              ind = ind1(idi,ida)
              C1(ind) = sign_tmp(1) / AllNoatHf(1)
              T1(ind) = C1(ind)
              if(abs(c1(ind)) > 1.0d-3) &
               write(6,'(3I5,F20.12)') ex(1,1),ex(2,1),ind, C1(ind)
           end if
        end do
! look for C2, and obtain T2
        write(6,*)'dongxia testing double excitations'
        do idet=1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if (IC==2) then
              call extract_sign(CurrentDets(:,idet),sign_tmp)
              idi = Ex(1,1)
              ida = Ex(2,1)-Nel
              idj = Ex(1,2)
              idb = Ex(2,2)-Nel
              ind = ind2(idi,idj,ida,idb)
              C2(ind) = sign_tmp(1)/AllNoatHF(1)
              T2(ind) = C2(ind)+T1(ind1(idj,ida))*T1(ind1(idi,idb)) &
                        -T1(ind1(idi,ida))*T1(ind1(idj,idb))
              if(abs(c2(ind)) > 1.0d-3) &
               write(6,'(4I10,2F20.12)') idi,idj,ida,idb,C2(ind),T2(ind)
           end if
        end do
! look for C3, and compare it with T3 (from T1^3, T1T2)
        write(6,*)'dongxia testing triple excitations'
        do idet = 1, int(totwalkers)
           IC = 4
           call get_bit_excitmat(iLutRef,CurrentDets(:,idet),ex,IC)
           if(IC==3)then
             c3 = 0.0_dp
             t3 = 0.0_dp
             call extract_sign(CurrentDets(:,idet),sign_tmp)
             idi = Ex(1,1)
             idj = Ex(1,2)
             idk = Ex(1,3)
             ida = Ex(2,1) - Nel
             idb = Ex(2,2) - Nel
             idc = Ex(2,3) - Nel
             C3 = sign_tmp(1)/AllNoatHF(1)
! calculate t3 (by t1^3, t1^t2 and t2^t1)
! contribution from t1^3/6.
             t3 = t3 &
                 +t1(ind1(idi,ida))*t1(ind1(idj,idb))*t1(ind1(idk,idc)) &
                 +t1(ind1(idi,idb))*t1(ind1(idj,idc))*t1(ind1(idk,ida)) &
                 +t1(ind1(idi,idc))*t1(ind1(idj,ida))*t1(ind1(idk,idb)) &
                 -t1(ind1(idi,ida))*t1(ind1(idj,idc))*t1(ind1(idk,idb)) &
                 -t1(ind1(idi,idc))*t1(ind1(idj,idb))*t1(ind1(idk,ida)) &
                 -t1(ind1(idi,idb))*t1(ind1(idj,ida))*t1(ind1(idk,idc)) 
! contribution from t1^t2 (and t2^t1, which is the same)
             t3 = t3 &
                 +t1(ind1(idi,ida))*t2(ind2(idj,idk,idb,idc)) &
                 -t1(ind1(idi,idb))*t2(ind2(idj,idk,ida,idc)) &
                 +t1(ind1(idi,idc))*t2(ind2(idj,idk,ida,idb)) &
                 -t1(ind1(idj,ida))*t2(ind2(idi,idk,idb,idc)) &
                 +t1(ind1(idj,idb))*t2(ind2(idi,idk,ida,idc)) &
                 -t1(ind1(idj,idc))*t2(ind2(idi,idk,ida,idb)) &
                 +t1(ind1(idk,ida))*t2(ind2(idi,idj,idb,idc)) &
                 -t1(ind1(idk,idb))*t2(ind2(idi,idj,ida,idc)) &
                 +t1(ind1(idk,idc))*t2(ind2(idi,idj,ida,idb))  

                      
             if(abs(c3)+abs(t3) > 1.0d-3) &
              write(6,'(3I5,4X,3I5,3F20.12)')idi,idj,idk,ida,idb,idc,c3,t3,t3/c3-1
           end if
        end do

        deallocate(T2,stat=ierr)
        deallocate(T1,stat=ierr)
        deallocate(C2,stat=ierr)
        deallocate(C1,stat=ierr)
! end
    end subroutine dongxia_amplitudes

      pure integer function ind1(i,a)

      implicit none
      integer, intent(in) :: i, a

      ind1 = (i-1)*(nbasis-nel)+a
      
      return
      end function ind1

      pure integer function ind2(i,j,a,b)

      implicit none
      integer, intent(in) :: i,j,a,b
      integer :: ij, ab, nvirt,nab
      
      ij = (j-1)*j/2 + i
      ab = (b-1)*b/2 + a
      nvirt = nbasis - nel
      nab = nvirt*(nvirt+1)/2
       
      ind2 = (ij-1)*nab + ab

      return
      end function ind2      

end module cc_amplitudes
