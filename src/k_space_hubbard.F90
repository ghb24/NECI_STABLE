#include "macros.h" 

! finally write a specific routine for the k-space hubbard, where every 
! necessary functionality for the k-space/momentum space hubbard is brought 
! together. since i want to implement a transcorrelation factor for the 
! k-space hubbard, which would necessitate a cumulative sum exciation creation 
! i decided it is better to start a new branch, instead of hacking into the 
! old branches of the code. 
! in the end i also want to combine this fully with the lattice_mod implementation
! so i can easily decide which lattice to choose and then decide if we want to 
! use a real or momentum space basis (and in the future maybe even wavelets) 
module k_space_hubbard 
    use SystemData, only: t_lattice_model, t_k_space_hubbard, t_trans_corr, & 
                    trans_corr_param, t_trans_corr_2body, trans_corr_param_2body, & 
                    nel, tHPHF, nOccBeta, nOccAlpha, nbasis, tLatticeGens, tHub, &
                    omega, bhub, nBasisMax, G1, BasisFN, NullBasisFn, TSPINPOLAR, & 
                    treal, ttilt, tExch, ElecPairs, MaxABPairs, Symmetry, SymEq, &
                    t_new_real_space_hubbard, SymmetrySize, tNoBrillouin, tUseBrillouin, &
                    excit_cache, t_uniform_excits, brr, uhub, lms

    use lattice_mod, only: get_helement_lattice_ex_mat, get_helement_lattice_general, &
                           determine_optimal_time_step, lattice, sort_unique, lat, &
                           dispersion_rel_cached, init_dispersion_rel_cache, &
                           epsilon_kvec

    use procedure_pointers, only: get_umat_el, generate_excitation

    use constants, only: n_int, dp, EPS, bits_n_int, int64

    use bit_rep_data, only: NIfTot, nifd

    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, ilut_lt, ilut_gt

    use real_space_hubbard, only: lat_tau_factor

    use fcimcdata, only: tsearchtau, tsearchtauoption, pDoubles, pParallel, &
                         excit_gen_store_type, pSingles

    use CalcData, only: tau, t_hist_tau_search, t_hist_tau_search_option, &
                        p_doubles_input, p_parallel_input, t_fill_frequency_hists

    use dsfmt_interface, only: genrand_real2_dsfmt

    use util_mod, only: binary_search_first_ge, binary_search

    use get_excit, only: make_double

    use OneEInts, only: GetTMatEl, tmat2d

    use sltcnd_mod, only: sltcnd_0

    use sym_mod, only: RoundSym, AddElecSym, SetupSym, lChkSym, mompbcsym, & 
                       TotSymRep, GenMolpSymTable, SymProd, gensymstatepairs

    use SymExcitDataMod, only: KPointToBasisFn, ScratchSize, SpinOrbSymLabel, &
                               SymTableLabels, SymInvLabel, SymLabelList2, SymLabelCounts2, & 
                               OrbClassCount, ktotal

    use sym_general_mod, only: ClassCountInd

    use sort_mod, only: sort 

    use IntegralsData, only: UMat

    use bit_reps, only: decode_bit_det

    use global_utilities, only: LogMemDealloc, LogMemAlloc

    use SymData, only: nSymLabels, SymClasses, Symlabels

    use SymData, only: nSym,SymConjTab

    use SymData, only: tAbelian,SymTable

    use SymData, only: tagSymConjTab,tagSymClasses,tagSymLabels

    use SymData, only: tagSymTable

    use ParallelHelper, only: iProcIndex, root

    use lattice_models_utils, only: make_ilutJ, get_ispn, get_orb_from_kpoints, &
                                    create_all_dets, find_minority_spin, & 
                                    get_orb_from_kpoints_three, swap_excitations, &
                                    pick_spin_opp_elecs, pick_from_cum_list, & 
                                    pick_spin_par_elecs, pick_three_opp_elecs

    use unit_test_helpers, only: print_matrix
                                    

    implicit none 

    integer, parameter :: ABORT_EXCITATION = 0
    integer, parameter :: N_DIM = 3

    real(dp) :: three_body_prefac = 0.0_dp
    real(dp) :: prefac_test = 2.0_dp
    real(dp), allocatable :: umat_cache_kspace(:,:)
    real(dp) :: n_opp(-1:1) = 0.0_dp

    ! temporary flag for the j optimization
    logical :: t_symmetric = .true.

    real(dp), allocatable :: two_body_transcorr_factor_matrix(:,:), &
                             three_body_const_mat(:,:,:)
    ! i especially need an interface for the matrix element calculation to 
    ! implement the transcorrelated hamiltonian 
    interface get_helement_k_space_hub
        module procedure get_helement_k_space_hub_ex_mat
        module procedure get_helement_k_space_hub_general
    end interface get_helement_k_space_hub

    interface get_one_body_diag
        module procedure get_one_body_diag_sym
        module procedure get_one_body_diag_kvec
    end interface get_one_body_diag

    interface same_spin_transcorr_factor
        module procedure same_spin_transcorr_factor_kvec
        module procedure same_spin_transcorr_factor_ksym
    end interface same_spin_transcorr_factor

    interface three_body_transcorr_fac
        module procedure three_body_transcorr_fac_kvec
        module procedure three_body_transcorr_fac_ksym
    end interface three_body_transcorr_fac

    interface two_body_transcorr_factor
        module procedure two_body_transcorr_factor_kvec
        module procedure two_body_transcorr_factor_ksym
    end interface two_body_transcorr_factor

    interface three_body_rpa_contrib
        module procedure rpa_contrib_ksym
        module procedure rpa_contrib_kvec
    end interface three_body_rpa_contrib

    interface two_body_contrib
        module procedure two_body_contrib_ksym
        module procedure two_body_contrib_kvec
    end interface two_body_contrib

    interface three_body_exchange_contrib
        module procedure exchange_contrib_ksym
        module procedure exchange_contrib_kvec
    end interface three_body_exchange_contrib

contains 

    subroutine setup_symmetry_table() 
        ! implement a new symmetry setup to decouple it from the 
        ! old hubbard.F code.. 

        character(*), parameter :: this_routine = "setup_symmetry_table"

        integer :: i, j, k, l, k_i(3), k_inv(3), k_j(3), ind, kmin(3), kmax(3)

        ! the only problem could be that we reorderd the orbitals already or? 
        ! so G1 has a different ordering than just 1, nBasis/2... 
        ! yes it really is reordered already! hm.. 
        ASSERT(associated(lat))
        ASSERT(associated(G1))

        nsym = nBasis/2
        nSymLabels = nsym

        ! copy the output from the old hubbard code: 
         WRITE(6,"(A,I3,A)") "Generating abelian symmetry table with", &
         nsym , " generators for Hubbard momentum" 
         if (allocated(SymLabels)) then
             write (6,'(a/a)') &
                     'Warning: symmetry info already allocated.', &
                     'Deallocating and reallocating.' 
             deallocate(SymLabels)
             call LogMemDealloc(this_routine,tagSymLabels)
         end if
         allocate(SymLabels(nSym))
         call LogMemAlloc('SymLabels',nSym,SymmetrySize,this_routine, tagSymLabels)
         if (associated(SymClasses)) then
             deallocate(SymClasses)
             call LogMemDealloc(this_routine,tagSymClasses)
         end if
         ! [W.D] 
         ! why is symclasses allocated to nBasis? it is actually only used 
         ! and filled up to nBasis/2, also in the rest of the code..
         allocate(SymClasses(nBasis/2))
         call LogMemAlloc('SymClasses',nBasis,4,this_routine,tagSymClasses)

        ! for some strange reason the (0,0,0) k-vector is treated 
        ! special in the old implementation.. since it is its own 
        ! symmetry inverse.. but this might be not the case in the 
        ! general lattices or? there can be other states which are 
        ! also its own inverse or? 

        ! so try to change that here.. 
        ! ok i should still treat the gamma point special and set its 
        ! sym label to 1 and maybe also order the symmetries by energy? 
       
        if (allocated(SymTable)) then
            deallocate(SymTable)
            call LogMemDealloc(this_routine,tagSymTable)
        end if
        allocate(SymTable(nSym,nSym))
        call LogMemAlloc('SymTable',nSym**2,SymmetrySize,this_routine,tagSymTable)
        if (allocated(SymConjTab)) then
            deallocate(SymConjTab)
            call LogMemDealloc(this_routine,tagSymConjTab)
        end if
        allocate(SymConjTab(nSym))
        call LogMemAlloc('SymConjTable',nSym,4,this_routine, tagSymConjTab)

        SYMTABLE = Symmetry(0)

        tAbelian = .false.

        ! also setup the stuff contained in the lattice class. 
        ASSERT(associated(lat))

        if (allocated(lat%k_to_sym))   deallocate(lat%k_to_sym)
        if (allocated(lat%sym_to_k))   deallocate(lat%sym_to_k)
        if (allocated(lat%mult_table)) deallocate(lat%mult_table)
        if (allocated(lat%inv_table))  deallocate(lat%inv_table)

        allocate(lat%sym_to_k(lat%get_nsites(),3))
        allocate(lat%mult_table(lat%get_nsites(), lat%get_nsites()))
        allocate(lat%inv_table(lat%get_nsites()))

        ! i have to setup the symlabels first ofc.. 
        do i = 1, lat%get_nsites() 
            ind = get_spatial(brr(2*i))
            SymClasses(ind) = i
!             SymClasses = [( i, i = 1, lat%get_nsites())]
            ! and also just encode the symmetry labels as integers, instead of 
           ! 2^(k-1), to be able to treat more than 64 orbitals (in the old 
            ! implementation, an integer overflow happened in this case!)
            SymLabels(ind)%s = i
            ! i also need it in G1: 
            ! is G1 already ordered by energy?? 

            call lat%set_sym(ind,i)

        end do

        kmin = 0 
        kmax = 0 
        do i = 1, lat%get_nsites() 
            k_i = lat%get_k_vec(i)
            do j = 1, lat%get_ndim() 
                if (k_i(j) < kmin(j)) kmin(j) = k_i(j)
                if (k_i(j) > kmax(j)) kmax(j) = k_i(j) 
            end do
        end do

        allocate(lat%k_to_sym(kmin(1):kmax(1),kmin(2):kmax(2),kmin(3):kmax(3)))

        lat%k_to_sym = 0

        ! now find the inverses: 
        do i = 1, lat%get_nsites() 
            ! and also just encode the symmetry labels as integers, instead of 
!            ! 2^(k-1), to be able to treat more than 64 orbitals (in the old 
!             ! implementation, an integer overflow happened in this case!)
!             SymLabels(i)%s = i
!             ! i also need it in G1: 
            G1(2*i-1)%Sym = SymLabels(i)
            G1(2*i)%Sym = SymLabels(i)

!             k_i = lat%get_k_vec(i) 
! 
!             k_inv = lat%map_k_vec(-k_i)
! 
!             print *, "i, k_i, k_inv: ", i, k_i, k_inv
! 
!             j = lat%get_orb_from_k_vec(k_inv)
! 
!             print *, j
            ! find the orbital of -k 
            j = lat%get_orb_from_k_vec(-lat%get_k_vec(i))

!             print *, "i, k(i): ", i, lat%get_k_vec(i)
!             print *, "j, k(j): ", j, lat%get_k_vec(j)

            ! since i have a linear encoding of the symmetries i do not need 
            ! to use SymClasses here.. 
            SymConjTab(SymClasses(i)) = SymClasses(j)

            lat%inv_table(SymClasses(i)) = SymClasses(j) 
            ! maybe also store the k_inverse of a orbital.. 
            
            lat%sym_to_k(SymClasses(i),:) = lat%get_k_vec(i)

            k_i = lat%get_k_vec(i) 
            lat%k_to_sym(k_i(1),k_i(2),k_i(3)) = SymClasses(i)

            ! and create the symmetry product of (i) with every other symmetry
            do k = 1, lat%get_nsites()
                ! i just have to add the momenta and map it to the first BZ 

!                 k_j = lat%get_k_vec(k) 

                l = lat%get_orb_from_k_vec(lat%get_k_vec(i) + lat%get_k_vec(k))
!                 l = lat%get_orb_from_k_vec(k_i + k_j)

                SymTable(SymClasses(i),SymClasses(k)) = SymLabels(l)

                lat%mult_table(SymClasses(i), SymClasses(k)) = SymClasses(l)

            end do
        end do
! 
#ifdef __DEBUG
        WRITE(6,*) "Symmetry, Symmetry Conjugate"
        do i = 1, lat%get_nsites()
            print *, i, SymConjTab(i)
        end do

        print *, "lat%sym_to_k: "
        do i = 1, lat%get_nsites() 
            print *, i, "|", SymClasses(i), "|", lat%sym_to_k(i,:) 
        end do

        print *, "lat%inv_table: "
        do i = 1, lat%get_nsites() 
            print *, i, "|", SymClasses(i), "|", lat%inv_table(i)
        end do

        print *, "lat%mult_table: "
        do i = 1, lat%get_nsites() 
            print *, lat%mult_table(i,:)
        end do

        print *, "lat%k_to_sym: "
        do i = 1, lat%get_nsites() 
            k_i = lat%get_k_vec(i) 
            print *, k_i, "|", lat%k_to_sym(k_i(1),k_i(2),k_i(3))
        end do
            
#endif

    end subroutine setup_symmetry_table

    subroutine gen_symreps()
        ! i have to figure out what exactly those symreps do and how 
        ! i should set them up.. 
        use SystemData, only: arr, brr
        use SymData, only: symreps

        integer :: i, j
        if (allocated(symreps)) deallocate(symreps) 
        allocate(symreps(2,nbasis))
        symreps = 0

        j = 0 

        print *, "arr: "
        do i = 1, nbasis
            print *, arr(i,:) 
        end do

        print *, "brr: " 
        do i = 1, nbasis 
            print *, brr(i)
        end do
    
        do i = 1, 2*lat%get_nsites()


        end do

    end subroutine gen_symreps

    function get_umat_kspace(i, j, k, l) result(hel) 
        ! simplify this get_umat function for the k-space hubbard.. 
        ! since there was a lot of unnecessary stuff going on in the other 
        ! essentially we only have to check if the momenta involved 
        ! fullfil k_k + k_l = k_i + k_j 
        integer, intent(in) :: i, j, k, l 
        HElement_t(dp) :: hel, hel2
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_umat_kspace"
#endif

        ! or just use the symtable information? 
        ! do i have to access this with the symmetry label or the orbital 
        ! number?? i think i need the symmetry labels.. 
        ! i could set this on the site level! 
        if (SymTable(lat%get_sym(i),lat%get_sym(j))%S == & 
            SymTable(lat%get_sym(k),lat%get_sym(l))%S) then 
!         if (SymTable(G1(2*i)%sym%s,G1(2*j)%sym%s)%S ==  & 
!             SymTable(G1(2*k)%sym%s,G1(2*l)%sym%s)%S) then 
            hel = Umat(1)
        else 
            hel = 0.0_dp
        end if
        
        ! old implo: 
!         if (all(lat%map_k_vec(lat%get_k_vec(i)+lat%get_k_vec(j)) == & 
!                 lat%map_k_vec(lat%get_k_vec(k)+lat%get_k_vec(l)))) then 
!             hel = UMat(1) 
!         else
!             hel = 0.0_dp
!         end if

    end function get_umat_kspace

    subroutine init_k_space_hubbard() 

        character(*), parameter :: this_routine = "init_k_space_hubbard"
        real(dp) :: tau_opt


        if (iProcIndex == root) then 
            print *, " new k-space hubbard implementation init:" 
        end if

        ! i have to set the incorrect excitaiton generator flags to false 
        tLatticeGens = .false.
        ! maybe i also need to turn off the hubbard keyword.. at this 
        ! point
        thub = .false. 

        tExch = .true. 
        tNoBrillouin = .false. 
        tUseBrillouin = .true.

        call check_k_space_hubbard_input()

        get_umat_el => get_umat_kspace

        call init_get_helement_k_space_hub()

        if (.not. tHPHF .and. .not. t_uniform_excits) then 
             generate_excitation => gen_excit_k_space_hub
         end if
        ! for more efficiency, use the uniform excitation generation
        if(t_uniform_excits) then 
            generate_excitation => gen_excit_uniform_k_space_hub
        end if 

        tau_opt = determine_optimal_time_step() 

        if (tau < EPS) then 
            if (iProcIndex == root) then 
                print *, "setting time-step to optimally determined time-step: ", tau_opt
                print *, "times: ", lat_tau_factor
            end if
            tau = lat_tau_factor * tau_opt

        else 
            if (iProcIndex == root) then 
                print *, "optimal time-step would be: ", tau_opt
                print *, "but tau specified in input!"
            end if
        end if

        ! [W.D. 7.3.2018]
        ! re-enable the tau-search and histogramming for the k-space 
        ! hubbard model with transcorrelation 
        ! since there the matrix elements are not equal for every 
        ! excitation anymore and for the Gutzwiller correlation factor 
        ! we also have additional pTriples and pParallel variables, which 
        ! we might want to adapt on-the-fly as in the ab-initio calculations
        ! and I also want to check if i got the excitation weighting correct 
        ! in the transcorrelated case or if i messed it up due to the 
        ! non-hermitian character of the hamiltonian 
        if (.not. (t_trans_corr_2body .or. t_trans_corr)) then 
            ! but in the "normal" hubbard model, still turn it off as it is 
            ! unnecessary! 
            tsearchtau = .false. 
            tsearchtauoption = .true.

            t_hist_tau_search = .false. 
            t_hist_tau_search_option = .false.

            t_fill_frequency_hists = .false.
        end if

        if (associated(lat)) then 
            call setup_nbasismax(lat) 
            call setup_g1(lat) 
            call setup_tmat_k_space(lat) 
            call setup_kPointToBasisFn(lat)
        end if

        if (t_trans_corr_2body) then 
            if (tHPHF) then
                call Stop_All(this_routine, "not yet implemented with HPHF")
            end if

            three_body_prefac = real(bhub,dp) * prefac_test * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
            ! i also have to set some generation probability parameters.. 

            pDoubles = p_doubles_input 
            ! use pSingles for triples! 
            ! BE CAREFUL and dont get confused! 
            pSingles = 1.0_dp - pDoubles
            pParallel = p_parallel_input

        end if

        if (.not. (t_trans_corr_2body .or. t_trans_corr)) then
            call initialize_excit_table()
        end if


    end subroutine init_k_space_hubbard

    subroutine initialize_excit_table()
      implicit none
      ! This cannot be a member of the lattice class because that would introduce
      ! a circulat dependency on get_offdiag_helement_k_sp_hub
      integer :: a, b, i, j, ex(2,2)
      ! nI is not needed anywhere, it is a redundant argument in the k_space_hubbard module
      integer :: nI(2)
      real(dp) :: matEL

      ! the buffer for excitations is (number of states)^3
      ! strictly speaking, we do not need to distinguish spin orbs here for simple hubbard
      ! but for transcorrelated, it might matter
      if (allocated(excit_cache)) deallocate(excit_cache)

      allocate(excit_cache(nbasis,nbasis,nbasis))
      ! loop over all pairs of orbitals
      do i = 1, nbasis
         do j = 1, nbasis
            nI = [i,j]
            ex(1,:) = [i,j]
            ! and for each, check all possible excitations
            do a = 1, nbasis
               matEl = 0.0_dp
               if(a /= i .and. a/=j) then
                  b = get_orb_from_kpoints(i,j,a)
                  if(b /= a .and. b /= i .and. b /= j) then
                     ex(2,:) = [a,b]    
                     matEl = abs(get_offdiag_helement_k_sp_hub(nI,ex,.false.))
                  endif
               endif
               ! most excitations will have a finite amplitude, so store all
               excit_cache(i,j,a) = matEl
            end do
         end do
      end do
      
    end subroutine initialize_excit_table

    subroutine check_k_space_hubbard_input()
        character(*), parameter :: this_routine = "check_k_space_hubbard_input"

        if (iProcIndex == root) then 
            print *, "checking input for k-space hubbard:" 
            !todo: find the incompatible input and abort here!
            print *, "input is fine!"
        end if

    end subroutine check_k_space_hubbard_input

    subroutine gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_k_space_hub"
#endif
        real(dp) :: p_elec, p_orb
        integer :: elecs(2), orbs(2), src(2)
        logical :: isvaliddet

        ! i first have to choose an electron pair (ij) at random 
        ! but with the condition that they have to have opposite spin! 
        call pick_spin_opp_elecs(nI, elecs, p_elec) 

        src = nI(elecs)

        call pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orb)

        if (orbs(1) == ABORT_EXCITATION) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ic = 2 

        ! and make the excitation 
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2) 

        ! i think in both the electrons and the orbitals i have twice the 
        ! probability to pick them
        ! already modified in the orbital picker.. 
        pgen = p_elec * p_orb

#ifdef __DEBUG
        if (.not. isvaliddet(nI,nel)) then
            if (nJ(1) /= 0) then 
                print *, "nI: ", nI
                print *, "nJ: ", nJ 
                print *, "src: ", ex(1,:)
                print *, "tgt: ", ex(2,:)
            end if
        end if
#endif

!         if (abs(pgen -0.5) > 1e-3) print *, "pgen: ", pgen 

    end subroutine gen_excit_k_space_hub

    subroutine gen_excit_uniform_k_space_hub(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, &
         tParity, pGen, hel, store, run)
      implicit none
      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ic, ex(2,2)
      real(dp), intent(out) :: pGen
      logical, intent(out) :: tParity
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: run

      ! not used
      HElement_t(dp), intent(out) :: hel
      integer(n_int), intent(out) :: ilutJ(0:NifTot)

      real(dp) :: p_elec, r
      integer :: i, a, b, ki(N_DIM), kj(N_DIM), ka(N_DIM), kb(N_DIM), elecs(2)
      integer, parameter :: maxTrials = 1000

      hel = 0.0_dp
      ilutJ = 0
      ic = 0

      nJ = nI

      ! first, get two electrons
      call pick_spin_opp_elecs(nI,elecs,p_elec)

      ! uniform random excit gen probability
      pGen = 1.0_dp/(nbasis-nel)*2.0_dp/(nOccAlpha*nOccBeta)

      ! try finding an allowed excitation
      do i = 1, maxTrials
         ! we do this by picking random momenta and checking if the targets are 
         ! empty
         r = genrand_real2_dSFMT()
         ! this is our random orb
         a = INT(nBasis*r)+1
         ! only empty targets are of interest
         if(IsOcc(ilutI,a)) cycle
         
         ! now, get the missing momentum
         b = get_orb_from_kpoints(nI(elecs(1)),nI(elecs(2)),a)
         ! and check if its empty and differs from a
         if(IsOcc(ilutI,b) .or. a==b) then
            ! if not, the excitation is rejected (!)
            nJ(1) = 0
            return
         endif
         
         call make_double(nI,nJ,elecs(1),elecs(2),a,b,ex,tParity)
         ic = 2
         exit
      end do
      
    end subroutine gen_excit_uniform_k_space_hub

    subroutine gen_excit_uniform_k_space_hub_transcorr(nI, ilutI, nJ, ilutJ, & 
            exFlag, ic, ex, tParity, pgen, hel, store, run) 
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot) 
        integer, intent(out) :: nJ(nel), ic, ex(2,3) 
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        real(dp), intent(out) :: pgen 
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel 
        type(excit_gen_store_type), intent(inout), target :: store 
        integer, intent(in), optional :: run 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "gen_excit_uniform_k_space_hub_transcorr" 
#endif
        integer :: temp_ex(2,3) , elecs(3), ispn, i, a, b, c, src(3), sum_ms, spin
        real(dp) :: p_elec, p_orb, p_orb_a
        integer, parameter :: max_trials = 1000

        hel = 0.0_dp 
        ilutJ = 0_n_int 
        ic = 0
        nJ(1) = 0
        elecs = 0

        if (genrand_real2_dsfmt() < pDoubles) then 

            ! double excitation
            ic = 2

            if (genrand_real2_dsfmt() < pParallel) then 

                ! pick two spin-parallel electrons 
                call pick_spin_par_elecs(nI, elecs(1:2), p_elec, ispn)

                ! i have to figure out the probabilty of two spin-parallel
                if (ispn == 1) then
                    spin = 1
                    p_orb = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 
                else if (ispn == 3) then 
                    spin = 2
                    p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)
#ifdef __DEBUG
                else 
                    call stop_all(this_routine, "no parallel spin!")
#endif
                end if

                ! pick 2 holes now 
                do i = 1, max_trials

                    a = 2*int(nbasis/2 * genrand_real2_dsfmt()) + spin

                    if (IsOcc(ilutI, a)) cycle 

                    b = get_orb_from_kpoints(nI(elecs(1)), nI(elecs(2)), a)

                    ! do we have to reject or can we cycle if not fitting?
                    ! a == b test has to be here for the spin-parallel 
                    ! excitations!
                    if (IsOcc(ilutI,b) .or. a == b) then  
                        nJ(1) = 0
                        return 
                    end if

                    call make_double(nI, nJ, elecs(1), elecs(2), a, b, ex, tParity)

                    ilutJ = make_ilutJ(ilutI, ex, 2)
                    exit 
                end do

                ! times 2, since both ab, ba orders are possible
                pgen = p_elec * p_orb * pDoubles * pParallel * 2.0_dp

            else 

                call pick_spin_opp_elecs(nI, elecs(1:2), p_elec) 

                p_orb = 2.0_dp / real(nbasis - nel, dp)

                pgen = p_elec * p_orb * pDoubles * (1.0_dp - pParallel)

                ! pick 2 holes now 
                do i = 1, max_trials

                    a = int(nbasis * genrand_real2_dsfmt()) + 1 

                    if (IsOcc(ilutI, a)) cycle 

                    b = get_orb_from_kpoints(nI(elecs(1)), nI(elecs(2)), a)

                    ! do we have to reject or can we cycle if not fitting?
                    ! a == b test has to be here for the spin-parallel 
                    ! excitations!
                    if (IsOcc(ilutI,b) .or. a == b) then  
                        nJ(1) = 0
                        return 
                    end if

                    call make_double(nI, nJ, elecs(1), elecs(2), a, b, ex, tParity)

                    ilutJ = make_ilutJ(ilutI, ex, 2)
                    exit 
                end do
            end if

        else
            ! triple excitation 
            ic = 3 

            call pick_three_opp_elecs(nI, elecs, p_elec, sum_ms) 
            src = nI(elecs)

            ASSERT(sum_ms == -1 .or. sum_ms == 1)

            call pick_a_orbital_hubbard(ilutI, a, p_orb_a, sum_ms)

            ! and now i have to pick orbital b and fitting c in a uniform 
            ! way.. i hope this still works with the probabilities 
            ! if A is beta, we need to pick a alpha B uniformly and vv.
            if (is_beta(a)) then 
                p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp) 
                ! also use a spin to specify the spin-orbital
                ! is a is beta we want an alpha -> so add +1
                ispn = 0
            else 
                p_orb = 1.0_dp / real(nBasis/2 - nOccBeta, dp)
                ispn = 1
            end if

            ! times 2 since BC <> CB is both possible 
            pgen = p_elec * p_orb * p_orb_a * (1.0_dp - pDoubles) * 2.0_dp
            do i = 1, max_trials

                b = 2 * (1 + int(genrand_real2_dsfmt() * nbasis/2)) - ispn

                if (IsOcc(ilutI,b)) cycle 

                c = get_orb_from_kpoints_three(src, a, b) 

                if (IsOcc(ilutI,c) .or. b == c) then 
                    nJ(1) = 0
                    return
                end if

                call make_triple(nI, nJ, elecs, [a,b,c], ex, tParity)

                ilutJ = make_ilutJ(ilutI, ex, 3) 
                exit 
            end do
        end if

#ifdef __DEBUG
        if (abs(pgen - calc_pgen_k_space_hubbard_uniform_transcorr(nI, ilutI, ex, ic))>EPS) then
            print *, "nI: ", nI
            print *, "nJ: ", nJ
            print *, "ic: ", ic
            print *, "calc. pgen: ",calc_pgen_k_space_hubbard_uniform_transcorr(nI, ilutI, ex, ic)
            print *, "prd. pgen: ", pgen
            call stop_all(this_routine, "pgens wrong!")
        end if
#endif

    end subroutine gen_excit_uniform_k_space_hub_transcorr 

    subroutine gen_excit_k_space_hub_transcorr (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic
        integer, intent(out) :: ex(2,3)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_k_space_hub_transcorr"
#endif
        integer :: temp_ex(2,3) 

        if (genrand_real2_dsfmt() < pDoubles) then 
            if (genrand_real2_dsfmt() < pParallel) then 
                ! do a parallel triple excitation, coming from the triples..
                call gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
                ic = 2
                pgen = pgen * pDoubles * pParallel
!                 if (nJ(1) /= 0) then 
!                     print *, "parallel: (",nI,") -> (", nJ, ")"
!                 end if
            else 
                ! do a "normal" hubbard k-space excitation 
                call gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

                pgen = pgen * pDoubles * (1.0_dp - pParallel)

            end if 
        else 
            ! otherwise to a triple.. 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen) 
            ic = 3 
            pgen = pgen * (1.0_dp - pDoubles)

!                 if (nJ(1) /= 0) then 
!                     print *, "triple: (",nI,") -> (", nJ, ")"
!                 end if

        end if

    end subroutine gen_excit_k_space_hub_transcorr

    ! make an exact copy of the transcorrelation excitation generator to 
    ! run the stochastic test driver on it! so it must have the same 
    ! interface as the other excitation generators!
    subroutine gen_excit_uniform_k_space_hub_test(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

        implicit none

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic
        integer, intent(out) :: ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_uniform_k_space_hub_test"
#endif
        integer :: temp_ex(2,3) , elecs(3), ispn, i, a, b, c, src(3), sum_ms, spin
        real(dp) :: p_elec, p_orb, p_orb_a
        integer, parameter :: max_trials = 1000

        hel = 0.0_dp 
        ilutJ = 0_n_int 
        ic = 0
        nJ(1) = 0
        elecs = 0

        if (genrand_real2_dsfmt() < pDoubles) then 

            ! double excitation
            ic = 2

            if (genrand_real2_dsfmt() < pParallel) then 

                ! pick two spin-parallel electrons 
                call pick_spin_par_elecs(nI, elecs(1:2), p_elec, ispn)

                ! i have to figure out the probabilty of two spin-parallel
                if (ispn == 1) then
                    spin = 1
                    p_orb = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 
                else if (ispn == 3) then 
                    spin = 2
                    p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)
#ifdef __DEBUG
                else 
                    call stop_all(this_routine, "no parallel spin!")
#endif
                end if

                do i = 1, max_trials
                    ! ispn == 1 means both beta 
                    a = 2*int(nbasis/2 * genrand_real2_dsfmt()) + spin 

                    if (IsOcc(ilutI, a)) cycle 

                    b = get_orb_from_kpoints(nI(elecs(1)), nI(elecs(2)), a)

                    ! do we have to reject or can we cycle if not fitting?
                    ! a == b test has to be here for the spin-parallel 
                    ! excitations!
                    if (IsOcc(ilutI,b) .or. a == b) then  
                        nJ(1) = 0
                        return 
                    end if

                    call make_double(nI, nJ, elecs(1), elecs(2), a, b, ex, tParity)

                    ilutJ = make_ilutJ(ilutI, ex, 2)
                    exit 
                end do

                ! times 2, since both ab, ba orders are possible
                pgen = p_elec * p_orb * pDoubles * pParallel * 2.0_dp

            else 

                call pick_spin_opp_elecs(nI, elecs(1:2), p_elec) 

                p_orb = 2.0_dp / real(nbasis - nel, dp)

                pgen = p_elec * p_orb * pDoubles * (1.0_dp - pParallel)

                ! pick 2 holes now 
                do i = 1, max_trials

                    a = int(nbasis * genrand_real2_dsfmt()) + 1 

                    if (IsOcc(ilutI, a)) cycle 

                    b = get_orb_from_kpoints(nI(elecs(1)), nI(elecs(2)), a)

                    ! do we have to reject or can we cycle if not fitting?
                    ! a == b test has to be here for the spin-parallel 
                    ! excitations!
                    if (IsOcc(ilutI,b) .or. a == b) then  
                        nJ(1) = 0
                        return 
                    end if

                    call make_double(nI, nJ, elecs(1), elecs(2), a, b, ex, tParity)

                    ilutJ = make_ilutJ(ilutI, ex, 2)
                    exit 
                end do
            end if
        else
            ! triple excitation 
            ic = 3 

            call pick_three_opp_elecs(nI, elecs, p_elec, sum_ms) 
            src = nI(elecs)

            ASSERT(sum_ms == -1 .or. sum_ms == 1)

            call pick_a_orbital_hubbard(ilutI, a, p_orb_a, sum_ms)

            ! and now i have to pick orbital b and fitting c in a uniform 
            ! way.. i hope this still works with the probabilities 
            ! if A is beta, we need to pick a alpha B uniformly and vv.
            if (is_beta(a)) then 
                p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp) 
                ! also use a spin to specify the spin-orbital
                ! is a is beta we want an alpha -> so add +1
                ispn = 0
            else 
                p_orb = 1.0_dp / real(nBasis/2 - nOccBeta, dp)
                ispn = 1
            end if

            ! times 2 since BC <> CB is both possible 
            pgen = p_elec * p_orb * p_orb_a * (1.0_dp - pDoubles) * 2.0_dp
            do i = 1, max_trials

                b = 2 * (1 + int(genrand_real2_dsfmt() * nbasis/2)) - ispn

                if (IsOcc(ilutI,b)) cycle 

                c = get_orb_from_kpoints_three(src, a, b) 

                if (IsOcc(ilutI,c) .or. b == c) then 
                    nJ(1) = 0
                    return
                end if

                call make_triple(nI, nJ, elecs, [a,b,c], temp_ex, tParity)

                ilutJ = make_ilutJ(ilutI, temp_ex, 3) 
                exit 
            end do
        end if

    end subroutine gen_excit_uniform_k_space_hub_test

    ! make an exact copy of the transcorrelation excitation generator to 
    ! run the stochastic test driver on it! so it must have the same 
    ! interface as the other excitation generators!
    subroutine gen_excit_k_space_hub_transcorr_test (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)


        implicit none

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic
        integer, intent(out) :: ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_k_space_hub_transcorr_test"
#endif
        integer :: temp_ex(2,3) 

        if (genrand_real2_dsfmt() < pDoubles) then 
            if (genrand_real2_dsfmt() < pParallel) then 
                ! do a parallel triple excitation, coming from the triples..
                call gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
                ic = 2
                pgen = pgen * pDoubles * pParallel
!                 if (nJ(1) /= 0) then 
!                     print *, "parallel: ", nJ
!                 end if
            else 
                ! do a "normal" hubbard k-space excitation 
                call gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

                pgen = pgen * pDoubles * (1.0_dp - pParallel)
            end if 
        else 
            ! otherwise to a triple.. 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, temp_ex, tParity, pgen) 
            ic = 3 
            pgen = pgen * (1.0_dp - pDoubles)

        end if

    end subroutine gen_excit_k_space_hub_transcorr_test


    subroutine gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen) 
        ! i think i should calculat the matrix element in here already! 
        ! in this case.. otherwise i have to carry the tParity and ex 
        ! all the way through the rest of the code and this makes problems 
        ! i guess, since triples are actually never considered .. todo
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,3) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_triple_hubbard" 
#endif
        integer :: elecs(3), orbs(3), src(3), sum_ms
        real(dp) :: p_elec, p_orb(2)

        ! first we pick 3 electrons in this case ofc.
        ! with the restriction, that they must not be all the same spin! 
        call pick_three_opp_elecs(nI, elecs, p_elec, sum_ms)

        src = nI(elecs) 
        ! then i pick 1 orbital? maybe.. 
        ! if this orbital is of the minority spin type, the other 2 orbitals 
        ! MUST be of parallel spin.. 
        ! if the orbital is of the majority spin type the remaining 
        ! orbitals must be of opposite spin type.. 
        ! can i take that into account with a tailored get_orb_from_kpoints() 
        ! function for 3 electrons or should i decide here if we always 
        ! want to do a specific picking order (restricting pgens, but making 
        ! it easier to handle algorithmically) or if we want full flexibility 
        ! (increasing pgens, but kind of making it a bit more difficult..) 
        ! NO: we decide to always pick the minority spin first! 
        call pick_a_orbital_hubbard(ilutI, orbs(1), p_orb(1), sum_ms) 

        ! and pick the remaining two orbitals (essentially it is only 
        ! one, since the last is restricted due to momentum conservation!)
        call pick_bc_orbitals_hubbard(nI, ilutI, src, orbs(1), orbs(2:3), p_orb(2))
        
        if (orbs(2) == 0) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ! so now.. did Robert make a routine like: 
        call make_triple(nI, nJ, elecs, orbs, ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 3) 

        ! i am not sure about the factor of 2 here.. maybe this is already 
        ! taken into account in the orbital creation
        pgen = p_elec * product(p_orb) 
        
    end subroutine gen_triple_hubbard

    subroutine pick_bc_orbitals_hubbard(nI, ilutI, src, orb_a, orbs, p_orb)
        ! this is the main routine, which picks the final (b,c) orbital for 
        ! the 3-body excitation
        integer, intent(in) :: nI(nel), src(3), orb_a
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2)
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_bc_orbitals_hubbard"
        real(dp) :: test
#endif
        real(dp) :: cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2, 2), ind

        ! do it similar to the other routines.. 
        call create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            orbs(1) = ABORT_EXCITATION
            return 
        end if

        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb) 

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! the influence of orb_a is important in the pgen recalc!!
        call create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong: in ", this_routine
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ijk: ", src 
            print *, "a: ", orb_a
            print *, "orbs: ", orbs
            print *, "cum_arr: ", cum_arr
            print *, "orb_list(:,1): ", orb_list(:,1)
            print *, "orb_list(:,2): ", orb_list(:,2)
        end if
#endif

        p_orb = 2.0_dp * p_orb 

    end subroutine pick_bc_orbitals_hubbard

    subroutine create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, & 
                cum_sum, tgt, cpt) 
        integer, intent(in) :: nI(nel), src(3), orb_a
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        integer, intent(out) :: orb_list(nbasis/2, 2) 
        real(dp), intent(out) :: cum_arr(nbasis/2), cum_sum
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_bc_list_hubbard" 
#endif
        integer :: b, c, ex(2,3), spin, orb_b
        real(dp) :: elem
        integer :: nJ(nel)
        integer, allocatable :: ex2(:,:)

        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 
        ex(1,:) = src
        ex(2,1) = orb_a 

        ! we want to do a restriction! to make it easier to recalculate the 
        ! pgens and stuff! 
        ! the restriction is, that the first picked orbital (a) is of the 
        ! minority spin of the picked electrons. 
        ! so if we have to alpha and 1 beta electron picked, orbital (a) is 
        ! a beta orbital and the last 2 orbitals will be alpha and v.v.

        ! decide that the  first orbital orb_a, which is already picked is the 
        ! minority spin electron, so now pick two electrons of the opposite 
        ! spin 
        if (is_beta(orb_a)) then 
            ! then we want alpha orbitals
            spin = 0 
            ! and also be sure that we did the right thing until now, 
            ! otherwise it breaks 
            ASSERT(sum(get_spin_pn(src)) == 1)
        else 
            ! otherwise we want beta now
            spin = 1 
            ASSERT(sum(get_spin_pn(src)) == -1)
        end if

        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp 

            ! does the spin of tgt fit? 
            ! it must be opposite of orb_a!
            if (same_spin(orb_a, tgt)) return

            !TODO: we only need to consider one spin-type!!
            do b = 1, nbasis/2
                elem = 0.0_dp

                ! convert to the appropriate spin-orbitals 
                orb_b = 2 * b - spin

                if (IsNotOcc(ilutI,orb_b)) then 
                    ! get the appropriate momentum conserverd orbital c
                    c = get_orb_from_kpoints_three(src, orb_a, orb_b) 
                    
                    if (c /= orb_b .and. IsNotOcc(ilutI,c)) then 

                        ex(2,2:3) = [orb_b, c]

                        ! actually i messed up with the non-hermiticity 
                        ! i should actually switch the order of the 
                        ! determinants in matrix element calculation
                        ! old one: 
!                         elem = abs(get_3_body_helement_ks_hub(nI, ex, .false.))
                        call swap_excitations(nI, ex, nJ, ex2)
                        elem = abs(get_3_body_helement_ks_hub(nJ, ex2, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem 
                if (tgt == orb_b) then 
                    cpt = elem 
                end if
            end do

            if (cum_sum < EPS) then 
                cpt = 0.0_dp 
            else 
                cpt = cpt / cum_sum
            end if
        else 
            do b = 1, nbasis/2 
                orb_b = 2 * b - spin 

                elem = 0.0_dp 

                if (IsNotOcc(ilutI, orb_b)) then 
                    c = get_orb_from_kpoints_three(src, orb_a, orb_b) 

                    if (c /= orb_b .and. IsNotOcc(ilutI,c)) then 

                        ex(2,2:3) = [orb_b, c] 
                        call swap_excitations(nI, ex, nJ, ex2)
                        elem = abs(get_3_body_helement_ks_hub(nJ, ex2, .false.))
 
                    end if
                end if 
                cum_sum = cum_sum + elem 
                cum_arr(b) = cum_sum 
                orb_list(b,:) = [orb_b,c] 
            end do
        end if
                
    end subroutine create_bc_list_hubbard

    subroutine pick_a_orbital_hubbard(ilutI, orb, p_orb, sum_ms)
        ! hm... i think the first orbital can be picked totally random out 
        ! of the empty orbitals or? since every spin and every momentum 
        ! is allowed.. and i do not want to overdo the weighting i guess, 
        ! especially not for the beginning 
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        integer, intent(out) :: orb 
        real(dp), intent(out) :: p_orb
        integer, intent(in), optional :: sum_ms
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_a_orbital_hubbard"
#endif
        integer :: spin

        ! if sum_ms is present, we pick the first orbital from the minority 
        ! spins in the picked electrons -> so it is the opposite 
        if (present(sum_ms)) then 
            if (sum_ms == -1) then 
                ! there a are 2 beta and one alpha electron picked -> 
                ! so pick alpha here! 
                spin = 0 
                p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp) 
            else if (sum_ms == 1) then 
                spin = -1
                p_orb = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 
            end if

            do 
                orb = 2*(1 + int(genrand_real2_dsfmt() * nbasis/2)) + spin

                if (IsNotOcc(ilutI, orb)) exit 

            end do
        else 

            do 
                orb = 1 + int(genrand_real2_dsfmt() * nbasis) 

                if (IsNotOcc(ilutI, orb)) exit 

            end do

            p_orb = 1.0_dp/real(nbasis - nel, dp) 
        end if

    end subroutine pick_a_orbital_hubbard

    subroutine gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
        integer, intent(in) :: nI(nel) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_parallel_double_hubbard"
#endif
        real(dp) :: p_elec, p_orb
        integer :: elecs(2), orbs(2), src(2), ispn

        ! in the transcorrelated case we have to decide 
        ! i first have to choose an electron pair (ij) at random 
        ! but with the condition that they have to have opposite spin! 
        ! this is the only difference: i pick two spin-parallel electrons.. 
        ! the rest stays the same.. i just have to adjust the 
        ! get_orb_from_kpoints routine 
        ! and the matrix element calculation
        call pick_spin_par_elecs(nI, elecs, p_elec) 

        src = nI(elecs)

        ! i realised i could reuse the already implemented orbital picker, 
        ! but the other one does not use the fact that we know that we 
        ! have parallel spin in this case! so implement a parallel 
        ! spin-orbital picker!! (also better for matrix elements.. so we 
        ! must not check if the spins fit, if we only take the correct ones!)
        call pick_ab_orbitals_par_hubbard(nI, ilutI, src, orbs, p_orb)

        if (orbs(1) == ABORT_EXCITATION) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ! and make the excitation 
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2) 

        ! i think in both the electrons and the orbitals i have twice the 
        ! probability to pick them
        ! already modified in the orbital picker! 
        ! i am not super sure about a factor of 2 here.. 
        pgen = p_elec * p_orb 

    end subroutine gen_parallel_double_hubbard

    subroutine pick_ab_orbitals_par_hubbard(nI, ilutI, src, orbs, p_orb) 
        integer, intent(in) :: nI(nel), src(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_ab_orbitals_par_hubbard"
        real(dp) :: test, test_arr(nBasis/2)
        integer :: ex(2,2)
#endif
        real(dp) :: cum_arr(nbasis/2)
        real(dp) :: cum_sum
        integer :: orb_list(nbasis/2, 2)
        integer :: ind

        ! without transcorrelation factor this is uniform, but with a 
        ! transcorrelation factor the matrix element might change and so also 
        ! the pgen should change. 

        call create_ab_list_par_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            p_orb = 0.0_dp
            orbs(1) = ABORT_EXCITATION
            return
        end if

        ! this stuff is also written so often i should finally make a routine 
        ! out of that 
        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! check that the other way of picking the orbital has the same 
        ! probability.. 
        call create_ab_list_par_hubbard(nI, ilutI, src, orb_list, test_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong in ", this_routine
            print *, "nI: ", nI
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ij: ", src
            print *, "orbs: ", orbs
            print *, "cum_arr: ", cum_arr
            print *, "test_arr: ", test_arr
        end if

        !todo: also call the calc_pgen_k_space_hubbard here and check 
        ! pgens 
        ex(1,:) = src
        ex(2,:) = orbs 

#endif

        ! do i have to recalc. the pgen the other way around? yes! 
        ! effectively reuse the above functionality
        ! i am pretty sure i just have to find the position in the 
        ! list.. OR: since in the hubbard it is just twice the 
        ! probability or? i am pretty sure yes.. but for all of them.. 
        ! so in the end it shouldnt matter again..
        p_orb = 2.0_dp * p_orb

    end subroutine pick_ab_orbitals_par_hubbard

    subroutine pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orb) 
        ! depending on the already picked electrons (ij) pick an orbital 
        ! (a) and the connected orbital (b)
        integer, intent(in) :: nI(nel), src(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_ab_orbitals_hubbard"
        real(dp) :: test
        integer :: ex(2,2)
#endif
        real(dp) :: cum_arr(nbasis)
        real(dp) :: cum_sum
        integer :: orb_list(nbasis, 2)
        integer :: ind

        ! without transcorrelation factor this is uniform, but with a 
        ! transcorrelation factor the matrix element might change and so also 
        ! the pgen should change. 
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            orbs(1) = ABORT_EXCITATION
            return
        end if

        ! this stuff is also written so often i should finally make a routine 
        ! out of that 
        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! check that the other way of picking the orbital has the same 
        ! probability.. 
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "orbs: ", orbs
        end if

        !todo: also call the calc_pgen_k_space_hubbard here and check 
        ! pgens 
        ex(1,:) = src
        ex(2,:) = orbs 

#endif

        ! do i have to recalc. the pgen the other way around? yes! 
        ! effectively reuse the above functionality
        ! i am pretty sure i just have to find the position in the 
        ! list.. OR: since in the hubbard it is just twice the 
        ! probability or? i am pretty sure yes.. but for all of them.. 
        ! so in the end it shouldnt matter again..
        p_orb = 2.0_dp * p_orb

    end subroutine pick_ab_orbitals_hubbard

    subroutine create_ab_list_par_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            tgt, cpt) 
        integer, intent(in) :: nI(nel), src(2) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orb_list(nbasis/2, 2)
        real(dp), intent(out) :: cum_arr(nbasis/2)
        real(dp), intent(out) :: cum_sum 
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_ab_list_par_hubbard"
#endif
        integer :: a, b, ex(2,2), spin, orb_a
        real(dp) :: elem
        integer :: nJ(nel) 
        integer, allocatable :: ex2(:,:)
        ! do the cum_arr for the k-space hubbard 
        ! i think here i might really use flags.. and not just do the 
        ! influence over the matrix elements.. since without transcorrelation 
        ! i would waste alot of effort if i calculate the matrix elements 
        ! here all the time.. 
        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 

        ex(1,:) = src
        
        ! this routine only checks for parallel spins, depending on src
        ASSERT(same_spin(src(1),src(2)))

        ! make a spin factor for the orbital conversion
        ! 0...alpha
        ! 1...beta
        spin = get_spin(src(1)) - 1

        ! and only loop over the correct spin
        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp 

            ! if target does not have the same spin, do we abort or return?
            if (.not. same_spin(src(1),tgt)) return

            do a = 1, nbasis/2
                elem = 0.0_dp 

                orb_a = 2 * a - spin

                if (IsNotOcc(ilutI,orb_a)) then 
                    ! modify get_orb_from_kpoints to give spin-parallel
                    ! it already does i think! 
                    b = get_orb_from_kpoints(src(1),src(2),orb_a)

                    if (b /= orb_a .and. IsNotOcc(ilutI,b)) then 

                        ex(2,:) = [orb_a,b] 

                        call swap_excitations(nI, ex, nJ, ex2)
                        elem = abs(get_offdiag_helement_k_sp_hub(nJ, ex2, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem
                if (tgt == orb_a) then 
                    cpt = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                cpt = 0.0_dp
            else 
                cpt = cpt / cum_sum 
            end if
        else 
            do a = 1, nbasis/2
                orb_a = 2 * a - spin 

                elem = 0.0_dp 

                if (IsNotOcc(ilutI,orb_a)) then 
                    b = get_orb_from_kpoints(src(1),src(2), orb_a)

                    if (b /= orb_a .and. IsNotOcc(ilutI, b)) then 
                        ex(2,:) = [orb_a, b]
                        call swap_excitations(nI, ex, nJ, ex2)
                        elem = abs(get_offdiag_helement_k_sp_hub(nJ, ex2, .false.))
                    end if
                end if
                cum_sum = cum_sum + elem 
                cum_arr(a) = cum_sum 
                orb_list(a,:) = [orb_a,b]
            end do
        end if 
        
    end subroutine create_ab_list_par_hubbard
      
    subroutine create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            tgt, cpt) 
        integer, intent(in) :: nI(nel), src(2) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orb_list(nbasis, 2)
        real(dp), intent(out) :: cum_arr(nbasis)
        real(dp), intent(out) :: cum_sum 
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_ab_list_hubbard"
#endif
        integer :: a, b, ex(2,2)
        real(dp) :: elem
        integer :: nJ(nel)
        integer, allocatable :: ex2(:,:)

        ! do the cum_arr for the k-space hubbard 
        ! i think here i might really use flags.. and not just do the 
        ! influence over the matrix elements.. since without transcorrelation 
        ! i would waste alot of effort if i calculate the matrix elements 
        ! here all the time.. 
        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 

        ex(1,:) = src

        ! we are also using this routine for the parallel excitations due to 
        ! the transcorrelation factor.. this is nice, since it is easily 
        ! reusable, but, loses alot of efficiency, since the we are looping 
        ! over all spin-orbital, although we know we only want to loop over a 
        ! certain spin! 
        ! todo
        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp

            ! OPTIMIZATION: Do not loop over nbasis here, but over a pre-computed
            ! lookup table of excitations for src (if possible, is not an option
            ! if the matrix element depends on nI)
            do a = 1, nbasis
                elem = 0.0_dp
                ! if a is empty
                if (IsNotOcc(ilutI, a)) then 
                    ! i have to rewrite get_orb, so it gives me the same 
                    ! spin if src has the same spin! todo
                    ! to take into account spin-parallel double 
                    ! excitations!

                   ! get the excitation
                   b = get_orb_from_kpoints(src(1), src(2), a)
                   ! we have yet to check if b is unoccupied
                   if(b /= a .and. IsNotOcc(ilutI,b)) then
                      ! assert that we hit opposite spins
                      if(.not. t_trans_corr_2body) then 
                         ASSERT(.not. same_spin(a,b))
                      endif
                      ! get the matrix element (from storage)
                      if (.not. (t_trans_corr .or. t_trans_corr_2body)) then 
                          elem = excit_cache(src(1),src(2),a)
                      else 
                          ex(2,:) = [a,b]
                          call swap_excitations(nI, ex, nJ, ex2)
                          elem = abs(get_offdiag_helement_k_sp_hub(nJ, ex2, .false.))
                      end if
                   endif
                end if
                cum_sum = cum_sum + elem 

                if (tgt == a)  then 
                    cpt = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                cpt = 0.0_dp 
            else 
                ! todo: maybe i have to multiply by 2 here.. since both 
                ! direction possible.. 
                cpt = cpt / cum_sum 
            end if
        else 
            do a = 1, nbasis
                elem = 0.0_dp 
                b = -1

                if (IsNotOcc(ilutI, a)) then 
                    b = get_orb_from_kpoints(src(1), src(2), a)

                    if (b /= a .and. IsNotOcc(ilutI, b)) then 
                        ! get the matrix element (from storage)
                        if (.not. (t_trans_corr .or. t_trans_corr_2body)) then 
                            elem = excit_cache(src(1),src(2),a)
                        else 
                            ex(2,:) = [a,b]
                            call swap_excitations(nI, ex, nJ, ex2)
                            elem = abs(get_offdiag_helement_k_sp_hub(nJ, ex2, .false.))
                        end if
                    end if
                end if
                cum_sum = cum_sum + elem 
                cum_arr(a) = cum_sum 
                orb_list(a,:) = [a,b] 
            end do
        end if

    end subroutine create_ab_list_hubbard

    function calc_pgen_k_space_hubbard_uniform_transcorr(nI, ilutI, ex, ic) result(pgen)
        ! need a calc pgen functionality for the uniform transcorrelated 
        ! excitation generator
        integer, intent(in) :: nI(nel), ex(:,:), ic
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_uniform_transcorr"
#endif
        real(dp) :: p_elec, p_orb
        integer :: sum_ms

        pgen = 0.0_dp

        if (ic == 2) then 

            pgen = pDoubles

            if (same_spin(ex(1,1),ex(1,2))) then 
                pgen = pgen * pParallel

                if (is_beta(ex(1,1))) then 
                    p_elec = 1.0_dp / real(nOccBeta*(nOccBeta-1),dp)
                    p_orb = 2.0_dp / real(nbasis/2 - nOccBeta, dp)
                else 
                    p_elec = 1.0_dp / real(nOccAlpha*(nOccAlpha-1),dp)
                    p_orb = 2.0_dp / real(nbasis/2 - nOccAlpha, dp)
                end if

            else 
                pgen = pgen * (1.0_dp - pParallel)
                p_elec = 1.0_dp / real(nOccBeta * nOccAlpha,dp)

                p_orb = 2.0_dp / real(nbasis - nel, dp)

            end if

        else 
            pgen = 1.0_dp - pDoubles

            sum_ms = sum(get_spin_pn(ex(1,:)))

            ASSERT(sum_ms == 1 .or. sum_ms == -1)

            if (sum_ms == 1) then 
                p_elec = 2.0_dp/real(nel*(nel-1),dp) * &
                    (1.0_dp/real(nOccBeta,dp) + 2.0_dp/real(nel-2,dp))

                p_orb = 1.0_dp / real(nbasis/2 - nOccBeta, dp) * & 
                        2.0_dp / real(nbasis/2 - nOccAlpha, dp) 

            else if (sum_ms == -1) then 
                p_elec = 2.0_dp/real(nel*(nel-1),dp) * & 
                    (1.0_dp/real(nOccAlpha,dp) + 2.0_dp/real(nel-2,dp))

                p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp) * & 
                        2.0_dp / real(nBasis/2 - nOccBeta, dp)

            end if
        end if

        pgen = pgen * p_elec * p_orb

    end function calc_pgen_k_space_hubbard_uniform_transcorr

    function calc_pgen_k_space_hubbard_transcorr(nI, ilutI, ex, ic) result(pgen)
        ! this function i have to rewrite for the transcorrelated to take 
        ! the same-spin doubles and triples into account! 
        ! NOTE: ex could be of form (2,3) in the case of triples!
        integer, intent(in) :: nI(nel), ex(:,:), ic 
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_transcorr" 
#endif 
    
        pgen = 0.0_dp

        if (ic == 2) then 
            if (same_spin(ex(1,1),ex(1,2))) then 
                ! parallel double excitation
                ! the spins are checked within the function:
                pgen = calc_pgen_k_space_hubbard_par(nI,ilutI,ex,ic)

                pgen = pgen * pDoubles * pParallel

            else 
                ! "normal" opposite spin excitation
                ! the spins are checked within the function: 
                pgen = calc_pgen_k_space_hubbard(nI, ilutI, ex, ic) 

                pgen = pgen * pDoubles * (1.0_dp - pParallel) 
            end if
        else if (ic == 3) then 
            pgen = calc_pgen_k_space_hubbard_triples(nI, ilutI, ex, ic)

            pgen = pgen * (1.0_dp - pDoubles) 

        end if

    end function calc_pgen_k_space_hubbard_transcorr

    function calc_pgen_k_space_hubbard_triples(nI, ilutI, ex, ic) result(pgen)
        integer, intent(in) :: nI(nel), ex(:,:), ic
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_triples"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb(2), cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2,2), sum_ms, orb_a, orbs(2)

        if (ic /= 3) then 
            pgen = 0.0_dp
            return
        end if

        sum_ms = sum(get_spin_pn(ex(1,:)))

        ! check spins
        if (.not. (sum_ms == 1 .or. sum_ms == -1) .or. sum_ms /= sum(get_spin_pn(ex(2,:)))) then
            pgen = 0.0_dp
            return
        end if

        ! get the probabilites for the electrons and orbital (a)
        if (sum_ms == 1) then 
            p_elec = 1.0_dp / real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccBeta,dp) + 2.0_dp / real(nel-2,dp))

            p_orb(1) = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 

        else 
            p_elec = 1.0_dp / real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccAlpha,dp) + 2.0_dp / real(nel-2,dp))

            p_orb(1) = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)

        end if

        ! for this i need the minority spin orbital (a)
        orb_a = find_minority_spin(ex(2,:))

        orbs = pack(ex(2,:), ex(2,:) /= orb_a) 

        call create_bc_list_hubbard(nI, ilutI, ex(1,:), orb_a, orb_list, cum_arr, & 
            cum_sum, orbs(1), p_orb(2))

        pgen = p_elec * product(p_orb) * 2.0_dp 

#ifdef __DEBUG 
        call create_bc_list_hubbard(nI, ilutI, ex(1,:), orb_a, orb_list, cum_arr, & 
            cum_sum, orbs(2), test)

        if (abs(test - p_orb(2)) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb(2)
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if
#endif
        
    end function calc_pgen_k_space_hubbard_triples

    function calc_pgen_k_space_hubbard_par(nI, ilutI, ex, ic) result(pgen) 
        integer, intent(in) :: nI(nel), ex(:,:), ic
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_par"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb, cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2,2)

        ! check ic:
        if (ic /= 2) then 
            pgen = 0.0_dp
            return
        end if

        ! check spin:
        if (.not.(same_spin(ex(1,1),ex(1,2)) .and. same_spin(ex(2,1),ex(2,2)) .and. & 
            same_spin(ex(1,1),ex(2,1)))) then 
            pgen = 0.0_dp
            return
        end if

        if (get_ispn(ex(1,1:2)) == 1) then 
            p_elec = 1.0_dp / real(nbasis/2 - nOccBeta, dp)
        else 
            p_elec = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)
        end if

        call create_ab_list_par_hubbard(nI, ilutI, ex(1,1:2), orb_list, cum_arr, & 
            cum_sum, ex(2,1), p_orb)

        pgen = p_elec * p_orb * 2.0_dp

#ifdef __DEBUG 
        call create_ab_list_par_hubbard(nI, ilutI, ex(1,1:2), orb_list, cum_arr, & 
            cum_sum, ex(2,1), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if

#endif

    end function calc_pgen_k_space_hubbard_par

    function calc_pgen_k_space_hubbard(nI, ilutI, ex, ic) result(pgen) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: nI(nel), ex(2,2), ic
        real(dp) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb, cum_arr(nbasis), cum_sum
        integer :: orb_list(nbasis,2), src(2)

        if (ic /= 2) then 
            pgen = 0.0_dp 
            return 
        end if

        if (same_spin(ex(1,1),ex(1,2)) .or. same_spin(ex(2,1),ex(2,2))) then 
            pgen = 0.0_dp 
            return
        end if

        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp) 

        src = get_src(ex)

        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
                ex(2,1), p_orb) 

#ifdef __DEBUG
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
                ex(2,2), test) 

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if

#endif

        ! i do not need to recalc, the p(b|ij) since it is the same.. 
        ! but i need a factor of 2 somewhere.. figure that out!
        pgen = p_elec * p_orb * 2.0_dp
 
    end function calc_pgen_k_space_hubbard 

    subroutine init_get_helement_k_space_hub

        if (iProcIndex == root) then 
            print *, "initialize k-space get_helement pointer"
        end if

        if (t_trans_corr_2body) then 
            three_body_prefac = real(bhub,dp) * prefac_test * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
        end if

        call init_dispersion_rel_cache()
        call init_tmat_kspace(lat)
        call init_two_body_trancorr_fac_matrix()
        n_opp(-1) = real(nel/2 + lms,dp)
        n_opp(1) = real(nel/2 - lms,dp)
        call init_three_body_const_mat()

        get_umat_el => get_umat_kspace
        ! i guess i should also set the transcorr factor here or?? 
        get_helement_lattice_ex_mat => get_helement_k_space_hub_ex_mat
        get_helement_lattice_general => get_helement_k_space_hub_general
        ! maybe i have to initialize more here, especially if we are using the 
        ! HPHF keyword I guess.. 

    end subroutine init_get_helement_k_space_hub

    function get_helement_k_space_hub_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ic, ex(2,ic)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel
#ifdef __DEBUG 
        character(*), parameter :: this_routine ="get_helement_k_space_hub_ex_mat"
#endif

        !todo: if 2-body-transcorrelation, we can have triple excitations now..
        ! fix that here.. (and also in a lot of other parts in the code..)

        if (ic == 0) then 
            ! the diagonal is just the sum of the occupied one-particle 
            ! basis states 
            hel = get_diag_helement_k_sp_hub(nI)

        else if (ic == 2) then 

            hel = get_offdiag_helement_k_sp_hub(nI, ex, tpar) 

        else if (ic == 3 .and. t_trans_corr_2body) then 

            hel = get_3_body_helement_ks_hub(nI, ex, tpar)

        else 

            hel = h_cast(0.0_dp) 

        end if

    end function get_helement_k_space_hub_ex_mat


    function get_helement_k_space_hub_general(nI, nJ, ic_ret) result(hel) 
        integer, intent(in) :: nI(nel), nJ(nel) 
        integer, intent(inout), optional :: ic_ret
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_helement_k_space_hub_general"
#endif
        integer :: ic, ex(2,3), ex_2(2,2)
        logical :: tpar 
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:niftot)

        !todo: if 2-body-transcorrelation, we can have triple excitations now..
        ! fix that here.. (and also in a lot of other parts in the code..)
        if (present(ic_ret)) then 
            if (ic_ret == 0) then 
                hel = get_diag_helement_k_sp_hub(nI)

            else if (ic_ret == 2) then 
                ex_2(1,1) = 2

                call GetExcitation(nI, nJ, nel, ex_2, tpar) 
                hel = get_offdiag_helement_k_sp_hub(nI, ex_2, tpar) 

            else if (ic_ret == 3 .and. t_trans_corr_2body) then 
                ex(1,1) = 3 
                call GetExcitation(nI, nJ, nel, ex, tpar) 
                hel = get_3_body_helement_ks_hub(nI, ex, tpar)

            else if (ic_ret == -1) then 
                call EncodeBitDet(nI, ilutI) 
                call EncodeBitDet(nJ, ilutJ) 

                ic_ret = FindBitExcitLevel(ilutI, ilutJ) 

                if (ic_ret == 0) then 
                    hel = get_diag_helement_k_sp_hub(nI)

                else if (ic_ret == 2) then 
                    ex_2(1,1) = 2 
                    call GetBitExcitation(ilutI, ilutJ, ex_2, tpar) 

                    hel = get_offdiag_helement_k_sp_hub(nI, ex_2, tpar) 

                else if (ic_ret == 3 .and. t_trans_corr_2body) then 
                    ex(1,1) = 3 
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar) 

                    hel = get_3_body_helement_ks_hub(nI, ex, tpar) 

                else 
                    hel = h_cast(0.0_dp) 
                end if
            else 
                hel = h_cast(0.0_dp) 
            end if
        else 
            call EncodeBitDet(nI, ilutI) 
            call EncodeBitDet(nJ, ilutJ) 

            ic = FindBitExcitLevel(ilutI, ilutJ) 

            if (ic == 0) then 
                hel = get_diag_helement_k_sp_hub(nI)
            else if (ic == 2) then
                ex_2(1,1) = 2 
                call GetBitExcitation(ilutI, ilutJ, ex_2, tpar) 

                hel = get_offdiag_helement_k_sp_hub(nI, ex_2, tpar) 

            else if (ic == 3 .and. t_trans_corr_2body) then 
                ex(1,1) = 3 
                call GetBitExcitation(ilutI, ilutJ, ex, tpar) 

                hel = get_3_body_helement_ks_hub(nI, ex, tpar) 

            else 
                hel = h_cast(0.0_dp) 
            end if 
        end if

    end function get_helement_k_space_hub_general

    ! i have not switched to this diag routine yet?!
    function get_diag_helement_k_sp_hub(nI) result(hel) 
        integer, intent(in) :: nI(nel) 
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_diag_helement_k_sp_hub" 
#endif
        integer :: i, j, id(nel), idX, idN, spin, k, k_vec(3), p_vec(3)
        HElement_t(dp) :: hel_sing, hel_doub, hel_par, hel_opp, hel_one, hel_three
        HElement_t(dp) :: temp_hel, temp_hel2
        type(symmetry) :: p_sym, k_sym

        ! todo: in the case of 2-body-transcorrelation, there are more 
        ! contributions.. 
        hel = h_cast(0.0_dp)

        if (t_trans_corr_2body) then 
            ! this is the 
            hel_sing = sum(GetTMatEl(nI,nI))

            id = get_spatial(nI)

            hel_doub = h_cast(0.0_dp) 
            hel_one = h_cast(0.0_dp)
            hel_three = h_cast(0.0_dp)

            temp_hel = 0.0_dp

            ! redo this whole shabang.. the formulas are actually much easier: 
            ! but just to be sure for now, do i explicetly without any use of 
            ! symmetry 
            do i = 1, nel 
                do j = 1, nel 
                    ! the restriction is, that i and j must have opposite 
                    ! spin! this also excludes i == j 
                    if (.not. same_spin(nI(i),nI(j))) then 

                        idX = max(id(i),id(j))
                        idN = min(id(i),id(j))

                        ! now we need 1/2, since we loop over all electrons
                        hel_doub = hel_doub + 0.5_dp * get_umat_kspace(idN,idX,idN,idX)

                        ! then we need the factor of the one-body transcorr influence
                        ! t is defined as -t in our code!, so bhub is usually -1
                        ! and look in the formulas it is actually -2t*cos(k)*2(cosh J - 1)
                        ! (with the k-vector of orbial i!
                        hel_one = hel_one + epsilon_kvec(G1(nI(i))%Sym) & 
                                * omega * three_body_prefac

                        ! and the next part is the three-body with the direct 
                        ! and the exchange parts 
                        do k = 1, nel 
                            ! and the convention is that j and k have same spin! 
                            ! and j == k is also allowed and part of it! 
                            if (j == k) cycle
                            if (same_spin(ni(j),nI(k))) then 
                                ! the k vector is of i and i + j - k
                                ! i need the electrons here ofc.. 
                                ! even better then the correct k-vector addition 
                                ! would be to store an epsilon-k in terms of 
                                ! the symmetry symbols! 
                                ! something like this but nicer! 
                                p_sym = G1(nI(i))%sym
                                k_sym = SymTable(G1(nI(j))%sym%s, SymConjTab(G1(nI(k))%sym%s))


                                hel_three = hel_three - three_body_prefac * (& 
                                    epsilon_kvec(p_sym) -  & 
                                    (epsilon_kvec(SymTable(p_sym%s,k_sym%s))))


                            end if
                        end do
                    end if
                end do
            end do

            hel = hel_sing + hel_doub + hel_one + hel_three

        else 
            hel = sltcnd_0(nI)
        end if

    end function get_diag_helement_k_sp_hub

    real(dp) function get_j_opt(nI, corr_J)
        ! routine to evaluate Hongjuns J-optimization formulas 
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: corr_J

        integer :: i, j, a, spin_p, spin_q, b
        integer(n_int) :: ilut(0:niftot)
        type(symmetry) :: p_sym, q_sym, a_sym, b_sym, k_sym

        integer :: src(2), tgt(2), ex(2,2), nJ(nel)
        real(dp) :: sgn
        real(dp) :: two, rpa, exchange, sum_3, tmp_hel, sum_hel
        logical :: tsign

        call EncodeBitDet(nI, ilut)

        get_j_opt = 0.0_dp

        two = 0.0_dp
        rpa = 0.0_dp
        exchange = 0.0_dp
        sum_3 = 0.0_dp
        sum_hel = 0.0_dp
        if (.not. t_symmetric) then
            do i = 1, nel! -1
                do j = 1, nel 
                    ! i only have a contribution if the spins of nI and nJ 
                    ! are not the same! 
                    if (.not. same_spin(nI(i),nI(j))) then 
                        ! and then I need to loop over the holes, but due to 
                        ! momentum conservation, only once! 
                        do a = 1, nBasis
                            ! if a is empty 
                            if (IsNotOcc(ilut,a)) then 
                                b = get_orb_from_kpoints(nI(i),nI(j),a)
                                if (IsNotOcc(ilut,b) .and. .not. same_spin(a,b)) then

                                    p_sym = G1(nI(i))%sym
                                    q_sym = G1(nI(j))%sym 
                                    spin_p = get_spin_pn(nI(i))
                                    spin_q = get_spin_pn(nI(j))

                                    ! and now i have to think how to correly 
                                    ! choose the momenta involved
                                    if (same_spin(nI(i), a)) then
                                        a_sym = G1(a)%sym
                                        b_sym = G1(b)%sym 
                                    else 
                                        a_sym = G1(b)%sym
                                        b_sym = G1(a)%sym 
                                    end if
                                    k_sym = SymTable(p_sym%s, SymConjTab(a_sym%s))

                                    ! since i loop over all possible i,j i do not need 
                                    ! the sum like below i think 
                                    get_j_opt = get_j_opt + & 
                                        two_body_contrib(corr_J, p_sym, a_sym) + & 
                                        three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) + & 
                                        three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_q)
                                    
                                    two = two + two_body_contrib(corr_J, p_sym, a_sym)
                                    rpa = rpa +  three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p)
                                    exchange = exchange + three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p)

                                    ! and i also need to think how to sum over 
                                    ! the spins correctly!
    !                                 get_j_opt = get_j_opt + & 
    !                                     (two_body_contrib(corr_J, p_sym, a_sym) + & 
    !                                      two_body_contrib(corr_J, q_sym, b_sym))/2.0_dp + &
    !                                     (three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) +  &
    !                                      three_body_rpa_contrib(corr_J, q_sym, a_sym, spin_q))/2.0_dp + &
    !                                     (three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p) + &
    !                                      three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, b_sym, spin_q))/2.0_dp
    ! 
                                     ! especially the three-body term, i am not sure 

                                end if
                            end if
                        end do
                    end if
                end do
            end do
        else 
            do i = 1, nel! - 1
                do j = 1, nel 
                    ! i only have a contribution if the spins of nI and nJ 
                    ! are not the same! 
                    if (.not. same_spin(nI(i),nI(j))) then 
                        ! and then I need to loop over the holes, but due to 
                        ! momentum conservation, only once! 
                        do a = 1, nBasis
                            ! if a is empty 
                            if (IsNotOcc(ilut,a)) then 
                                b = get_orb_from_kpoints(nI(i),nI(j),a)
                                if (IsNotOcc(ilut,b) .and. .not. same_spin(a,b))then! .and. a < b) then

!                                     print *, "---------"
!                                     print *, "i,j, a,b:", nI(i),nI(j),a,b

                                    src = [min(nI(i),nI(j)),max(nI(i),nI(j))]
                                    tgt = [min(a,b),max(a,b)]
! 
!                                     ex(1,:) = [i,j]
!                                     ex(2,:) = tgt
!                                     nJ = nI
! 
!                                     call findexcitdet(ex,nJ,2,tsign)
! 
!                                     tmp_hel = get_helement_k_space_hub_general(nI,nJ)
!                                     print *, "hel: ", tmp_hel
                                    
                                    if (is_beta(src(1))) then 
                                        p_sym = G1(src(1))%sym
                                        q_sym = G1(src(2))%sym 
                                    else
                                        p_sym = G1(src(2))%sym
                                        q_sym = G1(src(1))%sym 
                                    end if

                                    spin_p = get_spin_pn(src(1))
                                    spin_q = get_spin_pn(src(2))

                                    if (same_spin(src(1),tgt(1))) then 
                                        a_sym = G1(tgt(1))%sym
                                        b_sym = G1(tgt(2))%sym 
                                        sgn = 1.0_dp
                                    else 
                                        a_sym = G1(tgt(2))%sym
                                        b_sym = G1(tgt(1))%sym 
                                        sgn = -1.0_dp
                                    end if

!                                     print *, "j-sign:", sgn
!                                     if (tsign) sgn = -sgn
!                                     print *, "j-tpar:", tsign

!                                     print *, "tmp_hel*sgn: ", tmp_hel*sgn
!                                     sum_hel = sum_hel + tmp_hel*sgn
                                    k_sym = SymTable(p_sym%s, SymConjTab(a_sym%s))

                                    ! since i loop over all possible i,j i do not need 
                                    ! the sum like below i think 
! 				                    get_j_opt = get_j_opt + sgn*( & 
!                                         two_body_contrib(corr_J, p_sym, a_sym) + & 
!                                         three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) + & 
!                                         three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p))

				                    get_j_opt = get_j_opt + ( & 
                                        two_body_contrib(corr_J, p_sym, a_sym) + & 
                                        two_body_contrib(corr_J, q_sym, b_sym) + & 
                                        three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) + & 
                                        three_body_rpa_contrib(corr_J, q_sym, b_sym, spin_q) + & 
                                        three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p)+&
                                        three_body_exchange_contrib(nI, corr_J, q_sym, p_sym, b_sym, spin_q))
! 
!                                     two = two_body_contrib(corr_J, p_sym, a_sym) + &
!                                                 two_body_contrib(corr_J, q_sym, b_sym)
!                                     rpa = three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) + &
!                                                  three_body_rpa_contrib(corr_J, q_sym,b_sym,spin_q)
!                                     exchange = three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p) + &
!                                         three_body_exchange_contrib(nI, corr_J, q_sym,p_sym,b_sym,spin_q)
! 
! !                                     print *, "j_opt: ", sgn*(two+rpa+exchange)/real(omega,dp)
! 
!                                     sum_3 = sum_3 + rpa + exchange

                                    ! and i also need to think how to sum over 
                                    ! the spins correctly!
    !                                 get_j_opt = get_j_opt + & 
    !                                     (two_body_contrib(corr_J, p_sym, a_sym) + & 
    !                                      two_body_contrib(corr_J, q_sym, b_sym))/2.0_dp + &
    !                                     (three_body_rpa_contrib(corr_J, p_sym, a_sym, spin_p) +  &
    !                                      three_body_rpa_contrib(corr_J, q_sym, a_sym, spin_q))/2.0_dp + &
    !                                     (three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, a_sym, spin_p) + &
    !                                      three_body_exchange_contrib(nI, corr_J, p_sym, q_sym, b_sym, spin_q))/2.0_dp
    ! 
                                     ! especially the three-body term, i am not sure 

                                end if
                            end if
                        end do
                    end if
                end do
            end do
        end if
! 
!         print *, "-----"
! 
!         print *, "2body: ", two
!         print *, "rpa: ", rpa
!         print *, "exchange:", exchange
!         print *, "sum_3: ", sum_3
        get_j_opt = get_j_opt/real(omega**2,dp)
!         print *, "sum_j_opt: ", get_j_opt
!         print *, "sum_hel: ", sum_hel


    end function get_j_opt

    function get_one_body_diag_sym(nI, spin, k_sym,t_sign) result(hel)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: spin 
        type(symmetry), intent(in) :: k_sym
        logical, intent(in), optional :: t_sign
        HElement_t(dp) :: hel 

        integer :: i, sgn, k(3)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_one_body_diag_sym"
#endif 
        ! change this routine to also use just the symmetry symbols
        type(symmetry) :: sym

        ! the spin input: -1 is beta, +1 is alpha, 0 is both!
        ! if spin is not present, default is both!
        hel = h_cast(0.0_dp)

        ! k_sym is actually always present.. 
!         if (present(k_sym)) then 
        ! work on the newest, hopefully correct way to do this.. 
        ! i need -s k vector for the triples contribution to the doubles.. 
        if (present(t_sign) .and. t_sign) then 
            sgn = -1 
        else 
            sgn = 1 
        end if

        if (sgn == 1) then 
            ASSERT(spin == -1 .or. spin == 1) 
            if (spin == -1) then
                do i = 1, nel
                    if (is_beta(nI(i))) then 
                        sym = SymTable(G1(nI(i))%sym%s, k_sym%s)
                        hel = hel + epsilon_kvec(sym)
                    end if 
                end do
            else if (spin == 1) then 
                do i = 1, nel
                    if (is_alpha(nI(i))) then 
                        sym = SymTable(G1(nI(i))%sym%s, k_sym%s)
                        hel = hel + epsilon_kvec(sym)
                    end if
                end do
            end if
        else if (sgn == -1) then 
            ASSERT(spin == -1 .or. spin == 1) 
            if (spin == -1) then
                do i = 1, nel
                    if (is_beta(nI(i))) then 
                        sym = SymTable(k_sym%s, SymConjTab(G1(nI(i))%sym%s))
                        hel = hel + epsilon_kvec(sym)
                    end if 
                end do
            else if (spin == 1) then 
                do i = 1, nel
                    if (is_alpha(nI(i))) then 
                        sym = SymTable(k_sym%s, SymConjTab(G1(nI(i))%sym%s))
                        hel = hel + epsilon_kvec(sym)
                    end if
                end do
            end if
        end if

     end function get_one_body_diag_sym

    function get_one_body_diag_kvec(nI, spin, k_shift,t_sign) result(hel)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: spin, k_shift(3)
        logical, intent(in), optional :: t_sign
        HElement_t(dp) :: hel 

        integer :: i, sgn, k(3)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_one_body_diag_kvec"
#endif 
        ! change this routine to also use just the symmetry symbols
        integer :: sym_shift
        type(symmetry) :: sym

        ! the spin input: -1 is beta, +1 is alpha, 0 is both!
        ! if spin is not present, default is both!
        hel = h_cast(0.0_dp)

        ! k_shift is actually always present.. 
!         if (present(k_shift)) then 
        ! work on the newest, hopefully correct way to do this.. 
        ! i need -s k vector for the triples contribution to the doubles.. 
        if (present(t_sign) .and. t_sign) then 
            sgn = -1 
        else 
            sgn = 1 
        end if

        sym_shift = lat%k_to_sym(k_shift(1),k_shift(2),k_shift(3))

        if (sgn == 1) then 
            ASSERT(spin == -1 .or. spin == 1) 
            if (spin == -1) then
                do i = 1, nel
                    if (is_beta(nI(i))) then 
!                                 k = lat%add_k_vec(G1(nI(i))%k, k_shift)
!                                 hel = hel + epsilon_kvec(k)
                        sym = SymTable(G1(nI(i))%sym%s, sym_shift)
                        hel = hel + epsilon_kvec(sym)
                    end if 
                end do
            else if (spin == 1) then 
                do i = 1, nel
                    if (is_alpha(nI(i))) then 
!                                 k = lat%add_k_vec(G1(nI(i))%k, k_shift)
!                                 hel = hel + epsilon_kvec(k)
                        sym = SymTable(G1(nI(i))%sym%s, sym_shift)
                        hel = hel + epsilon_kvec(sym)
                    end if
                end do
            end if
        else if (sgn == -1) then 
            ASSERT(spin == -1 .or. spin == 1) 
            if (spin == -1) then
                do i = 1, nel
                    if (is_beta(nI(i))) then 
!                                 k = lat%subtract_k_vec(k_shift, G1(nI(i))%k)
!                                 hel = hel + epsilon_kvec(k)
                        sym = SymTable(sym_shift, SymConjTab(G1(nI(i))%sym%s))
                        hel = hel + epsilon_kvec(sym)
                    end if 
                end do
            else if (spin == 1) then 
                do i = 1, nel
                    if (is_alpha(nI(i))) then 
!                                 k = lat%subtract_k_vec(k_shift, G1(nI(i))%k)
!                                 hel = hel + epsilon_kvec(k)
                        sym = SymTable(sym_shift, SymConjTab(G1(nI(i))%sym%s))
                        hel = hel + epsilon_kvec(sym)
                    end if
                end do
            end if
        end if

    end function get_one_body_diag_kvec

    function get_offdiag_helement_k_sp_hub(nI, ex, tpar) result(hel) 
        ! this routine is called for the double excitations in the 
        ! k-space hubbard. in case of transcorrelation, this can also be 
        ! spin-parallel excitations now. the triple excitation have a 
        ! seperate routine!
        ! The result does not depend on nI! 
        integer, intent(in) :: nI(nel), ex(2,2)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_offdiag_helement_k_sp_hub"
#endif
        integer :: src(2), tgt(2), ij(2), ab(2), k_vec_a(3), spin, k_vec_b(3)
        integer :: k_vec_c(3), k_vec_d(3)
        type(symmetry) :: k_sym_a, k_sym_b, k_sym_c, k_sym_d
        HElement_t(dp) :: temp
        real(dp) :: sgn

        src = get_src(ex)
        tgt = get_tgt(ex)

        if (.not. t_trans_corr_2body) then 
            if (same_spin(src(1),src(2)) .or. same_spin(tgt(1),tgt(2))) then 
                hel = h_cast(0.0_dp)
                return 
            end if
        else 
            ! if src has same spin but tgt has opposite spin -> 0 mat ele
            if (same_spin(src(1),src(2)) .and. (.not. same_spin(tgt(1),tgt(2)) &
                .or. .not. same_spin(src(1),tgt(1)))) then
                
                hel = h_cast(0.0_dp)
                return
            end if
        end if

        ij = get_spatial(src)
        ab = get_spatial(tgt)
!         ij = gtid(src)
!         ab = gtid(tgt) 
        ! that about the spin?? must spin(a) be equal spin(i) and same for 
        ! b and j? does this have an effect on the sign of the matrix element? 

        ! in the case of 2-body transcorrelation, parallel spin double exciattions 
        ! are possible todo: check if we get the coulomb and exchange contributions
        ! correct..

        ! the U part is still just the the spin-opposite part
        ! damn... i need a sign convention here too..
        if (same_spin(src(1),tgt(1)) .and. same_spin(src(2),tgt(2))) then 
            hel = get_umat_kspace(ij(1),ij(2),ab(1),ab(2))
        else if (same_spin(src(1),tgt(2)) .and. same_spin(src(2),tgt(1))) then 
            hel = -get_umat_kspace(ij(1),ij(2),ab(1),ab(2))
        end if

        ! if hel == 0, due to momentum conservation violation we can already 
        ! exit here, since this means this excitation is just no possible! 
        ! is hel only 0 due to momentum conservation? 
        if (abs(hel) < EPS) return

        if (t_trans_corr) then 
            ! do something 
            ! here the one-body term with out (-t) is necessary
!             hel = hel * exp(trans_corr_param/2.0_dp * & 
!                 (epsilon_kvec(G1(src(1))%k) + epsilon_kvec(G1(src(2))%k) & 
!                 - epsilon_kvec(G1(tgt(1))%k) - epsilon_kvec(G1(tgt(2))%k)))

            ! optimized version: 
            hel = hel * exp(trans_corr_param/2.0_dp * & 
                (epsilon_kvec(G1(src(1))%Sym) + epsilon_kvec(G1(src(2))%Sym) & 
                - epsilon_kvec(G1(tgt(1))%Sym) - epsilon_kvec(G1(tgt(2))%Sym)))
        end if

        if (t_trans_corr_2body) then 
            ! i need the k-vector of the transferred momentum.. 
            ! i am not sure if the orbitals involved in ex() are every 
            ! re-shuffled.. if yes, it is not so easy in the spin-parallel 
            ! case to reobtain the transferred momentum. although it must be 
            ! possible. for now just assume (ex(2,2)) is the final orbital b 
            ! with momentum k_i + k_j - k_a and we need the 
            ! k_j - k_a momentum
            
            if (same_spin(src(1),src(2))) then
                spin = get_spin_pn(src(1))
                ! we need the spin of the excitation here if it is parallel

                ! in the same-spin case, this is the only contribution to the 
                ! matrix element
                ! and maybe i have to take the sign additionally into 
                ! account here?? or is this taken care of with tpar??

                ! thanks to Manu i have figured it out. we have to take 
                ! the momentum between the to equally possible excitations: 
                ! c^+_b c^+_a c_q c_p with W(q-a) 
                ! and 
                !-c^+_b c^+_a c_p c_q with W(p-a) 
                ! with one of the orbital spins. 
                ! i think it doesnt matter, which one. 
                ! although for the sign it maybe does.. check thate
!                 hel = same_spin_transcorr_factor(nI, G1(ex(1,1))%k - G1(ex(2,1))%k, spin) &
!                     - same_spin_transcorr_factor(nI, G1(ex(1,2))%k - G1(ex(2,1))%k, spin)
                ! TODO: i am not sure about the sign here... 
                ! with a + i get nice symmetric results.. but i am really
                ! not sure damn.. ask ALI!
                ! i have to define an order of the input! 
                ! maybe only look at i < j and a < b, as in the rest of the 
                ! code! and then take the symmetrized matrix element 

                src = [minval(src),maxval(src)]
                tgt = [minval(tgt),maxval(tgt)]

!                 k_vec_a = lat%subtract_k_vec(G1(src(1))%k, G1(tgt(1))%k)
!                 k_vec_b = lat%subtract_k_vec(G1(src(1))%k, G1(tgt(2))%k)

                k_sym_a = SymTable(G1(src(1))%sym%s, SymConjTab(G1(tgt(1))%sym%s))
                k_sym_b = SymTable(G1(src(1))%sym%s, SymConjTab(G1(tgt(2))%sym%s))
                ! old:
!                 k_vec_a = G1(src(1))%k - G1(tgt(1))%k 
!                 k_vec_b = (G1(src(1))%k - G1(tgt(2))%k)
! 
!                 call mompbcsym(k_vec_a, nBasisMax)
!                 call mompbcsym(k_vec_b, nBasisMax)
!                 call mompbcsym(k_vec_c, nBasisMax)
!                 call mompbcsym(k_vec_d, nBasisMax)
!                 print *, "ka: ", k_vec_a(1)
!                 print *, "kb: ", k_vec_b(1)
!                 print *, "kc: ", k_vec_c(1)
!                 print *, "kd: ", k_vec_d(1)

                ! fuck.. i am really not sure how to deal with that.. 
                ! yes this is it below! i just have to be sure that src and 
                ! tgt are ordered.. we need a convention for these matrix 
                ! elements!
                spin = get_spin_pn(src(1))
!                 hel = (same_spin_transcorr_factor(nI, k_vec_a, spin) & 
!                     - same_spin_transcorr_factor(nI, k_vec_b, spin))

                hel = (same_spin_transcorr_factor(nI, k_sym_a, spin) & 
                    - same_spin_transcorr_factor(nI, k_sym_b, spin))! &

            else 
                ! else we need the opposite spin contribution
                
                ! the two-body contribution needs two k-vector inputs. 
                ! figure out what momentum is necessary there! 
                ! i need the transferred momentum 
                ! and the momentum of other involved electron 
                ! which by definition of k, is always ex(1,2) todo: 
                ! check if this works as intented
                ! TODO no it is not!! I have to get the signed contribution 
                ! here correct.. order in EX is not ensured!
                ! see above for same-spin excitations! 
                ! what is k-vec now?? 
                ! this seems to have the correct symmetry.. 
                ! todo.. still the check if i need 1/2 factor or smth..
                ! and not sure about the sign between those two.. 
                ! here i am still not sure why i need two factors.. 
                ! i think i could get away with a convention, which momentum 
                ! to take depending on the spin.. or i just symmetrize it.. 
                ! which hopefully is ok.. 
                ! because if i put it like that with k and -k it apparently 
                ! cancels.. 
                ! maybe i also need a convention of an ordered input of ex.. 

                sgn = 1.0_dp
                ! also adapt this two body factor.. i hope this is correct now
                if (same_spin(src(1),tgt(1))) then 
                    ! i need the right hole-momenta
                    k_sym_c = G1(tgt(1))%sym
                    k_sym_d = G1(tgt(2))%sym
                    sgn = 1.0_dp
                else 
                    k_sym_c = G1(tgt(2))%sym
                    k_sym_d = G1(tgt(1))%sym
                    sgn = -1.0_dp
                end if

!                 print *, "sgn hel:", sgn

                hel = hel + sgn*(two_body_transcorr_factor(G1(src(1))%sym, k_sym_c) & 
                          + two_body_transcorr_factor(G1(src(2))%sym, k_sym_d))

!                 print *, "two-body hel: ", sgn*(two_body_transcorr_factor(G1(src(1))%sym, k_sym_c) & 
!                           + two_body_transcorr_factor(G1(src(2))%sym, k_sym_d))

                ! and now the 3-body contribution: 
                ! which also needs the third involved mometum, which then 
                ! again is ex(1,1)
                ! todo.. figure out spins! 
                ! also check that! which electron momentum one has to take! 
                ! maybe this cancels in the end.. who knows.. 

                ! what should i take as the spin here?? electron 1 or 2? 
                ! i have to account for the sum of both possible spin 
                ! influences!! damn.. todo! 
                ! and this then determines which momentum i have to take.. or?
                
                hel = hel + sgn*(three_body_transcorr_fac(nI, G1(src(1))%sym, & 
                                G1(src(2))%sym, k_sym_c, get_spin_pn(src(1))) & 
                          + three_body_transcorr_fac(nI, G1(src(2))%sym, & 
                                G1(src(1))%sym, k_sym_d, get_spin_pn(src(2))))

!                 print *, "3-body hel:", sgn*(three_body_transcorr_fac(nI, G1(src(1))%sym, & 
!                                 G1(src(2))%sym, k_sym_c, get_spin_pn(src(1))) & 
!                           + three_body_transcorr_fac(nI, G1(src(2))%sym, & 
!                                 G1(src(1))%sym, k_sym_d, get_spin_pn(src(2))))
            end if
        end if

        if (tpar) hel = -hel 

!         if (abs(hel - 1.0_dp/3.0_dp) > 1e-8) print *, "hel: ", hel

    end function get_offdiag_helement_k_sp_hub

    subroutine get_transferred_momenta(ex, k_vec_a, k_vec_b)
        ! routine to reobtain transferred momentum from a given excitation 
        ! for spin-opposite double excitations i am pretty sure how, but 
        ! for triple excitations and spin-parallel doubles not so much.. todo
        integer, intent(in) :: ex(:,:) 
        integer, intent(out) :: k_vec_a(3), k_vec_b(3)
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_transferred_momenta"
#endif
        integer :: sort_ex(2,size(ex,2))
        
        ! just to be sure, sort ex again.. 
        sort_ex(1,:) = [minval(ex(1,:)),maxval(ex(1,:))]
        sort_ex(2,:) = [minval(ex(2,:)),maxval(ex(2,:))]

        ASSERT(size(ex,1) == 2)
        ASSERT(size(ex,2) == 2 .or. size(ex,2) == 3)

        if (size(sort_ex,2) == 2) then 
            ! double excitation
            if (same_spin(sort_ex(1,1),sort_ex(1,2))) then 
                ! spin-parallel excitation
                ASSERT(same_spin(sort_ex(2,1),sort_ex(2,2)))
                ASSERT(same_spin(sort_ex(1,1),sort_ex(2,1)))

                ! for now just take the momentum of ex(1,2) - ex(2,1) 
                ! and ex(2,1) - ex(1,1)
!                 k_vec_a = G1(sort_ex(1,1))%k - G1(sort_ex(2,1))%k 
!                 k_vec_b = G1(sort_ex(1,2))%k - G1(sort_ex(2,1))%k
                k_vec_a = lat%subtract_k_vec(G1(sort_ex(1,1))%k, G1(sort_ex(2,1))%k)
                k_vec_b = lat%subtract_k_vec(G1(sort_ex(1,2))%k, G1(sort_ex(2,1))%k)


                if (t_k_space_hubbard) then
!                     k_vec_a = lat%map_k_vec(k_vec_a)
!                     k_vec_b = lat%map_k_vec(k_vec_b)
                else 
                    call mompbcsym(k_vec_a, nBasisMax)
                    call mompbcsym(k_vec_b, nBasisMax)
                end if

            else 
                ! "normal" hubbard spin-opposite excitation
                ASSERT(.not. same_spin(ex(2,1),ex(2,2)))
                ! here it is easier, we need the momentum difference of the 
                ! same spin-electrons 
                ! the sign of k should be irrelevant or? todo!
                if (same_spin(ex(1,1),ex(2,1))) then 

!                     k_vec_a = G1(ex(1,1))%k - G1(ex(2,1))%k 
!                     k_vec_b = G1(ex(1,2))%k - G1(ex(2,2))%k 
                    k_vec_a = lat%subtract_k_vec(G1(ex(1,1))%k, G1(ex(2,1))%k)
                    k_vec_a = lat%subtract_k_vec(G1(ex(1,2))%k, G1(ex(2,2))%k)
                    
                else 
!                     k_vec_a = G1(ex(1,1))%k - G1(ex(2,2))%k 
!                     k_vec_b = G1(ex(1,2))%k - G1(ex(2,1))%k 
                    k_vec_a = lat%subtract_k_vec(G1(ex(1,1))%k, G1(ex(2,2))%k)
                    k_vec_a = lat%subtract_k_vec(G1(ex(1,2))%k, G1(ex(2,1))%k)

                end if

                if (t_k_space_hubbard) then 
!                     k_vec_a = lat%map_k_vec(k_vec_a)
!                     k_vec_b = lat%map_k_vec(k_vec_b)
                else 
                    call mompbcsym(k_vec_a, nBasisMax)
                    call mompbcsym(k_vec_b, nBasisMax)
                end if
            end if
        else 
            ! triple excitations..
            ! i think i do not really need the triples.. 
            ASSERT(.false.)
        end if

    end subroutine get_transferred_momenta

    subroutine setup_k_total(nI) 
        integer, intent(in), optional :: nI(nel) 
        character(*), parameter :: this_routine = "setup_k_total"

        integer :: i

        if (present(nI)) then 

            ktotal = 0 

            do i = 1, nel 
                kTotal = lat%add_k_vec(kTotal, G1(nI(i))%k)
!                 kTotal = kTotal + G1(nI(i))%k 
!                 kTotal = lat%map_k_vec(kTotal)
            end do

            if (t_k_space_hubbard) then 
!                 ktotal = lat%map_k_vec(kTotal)
            else 
                call MomPbcSym(kTotal, nBasisMax)
            end if

        else 
            ! do i based on the HF det! 
            call stop_all(this_routine, "not yet implemented")
        end if

    end subroutine setup_k_total


    ! finally write the functions to setup up the pesky G1 and nBasisMax 
    ! quantities to be consistent with the rest of the old code 
    subroutine setup_k_space_hub_sym(in_lat) 
        class(lattice), intent(in), optional :: in_lat 
        character(*), parameter :: this_routine = "setup_k_space_hub_sym"

        INTEGER :: i,j,SymInd,Lab, spin, sym0
        INTEGER , ALLOCATABLE :: Temp(:)
        ! These are for the hubbard and UEG model look-up table
        type(Symmetry) :: SymProduct, SymI, SymJ

        if (present(in_lat)) then 
            if (.not. associated(G1)) then 
                call setup_g1(in_lat)
            end if
            if (all(nBasisMax == 0)) then 
                call setup_nbasismax(in_lat)
            end if
            
            ! although this is already setup: 
!             call GenHubMomIrrepsSymTable(G1, in_lat%get_nsites()*2, nBasisMax)

            ! do only the necessary setup here! 
            ! this is essentially from symrandexcit2.F90 SpinOrbSymSetup()

            ElecPairs=(NEl*(NEl-1))/2
            MaxABPairs=(nBasis*(nBasis-1)/2)
            
            ScratchSize = 2 * nSymLabels

            if(allocated(SpinOrbSymLabel)) deallocate(SpinOrbSymLabel)

            Allocate(SpinOrbSymLabel(nBasis))
            do i=1,nBasis
            !This ensures that the symmetry labels go from 0 -> nSymLabels-1
                SpinOrbSymLabel(i)=SymClasses(((i+1)/2))-1        
            end do
#ifdef __DEBUG
            WRITE(6,*) "SpinOrbSymLabel: "
            do i=1,nBasis
                WRITE(6,*) i,SpinOrbSymLabel(i)
            enddo
#endif
            if (allocated(SymTableLabels)) deallocate(SymTableLabels)

            Allocate(SymTableLabels(0:nSymLabels-1,0:nSymLabels-1))
            SymTableLabels(:,:)=-9000    !To make it easier to track bugs
            do i=0,nSymLabels-1
                do j=0,i
                    SymI=SymLabels(i+1)        !Convert to the other symlabel convention to use SymLabels - 
                                            !TODO: I will fix this to make them consistent when working (ghb24)!
                    SymJ=SymLabels(j+1)
                    SymProduct=SymProd(SymI,SymJ)
                    !Now, we need to find the label according to this symmetry!
                    !Run through all symmetries to make working (this could be far more efficient, but its only once, so sod it...
                    do Lab=1,nSymLabels
                        if(SymLabels(Lab)%S.eq.SymProduct%S) then
                            EXIT
                        endif
                    enddo
                    if(Lab.eq.nSymLabels+1) then
                        call stop_all("SpinOrbSymSetup","Cannot find symmetry label")
                    endif
                    SymTableLabels(i,j)=Lab-1
                    SymTableLabels(j,i)=Lab-1
                enddo
            enddo
#ifdef __DEBUG
            WRITE(6,*) "SymTable:"
            do i=0,nSymLabels-1
                do j=0,nSymLabels-1
                    WRITE(6,"(I6)",advance='no') SymTableLabels(i,j)
                enddo
                WRITE(6,*) ""
            enddo
#endif        

            !SymInvLabel takes the label (0 -> nSymLabels-1) of a spin orbital, and returns the inverse symmetry label, suitable for
            !use in ClassCountInd.
            if(allocated(SymInvLabel)) deallocate(SymInvLabel)
            Allocate(SymInvLabel(0:nSymLabels-1))
            SymInvLabel=-999

            ! Dongxia changes the gamma point away from center.
            ! SDS: Provide a default sym0 for cases where this doesn't apply
            sym0 = 0
            do i = 1, nsymlabels
                if (symlabels(i)%s == 0) sym0 = i - 1
            end do

            do i = 0, nSymLabels - 1
                ! Change the sym label back to the representation used by the rest
                ! of the code, use SymConjTab, then change back to other rep of
                ! labels SymConjTab only works when all irreps are self-inverse.
                ! Therefore, instead, we will calculate the inverses by just
                ! finding the symmetry which will give A1.
                do j = 0, nSymLabels - 1
                    ! Run through all labels to find what gives totally symmetric
                    ! rep
                    if(SymTableLabels(i,j) == sym0) then
                        if(SymInvLabel(i) /= -999) then
                            write(6,*) "SymLabel: ", i
                            call stop_all(this_routine, &
                                           "Multiple inverse irreps found - error")
                        endif
                        ! This is the inverse
                        SymInvLabel(i) = j
                    endif
                enddo
                if (SymInvLabel(i) == -999) then
                    write(6,*) "SymLabel: ", i
                    call stop_all(this_routine,"No inverse symmetry found - error")
                endif
            enddo
#ifdef __DEBUG
            write(6,*) "SymInvLabel: "
            do i = 0, nSymLabels - 1
                write(6,*) i, SymInvLabel(i)
            enddo
#endif

        !SymLabelList2 and SymLabelCounts2 are now organised differently, so that it is more efficient, and easier to add new symmetries.
        !SymLabelCounts is of size (2,ScratchSize), where 1,x gives the index in SymlabelList2 where the orbitals of symmetry x start.
        !SymLabelCounts(2,x) tells us the number of orbitals of spin & symmetry x there are.

        !Therefore, if you want to run over all orbitals of a specific symmetry, you want to run over 
        !SymLabelList from SymLabelCounts(1,sym) to SymLabelCounts(1,sym)+SymLabelCounts(2,sym)-1

            if (allocated(SymLabelList2)) deallocate(SymLabelList2)
            if (allocated(SymLabelCounts2)) deallocate(SymLabelCounts2)

            Allocate(SymLabelList2(nBasis))
            Allocate(SymLabelCounts2(2,ScratchSize))
            SymLabelList2(:)=0          !Indices:   spin-orbital number
            SymLabelCounts2(:,:)=0      !Indices:   index/Number , symmetry(inc. spin)
            Allocate(Temp(ScratchSize))
            
            do j=1,nBasis
                IF(G1(j)%Ms.eq.1) THEN
                    Spin=1
                ELSE
                    Spin=2
                ENDIF
        !        WRITE(6,*) "BASIS FN ",j,G1(j)%Sym,SymClasses((j+1)/2)
                SymInd=ClassCountInd(Spin,SpinOrbSymLabel(j),G1(j)%Ml)
                SymLabelCounts2(2,SymInd)=SymLabelCounts2(2,SymInd)+1
            enddo
            SymLabelCounts2(1,1)=1
            do j=2,ScratchSize
                SymLabelCounts2(1,j)=SymLabelCounts2(1,j-1)+SymLabelCounts2(2,j-1)
            enddo
            Temp(:)=SymLabelCounts2(1,:)
            do j=1,nBasis
                IF(G1(j)%Ms.eq.1) THEN
                    Spin=1
                ELSE
                    Spin=2
                ENDIF
                SymInd=ClassCountInd(Spin,SpinOrbSymLabel(j),G1(j)%Ml)
                SymLabelList2(Temp(SymInd))=j
                Temp(SymInd)=Temp(SymInd)+1
            enddo

        !    write(6,*) "SymLabelCounts2: ",SymLabelCounts2(1,:)
        !    write(6,*) "SymLabelCounts2: ",SymLabelCounts2(2,:)
            Deallocate(Temp)

            if (allocated(OrbClassCount)) deallocate(OrbClassCount)
            ALLOCATE(OrbClassCount(ScratchSize))
            OrbClassCount(:)=0
            do i=1,nBasis
                IF(G1(i)%Ms.eq.1) THEN
    !                WRITE(6,*) "Index: ",ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml)
    !                WRITE(6,*) i,"SpinOrbSymLabel: ",SpinOrbSymLabel(i)
                    OrbClassCount(ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml))= &
                    & OrbClassCount(ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml))+1
                ELSE
    !                WRITE(6,*) "Index: ",ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml)
    !                WRITE(6,*) i,"SpinOrbSymLabel: ",SpinOrbSymLabel(i)
                    OrbClassCount(ClassCountInd(2,SpinOrbSymLabel(i),G1(i)%Ml))= &
                    & OrbClassCount(ClassCountInd(2,SpinOrbSymLabel(i),G1(i)%Ml))+1
                ENDIF
            enddo

            call setup_kPointToBasisFn(in_lat) 

            call gensymstatepairs(nbasis/2,.false.)

        else 
            call Stop_All(this_routine, "not yet implemented")
        end if

    end subroutine setup_k_space_hub_sym

    subroutine setup_g1(in_lat) 
        use SystemData, only: G1
        class(lattice), intent(in), optional :: in_lat
        character(*), parameter :: this_routine = "setup_g1"

        type(BasisFN) :: temp_g
        integer :: i,j,k,l,ind
        logical :: kallowed

        ! i think everything is in the System_neci file
        if (present(in_lat)) then 
            ! only do it if G1 has not been setup yet!
            if (.not. associated(G1)) then
                ! i need number of spin-orbitals
                allocate(G1(in_lat%get_nsites()*2))
                G1 = NullBasisFn

                ! should i rely on the already setup nBasisMax?
                if (all(nBasisMax == 0)) then 
                    call setup_nbasismax(in_lat)
                end if

                ! do a new setup of the G1.. 
                do i = 1, in_lat%get_nsites() 
                    G1(2*i-1)%k = in_lat%get_k_vec(i)
                    G1(2*i-1)%ms = -1 
                    ! can i already write the symmetry representation in here?
                    ! i guess so.. 
                    G1(2*i-1)%Sym = Symmetry(i)

                    G1(2*i)%k = in_lat%get_k_vec(i)
                    G1(2*i)%ms = 1
                    G1(2*i)%Sym = Symmetry(i) 
                end do



!                 ind = 0
!                 do i = nBasisMax(1,1), nBasisMax(1,2)
!                     do j = nBasisMax(2,1), nBasisMax(2,2)
!                         do k = nBasisMax(3,1), nBasisMax(3,2)
!                             do l  = nBasisMax(4,1), nBasisMax(4,2), 2
!                                
!                                 temp_g%k = [i,j,k]
!                                 temp_g%ms = l 
!                                 if ((treal .and. .not. ttilt) .or. kallowed(temp_g, nBasisMax)) then
!                                     ind = ind + 1 
!                                     G1(ind)%k = [i,j,k] 
!                                     G1(ind)%ms = l
!                                     G1(ind)%Sym = TotSymRep()
!                                     if (.not. in_lat%is_k_space()) then 
!                                         ! turn off- symmetry in the hubbard case
!                                         G1(ind)%sym%s = 0
!                                     end if
!                                 end if
!                             end do
!                         end do
!                     end do
!                 end do
                if (in_lat%is_k_space()) then 
                    call setup_symmetry_table
!                     call GenHubMomIrrepsSymTable(G1, in_lat%get_nsites()*2, nbasismax)

                else 
                    ! also to the rest of the symmetry stuff here: 
                    ! in case of real-space turn off symmetry completely: 
                    call GenMolpSymTable(1, G1, in_lat%get_nsites()*2)
                    ! and i have to redo the symmetry setting to 0 
                    do i = 1, in_lat%get_nsites()*2
                        G1(i)%sym%s = 0
                    end do
                end if

            end if
        else 
            ! not yet implemented!
            call Stop_All(this_routine, "not yet implemented")
        end if

    end subroutine setup_g1

    subroutine setup_nbasismax(in_lat) 
        class(lattice), intent(in), optional :: in_lat
        character(*), parameter :: this_routine =  "setup_nbasismax"

        integer :: dummy_size
        ! thats a fucking pain in the ass.. i do not want to do that now!
        if (present(in_lat)) then 
            if (all(nBasisMax == 0)) then 
                ! only do smth if nbasismax was not changed yet

                ! whatever spin-polar means: 
                if (TSPINPOLAR) then 
                    nBasisMax(4,1) = 1 
                    nBasisMax(2,3) = 1 
                else 
                    nBasisMax(4,1) = -1
                    if (nBasisMax(2,3) == 0) nBasisMax(2,3) = 2 
                end if

                if (t_k_space_hubbard) then 
                    nBasisMax(1,3) = 0
                else if (t_new_real_space_hubbard) then 
                    nBasisMax(1,3) = 4
                    nBasisMax(3,3) = 0
                    nBasisMax(2,3) = 1
                end if

                ! this is never explained: 
                nBasisMax(4,2) = 1

                return

                ! i should give lattice also a member type and a k-space flag..
                if (trim(in_lat%get_name()) == 'tilted') then 
                    ! how do i get nmaxx and the rest effectively?? 
!                     call SETBASISLIM_HUBTILT(nBasisMax, nmaxx, nmaxy, nmaxz, & 
!                         in_lat%gen_nsites()*2, in_lat%is_periodic(), itiltx, itilty))
                    ! if it is tilted the nmax stuff is usualy 1 or?? 
                    call SETBASISLIM_HUBTILT(nBasisMax, 1,1,1, dummy_size, & 
                        in_lat%is_periodic(), in_lat%get_length(1), in_lat%get_length(2))
                else 
                    call SETBASISLIM_HUB(nBasisMax, in_lat%get_length(1), & 
                        in_lat%get_length(2), in_lat%get_length(3), dummy_size, & 
                        in_lat%is_periodic(), .not. in_lat%is_k_space())
                end if
                if (thub .and. treal) then
                    ! apparently this allows integrals between different 
                    ! spins: so in the transcorrelated hubbard this should 
                    ! be changed also maybe? 
                    nBasisMax(2,3) = 1
                end if
                ASSERT(dummy_size == in_lat%get_nsites()*2)
            end if
        else 
            call Stop_All(this_routine, "not yet implemented")
        end if

    end subroutine setup_nbasismax

    subroutine setup_kPointToBasisFn(in_lat) 
        class(lattice), intent(in), optional :: in_lat 
        character(*), parameter :: this_routine = "setup_kPointToBasisFn"

        integer :: i, kmaxX, kminX, kmaxY, kminY, kminZ, kmaxZ, iSpinIndex

        if (present(in_lat)) then 
            
            if (allocated(kPointToBasisFn)) return

            if (all(nBasisMax == 0)) then 
                call setup_nbasismax(in_lat)
                call setup_g1(in_lat)
            end if

            if (.not. associated(G1)) then 
                call setup_g1(in_lat)
            end if

            kmaxX=0
            kminX=0
            kmaxY=0
            kminY=0
            ! [W.D:] can we make this the same as in the UEG: 
            kminZ = 0
            kmaxZ = 0
            do i=1, 2 * in_lat%get_nsites() 
                ! In the hubbard model with tilted lattice boundary conditions, 
                ! it's unobvious what the maximum values of
                ! kx and ky are, so this should be found
                IF(G1(i)%k(1).gt.kmaxX) kmaxX=G1(i)%k(1)
                IF(G1(i)%k(1).lt.kminX) kminX=G1(i)%k(1)
                IF(G1(i)%k(2).gt.kmaxY) kmaxY=G1(i)%k(2)
                IF(G1(i)%k(2).lt.kminY) kminY=G1(i)%k(2)
                if (G1(i)%k(3) > kmaxz) kmaxz = g1(i)%k(3)
                if (G1(i)%k(3) < kminz) kminz = g1(i)%k(3)
            enddo
            ALLOCATE(kPointToBasisFn(kminX:kmaxX,kminY:kmaxY,kminz:kmaxz,2))
            !Init to invalid
            kPointToBasisFn=-1 
            do i=1, 2*in_lat%get_nsites()
                ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
                iSpinIndex=(G1(i)%Ms+1)/2+1 
                kPointToBasisFn(G1(i)%k(1),G1(i)%k(2),G1(i)%k(3),iSpinIndex)=i
            enddo
        else 
            call Stop_All(this_routine, "not yet implemented!")
        end if

    end subroutine setup_kPointToBasisFn

    ! create the necessary routines for the triple excitation in the 
    ! 2-body transcorrelated k-space hamiltonian 
    HElement_t(dp) function same_spin_transcorr_factor_kvec(nI, k_vec, spin)
        ! this is the term coming appearing in the spin-parallel 
        ! excitations coming from the k = 0 triple excitation
        integer, intent(in) :: nI(nel), k_vec(N_DIM), spin

        same_spin_transcorr_factor_kvec = -three_body_prefac * ( & 
            get_one_body_diag(nI,-spin,k_vec) + get_one_body_diag(nI,-spin,k_vec,.true.))

    end function same_spin_transcorr_factor_kvec

    ! create the necessary routines for the triple excitation in the 
    ! 2-body transcorrelated k-space hamiltonian 
    HElement_t(dp) function same_spin_transcorr_factor_ksym(nI, k_sym, spin)
        ! this is the term coming appearing in the spin-parallel 
        ! excitations coming from the k = 0 triple excitation
        integer, intent(in) :: nI(nel), spin
        type(symmetry), intent(in) :: k_sym

        same_spin_transcorr_factor_ksym = -three_body_prefac * ( & 
            get_one_body_diag(nI,-spin,k_sym) + get_one_body_diag(nI,-spin,k_sym,.true.))

    end function same_spin_transcorr_factor_ksym

    HElement_t(dp) function rpa_contrib_kvec(J, p, k, spin) 
        ! gives the rpa contribution in the J-optimization
        real(dp), intent(in) :: J
        integer, intent(in) :: p(N_DIM), k(N_DIM), spin
#ifdef __DEBUG
        character(*), parameter :: this_routine = "rpa_contrib_kvec"
#endif
        integer :: q(N_DIM)

        q = lat%subtract_k_vec(p,k)

        ASSERT(spin == 1 .or. spin == -1)

!         if (spin == -1) then
!             n_opp_loc = real(nOccAlpha,dp)
!         else
!             n_opp_loc = real(nOccBeta,dp)
!         end if

        rpa_contrib_kvec = real(bhub,dp) * (cosh(J) - 1.0_dp) / real(omega,dp) * &
            (n_opp(spin) - 1.0_dp) * (epsilon_kvec(p) + epsilon_kvec(q))

    end function rpa_contrib_kvec

    HElement_t(dp) function rpa_contrib_ksym(J, p, a, spin)
        ! same as above just with symmetry symbols instead of vectors 
        ! BUT here i have to be careful to determine the substraction p - k
        ! already before calling this function! it is the k-symbol of the 
        ! hole correspinding to p!
        real(dp), intent(in) :: J
        type(symmetry), intent(in) :: p, a
        integer, intent(in) :: spin 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "rpa_contrib_ksym"
#endif 
        real(dp) :: n_opp_loc

        ASSERT(spin == -1 .or. spin == 1)

        if (spin == -1) then
            n_opp_loc = real(nOccAlpha,dp)
        else
            n_opp_loc = real(nOccBeta,dp)
        end if

        rpa_contrib_ksym = -2.0_dp * real(bhub,dp)*(cosh(J) - 1.0_dp) / real(omega, dp) * &
            n_opp_loc * (epsilon_kvec(p) + epsilon_kvec(a))

!         rpa_contrib_ksym = 2.0_dp * real(bhub,dp)*(cosh(J) - 1.0_dp) / real(omega, dp) * &
!             n_opp(spin) * (epsilon_kvec(p) + epsilon_kvec(a))

    end function rpa_contrib_ksym

    HElement_t(dp) function two_body_contrib_kvec(J, p, k) 
        ! two body contribution for the J optimization!
        real(dp), intent(in) :: J
        integer, intent(in) :: p(N_DIM), k(N_DIM)
        
        integer :: q(N_DIM)

        q = lat%subtract_k_vec(p,k)

        ! i still have to decide how to loop over the HF.. maybe i dont 
        ! double count and then i dont need the /2 here! 
        two_body_contrib_kvec = real(uhub, dp) / 2.0_dp - real(bhub,dp) * & 
            ((exp(J) - 1.0_dp) * epsilon_kvec(p) + (exp(-J) - 1.0) * epsilon_kvec(q))

    end function two_body_contrib_kvec

    HElement_t(dp) function two_body_contrib_ksym(J, p, a) 
        ! same as above just with symmetry symbols
        ! AND: we have  to do the p - k before calling this function! 
        real(dp), intent(in) :: J
        type(symmetry), intent(in) :: p, a

        if (.not. t_symmetric) then
            two_body_contrib_ksym = real(uhub,dp) / 2.0_dp + real(bhub,dp) * & 
                ((exp(J) - 1.0_dp) * epsilon_kvec(a) + (exp(-J) - 1.0) * epsilon_kvec(p))
        else 
            two_body_contrib_ksym = real(uhub,dp)/2.0_dp + real(bhub,dp) * & 
                ((exp(J) - 1.0_dp) * epsilon_kvec(a) + (exp(-J) - 1.0) * epsilon_kvec(p))
        end if


    end function two_body_contrib_ksym

    HElement_t(dp) function exchange_contrib_kvec(nI, J, p, q, k, spin)
        ! the 3-body exchange contribution for the J optimization
        integer, intent(in) :: nI(:), p(N_DIM), q(N_DIM), k(N_DIM), spin
        real(dp), intent(in) :: J
#ifdef __DEBUG
        character(*), parameter :: this_routine = "exchange_contrib_kvec"
#endif 
        integer :: k1(N_DIM), k2(N_DIM) 

        ASSERT(spin == -1 .or. spin == 1) 

        k1 = lat%add_k_vec(p,q)
        k2 = lat%subtract_k_vec(p,q)
        k2 = lat%subtract_k_vec(k2,k)

        exchange_contrib_kvec = -2.0_dp * real(bhub,dp)*(cosh(J)-1.0_dp)/real(omega,dp) & 
            * (get_one_body_diag(nI,-spin,k1,.true.) + get_one_body_diag(nI,-spin,k2,.true.))

    end function exchange_contrib_kvec

    HElement_t(dp) function exchange_contrib_ksym(nI, J, p, q, a, spin) 
        ! sym-symbol version of above! 
        ! BUT here: p and q are the symbols of the electrons and q is the 
        ! symbol of 1 hole! so i have to call this function also for the 
        ! exchanged version! 
        integer, intent(in) :: nI(:), spin
        real(dp), intent(in) :: J
        type(symmetry), intent(in) :: p, q, a
#ifdef __DEBUG
        character(*), parameter :: this_routine = "exchange_contrib_ksym"
#endif
        type(symmetry) :: k1, k2    

        ASSERT(spin == -1 .or. spin == 1)
        
        k1 = SymTable(p%s, q%s)
        ! the subtraction has something to do with the inputted spin!! 
        ! todo
        ! spin is chosen from p momentum! so a is the p-k = a hole 
        ! and we need the dipersion of p-k - q = a - q
        k2 = SymTable(a%s, SymConjTab(q%s))

        exchange_contrib_ksym = 2.0_dp * real(bhub,dp)*(cosh(J)-1.0_dp)/real(omega,dp) & 
            * (get_one_body_diag(nI,-spin,k1,.true.) + get_one_body_diag(nI,-spin,k2))

    end function exchange_contrib_ksym

    HElement_t(dp) function two_body_transcorr_factor_kvec(p,k) 
        integer, intent(in) :: p(N_DIM), k(N_DIM) 
        
        ! take out the part with U/2 since this is already covered in the 
        ! "normal" matrix elements
!         two_body_transcorr_factor_kvec = real(bhub,dp)/real(omega,dp)*(&
!             (exp(trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(p - k) + &
!             (exp(-trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(p))

        two_body_transcorr_factor_kvec = real(bhub,dp)/real(omega,dp) * ( & 
            (exp(trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(k) + & 
            (exp(-trans_corr_param_2body) -1.0_dp) * epsilon_kvec(p))

    end function two_body_transcorr_factor_kvec

    subroutine init_two_body_trancorr_fac_matrix() 
        integer :: i, j
        type(symmetry) :: sym_i, sym_j

        ! for more efficiency, precompute the two-body factor for all possible 
        ! symmetry symbols 
        if (allocated(two_body_transcorr_factor_matrix)) deallocate(two_body_transcorr_factor_matrix)

        allocate(two_body_transcorr_factor_matrix(nBasis/2,nBasis/2), source = 0.0_dp)

        ! loop over spatial orbitals 
        do i = 1, nBasis/2
            sym_i = G1(2*i)%sym
            do j = 1, nBasis/2
                sym_j = G1(2*j)%Sym

                two_body_transcorr_factor_matrix(sym_j%s,sym_i%s) = & 
                    real(bhub,dp)/real(omega,dp) * &
                    ((exp(trans_corr_param_2body) - 1.0_dp)*epsilon_kvec(sym_i) + & 
                     (exp(-trans_corr_param_2body) - 1.0_dp)*epsilon_kvec(sym_j))

            end do
        end do

    end subroutine init_two_body_trancorr_fac_matrix

    HElement_t(dp) function two_body_transcorr_factor_ksym(p,k) 
        type(symmetry), intent(in) :: p, k
        
        ! take out the part with U/2 since this is already covered in the 
        ! "normal" matrix elements

        ! optimize this better and precompute more stuff! 
!         two_body_transcorr_factor_ksym = real(bhub,dp)/real(omega,dp) * ( & 
!             (exp(trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(k) + & 
!             (exp(-trans_corr_param_2body) -1.0_dp) * epsilon_kvec(p))

        two_body_transcorr_factor_ksym = two_body_transcorr_factor_matrix(p%s,k%s)

    end function two_body_transcorr_factor_ksym

    HElement_t(dp) function three_body_transcorr_fac_kvec(nI, p, q, k, spin) 
        integer, intent(in) :: nI(nel), p(N_DIM), q(N_DIM), k(N_DIM), spin 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "three_body_transcorr_fac_kvec"
#endif
        integer :: k1(3), k2(3)
        real(dp) :: n_opp_loc

        ASSERT(spin == 1 .or. spin == -1)

!         if (spin == -1) then
!             n_opp_loc = real(nOccAlpha,dp)
!         else
!             n_opp_loc = real(nOccBeta,dp)
!         end if
        ! update: add k-vec to not leave first BZ 
        k1 = lat%subtract_k_vec(k,q)
        k2 = lat%add_k_vec(p,q)
        three_body_transcorr_fac_kvec = -three_body_prefac * (& 
            n_opp(spin) * (epsilon_kvec(p) + epsilon_kvec(k)) - (& 
            get_one_body_diag(nI,-spin,k1) + get_one_body_diag(nI,-spin,k2,.true.)))

    end function three_body_transcorr_fac_kvec

    HElement_t(dp) function three_body_transcorr_fac_ksym(nI, p, q, k, spin) 
        integer, intent(in) :: nI(nel), spin
        type(symmetry) :: p, q, k
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "three_body_transcorr_fac_ksym"
#endif
        type(symmetry) :: k1, k2
        real(dp) :: n_opp_loc

        ASSERT(spin == 1 .or. spin == -1)

        if (spin == -1) then
            n_opp_loc = real(nOccAlpha,dp)
        else
            n_opp_loc = real(nOccBeta,dp)
        end if

        k1 = SymTable(k%s, SymConjTab(q%s))
        k2 = SymTable(p%s, q%s)

!         three_body_transcorr_fac_ksym = -three_body_prefac * (& 
!             n_opp_loc * (epsilon_kvec(p) + epsilon_kvec(k)) - (& 
!             get_one_body_diag(nI,-spin,k1) + get_one_body_diag(nI,-spin,k2,.true.)))

!         three_body_transcorr_fac_ksym = -three_body_prefac * (& 
!             n_opp(spin) * (epsilon_kvec(p) + epsilon_kvec(k)) - (& 
!             get_one_body_diag(nI,-spin,k1) + get_one_body_diag(nI,-spin,k2,.true.)))

        three_body_transcorr_fac_ksym = three_body_const_mat(p%s,k%s,spin) + & 
            three_body_prefac * (get_one_body_diag(nI,-spin,k1) + get_one_body_diag(nI,-spin,k2,.true.))

    end function three_body_transcorr_fac_ksym

    subroutine init_three_body_const_mat()
        integer :: i, j
        type(symmetry) :: sym_i, sym_j

        if (allocated(three_body_const_mat)) deallocate(three_body_const_mat)
        allocate(three_body_const_mat(nBasis/2,nBasis/2,-1:1), source = 0.0_dp)

        do i = 1, nBasis/2
            sym_i = G1(2*i)%Sym
            do j = 1, nBasis/2
                sym_j = G1(2*j)%Sym

                three_body_const_mat(sym_i%s,sym_j%s,-1) = -three_body_prefac * &
                    n_opp(-1)*(epsilon_kvec(sym_i) + epsilon_kvec(sym_j))

                three_body_const_mat(sym_i%s,sym_j%s,1) = -three_body_prefac * &
                    n_opp(1)*(epsilon_kvec(sym_i) + epsilon_kvec(sym_j))
                
            end do
        end do

    end subroutine init_three_body_const_mat

    function get_3_body_helement_ks_hub(nI, ex, tpar) result(hel)
        ! the 3-body matrix element.. here i have to be careful about 
        ! the sign and stuff.. and also if momentum conservation is 
        ! fullfilled .. 
        integer, intent(in) :: nI(nel), ex(2,3)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_3_body_helement_ks_hub"
#endif
        integer :: ms_elec, ms_orbs, opp_elec, opp_orb, par_elecs(2), par_orbs(2)
        integer :: p_vec(3), k1(3), k2(3), k_vec(3), hole_k(3), ka(3), kb(3), kc(3), kd(3)
        logical :: sgn
        type(symmetry) :: p_sym, hole_sym, k_sym, k1_sym, k2_sym
        type(symmetry) :: ka_sym, kb_sym, kc_sym, kd_sym

        hel = h_cast(0.0_dp)

        ms_elec = sum(get_spin_pn(ex(1,:)))
        ms_orbs = sum(get_spin_pn(ex(2,:)))

        ! check spin:
        if (.not.(ms_elec == ms_orbs)) return 
        if (.not.(ms_elec == 1 .or. ms_elec == -1)) return 

        ! check momentum conservation: 
        if (.not. check_momentum_sym(ex(1,:),ex(2,:))) return

        ! i have to get the correct momenta for the epsilon contribution 
        ! see the sheets for the k-vec relations. 
        ! we need the momentum p + (s-b) or variations thereof.. 
        ! p, we know, since it is the momentum of the minority spin-electron
        ! and (s-b) is the difference of an majority spin electron and the 
        ! fitting hole, so the total momentum conservation a + b + c = s + p + q
        ! is fulfilled
        ! the same ofc is a + (c-q) 
        ! is the minority hole always in ex(2,1)? otherwise we have to find it
        ! it is not! 

        ! i think i have figured it out with the help of Manu 
        ! the k-vector of the minority spin is always involved 
        ! but of the electron.. or can we transform it? 
        ! any way we have to calculate 
        ! W(k_p + k_s - k_b) - W(k_p + k_q - k_b) 
        ! for this we have to figure out what the minority and majority 
        ! electrons are! 
        opp_elec = find_minority_spin(ex(1,:))

        ! although i really can't be sure about the minority whole always 
        ! being at the first position in ex(2,:).. 
        opp_orb = find_minority_spin(ex(2,:)) 

        p_sym = G1(opp_elec)%sym
        hole_sym = G1(opp_orb)%sym

        par_elecs = pack(ex(1,:),ex(1,:) /= opp_elec)
        par_orbs = pack(ex(2,:),ex(2,:) /= opp_orb)

        k_sym = SymTable(p_sym%s, hole_sym%s)


        ! we have to define an order here too 
        par_elecs = [minval(par_elecs), maxval(par_elecs)]
        par_orbs = [minval(par_orbs), maxval(par_orbs)] 

        ! BZ conserving addition: 
        k1_sym = SymTable(G1(par_orbs(1))%sym%s, SymConjTab(G1(par_elecs(1))%sym%s))
        k2_sym = SymTable(G1(par_orbs(1))%sym%s, SymConjTab(G1(par_elecs(2))%sym%s))

        ! need to do the correct k additions! 
        ! for some reason the compiler does not recognize the output of 
        ! add_k_vec and subtract_k_vec as a vector... 
        ! so do it intermediately 
        ka_sym = SymTable(hole_sym%s, k1_sym%s)
        kb_sym = SymTable(hole_sym%s, k2_sym%s)
        kc_sym = SymTable(k1_sym%s, SymConjTab(p_sym%s))
        kd_sym = SymTable(k2_sym%s, SymConjTab(p_sym%s))

        hel = three_body_prefac * ( & 
            epsilon_kvec(ka_sym) - epsilon_kvec(kb_sym) +  & 
            epsilon_kvec(kc_sym) - epsilon_kvec(kd_sym))

        ! i have to decide on a sign here depending on the order of the 
        ! operators.. todo! 
        sgn = get_3body_sign(ex)

        if (.not.sgn) hel = -hel

        if (tpar) hel = -hel

    end function get_3_body_helement_ks_hub

    logical function get_3body_sign(ex)
        ! i need to find some sign convention on the 3-body term, depending 
        ! on the spin of the involved orbitals 
        integer, intent(in) :: ex(2,3)

        integer :: src(3), tgt(3),i, elec_pos, orb_pos

        ! we also have to define an order of the parallel spins.. 
        ! or is this ensured? i guess it should.. 
        src = get_src(ex)
        tgt = get_tgt(ex)

        if (sum(get_spin_pn(src)) == -1) then 
            ! then alpha is the opposite spin 
            do i = 1, 3
                if (is_alpha(src(i))) elec_pos = i
                if (is_alpha(tgt(i))) orb_pos = i 
            end do

        else 
            ! otherwise beta is minority 
            do i = 1, 3 
                if (is_beta(src(i))) elec_pos = i 
                if (is_beta(tgt(i))) orb_pos = i 
            end do

        end if

        if (elec_pos == orb_pos .or. abs(elec_pos - orb_pos) == 2) then 
            get_3body_sign = .false.
        else 
            get_3body_sign = .true.
        end if

    end function get_3body_sign

    logical function check_momentum_sym(elecs, orbs) 
        ! routine to check the momentum conservation for double and triple 
        ! spawns 
        ! although this could in fact be used for a general check of 
        ! symmetry adaptability
        integer, intent(in) :: elecs(:), orbs(:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "check_momentum_sym"
#endif 
!         type(BasisFN) :: ka, kb 
        integer :: i
        type(Symmetry) :: sym_1, sym_2

        ASSERT(size(elecs) == size(orbs))

        ! make this more efficient here! and do not use all the old 
        ! functionality 

        check_momentum_sym = .true. 

        if (.not. sum(get_spin_pn(elecs)) == sum(get_spin_pn(orbs))) then 
            check_momentum_sym = .false.
            return 
        end if

        ! get the symmetry symbol
        sym_1 = G1(elecs(1))%Sym
        sym_2 = G1(orbs(1))%Sym
        do i = 2, size(elecs) 
            ! i could ofc also use G1 here again
            sym_1 = SymTable(sym_1%s, G1(elecs(i))%sym%s)
            sym_2 = SymTable(sym_2%s, G1(orbs(i))%sym%s)
        end do

        if (sym_1%s /= sym_2%s) then 
            check_momentum_sym = .false. 
        end if

        ! old implo: 
!         call SetupSym(ka)
!         call SetupSym(kb) 
! 
!         do i = 1, size(elecs)
!             call AddElecSym(elecs(i), G1, nBasisMax, ka)
!         end do
!         do i = 1, size(orbs)
!             call AddElecSym(orbs(i), G1, nBasisMax, kb)
!         end do
!         
!         ! apply periodic BC:
!         call RoundSym(ka, nBasisMax)
!         call RoundSym(kb, nBasisMax)
! 
!         ! and check sym:
!         ! i want to switch from this old functionality.. 
!         ! since this works with these weird symconj functionality.. 
! !         check_momentum_sym = (lChkSym(ka, kb)) 
!         check_momentum_sym = sym_equal(ka,kb)

    end function check_momentum_sym

    logical function sym_equal(sym_1, sym_2) 
        type(BasisFN), intent(in) :: sym_1, sym_2 

        sym_equal = .true. 

        ! just check if every entries are the same! 
        if (.not. all(sym_1%k == sym_2%k)) sym_equal = .false. 
        if (sym_1%ms /= sym_2%ms) sym_equal = .false. 
        if (sym_1%ml /= sym_2%ml) sym_equal = .false. 
        if (sym_1%Sym%s /= sym_2%Sym%s) sym_equal = .false. 

    end function sym_equal 

    subroutine make_triple(nI, nJ, elecs, orbs, ex, tPar)
        integer, intent(in) :: nI(nel), elecs(3), orbs(3)
        integer, intent(out) :: nJ(nel), ex(2,3)
        logical, intent(out) :: tPar
#ifdef __DEBUG
        character(*), parameter :: this_routine = "make_triple"
#endif
        integer :: sort_elecs(3), sort_orbs(3), src(3), pos_moved, k, i

        ! i should also check if this excitation is actually possible! 
        ! and which spin i move to where or?? 

        ! figure out how to do triples efficiently.. 
        ! can we do a single and then a double? 

        ! NO: talk to Manu and do this with integer representation! not
        ! with nI and nJ, since this can be done way more effective as 
        ! via the occupied orbitals.. 

        ! TODO: thats a super strange convention, .. talk with Ali and 
        ! Simon about that.. but for now accept it as it is.. 

        sort_elecs = sort_unique(elecs)
        sort_orbs = sort_unique(orbs)
        
        src = nI(sort_elecs)

        ex(1,:) = src
        ex(2,:) = sort_orbs

        nJ = nI 

        ASSERT(sum(get_spin_pn(src)) == sum(get_spin_pn(orbs)))

        ! i should do some stuff depending on the order of src and tgt 
        ! we have to check if electrons hop over other electrons, so we 
        ! might have to change the indexing to adapt to the change in nJ! 

        ! or check it individually: 
        if (src(2) < sort_orbs(1)) then 
            ! then i hops over j: 
            sort_elecs(2) = sort_elecs(2) - 1 
        end if
        if (src(3) < sort_orbs(1)) then 
            ! then i hops over k 
            ! (note: this also implies that j hops over k, but treat that 
            ! seperately below, to also cover the case, where this if here 
            ! is not fullfilled!) 
            sort_elecs(3) = sort_elecs(3) - 1 
        end if 
        if (src(3) < sort_orbs(2)) then 
            ! then j hops over k 
            sort_elecs(3) = sort_elecs(3) - 1
        end if

        pos_moved = 0 

        do k = 1, 3 
            if (src(k) < sort_orbs(k)) then 
                if (sort_elecs(k) == nel) then 
                    ! this can only happen for k == 3 
                    i = nel + 1
                    nJ(nel) = sort_orbs(k) 
                else 
                    do i = sort_elecs(k) + 1, nel 
                        if (sort_orbs(k) < nJ(i)) then 
                            nJ(i-1) = sort_orbs(k)
                            exit 
                        else 
                            nJ(i-1) = nJ(i)
                        end if
                    end do
                    if (i == nel + 1) then 
                        nJ(nel) = sort_orbs(k)
                    end if
                end if
            else 
                if (sort_elecs(k) == 1) then 
                    i = 0
                    nJ(1) = sort_orbs(k)
                else 
                    do i = sort_elecs(k)-1, 1, -1
                        if (sort_orbs(k) > nJ(i)) then 
                            nJ(i+1) = sort_orbs(k)
                            exit 
                        else 
                            nJ(i+1) = nJ(i)
                        end if
                    end do
                    if (i == 0) then 
                        nJ(1) = sort_orbs(k)
                    end if
                end if
            end if

            pos_moved = pos_moved + sort_elecs(k) - i + 1

        end do

        tPar = btest(pos_moved, 0)
        
    end subroutine make_triple

    subroutine init_tmat_kspace(in_lat) 
        ! similar to the real-space tmat setup also do this here based on 
        ! the inputted lattice! 
        class(lattice), optional :: in_lat  
        character(*), parameter :: this_routine = "init_tmat_kspace"

        integer :: i 

        if (present(in_lat)) then 
            if (associated(tmat2d)) deallocate(tmat2d)

            allocate(tmat2d(nbasis, nbasis))
            tmat2d = 0.0_dp 

            do i = 1, in_lat%get_nsites() 
                tmat2d(2*i-1,2*i-1) = bhub * in_lat%dispersion_rel_orb(i) 
                tmat2d(2*i,2*i) = bhub * in_lat%dispersion_rel_orb(i)
            end do

        else 
            call Stop_All(this_routine, "not yet implemented!")
        end if

    end subroutine init_tmat_kspace

    subroutine setup_tmat_k_space(in_lat)
        ! routine which sets up the (diagonal) t-matrix in the k-space 
        ! the dimensionality and connectivity of nearest and next-nearest 
        ! neighbors influences that! 
        class(lattice), intent(in), optional :: in_lat 
        character(*), parameter :: this_routine = "setup_tmat_k_space" 

        if (present(in_lat)) then 
            if (all(nBasisMax == 0)) then 
                call setup_nbasismax(in_lat)
            end if

            if (.not. associated(G1)) then 
                call setup_g1(in_lat) 
            end if
            ! else assume it is already setup correctly 

            if (.not. associated(tmat2d)) then 
                ! call the already implemented hubbard tmat calculator.. 
                ! for now only.. in the future we should do this standalone
                call CALCTMATHUB(in_lat%get_nsites()*2, nBasisMax, bhub, & 
                    ttilt,G1,.not. in_lat%is_k_space(), in_lat%is_periodic())

            end if

        else 
            call Stop_All(this_routine, "not yet implemented!")
        end if

    end subroutine setup_tmat_k_space

end module k_space_hubbard
