#include<macros.h>
! robert.anderson@kcl.ac.uk
! April 2016
module SymExcit4

    use SystemData, only : nEl, nBasis, tkPntSym, tReltvy, G1, BasisFn
    use SymExcitDataMod, only : SpinOrbSymLabel
    use GenRandSymExcitNUMod, only : RandExcitSymLabelProd

    implicit none

    interface GenExcitations4
        module procedure GenExcitations4_non_initd
        module procedure GenExcitations4_initd
        module procedure GenExcitations4_compat_non_initd
    end interface

    type ExcitGenSessionType
        ! the (inclusive) lower and upper bounds on the excitation rank
        ! where rank is the number of elecs moved
        integer :: minRank, maxRank, rank
        ! the (inclusive) lower and upper bounds on the overall absolute
        ! spin number of the excitation
        integer :: minSpinDiff, maxSpinDiff
        ! parent determinant
        integer, allocatable :: nI(:)
        ! unoccupied spin orbitals
        integer, allocatable :: holes(:)
        ! the indices of nI corresponding to the selected electrons
        integer, allocatable :: elecIndices(:)
        ! the indices of holes corresponding to the selected unocc orbs
        integer, allocatable :: holeIndices(:)
        ! selected occupied spin orbitals
        integer, allocatable :: elecSpinOrbs(:)
        ! selected unoccupied spin orbitals
        integer, allocatable :: holeSpinOrbs(:)
        ! the difference of the numbers of beta spin orbitals being vacated 
        ! and filled by the excitation
        integer :: spinDiff
        ! the overall ml of the selected occ orbs
        integer :: elecTotMl
        ! the overall ml of the selected unocc orbs
        integer :: holeTotMl
        ! the overall point group symmetry label of the selected occ orbs
        integer :: elecSymLabel
        ! the overall point group symmetry label of the selected unocc orbs
        integer :: holeSymLabel
        ! stored for computing k-point symmetry
        type(BasisFn) :: nISym
        ! indicates whether the InitExcitGenSession method has been called
        logical :: tInitialised = .false.
    end type

    type(ExcitGenSessionType) :: storedSession

    contains

    function InitExcitGenSession(nI, minRank, maxRank, minSpinDiff, maxSpinDiff) result (session)
        use sym_mod, only : getsym_wrapper
        implicit none
        integer, intent(in) :: nI(nEl)
        integer, intent(in) :: minRank, maxRank, minSpinDiff, maxSpinDiff
        type(ExcitGenSessionType) :: session
        integer :: i, orb

        allocate(session%nI(nEl))
        allocate(session%holes(nBasis-nEl))
        
        session%nI = nI
        session%minRank = minRank
        session%maxRank = maxRank
        session%minSpinDiff = minSpinDiff
        session%maxSpinDiff = maxSpinDiff
        session%rank = minRank
        call InitExcitVecs(session) 
        i = 1
        do orb = 1, nBasis
            if (.not. any(session%nI==orb)) then
                session%holes(i) = orb
                i = i+1
            endif
        enddo
        if (tkPntSym) then
            call getsym_wrapper(session%nI, session%nISym)
        endif
        session%tInitialised = .true.
    end function InitExcitGenSession

    subroutine ResetIndices(vec, maxIdx)
        implicit none
        integer, intent(inout) :: vec(:)
        integer, intent(in) :: maxIdx
        integer :: i
        do i = 1, maxIdx
            vec(i) = i
        enddo
    end subroutine ResetIndices

    subroutine InitExcitVecs(session)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        call DestructExcitVecs(session)
        allocate(session%elecIndices(session%rank))
        allocate(session%holeIndices(session%rank))
        ! to indicate initialised indices
        session%elecIndices(1) = 0
        session%holeIndices(1) = 0
        allocate(session%elecSpinOrbs(session%rank))
        allocate(session%holeSpinOrbs(session%rank))
    end subroutine InitExcitVecs

    subroutine DestructExcitVecs(session)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        if (allocated(session%elecIndices)) deallocate(session%elecIndices)
        if (allocated(session%holeIndices)) deallocate(session%holeIndices)
        if (allocated(session%elecSpinOrbs)) deallocate(session%elecSpinOrbs)
        if (allocated(session%holeSpinOrbs)) deallocate(session%holeSpinOrbs)
    end subroutine DestructExcitVecs

    subroutine DestructSession(session)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        if (allocated(session%nI)) deallocate(session%nI)
        if (allocated(session%holes)) deallocate(session%holes)
        call DestructExcitVecs(session)
    end subroutine DestructSession

    subroutine IncrementIndex(vec, rank, limit, tReachedLimit)
        implicit none
        integer, intent(inout) :: vec(:)
        integer, intent(in) :: rank, limit
        logical, intent(inout) :: tReachedLimit
        integer :: i

        tReachedLimit = .false.

        ! start from the left, increment any element which differs in value by at least 2
        ! with its right-hand neighbour

        ! if we are starting from empty indices:
        if (vec(1)==0) then
            call resetIndices(vec, rank)
            return
        endif

        do i=1,rank
            if (i==rank) then
                if (vec(rank)<limit) then
                    vec(rank) = vec(rank) + 1
                    call ResetIndices(vec, rank-1)
                    return
                else
                    tReachedLimit = .true.
                    return 
                endif
            else if (vec(i+1)-vec(i)>1) then
                vec(i) = vec(i) + 1
                call ResetIndices(vec, i-1)
                return
            endif
        enddo
    end subroutine IncrementIndex

    function GetPosIdentifier(vec, rank, limit) result (posId)
        implicit none
        ! this is not a triangular index, it is only computed for the
        ! purpose of evaluating ti_lt_a_only condition
        integer, intent(in) :: vec(:), rank, limit
        integer :: posId, i
        posId = 0
        do i = 1, rank
            posId = posId + vec(rank)+limit*(rank-1)
        enddo
    end function GetPosIdentifier

    function GetElecPosIdentifier(session) result (posId)
        implicit none
        type(ExcitGenSessionType), intent(in) :: session
        integer :: posId
        posId = GetPosIdentifier(session%elecIndices, session%rank, nEl)
    end function GetElecPosIdentifier

    function GetHolePosIdentifier(session) result (posId)
        implicit none
        type(excitGenSessionType), intent(in) :: session
        integer :: posId
        posId = GetPosIdentifier(session%holeIndices, session%rank, nBasis-nEl)
    end function GetHolePosIdentifier

    subroutine GoToNextRank(session, tReachedLimit)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        logical, intent(out) :: tReachedLimit
        tReachedLimit = .false.
        if (session%rank<session%maxRank) then
            session%rank = session%rank+1
            call initExcitVecs(session)
            ! reset bot occupied and unoccupied indices
            call ResetIndices(session%elecIndices, session%rank)
            call ResetIndices(session%holeIndices, session%rank)
        else 
            tReachedLimit = .true.
        endif
    end subroutine GoToNextRank

    subroutine GoToNextElecIndices(session, tReachedLimit)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        logical, intent(out) :: tReachedLimit
        call IncrementIndex(session%elecIndices, session%rank, nEl, tReachedLimit)
        ! reset unoccupied indices
        call ResetIndices(session%holeIndices, session%rank)
    end subroutine GoToNextElecIndices

    subroutine GoToNextHoleIndices(session, tReachedLimit)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        logical, intent(out) :: tReachedLimit
        call IncrementIndex(session%holeIndices, session%rank, nBasis-nEl, tReachedLimit)
    end subroutine GoToNextHoleIndices

    subroutine SetSpinOrbs(session)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        integer :: i, elecBetaOrbCount, holeBetaOrbCount
        ! initialise totals
        elecBetaOrbCount = 0
        session%elecTotMl = 0
        session%elecSymLabel = 1

        holeBetaOrbCount = 0
        session%holeTotMl = 0
        session%holeSymLabel = 1

        do i = 1, session%rank
            ! fill selected spin orbs from parent determinant
            session%elecSpinOrbs(i) = session%nI(session%elecIndices(i))
            session%holeSpinOrbs(i) = session%holes(session%holeIndices(i))
            ! set total spin value
            if (is_beta(session%elecSpinOrbs(i))) then
                ! odd => beta
                elecBetaOrbCount = elecBetaOrbCount + 1
            endif
            if (is_beta(session%holeSpinOrbs(i))) then
                ! odd => beta
                holeBetaOrbCount = holeBetaOrbCount + 1
            endif
            session%spinDiff = abs(elecBetaOrbCount - holeBetaOrbCount)
            ! set total ml value
            session%elecTotMl = session%elecTotMl + G1(session%elecSpinOrbs(i))%Ml
            session%holeTotMl = session%holeTotMl + G1(session%holeSpinOrbs(i))%Ml
            ! accumulate total point group symmetry label
            session%elecSymLabel = RandExcitSymLabelProd( &
                session%elecSymLabel, SpinOrbSymLabel(session%elecSpinOrbs(i)))
            session%holeSymLabel = RandExcitSymLabelProd( &
                session%holeSymLabel, SpinOrbSymLabel(session%holeSpinOrbs(i)))
        enddo

    end subroutine SetSpinOrbs

    subroutine FindNewDet(session, nJ, tParity)
        implicit none
        type(ExcitGenSessionType), intent(in) :: session
        integer, intent(out) :: nJ(nEl)
        logical, intent(out) :: tParity
        integer :: nJtmp(nEL)
        integer :: minOrbBound, maxOrbBound
        integer :: holes(nBasis-nEl)
        integer :: i, j, switchedElecs
        integer :: elecIdx, elecOrb, holeOrb
        nJ = session%nI
        holes = session%holes
        
        tParity = .false. ! even parity

        do i = 1, session%rank
            ! loop over electron hole pairs
            ! spin orb arrays in session do the same job as an excitmat

            elecOrb = session%elecSpinOrbs(i)
            holeOrb = session%holeSpinOrbs(i)

            minOrbBound = min(elecOrb, holeOrb)
            maxOrbBound = max(elecOrb, holeOrb)
            switchedElecs = 0
            do j = 1, nEl
                ! the number of electrons inbetween determines parity
                ! and also determines insertion position
                if (nJ(j)>minOrbBound .and. nJ(j)<maxOrbBound) then
                    switchedElecs = switchedElecs + 1
                    tParity = .not. tParity
                endif

            enddo

            ! Suppose we have the following excitation
            !
            !   1  2  3  4  5  6  7  8   (spinorbs)
            !
            ! [ 1, 0, 1, 1, 1, 0, 0, 0 ]
            !              |
            !              V
            ! [ 1, 1, 0, 0, 1, 0, 1, 0 ]
            !
            ! this is a rank 2 excitation via:
            ! elecIndices = (2, 3)
            ! holeIndices = (3, 1)
            !
            ! which gives:
            ! elecSpinOrbs = (3, 4)
            ! holeSpinOrbs = (7, 2)
            !
            ! the first excitation is 3 -> 7
            !
            ! nJ was originally:
            ! ( 1, 3, 4, 5 )
            ! first iteration we should get:
            ! ( 1, 4, 5, 7 )
            ! this routine must yield:
            ! ( 1, 2, 5, 7 )
            !
            ! the loop above has given us the number of electrons
            ! interchanged in the excitation
            !
            ! we have all the information we need to manipulate nJ
            ! such that it is still ordered after the first excitation
            !
            ! for a hole orbital greater than the electron orbital
            ! the insertion point is at elecIdx + switchedElecs
            ! 
            ! 2 switched elecs
            ! use a temporary array nJtmp
            !
            ! business as usual for the part of the array before the moving electron
            ! nJtmp(1:elecIdx-1) = nJ(1:elecIdx-1)
            !    ( 1, 0, 0, 0 )
            !
            ! move down the elements between origin and destination indices
            ! nJtmp(elecIdx:elecIdx+switchedElecs-1) = nJ(elecIdx+1:elecIdx+switchedElecs)
            !    ( 1, 4, 5, 0 )
            !
            ! insert the electron at its destination orbital
            ! nJtmp(elecIdx+switchedElecs) = holeOrb
            !    ( 1, 4, 5, 7 )
            !
            ! finally, simply copy the rest of the array
            ! nJtmp(elecIdx+switchedElecs+1:nEl) = nJ(elecIdx+switchedElecs+1:nEl)
            !    ( 1, 4, 5, 7 )
            !
            ! 
            ! the second excitation is 4 -> 2
            !
            ! likewise, for an electron orbital greater than the hole
            ! orbital the insertion point is at elecIdx - switchedElecs
            !
            ! but there are 0 switched elements, so simply do
            ! nJ(elecIdx) = holeOrb

            
            ! we can't just take elecIdx from session%elecIndices, because of the re-ordering
            ! that occurs in this loop, so first we have to find it

            do j = 1, nEl
                if (nJ(j)==elecOrb) then
                    elecIdx = j
                endif
            enddo

            ! then construct the arrays
            if (switchedElecs == 0) then
                ! this is easy
                nJ(elecIdx) = holeOrb
            else
                if (holeOrb-elecOrb>0) then
                    nJtmp(1:elecIdx-1) = nJ(1:elecIdx-1)
                    nJtmp(elecIdx:elecIdx+switchedElecs-1) = nJ(elecIdx+1:elecIdx+switchedElecs)
                    nJtmp(elecIdx+switchedElecs) = holeOrb
                    nJtmp(elecIdx+switchedElecs+1:nEl) = nJ(elecIdx+switchedElecs+1:nEl)
                else
                    ! same as above, but with elecIdx -> elecIdx - switchedElecs
                    nJtmp(1:elecIdx-switchedElecs-1) = nJ(1:elecIdx-switchedElecs-1)
                    nJtmp(elecIdx-switchedElecs:elecIdx-1) = nJ(elecIdx-switchedElecs+1:elecIdx)
                    nJtmp(elecIdx) = holeOrb
                    nJtmp(elecIdx+1:nEl) = nJ(elecIdx+1:nEl)
                endif
                ! copy back into target det array
                nJ(1:nEl) = nJtmp
            endif
        enddo
    end subroutine FindNewDet

    pure function InverseTriMap(i) result (xy)
        implicit none
        ! computes the matrix indices of a triangular map
        ! primary index x is ceiling of triangular root
        integer, intent(in) :: i
        integer :: x, y, xy(2)
        x = ceiling(((8*i+1)**5d-1-1)*5d-1)
        y = i - (x*(x-1))/2
        xy = (/ x, y /)
    end function InverseTriMap
   
    subroutine GenExcitations4_non_initd(session, nI, nJ, tParity, tAllExcitFound, ti_lt_a_only)
        ! this is the variant for use in the case that session has not been initialised
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: nJ(nEl)
        logical, intent(out) :: tParity, tAllExcitFound
        logical, intent(in) :: ti_lt_a_only
        logical :: tReachedLimit
        integer :: elecPosIdentifier, holePosIdentifier
        integer :: i, spinDiff
        integer :: minRankDefault, maxRankDefault, minSpinDiffDefault, maxSpinDiffDefault
        ! this is defined in symrandexcit2.F90, but it is external to the
        ! main module defined therein
        logical :: IsMomAllowedDetAnyParent
        
        tAllExcitFound = .false.

        if (.not. session%tInitialised) then
            minRankDefault = 1
            maxRankDefault = 2
            minSpinDiffDefault = 0
            if (tReltvy) then
               maxSpinDiffDefault = 2
           else
               maxSpinDiffDefault = 2
           endif
           session = InitExcitGenSession(nI, minRankDefault, maxRankDefault, minSpinDiffDefault, maxSpinDiffDefault)
        endif

        do while(.true.)
            ! exit loop when we find a legal excitation
            if (session%elecIndices(1)==0) then
                ! new session
                ! hole indices get reset in the following call:
                call GoToNextElecIndices(session, tReachedLimit)
            else
                call GoToNextHoleIndices(session, tReachedLimit)
                if (tReachedLimit) then
                    ! reached the end of the combinations of unocc spin orbs
                    call GoToNextElecIndices(session, tReachedLimit)
                    if (tReachedLimit) then
                        ! reached the end of the combinations of occ spin orbs
                        call GoToNextRank(session, tReachedLimit)
                        if (tReachedLimit) then
                            ! exhausted all excitation classes
                            tAllExcitFound = .true.
                            exit
                        endif
                    endif
                endif
            endif

            if (ti_lt_a_only) then
                elecPosIdentifier = getElecPosIdentifier(session)
                holePosIdentifier = getHolePosIdentifier(session)
                if (elecPosIdentifier>=holePosIdentifier) then
                    cycle
                endif
            endif

            ! first, dereference the indices
            call setSpinOrbs(session)               
            ! we now have a new set of session%rank electron-hole pairs
            ! to return an nJ, the selection must satisfy the following criteria
            ! 1. total spin difference between occupied and unoccupied spin orbs must
            !    lie within the bounds specified in the initialisation of the session
            ! 2. point group spatial symmetry must be conserved
            ! 3. mL must be conserved
            ! 4. k-point symmetry must be conserved

            ! 1. spin difference
            if (session%spinDiff > session%maxSpinDiff .or. session%spinDiff < session%minSpinDiff) then
                cycle
            endif

            ! 2. spatial symmetry
            if (session%elecSymLabel-session%holeSymLabel/=0) then
                cycle
            endif

            ! 3. mL number conservation
            if (session%elecTotMl-session%holeTotMl/=0) then
                cycle
            endif

            ! 4. k-point symmetry
            ! at this point we have to generate the target determinant
            call FindNewDet(session, nJ, tParity)
            if (tkPntSym) then
                if (.not. IsMomAllowedDetAnyParent(nJ, session%nISym%Sym)) then
                    cycle
                endif
            endif
            
            ! if we made it this far, the selected electron hole pairs have
            ! passed all the selection criteria
            return
        enddo

        if (tAllExcitFound) then
            ! tidy up
            call DestructSession(session)
        endif
    end subroutine GenExcitations4_non_initd

    subroutine GenExcitations4_initd(session, nJ, tParity, tAllExcitFound, ti_lt_a_only)
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        integer, intent(out) :: nJ(nEl)
        logical, intent(out) :: tParity, tAllExcitFound
        logical, intent(in) :: ti_lt_a_only
        call GenExcitations4_non_initd(session, session%nI, nJ, tParity, tAllExcitFound, ti_lt_a_only)
    end subroutine GenExcitations4_initd


    subroutine GenExcitations4_compat_non_initd(session, nI, nJ, exFlag, excitMat, tParity, tAllExcitFound, ti_lt_a_only)
        ! this routine is only included to provide an interface consistent with that of
        ! the existing GenExcitations3. Here we have to assume a max rank of 2
        implicit none
        type(ExcitGenSessionType), intent(inout) :: session
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: nJ(nEl), exFlag, excitMat(2,2)
        logical, intent(out) :: tParity, tAllExcitFound
        logical, intent(in) :: ti_lt_a_only
        integer :: i
        call GenExcitations4_non_initd(session, nI, nJ, tParity, tAllExcitFound, ti_lt_a_only)
        if (tAllExcitFound) return
        exFlag = session%rank
        ! fill the excitation matrix
        excitMat(:,:) = 0
        do i = 1, exFlag
            excitMat(:,i) = (/ session%elecSpinOrbs(i), session%holeSpinOrbs(i) /)
        enddo
        
    end subroutine

    subroutine CountExcitations4(nI, minRank, maxRank, minSpinDiff, maxSpinDiff, tot)
        implicit none
        integer, intent(in) :: nI(nEl), minRank, maxRank, minSpinDiff, maxSpinDiff
        integer, intent(out) :: tot
        logical :: tAllExcitFound, tParity
        integer :: nJ(nEl)
        type(excitGenSessionType) :: session

        session = InitExcitGenSession(nI, minRank, maxRank, minSpinDiff, maxSpinDiff)
        
        tot = 0
        tAllExcitFound = .false.
        do while(.true.)
            call GenExcitations4(session, nJ, tParity, tAllExcitFound, .false.)
            if (tAllExcitFound) then
                exit
            endif
            tot = tot+1
        enddo
    end subroutine CountExcitations4
end module
