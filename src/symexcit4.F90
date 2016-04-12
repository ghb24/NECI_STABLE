#include<macros.h>

module SymExcit4

    use SystemData, only : nEl, nBasis
    use SymExcitDataMod, only : SpinOrbSymLabel
    use GenRandSymExcitNUMod, only : RandExcitSymLabelProd

    implicit none

    type excitIndexVecsType
        integer, allocatable :: elecIndices
        integer, allocatable :: holeIndices
        integer, allocatable :: elecSpinOrbs
        integer, allocatable :: holeSpinOrbs
        integer :: elecTotSpin
        integer :: holeTotSpin
        integer :: elecSymLabel
        integer :: holeSymLabel
    end type

    type excitGenSessionType
        ! the (inclusive) lower and upper bounds on the excitation order
        ! where order is the number of elecs moved
        integer :: lOrder, uOrder, order
        ! the (inclusive) lower and upper bounds on the overall absolute
        ! spin number of the excitation
        integer :: lSpinDiff, uSpinDiff
        integer, allocatable :: counts(:)
        type(excitIndexVecsType) eivs
        ! parent determinant
        integer :: nI(NEl)
        ! unoccupied spin orbitals
        integer :: holes(nBasis-nEl)
        logical :: ti_lt_a_only

    end type


    contains

    function initExcitGenSession(nI, lOrder, uOrder, lSpinDiff, uSpinDiff, ti_lt_a_only) result (session)
        integer, intent(in) :: nI
        integer, intent(in) :: lOrder, uOrder, lSpinDiff, uSpinDiff
        logical, intent(in) :: ti_lt_a_only
        type(ExcitGenSessionType) :: session
        integer :: i, orb

        session%nI = nI
        session%lOrder = lOrder
        session%uOrder = uOrder
        session%lSpinDiff = lSpinDiff
        session%uSpinDiff = uSpinDiff
        session%ti_lt_a_only = ti_lt_a_only
        session%order = lOrder
        session%eivs = initExcitIndexVecs(lOrder)
        i = 1
        do orb = 1, nBasis
            if (.not. any(session%nI==orb)) then
                session%holes(j) = orb
                i = i+1
            endif
        enddo
    end function

    subroutine goToNextOrder(session, tReachedLimit)
        type(excitGenSessionType), intent(inout) :: session
        logical, intent(out) :: tReachedLimt
        tReachedLimit = .false.
        if (session%order<session%uOrder) then
            session%order = session%order+1
            session%eivs = initExcitIndexVecs(session%order)
        else 
            tReachedLimit = .true.
        endif
    end subroutine

    subroutine resetIndices(vec, maxIdx)
        integer, intent(inout) :: vec(:)
        integer, intent(in) :: maxIdx
        integer :: i
        do i = 1, maxIdx
            vec(i) = i
        enddo
    end subroutine

    function initExcitIndexVecs(order) result (eiv)
        integer, intent(in) :: order
        type(excitIndexVecsType) :: eivs
        allocate(eivs%elecIndices(order))
        allocate(eivs%holeIndices(order))
        call resetIndices(eivs%elecIndices, order)
        call resetIndices(eivs%holeIndices, order)
        allocate(eivs%elecSpinOrbs(order))
        allocate(eivs%holeSpinOrbs(order))
    end function

    subroutine incrementIndex(eiv, order, limit, tReachedLimit)
        integer, intent(inout) :: eiv(:)
        integer, intent(in) :: order, limit
        logical, intent(inout) :: tReachedLimit
        integer :: i

        tReachedLimit = .false.

        ! start from the left, increment any element which differs in value by at least 2
        ! with its right-hand neighbour

        do i=1,order
            if (i==order) then
                if (eiv(order)<limit) then
                    eiv(order) = eiv(order) + 1
                    call resetIndices(eiv, order-1)
                    return
                else
                    tReachedLimit = .true.
                    return 
                endif
            else if (eiv(i+1)-eiv(i)>1) then
                eiv(i) = eiv(i) + 1
                call resetIndices(eiv, i-1)
                return
            endif
        enddo
    end subroutine

    function getPosIdentifier(vec, order, limit) result (posId)
        ! this is not a triangular index, it is only computed for the
        ! purpose of evaluating ti_lt_a_only condition
        integer, intent(in) :: vec(:)
        integer :: posId, i
        posId = 0
        do i = 1, order
            posId = posId + vec(order)+limit*(order-1)
        enddo
    end function

    function getElecPosIdentifier(session) result (posId)
        type(excitGenSessionType), intent(in) :: session
        integer :: posId
        posId = getPosIdentifier(session%eivs%elecIndices, session%order, nEl)
    end function

    function getHolePosIdentifier(session) result (posId)
        type(excitGenSessionType), intent(in) :: session
        integer :: posId
        posId = getPosIdentifier(session%eivs%holeIndices, session%order, nBasis-nEl)
    end function

    subroutine goToNextElecIndices(session, tReachedLimit)
        integer, intent(inout) :: session
        logical, intent(out) tReachedLimit
        call incrementIndex(session%eivs%elecIndices, session%order, nEl, tReachedLimit)
        ! reset unoccupied
        call resetIndices(session%eivs%holeIndices, session%order)
    end subroutine

    subroutine goToNextHoleIndices(session, tReachedLimit)
        integer, intent(inout) :: session
        logical, intent(out) tReachedLimit
        call incrementIndex(session%eivs%holeIndices, session%order, nBasis-nEl, tReachedLimit)
    end subroutine

    subroutine setSpinOrbs(session)
        type(excitGenSessionType), intent(inout) :: session
        integer :: i
        ! initialise total spin
        session%eivs%elecTotalSpin = 0
        session%eivs%elecSymLabel = 1
        do i = 1, nEl
            ! fill selected spin orbs from parent determinant
            session%eivs%elecSpinObs(i) = session%nI(session%eivs%elecIndices(i))
            ! set total spin value
            if (btest(session%eivs%elecSpinOrbs(i), 0) then
                ! odd => beta
                session%eivs%elecTotalSpin = session%eivs%elecTotalSpin - 1
            else
                session%eivs%elecTotalSpin = session%eivs%elecTotalSpin + 1
            endif
            ! accumulate total point group symmetry label
            session%eivs%elecSymLabel = RandExcitSymLabelProd( &
                session%eivs%elecSymLabel, SpinOrbSymLabel(session%eivs%elecSpinOrbs(i)))
        enddo
        session%eivs%holeTotalSpin = 0
        session%eivs%holeSymLabel = 1
        do i = nBasis-nEl
            ! fill selected spin orbs from unoccupied complement of parent determinant
            session%eivs%holeSpinObs(i) = session%holes(session%eivs%holeIndices(i))
            ! set total spin value
            if (btest(session%eivs%holeSpinOrbs(i), 0) then
                ! odd => beta
                session%eivs%holeTotalSpin = session%eivs%holeTotalSpin - 1
            else
                session%eivs%holeTotalSpin = session%eivs%holeTotalSpin + 1
            endif
            ! accumulate total point group symmetry label
            session%eivs%holeSymLabel = RandExcitSymLabelProd( &
                session%eivs%holeSymLabel, SpinOrbSymLabel(session%eivs%holeSpinOrbs(i)))
        enddo
    end subroutine

    pure function inverseTriMap(i) result (xy)
        ! computes the matrix indices of a triangular map
        ! primary index x is ceiling of triangular root
        integer, intent(in) :: i
        integer :: x, y, xy(2)
        x = ceiling(((8*i+1)**5d-1-1)*5d-1)
        y = i - (x*(x-1))/2
        xy = (/ x, y /)
    end function
    
    subroutine genExcitations4(session, nJ, tParity, tAllExcitFound, ti_lt_a_only)
        type(ExcitGenSessionType), intent(inout) :: session
        integer, intent(out) :: nJ(nEl)
        logical, intent(out) :: tParity, tAllExcitFound
        logical, intent(in) :: ti_lt_a_only
        logical :: tReachedLimit
        integer :: elecPosIdentifier, holePosIdentifier
        integer :: i, spinDiff

        tAllExcitFound = .false.

        do
            ! exit when we find a legal excitation
            call goToNextHoleIndices(session, tReachedLimit)
            if (tReachedLimit) then
                ! reached the end of the combinations of unocc spin orbs
                call goToNextElecIndices(session, tReachedLimit)
                if (tReachedLimit) then
                    ! reached the end of the combinations of occ spin orbs
                    call goToNextElecIndices(session, tReachedLimit)
                    if (tReachedLimit) then
                        ! exhausted all excitation classes
                        tAllExcitFound = .true.
                        return
                    endif
                endif
            endif

            if (ti_lt_a_only)
                elecPosIdentifier = getElecPosIdentifier(session)
                holePosIdentifier = getHolePosIdentifier(session)
                if (elecPosIdentifier>=holePosIdentifier) then
                    cycle
                endif
            endif

            ! first, dereference the indices
            call setSpinOrbs(session)               
            ! we now have a new set of session%order electron-hole pairs
            ! to return an nJ, the selection must satisfy the following criteria
            ! 1. total spin difference between occupied and unoccupied spin orbs must
            !    lie within the bounds specified in the initialisation of the session
            ! 2. point group spatial symmetry must be conserved
            ! 3. mL must be conserved
            ! 4. k-point symmetry must be conserved

            ! 1. spin difference
            spinDiff = abs(session%eivs%elecTotal - session%eivs%holeTotal)
            if (spinDiff > session%uSpinDiff .or. spinDiff < session%lSpinDiff) then
                cycle
            endif

            ! 2. spatial symmetry
            if (session%eivs%elecSymLabel-session%eivs%holeSymLabel/=0)
                cycle
            endif
            
            ! if we made it this far, the selected electron hole pairs have
            ! passed all the logical tests
            return
        enddo

    end subroutine

    subroutine countExcitations4(nI, lOrder, uOrder, lSpinDiff, uSpinDiff, ti_lt_a_only, tot)
        integer, intent(in) :: nI(nEl), lOrder, uOrder, lSpinDiff, uSpinDiff
        logical, intent(in) :: ti_lt_a_only
        integer, intent(out) :: tot
        logical :: tAllExcitFound, tParity
        integer :: nJ(nEl)
        type(excitGenSessionType) :: session

        session = initExcitGenSession(nI, lOrder, uOrder, lSpinDiff, uSpinDiff, ti_lt_a_only)
        
        tot = -1
        tAllExcitFound = .false.
        while(.not.tAllExcitFound) do
            call genExctitations4(session, nJ, tParity, tAllExcitFound, ti_lt_a_only)
            tot+1
        enddo
end module
