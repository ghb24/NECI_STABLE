#include "macros.h"
! GUGA Data module:
! containing all necessary declarations and basic function definitions for the
! GUGA approach.
module guga_data
    ! dependencies: be EXPLICIT about them!
    use SystemData, only: nBasis, tCSF, tSPN, tHPHF, lNoSymmetry, STOT, nEl, &
                          tNoBrillouin, tExactSizeSpace, tUHF, tGUGA
    use constants, only: dp, Root2, OverR2, n_int
    use MemoryManager, only: TagIntType

    implicit none

    ! ========================== type defs ===================================

    ! define types for the probabilistic weights functions used in the
    ! stochastic excitations generations
    type :: weight_data
        real(dp) :: F = 0.0_dp
        real(dp) :: G = 0.0_dp
        real(dp) :: minus = 0.0_dp
        real(dp) :: plus = 0.0_dp
        real(dp) :: zero = 0.0_dp
    end type weight_data

    ! define type structs for the probabilistic weights
    type :: branchWeight
        logical :: initialized = .false.
        real(dp) :: plusWeight = 0.0_dp
        real(dp) :: minusWeight = 0.0_dp
        real(dp) :: zeroWeight = 0.0_dp, &
                F = 0.0_dp, G = 0.0_dp, L = 0.0_dp
        procedure(dummyFunction), pointer, nopass :: minus => null(), plus => null(), &
                                                     zero => null()
    end type branchWeight

    ! define a type structure to store excitation information between two
    ! CSFs needed in the matrix element calculation between them
    ! this may also be used/needed for the excitation generation
    type :: ExcitationInformation_t
        integer :: typ = -1
        ! save type of excitation encoded as integer: all different possibs:
        ! 0 ... all kind of excitations which dont need much care
        ! 1 ... weight raising
        ! 2 ... weight lowering
        ! 3 ... non overlap
        ! 4 ... single overlap 2 lowering
        ! 5 ... single overlap 2 raising
        ! 6 ... single overlap lowering into raising
        ! 7 ... single overlap raising into lowering
        ! 8 ... normal double two lowering
        ! 9 ... normal double two raising
        ! 10 .. lowering into raising into lowering
        ! 11 .. raising into lowering into raising
        ! 12 .. lowering into raising double
        ! 13 .. raising into lowering double
        ! 14 .. full stop 2 lowering
        ! 15 .. full stop 2 raising
        ! 16 .. full stop lowering into raising
        ! 17 .. full stop raising into lowering
        ! 18 .. full start 2 lowering
        ! 19 .. full start 2 raising
        ! 20 .. full start lowering into raising
        ! 21 .. full start raising into lowering
        ! 22 .. full start into full stop alike
        ! 23 .. full start into full stop mixed

        ! need the involved indices of the excitation: list of integers
        ! for now the convention is, that they are given in an ordered form
        ! and is not related to the involved generators E_{ij} (E_{kl})
        ! directly, for single excitations ofc. only two entries of this
        ! vector needed.
        ! update:
        ! new convention store, original indiced and the ordered ones in
        ! the fullStart, etc. indices.
        integer :: i = -1
        integer :: j = -1
        integer :: k = -1
        integer :: l = -1
        integer :: fullStart = -1
        integer :: fullEnd = -1
        integer :: secondStart = -1
        integer :: firstEnd = -1
        integer :: weight = -1 ! can get rid of this in future!
        ! misuse secondstart firstend -> as weight as it is not used in the typ
        ! of excitations where weights is needed

        ! also need the overlap and nonOverlapRange of the excitations for
        ! efficiently identifying and calculatint excitations for a given CSF
        ! and indices i,j,k,l
!         integer, allocatable :: overlapRange(:), nonOverlapRange(:)
        ! update:
        ! dont need overlaprange, and nonoverlaprange anywhere, just need to
        ! indicate if its a non-overlap, single overlap or proper double!
        integer :: overlap = -1 ! in rework get rid of this too!
        ! not needed, %typ could be used to get same results
        ! since eg. calcRemainingSwitches is only needed in cases, where there
        ! is atleast one overlap site or in single excitations!
        ! where there is no overlap at all!
        ! 0 ... no overlap or single
        ! 1 ... single overlap
        ! >1 ... proper double overlap

        ! generator flags: necessary information on involved generators.
        ! could store it as flags (0,1) but that would mask the meaning and
        ! also ruin some matrix accessing functionality: so for now store
        ! lowering: -1
        ! raising : +1
        ! weight:    0
        ! maybe need firstGen, lastGen, (even secondGen) maybe?
        integer :: gen1 = -2
        integer :: gen2 = -2
        integer :: firstGen = -2
        integer :: lastGen = -2
        integer :: currentGen = -2
        ! also store excitation level(number of occupation number differences)
        ! in this type, to have everything at one place
        integer :: excitLvl = -1 ! definetly get rid of that! never used
        ! at all! -> UPDATE! with changing of relative probabilities of
        ! those excitations, i definetly need this type of information!
        ! misuse it in such a way, that i store 5 different types of double
        ! excitations!:
        ! 0.. (ii,jj) RR/LL
        ! 1.. (ii,jj) RL
        ! 2.. (ii,jk) RR/LL
        ! 3.. (ii,jk) RL
        ! 4.. (ij,kl) x

        ! additional flags:
        ! for a 4 index double excitation, there is an additional flag
        ! necessary to distiguish between certain kind of equally possible
        ! excitations:
!         logical :: fourFlag
        ! the order of generators is some excitations has an influence on the
        ! relative sign of the x_1 semi-stop matrix elements
        ! use a real here: 1.0 or -1.0 and just multiply x1 element
        real(dp) :: order = 0.0_dp
        real(dp) :: order1 = 0.0_dp
        ! maybe can get rid of order parameters.. since i could in general
        ! always choose such generators that the order parameter is +1
        ! but that did not work beforhand... ? hm
        !
        !TODO maybe more necessary info needed.
        ! add a flag to indicate if the excitation is even possible or other
        ! wise cancel the excitation
        logical :: valid = .false.
        ! for the exact calculation, to avoid calculating non-overlap
        ! contributions to matrix elements which are not possible, due to
        ! spin-coupling changes in the overlap range use a flag to
        ! indicate if a spin_change happened
        logical :: spin_change = .false.

    end type ExcitationInformation_t

    ! logical to indicate that GUGA and core space, like doubles and singles
    ! are used as the semi-stochastic core space
    logical :: tGUGACore = .false.

    ! also use a memory-tag for the exact excitations array to keep track
    ! how big those arrays get..  and also of the bigger temp_excits
    integer(TagIntType) :: tag_excitations = 0, tag_tmp_excits = 0
    ! also store the projected energy lists.. but only for 1 run not for
    ! all of them i guess.
    integer(TagIntType) :: tag_proje_list = 0

    ! ======================== end TYPE defs ================================

    ! ========================= INTERFACES ==================================
    ! create derived type to use array of procedure pointers to efficiently
    ! access the functions necessary for the matrix element calculation
    abstract interface
        function dummyFunction(bValue) result(ret)
            use constants, only: dp
            implicit none
            real(dp), intent(in) :: bValue
            real(dp) :: ret
        end function dummyFunction

        subroutine dummySubroutine(bValue, x0, x1)
            use constants, only: dp
            implicit none
            real(dp), intent(in) :: bValue
            real(dp), intent(out), optional :: x0, x1
        end subroutine dummySubroutine

    end interface

    type :: procedurePtrArray
        procedure(dummyFunction), pointer, nopass :: ptr => null()
        ! why the nopass flag is needed see:
        ! https://software.intel.com/en-us/forums/topic/508530
    end type procedurePtrArray

    type :: procPtrArrTwo
        procedure(dummySubroutine), pointer, nopass :: ptr => null()
    end type procPtrArrTwo

    ! all procedure pointers in an array have to use the same interface!
    ! no generic or differing ones for each element are allowed! -> this
    ! means that one has to use the most general one(3 input, 1 output) for
    ! all, even if everytime a 1 is given as output, as in many cases ->
    ! ask Simon if there is room for improvement here

    ! create an array of procedure pointers to all necessary functions for
    ! single excitation matrix element calculation
    type(procedurePtrArray) :: singleMatElesGUGA(15)


    ! define an array to give the correct indices to access matrix element
    ! function for given stepvectors, delta b, and generator type
    ! delta b value = {-1,1} and generator flags = {-1...lowering,1...raising}
    ! fill up with correct values, -1 corresponds to an impossible combination
    ! and will cause an error when tried to use as an index for matrix element
    ! calculation! are combined to a single index ranging from -1:2
    ! by accessing through: db + (G+1)/2
    integer, dimension(0:3,0:3,-1:2) :: indArrOne =  reshape( (/ &
    !& 00,10,20, 30,01,11, 21,31,02,12, 22,32, 03,13,23,33:
        2, 2, 2, -1, 2,10, 14, 5, 2, -1, 3, 7, -1, 7, 5, 3, & ! DeltaB = -1 & L
        2, 2, 2, -1, 2, 3, 12, 6, 2, -1,11, 8, -1, 5, 7, 3, & ! DeltaB = -1 & R
        2, 2, 2, -1, 2, 3, -1, 5, 2, 13,10, 7, -1, 7, 5, 3, & ! DeltaB = +1 & L
        2, 2, 2, -1, 2, 9, -1, 6, 2, 15, 3, 8, -1, 5, 7, 3  & ! DeltaB = +1 & R
        /), (/ 4, 4, 4 /))

    ! access matrix element terms just by
    !singleMatElesGUGA(indArrOne(d',d,dB,G,b))

    ! do the same for double-excitation matrix elements
!     type(procedurePtrArray) :: doubleMatElesGUGA(60)

    ! use two different arrays of procedure pointers for the x=0 and x=1
    ! double excitation matrix element
    ! x1 elements:
    type(procedurePtrArray) :: doubleMatEleX1GUGA(45)
    ! and make access index-matrix:
    ! probably an error here: element 45-> links to non-existent function
    ! funB(b,3,2) and may be in the wrong position too..
    integer, dimension(0:3,0:3,-7:7) :: indArrTwoX1 = reshape( (/ &
    !& !00,10, 20, 30, 01, 11, 21, 31, 02, 12, 22, 32, 03, 13, 23, 33:
        2, -1,  2, -1, -1, 34, 32, 18, -1, -1,  2, -1, -1, -1, -1,  2, &! db = -2 & LL
        2, -1,  2, -1,  2, 39, 43, 23, -1, -1, 39, -1, -1, -1, 15,  2, &! db = -2 & RL
        2, -1, -1, -1,  2,  2, 33, -1, -1, -1, 45, -1, -1, -1,  3,  2, &! db = -2 & RR
       -1, -1, -1, -1, 16, -1, -1, -1,  2, -1, -1, -1, -1,  3,  7, -1, &! db= -1 & LL
        1, 38,  7,  1, 12, 12, 40, 17, 38, 38,  7, 23,  1,  4, 20,  1, &! db = -1 & RL + full-start!
       -1,  2, 21, -1, -1, -1, -1, 12, -1, -1, -1, 22, -1, -1, -1, -1, &! db = -1 & RR->maybe store info here
        2, 21, 16,  1, -1, 31, 33, 28, -1, 32, 35, 26,  1, -1, -1,  2, &! db = 0 & LL
        2,  7, 12, -1, 28, 30, 44, 10, 26, 42, 34,  6, -1, 20, 17,  2, &! db = 0 & RL
        2, -1, -1,  1, 11, 31, 37, -1,  5, 29, 35, -1,  1,  7, 12,  2, &! db = 0 & RR
       -1, -1, -1, -1,  2, -1, -1, -1, 21, -1, -1, -1, -1, 12,  8, -1, &! db = +1 & LL
       -1, 12, 40, -1, 40, -1, -1, 25,  7, -1, -1, 20, -1, 17,  9, -1, &! db = +1 & RL
       -1, 16,  2, -1, -1, -1, -1, 24, -1, -1, -1,  7, -1, -1, -1, -1, &! db = +1 & RR
        2,  2, -1, -1, -1,  2, -1, -1, -1, 33, 30, 14, -1, -1, -1,  2, &! db = +2 & LL
        2,  2, -1, -1, -1, 41, -1, -1,  2, 43, 41, 25, -1, 19, -1,  2, &! db = +2 & RL
        2, -1, -1, -1, -1, 36, -1, -1,  2, 32,  2, -1, -1,  8, -1,  2  &! db = +2 & RR
        /), (/4, 4, 15/))

    ! x0 elemets:
    type(procedurePtrArray) :: doubleMatEleX0GUGA(17)

    ! build index matrix
    ! could make  some matrix elements indpendent of deltaB value, by filling
    ! every entry with the same matrix element, but for now choose to only
    ! make certain starts possible...
    ! also could point some impossible stepvector combinations to the
    ! zero function, and not making the index -1 -> decide on that!
    ! i could do smth like abort the matrix element calculation if its a -1
    ! index in the matrix, like abort the calcuation and output 0..

    ! REMEMBER: x0 matrix element only needed when Delta b = 0 all over excitation...
    ! maybe room for improvement there... -> finish that but think about that!
    integer, dimension(0:3,0:3,-7:7) :: indArrTwoX0 =  reshape( (/ &
    !& !00, 10,20, 30, 01,11,21, 31, 02, 12,22, 32, 03, 13, 23,33:
        1, -1,  1, -1, -1, 1, 1,  1, -1, -1, 1, -1, -1, -1, -1, 1, & ! db=-2 & LL -> this line always leads to 0 matEle...
        1, -1,  1, -1,  1, 1, 1,  1, -1, -1, 1, -1, -1, -1,  1, 1, & ! db=-2 & RL -> also this 0..
        1, -1, -1, -1,  1, 1, 1, -1,  1, -1, 1, -1, -1, -1,  1, 1, & ! db=-2 & RR -> this line is always 0...
        1, -1, -1, -1, 14, 1, 2, -1,  1, -1, 2, -1, -1,  1,  7, 2, & ! db=-1 & LL
        1,  1,  7, 10,  7, 7, 1, 15,  1,  1, 7,  1, 11,  1, 17, 5, & ! db=-1 & RL -> also store full-start elements here
        1,  1, 16, -1, -1, 2, 2,  7, -1, -1, 1,  1, -1, -1, -1, 2, & ! db=-1 & RR -> do include single overlap matrix elements!
        2, 16, 14,  4, -1, 3, 1,  7, -1,  1, 3,  7,  4, -1, -1, 2, & ! db= 0 & LL
        2,  6,  6, -1,  6, 2, 1, 13,  6,  1, 2, 12, -1, 16, 14, 2, & ! db= 0 & RL
        2, -1, -1,  4, 13, 3, 1, -1, 12,  1, 3, -1,  4,  7,  7, 2, & ! db= 0 & RR
        1, -1, -1, -1,  1, 2,-1, -1, 16,  2, 1, -1, -1,  7,  1, 2, & ! db=+1 & LL
       -1,  7,  1,  9,  1,-1,-1,  1,  7, -1,-1, 17,  8, 15,  1,-1, & ! db=+1 & RL
        1, 14,  1, -1, -1, 1,-1,  1, -1,  2, 2,  7, -1, -1, -1, 2, & ! db=+1 & RR
        1,  1, -1, -1, -1, 1,-1, -1, -1,  1, 1,  1, -1, -1, -1, 1, & ! db=+2 & LL -> always 0
        1,  1, -1, -1, -1, 1,-1, -1,  1,  1, 1,  1, -1,  1, -1, 1, & ! db=+2 & RL
        1, -1, -1, -1, -1, 1,-1, -1,  1,  1, 1, -1, -1,  1, -1, 1  & ! db=+2 & RR -> always 0
        /),(/ 4, 4, 15 /))
    ! to index third matrix dimension use: (db + 2)*3 + (G1 + G2)/2 + 2
    ! where G = +1 for R and -1 for L
    ! unfortunately for mixed generators, full-stop matrix element is needed
    ! seperately, because otherwise there would be ambiguities with RL
    ! intermediate matrix elements
    ! although procedure pointer arrays too much for this case.. as it is
    ! pretty special
    type(procPtrArrTwo) :: mixedGenFullStopMatEle(5)

    ! find good indexing converting function to efficiently access those
    ! functions. As only some combinations are needed, but the restriction on
    ! the delta b values has to be included here too.

    ! could write function or do matrix again, but almost only -1 in matrix...
    integer, dimension(0:3,0:3,-1:1) :: indArrEnd = reshape( [ &
        -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, & ! db = -2
         1, -1, -1, -1, -1,  2,-1, -1, -1, -1,  3, -1, -1, -1, -1,  4, & ! db = 0
        -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, -1  & ! db = +2
        ], [4,4,3])
    ! and create a indexing array to acces matrix elements through:
    ! doubleMatEle(indArrTwo(d',d',db,G1,G2,x0,x1,b))
    ! think about a good way to combine possible


    ! also need a similar procedure pointer array for the necessary calculation
    ! of the r_k terms, coming from the two-particle contributions to single
    ! excitations. they behave very similar to the usual single excitation
    ! product terms and hence, it is probably convenient to calculate them
    ! in a similar fashion. although it might be a overkill, since only 2
    ! stepvector combinations correspond to a b-dependent function and the
    ! rest are just constants... think ybout that efficiency and ask simon it
    ! this is too costly..
    type(procedurePtrArray) :: doubleContribution(7)


    ! write similar matrix indication table
    ! delete out the start and end values, as they are never needed!
    ! and should not be accessed!
    integer, dimension(0:3,0:3,-1:2) :: indContr = reshape( [ &
    !&! 00, 10, 20, 30, 01, 11, 21, 31, 02, 12, 22, 32, 03, 13, 23, 33
       1,  -1, -1, -1,  -1,  1,  7, -1, -1, -1,  3, -1, -1, -1, -1,  3, &! db = -1 & L
       1,  -1, -1, -1,  -1,  3,  5, -1, -1, -1,  1, -1, -1, -1, -1,  3, &! db = -1 & R
       1,  -1, -1, -1,  -1,  3, -1, -1, -1,  6,  1, -1, -1, -1, -1,  3, &! db = +1 & L
       1,  -1, -1, -1,  -1,  1, -1, -1, -1,  4,  3, -1, -1, -1, -1,  3 &! db = +1 & R
       ], [4,4,4])

    ! =========================== VARIABLES =================================
    ! b vector of the reference determinant should be kept as a saved variable
    ! as it is always needed in the H|D> calculation to calc. the refence energy
    ! not sure if still needed if a currentB_vector variable is used within
    ! the guga_excitation module... -> decide later how to implement
    ! ahh this is the bvecor of the reference determinant.. but make it real
    real(dp), allocatable :: bVectorRef_ilut(:), bVectorRef_nI(:)

    ! also need a list of determinants and matrix elements connceted to the
    ! reference determinant
    ! adapt that to multiple neci runs.. hope that works as intended..
    ! probably have to use it as a type, to store lists or different lists
    ! in it, and also be able to (de)allocate them individually
    type projE_type
        integer(n_int), allocatable :: projE_ilut_list(:,:)
        HElement_t(dp), allocatable :: projE_hel_list(:)
        ! also store the excitation level in the projected list, since otherwise
        ! it is really hard to determine it in the GUGA formalism
        integer, allocatable :: exlevel(:)
        ! also store the number of entries to correctly binary search
        integer :: num_entries
    end type projE_type

    type(projE_type), allocatable :: projE_replica(:)

    ! also make a global integer list of orbital indices, so i do not have to
    ! remake them in every random orbital picker!
    integer, allocatable :: orbitalIndex(:)

    ! also define a global variable nBasis/2 = nSpatOrbs, since otherwise
    ! integer division nBAsis/2 is done to often!
    ! but store that in SystemData module
!     integer :: nSpatOrbs

    ! in the end to make logic in excitation generation more efficient
    ! do create an allocatable integer array which stores the current
    ! stepvector for a CSF and do a select case() abfrage
!     integer, allocatable :: current_stepvector(:)

    ! use a global flag to indicate a switch to a new determinant in the
    ! main routine to avoid recalculating b vector occupation and
    ! stepvector
    logical :: tNewDet

    ! define a new flag to specify if we really want to run the full test
    ! suite -> You have to ensure that H = 1 is used otherwise it crashes!
    ! for now, still define that in SystemData, like the other guga stuff.
!      logical :: t_full_guga_tests

contains

    subroutine init_guga_data_procPtrs()
        ! this subroutine initializes the procedure pointer arrays needed for
        ! the matrix element calculation

        ! -------- single excitation procedure pointers:---------------------
        ! and point them to the corresponding functions -> for now write down all
        ! differing functions explicetly, as i can't think of somw fancy way to
        ! include that better
        ! store only procedure pointer to differing functions and convert to
        ! indices in such a way that correct ones get picked!
        ! also store weight-generator matrix element, as it may be of some use
        singleMatElesGUGA(1)%ptr => funZero
        singleMatElesGUGA(2)%ptr => funPlus1
        singleMatElesGUGA(3)%ptr => funMinus1
        singleMatElesGUGA(4)%ptr => funTwo
        singleMatElesGUGA(5)%ptr => funA_0_1
        singleMatElesGUGA(6)%ptr => funA_1_0
        singleMatElesGUGA(7)%ptr => funA_2_1
        singleMatElesGUGA(8)%ptr => funA_1_2
        singleMatElesGUGA(9)%ptr => funC_0
        singleMatElesGUGA(10)%ptr => funC_1
        singleMatElesGUGA(11)%ptr => funC_2
        singleMatElesGUGA(12)%ptr => funOverB_0
        singleMatElesGUGA(13)%ptr => funOverB_1
        singleMatElesGUGA(14)%ptr => minFunOverB_1
        singleMatElesGUGA(15)%ptr => minFunOverB_2

        ! ----------- double excitation procedure pointer --------------------

        ! ------ double excitation x0 matrix element part --------------------
        ! and fill it up with all unique values needed for x=1 matrix elements
        doubleMatEleX1GUGA(1)%ptr => funZero
        doubleMatEleX1GUGA(2)%ptr => funPlus1
        !     doubleMatEleX1GUGA(3)%ptr => funMinus1
        doubleMatEleX1GUGA(3)%ptr => funA_3_2
        doubleMatEleX1GUGA(4)%ptr => minFunA_3_2
        doubleMatEleX1GUGA(5)%ptr => funA_3_2overR2
        doubleMatEleX1GUGA(6)%ptr => minFunA_3_2overR2
        !     doubleMatEleX1GUGA(8)%ptr => funA_0_2overR2
        doubleMatEleX1GUGA(7)%ptr => minFunA_0_2_overR2
        doubleMatEleX1GUGA(8)%ptr => funA_m1_0
        doubleMatEleX1GUGA(9)%ptr => minFunA_m1_0
        doubleMatEleX1GUGA(10)%ptr => funA_m1_0_overR2
        doubleMatEleX1GUGA(11)%ptr => minFunA_m1_0_overR2
        doubleMatEleX1GUGA(12)%ptr => funA_2_0_overR2
        doubleMatEleX1GUGA(13)%ptr => minFunA_2_0_overR2
        doubleMatEleX1GUGA(14)%ptr => funA_2_1
        doubleMatEleX1GUGA(15)%ptr => minFunA_2_1
        doubleMatEleX1GUGA(16)%ptr => funA_2_1_overR2
        doubleMatEleX1GUGA(17)%ptr => minFunA_2_1_overR2
        doubleMatEleX1GUGA(18)%ptr => funA_0_1
        doubleMatEleX1GUGA(19)%ptr => minFunA_0_1
        doubleMatEleX1GUGA(20)%ptr => funA_0_1_overR2
        doubleMatEleX1GUGA(21)%ptr => minFunA_0_1_overR2
        doubleMatEleX1GUGA(22)%ptr => funA_1_2
        doubleMatEleX1GUGA(23)%ptr => minFunA_1_2
        doubleMatEleX1GUGA(24)%ptr => funA_1_0
        doubleMatEleX1GUGA(25)%ptr => minFunA_1_0
        doubleMatEleX1GUGA(26)%ptr => funA_3_1_overR2
        doubleMatEleX1GUGA(27)%ptr => minFunA_3_1_overR2
        !     doubleMatEleX1GUGA(30)%ptr => funA_m1_1_overR2
        doubleMatEleX1GUGA(28)%ptr => minFunA_m1_1_overR2
        doubleMatEleX1GUGA(29)%ptr => funB_2_3
        doubleMatEleX1GUGA(30)%ptr => funD_0
        doubleMatEleX1GUGA(31)%ptr => minFunD_0
        doubleMatEleX1GUGA(32)%ptr => funB_1_2
        doubleMatEleX1GUGA(33)%ptr => funB_0_1
        doubleMatEleX1GUGA(34)%ptr => funD_1
        doubleMatEleX1GUGA(35)%ptr => minFunD_1
        doubleMatEleX1GUGA(36)%ptr => funD_m1
        doubleMatEleX1GUGA(37)%ptr => funB_m1_0
        doubleMatEleX1GUGA(38)%ptr => funC_2
        doubleMatEleX1GUGA(39)%ptr => minFunC_2
        doubleMatEleX1GUGA(40)%ptr => funC_0
        doubleMatEleX1GUGA(41)%ptr => minFunC_0
        doubleMatEleX1GUGA(42)%ptr => minFunOverB_2_R2
        doubleMatEleX1GUGA(43)%ptr => minFunB_0_2
        doubleMatEleX1GUGA(44)%ptr => minFunOverB_0_R2
        doubleMatEleX1GUGA(45)%ptr => funD_2


        ! -------------- double excitation x0 matrix element part ------------
        doubleMatEleX0GUGA(1)%ptr => funZero
        doubleMatEleX0GUGA(2)%ptr => funPlus1
        doubleMatEleX0GUGA(3)%ptr => funMinus1
        doubleMatEleX0GUGA(4)%ptr => funSqrt2
        doubleMatEleX0GUGA(5)%ptr => minFunSqrt2
        doubleMatEleX0GUGA(6)%ptr => funOverRoot2
        doubleMatEleX0GUGA(7)%ptr => minFunOverR2
        doubleMatEleX0GUGA(8)%ptr => funA_0_1
        doubleMatEleX0GUGA(9)%ptr => funA_1_0
        doubleMatEleX0GUGA(10)%ptr => funA_1_2
        doubleMatEleX0GUGA(11)%ptr => funA_2_1
        doubleMatEleX0GUGA(12)%ptr => funA_1_2overR2
        doubleMatEleX0GUGA(13)%ptr => funA_1_0_overR2
        doubleMatEleX0GUGA(14)%ptr => funA_0_1_overR2
        doubleMatEleX0GUGA(15)%ptr => minFunA_0_1_overR2
        doubleMatEleX0GUGA(16)%ptr => funA_2_1_overR2
        doubleMatEleX0GUGA(17)%ptr => minFunA_2_1_overR2


        ! ------ mixed generator full-stop matrix elements--------------------
        mixedGenFullStopMatEle(1)%ptr => fullStop_00
        mixedGenFullStopMatEle(2)%ptr => fullStop_11
        mixedGenFullStopMatEle(3)%ptr => fullStop_22
        mixedGenFullStopMatEle(4)%ptr => fullStop_33
        mixedGenFullStopMatEle(5)%ptr => fullStop_12


        ! --------- double contributions to single excitaitons----------------
        doubleContribution(1)%ptr => funZero
        doubleContribution(2)%ptr => funPlus1
        doubleContribution(3)%ptr => funMinus1
        doubleContribution(4)%ptr => minFunBplus2
        doubleContribution(5)%ptr => funBplus0
        doubleContribution(6)%ptr => funBplus1
        doubleContribution(7)%ptr => minFunBplus1

    end subroutine init_guga_data_procPtrs

    subroutine nullify_guga_data_procPtrs()

        integer :: i

        do i = 1, 15
            nullify(singleMatElesGUGA(i)%ptr)
        end do

        do i = 1,45
            nullify(doubleMatEleX1GUGA(i)%ptr)
        end do

        do i = 1, 17
            nullify(doubleMatEleX0GUGA(i)%ptr)
        end do

        do i = 1, 5
            nullify(mixedGenFullStopMatEle(i)%ptr)
        end do

        do i = 1, 7
            nullify(doubleContribution(i)%ptr)
        end do

    end subroutine nullify_guga_data_procPtrs

    ! wrapper functions to access matrix element terms

    ! for consistency reasons also write an get.. function for the specific
    ! mixed gen full-stops.

    subroutine getMixedFullStop(step1, step2, deltaB, bValue, x0_element, &
            x1_element)
        ! function to access the special mixed generator full-stop elements
        ! which due to storage reasons are stored in a seperate func. pointer
        ! array
        integer, intent(in) :: step1, step2, deltaB
        real(dp), intent(in) :: bValue
        real(dp), intent(out), optional :: x0_element, x1_element
        integer :: ind

        ! get index:
        ind = indArrEnd(step1, step2, deltaB/2)

        ! with the optional output arguments can also just calc. x0 or x1
        call mixedGenFullStopMatEle(ind)%ptr(bValue, x0_element, x1_element)

    end subroutine getMixedFullStop

    function getDoubleContribution(step1,step2,deltaB,genFlag,bValue) &
            result (doubleContr)
        ! Access necessary two-particle contribution to single excitation
        ! matrix elements.
        !
        ! input:
        ! step1/2 ... stepvector values of CSFs, step1 is from <d'|
        ! deltaB  ... current delta b value of excitation
        ! genFlag ... generator flag of excitation
        ! bValue  ... current b value of CSF
        !
        ! output:
        ! doubleContr ... product term from shavitt graph rules used to calculate
        !               the matrix element between two given CSFs
        integer, intent(in) :: step1, step2, deltaB, genFlag
        real(dp), intent(in) :: bValue
        real(dp) :: doubleContr
        integer :: ind

        ! get index
        ind = indContr(step1, step2, deltaB + (genFlag + 1)/2)

        doubleContr = doubleContribution(ind)%ptr(bValue)
    end function getDoubleContribution

    function getSingleMatrixElement(step1,step2,deltaB,genFlag,bValue) &
            result(hElement)
        ! Access the necessary single excitation product terms for the H matrix
        ! element calculation.
        !
        ! input:
        ! step1/2 ... stepvector values of CSFs, step1 is from <d'|
        ! deltaB  ... current delta b value of excitation
        ! genFlag ... generator flag of excitation
        ! bValue  ... current b value of CSF
        !
        ! output:
        ! hElement ... product term from shavitt graph rules used to calculate
        !               the matrix element between two given CSFs
        integer, intent(in) :: step1, step2, deltaB, genFlag
        real(dp), intent(in) :: bValue
        real(dp) :: hElement
        integer :: ind

        ! only need for this function is to correctly access the procedure
        ! pointers to get correct function -> get index from index matrix:
        ind = indArrOne(step1, step2, deltaB + (genFlag + 1)/2)

        ! call correct function:
        hElement = singleMatElesGUGA(ind)%ptr(bValue)

    end function getSingleMatrixElement

    subroutine getDoubleMatrixElement(step1, step2, deltaB, genFlag1, genFlag2, &
            bValue, order, x0_element, x1_element)
        ! access the necessary double excitation product terms for the H matrix
        ! element calculation
        !
        ! input:
        ! step1/2   ...     stepvector values if given CSFs, step1 form <d'|
        ! deltaB    ...     current delta b value of excitation
        ! genFlag1/2...     generator types involved in excitation
        ! bValue    ...     current b value of CSF
        ! order     ...     order parameter determining the sign of x1 element
        !
        ! output:
        ! x0_element...     x0 product term (optional)
        ! x1_element...     x1 product term (optional)
        ! -> call this function by x0_element=var, x1_element=var2, if only
        ! specific matrix element is wanted. x0_element is default if only
        ! one output variable is present
        integer, intent(in) :: step1, step2, deltaB, genFlag1, genFlag2
        real(dp), intent(in) :: bValue, order
        real(dp), intent(out), optional :: x0_element, x1_element
        integer :: x0_ind, x1_ind

        if (present(x0_element)) then
            ! first get correct indices to access procedure pointer array
            x0_ind = indArrTwoX0(step1,step2,3*deltaB +(genFlag1 + genFlag2)/2)

            ! then call corresponding function
            x0_element = doubleMatEleX0GUGA(x0_ind)%ptr(bValue)
        end if

        ! same for x1 element
        if (present(x1_element)) then
            x1_ind = indArrTwoX1(step1,step2,3*deltaB+(genFlag1+genFlag2)/2)

            x1_element = order * doubleMatEleX1GUGA(x1_ind)%ptr(bValue)
        end if

    end subroutine getDoubleMatrixElement
    ! maybe write specific functions to only access x0/x1 matrix elements to
    ! avoid double conditioning...

    ! ===== special functions for mixed generator full-stop elements ========
    ! maybe create subroutines to output x0 and x1 matrix elements for all
    ! needed terms, and let the procedure pointers point to them
    subroutine fullStop_00(b, x0, x1)
        real(dp), intent(in) :: b
        real(dp), intent(out), optional :: x0, x1

        unused_variable(b)

        if (present(x0)) x0 = 0.0_dp
        if (present(x1)) x1 = 0.0_dp

    end subroutine fullStop_00

    subroutine fullStop_11(b, x0, x1)
        real(dp), intent(in) :: b
        real(dp), intent(out), optional :: x0, x1

        if (present(x0)) x0 = OverR2
        if (present(x1)) x1 = minFunA_m1_1_overR2(b)

    end subroutine fullStop_11

    subroutine fullStop_12(b, x0, x1)
        ! is the same for switched stepvector values
        real(dp), intent(in) :: b
        real(dp), intent(out), optional :: x0, x1

        unused_variable(b)

        if (present(x0)) x0 = 0.0_dp
        if (present(x1)) x1 = 1.0_dp

    end subroutine fullStop_12

    subroutine fullStop_22(b, x0, x1)
        real(dp), intent(in) :: b
        real(dp), intent(out), optional :: x0, x1

        if (present(x0)) x0 = OverR2
        if (present(x1)) x1 = funA_3_1_overR2(b)

    end subroutine fullStop_22

    subroutine fullStop_33(b, x0, x1)
        real(dp), intent(in) :: b
        real(dp), intent(out), optional :: x0, x1

        unused_variable(b)

        if (present(x0)) x0 = Root2
        if (present(x1)) x1 = 0.0_dp

    end  subroutine fullStop_33

    ! ===== special functions for the double contribution to single excitation
    ! matrix elements
    function minFunBplus2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -(b + 2.0_dp)
    end function minFunBplus2

    function funBplus0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = b
    end function funBplus0

    function minFunBplus1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -(b + 1.0_dp)
    end function minFunBplus1

    function funBplus1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = b + 1.0_dp
    end function funBplus1

    ! =========== additional double excitation matrix elements ===============
    function funMinusTwo(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = -2.0_dp
    end function funMinusTwo

    function funSqrt2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = Root2
    end function funSqrt2

    function minFunSqrt2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = -Root2
    end function minFunSqrt2

    function funOverRoot2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = OverR2
    end function funOverRoot2

    function minFunOverR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = -OverR2
    end function minFunOverR2

    function funA_1_2overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_1_2(b)
    end function funA_1_2overR2

    function funA_3_2overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_3_2(b)
    end function funA_3_2overR2

    function minFunA_3_2overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_3_2overR2(b)
    end function minFunA_3_2overR2

    function funA_0_2overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA(b, 0.0_dp, 2.0_dp)
    end function funA_0_2overR2

    function minFunA_0_2_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_0_2overR2(b)
    end function minFunA_0_2_overR2

    function funA_3_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, 3.0_dp, 2.0_dp)
    end function funA_3_2

    function minFunA_3_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_3_2(b)
    end function minFunA_3_2

    function funA_1_0_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_1_0(b)
    end function funA_1_0_overR2

    function funA_m1_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, -1.0_dp, 0.0_dp)
    end function funA_m1_0

    function minFunA_m1_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_m1_0(b)
    end function minFunA_m1_0

    function funA_m1_0_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_m1_0(b)
    end function funA_m1_0_overR2

    function minFunA_m1_0_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_m1_0_overR2(b)
    end function minFunA_m1_0_overR2

    function funA_2_0_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA(b, 2.0_dp, 0.0_dp)
    end function funA_2_0_overR2

    function minFunA_2_0_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_2_0_overR2(b)
    end function minFunA_2_0_overR2

    function funA_0_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_0_1(b)
    end function funA_0_1_overR2

    function minFunA_0_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_0_1_overR2(b)
    end function minFunA_0_1_overR2

    function funA_2_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA_2_1(b)
    end function funA_2_1_overR2

    function minFunA_2_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_2_1_overR2(b)
    end function minFunA_2_1_overR2

    function minFunA_1_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_1_2(b)
    end function minFunA_1_2

    function minFunA_1_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_1_0(b)
    end function minFunA_1_0

    function minFunA_0_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_0_1(b)
    end function minFunA_0_1

    function funA_3_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 * funA(b, 3.0_dp, 1.0_dp)
    end function funA_3_1_overR2

    function minFunA_3_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_3_1_overR2(b)
    end function minFunA_3_1_overR2

    function funA_m1_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = OverR2 *funA(b, -1.0_dp, 1.0_dp)
    end function funA_m1_1_overR2

    function minFunA_m1_1_overR2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_m1_1_overR2(b)
    end function minFunA_m1_1_overR2

    function funB_2_3(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funB(b, 2.0_dp, 3.0_dp)
    end function funB_2_3

    function funD_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funD(b, 2.0_dp)
    end function funD_2

    function funD_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funD(b, 0.0_dp)
    end function funD_0

    function minFunD_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funD_0(b)
    end function minFunD_0

    function funB_1_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funB(b, 1.0_dp, 2.0_dp)
    end function funB_1_2

    function funB_0_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funB(b, 0.0_dp, 1.0_dp)
    end function funB_0_1

    function funD_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funD(b, 1.0_dp)
    end function funD_1

    function minFunD_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funD_1(b)
    end function minFunD_1

    function funD_m1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funD(b, -1.0_dp)
    end function funD_m1

    function funB_m1_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funB(b, -1.0_dp, 0.0_dp)
    end function funB_m1_0

    ! mixed generator additional functions:
    function minFunA_2_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funA_2_1(b)
    end function minFunA_2_1

    function minFunC_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funC_2(b)
    end function minFunC_2

    function minFunC_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funC_0(b)
    end function minFunC_0

    function minFunOverB_2_R2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = Root2 * minFunOverB_2(b)
    end function minFunOverB_2_R2

    function minFunOverB_0_R2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -Root2 * funOverB_0(b)
    end function minFunOverB_0_R2

    function minFunB_0_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funB(b, 0.0_dp, 2.0_dp)
    end function minFunB_0_2






    !========= function for the single particle matrix calculation ===========
    ! ASSERT() probably not usable in "elemental" function, due to side-effects

    function funZero(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = 0.0_dp
    end function funZero

    function funPlus1(b) result(ret)
        ! think of a better way to include that -> wasted time
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = 1.0_dp
    end function funPlus1

    function funMinus1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = -1.0_dp
    end function funMinus1

    function funTwo(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        unused_variable(b)
        ret = 2.0_dp
    end function funTwo

    function funA_0_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, 0.0_dp, 1.0_dp)
    end function funA_0_1

    function funA_2_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, 2.0_dp, 1.0_dp)
    end function funA_2_1

    function funA_1_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, 1.0_dp, 0.0_dp)
    end function funA_1_0

    function funA_1_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funA(b, 1.0_dp, 2.0_dp)
    end function funA_1_2

    function funC_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funC(b, 0.0_dp)
    end function funC_0

    function funC_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funC(b, 1.0_dp)
    end function funC_1

    function funC_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funC(b, 2.0_dp)
    end function funC_2

    function minFunOverB_2(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funOverB(b, 2.0_dp)
    end function minFunOverB_2

    function funOverB_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = funOverB(b, 1.0_dp)
    end function funOverB_1

    function minFunOverB_1(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        ret = -funOverB(b, 1.0_dp)
    end function minFunOverB_1

    function funOverB_0(b) result(ret)
        real(dp), intent(in) :: b
        real(dp) :: ret
        !ASSERT( b > 0.0_dp)
        ret = funOverB(b, 0.0_dp)
    end function funOverB_0

    ! ==== end of specific single particle matrix element terms =============

    ! =========== generic necessary functions ===============================
    ! use of ASSERT() probably not possible here, as it causes side-effects!
    function funA(b, x, y) result(ret)
        real(dp), intent(in) :: b, x, y
        real(dp) :: ret
        !ASSERT( (b + y) >= 0.0_dp)
        ret = sqrt((b + x)/(b + y))
    end function funA

    function funB(b, x, y) result(ret)
        real(dp), intent(in) :: b, x, y
        real(dp) :: ret
        !ASSERT( (b + x) > 0.0_dp)
        !ASSERT( (b + y) > 0.0_dp)
        ret = sqrt(2.0_dp/((b + x)*(b + y)))
    end function funB

    function funC(b, x) result(ret)
        real(dp), intent(in) :: b, x
        real(dp) :: ret
        !ASSERT( (b + x) > 0.0_dp)
        ret = funA(b, x - 1.0_dp, x) * funA(b, x + 1.0_dp,x)
    end function funC

    function funD(b, x) result(ret)
        real(dp), intent(in) :: b, x
        real(dp) :: ret
        !ASSERT( (b + x) > 0.0_dp)
        ret = funA(b, x + 2.0_dp, x) * funA(b, x - 1.0_dp, x + 1.0_dp)
    end function funD

    function funOverB(b, x) result(ret)
        real(dp), intent(in) :: b, x
        real(dp) :: ret
        !ASSERT( (b + x) > 0.0_dp)
        ret = 1.0_dp/(b + x)
    end function funOverB

    ! =========== end of generic necessary functions ========================






end module
