#include "macros.h"
module gen_coul_mod
    use SystemData, only: BasisFN
    use UMatCache, only: GTID, UMatInd
    use global_utilities, only: timer, set_timer, halt_timer
    use util_mod, only: near_zero
    use constants, only: dp, stdout
    use fcoul_mod, only: SlatCoulFou
    better_implicit_none
    private
    public :: gen_coul

contains
    SUBROUTINE GEN_COUL(NBASISMAX, NHG, G1, NMSH, NMAX, FCK, UMAT, ZIA)
        INTEGER :: NMSH, NMAX, NHG, I, J, K, L, ID1, ID2, ID3, ID4, ISPINSKIP, II
        TYPE(BASISFN) G1(NHG)
        INTEGER nBasisMax(5, *)
        type(timer), save :: proc_timer
        real(dp) :: UMAT(*)
        !..Cube arrays
        complex(dp) FCK(NMSH, NMSH, NMSH)
        complex(dp) ZIA(-NMSH / 2:NMSH / 2, NMAX, NMAX)
        real(dp) :: SUM
        !..This routine generates *ALL* possible combinations of Coulomb integrals
        !..stored in the form (u1 u2 | U | u1' u2') = UMAT(n1 n2 n3 n4)
        !..
        !..This first call calculates the inner integral
        !..The call to SCOUL calculates the outer integral
        ! ==--------------------------------------------------------------------==
        proc_timer%timer_name = 'GEN_COUL  '
        call set_timer(proc_timer)
        ! ==--------------------------------------------------------------------==
        OPEN(10, FILE='UMAT', STATUS='UNKNOWN')
        II = 0
        ISPINSKIP = NBASISMAX(2, 3)
        DO I = 1, NHG, ISPINSKIP
            DO J = 1, NHG, ISPINSKIP
                DO K = 1, NHG, ISPINSKIP
                    DO L = 1, NHG, ISPINSKIP
                        SUM = 0.0_dp
                        !..Original call
                        !              CALL SLATCOULFOU(G1(1,I),G1(1,J),
                        !     &        G1(1,K),G1(1,L),NMSH,FCK,NMAX,SUM)
                        ID1 = GTID(I)
                        ID2 = GTID(J)
                        ID3 = GTID(K)
                        ID4 = GTID(L)
                        IF (near_zero(UMAT(UMatInd(ID1, ID2, ID3, ID4)))) then
                            CALL SLATCOULFOU(G1(I)%k, G1(J)%k, G1(K)%k, G1(L)%k, NMSH, FCK, NMAX, ZIA, SUM)
                            !              SUM=SUM*COUPLE
                            !..Original call
                            UMAT(UMatInd(ID1, ID2, ID3, ID4)) = SUM
                            !              UMAT(ID1,ID2,ID3,ID4)=SUM
                            !..Symmetries
                            !              UMAT(ID3,ID4,ID1,ID2)=SUM
                            !              UMAT(ID2,ID1,ID4,ID3)=SUM
                            !              UMAT(ID4,ID3,ID2,ID1)=SUM
                            IF (SUM > 1.0e-10_dp) WRITE(10, '(4I7,F19.9)') ID1, ID2, ID3, ID4, SUM
                        end if
                    END DO
                END DO
            END DO
        END DO
        CLOSE(10)
        WRITE(stdout, *) ' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
        ! ==--------------------------------------------------------------------==
        call halt_timer(proc_timer)
        ! ==--------------------------------------------------------------------==
    END subroutine
    ! ==------------------------------------------------------------------==
end module gen_coul_mod
