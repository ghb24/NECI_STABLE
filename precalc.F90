SUBROUTINE GETVARS(NI,BETA,I_P,IPATH,I,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,              &
     &         FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,     &
     &         DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,TLOGP)

     USE HElement
     IMPLICIT NONE
     include 'vmc.inc'
     include 'basis.inc'
     TYPE(BasisFN) G1(*)
     INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,I,L,LT,Q
     INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX),NBASIS
     INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),GIDHO,n
     INTEGER iters,r,s
     COMPLEX*16 FCK(*)
     REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM
     REAL*8 ax,bx,cx,minvar,xmin,originalc,originalimport
     REAL*8 ZEROVAR,p(2),pl(2),ph(2),xi(2,2),originala,originalb,originald
     REAL*8 fret,MCPATHSPRE,fa,fb,fc
     REAL*8 originalpolyb1,originalpolyb2,originalpolya
     REAL*8 polyp(3),polyxi(3,3),bestvalues(2,PREIV_MAX),bestvaluespoly(3,PREIV_MAX)
     REAL*8 bestvalpolyboth(4,PREIV_MAX),bestxipolyboth(4,4),polypboth(4)
     REAL*8 bestxipoly(3,3),bestxi(2,2),polyxiboth(4,4)
     REAL*8 origpolyboth(4),ORBENERGY,ENERGYLIMS(2)
     LOGICAL TSYM,PREVAR,NOTHING,TLOGP
     TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
     TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
     TYPE(HElement) HIJS(0:PREIV_MAX)
     TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
     EXTERNAL MCPATHSPRE
     
     ENERGYLIMS(1)=ORBENERGY(NMAX,1)
     ENERGYLIMS(2)=ORBENERGY(NMAX,NBASIS)
     
!     do b=1,nbasis
!        energy=ORBENERGY(NMAX,b)
!        WRITE(6,*) "ENERGIES ARE", b, energy
!     enddo
!     CALL FLUSH(6)
     
     bestvalues=0.D0
     bestvaluespoly=0.D0
     bestvalpolyboth=0.D0
     
     !Only open PRECALC file if logging option on
     IF (TLOGP) THEN
        OPEN(31,FILE="PRECALC",STATUS="UNKNOWN")
        WRITE(31,"(A,G25.16,A,G25.16)") "Energy limits on sigma given by:",ENERGYLIMS(1),", and ",ENERGYLIMS(2)
     ENDIF
      
     !Loop over vertex levels to look at
     DO Q=I,PREIV_MAX
         
         NOTHING=.FALSE.
     
         ! IF NONE IS SPECIFIED
         IF (PRE_TAY(3,Q).eq.0) THEN
             
             GIDHO=6
             VARSUM=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
     &                  FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,  &
     &                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
             
             CYCLE
         ENDIF
     
         WRITE(6,*) ""
         WRITE(6,"(A,I2,A)") "For a vertex level of", Q, ", PRECALC finds:"
         
        !IF NOTHING SPECIFIED
        IF ((.not.PRE_TAYLOG(1,Q)).and.(.not.PRE_TAYLOG(2,Q)).and.(.not.PRE_TAYLOG(3,Q)).and.(.not.PRE_TAYLOG(4,Q)).and.  &
     &      (PRE_TAYREAL(1,Q).eq.0.D0).and.(PRE_TAY(3,Q).ne.0).and.(.not.PRE_TAYLOG(5,Q)).and.(.not.PRE_TAYLOG(6,Q))) THEN
            NOTHING=.TRUE.
        ENDIF
        
        !If POLYEXCITBOTH is specified in the input, this means that both the i and j orbitals, as well as the kl orbitals are using a polynomial expansion as the weighting function
        IF(TPOLYEXCITBOTH.and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "SigmaA, nA, sigmaB and nB POLYEXCITWEIGHTING: (Vertex level, Iteration number, Parameter values, Expected Variance)"
            ENDIF

            !Initial bracketing
            IF((bestvalpolyboth(1,(Q-1)).eq.0.D0).and.(bestvalpolyboth(2,(Q-1)).eq.0.D0).and.(bestvalpolyboth(3,(Q-1)).eq.0.D0).and.(bestvalpolyboth(4,(Q-1)).eq.0.D0)) THEN

                polypboth=(/ 0,1,0,1 /)     !Initial sigmaA, nA, sigmaB and nB.
                polyxiboth=RESHAPE( (/ 1.D0,((0.D0,r=1,4),1.D0,s=1,3) /), (/ 4, 4 /) )!Initial directions - unit vectors
            ELSE
                !Choose values which the previous vertex level found as optimum
                polypboth=bestvalpolyboth(:,(Q-1))
                polyxiboth=bestxipolyboth
            ENDIF
            n=4                 !Dimensions
            GIDHO=7
            
            origpolyboth(1)=g_VMC_PolyExcitFromWeight1
            origpolyboth(2)=g_VMC_PolyExcitFromWeight2
            origpolyboth(3)=g_VMC_PolyExcitToWeight1
            origpolyboth(4)=g_VMC_PolyExcitToWeight2

            CALL POWELL(polypboth,polyxiboth,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,    & 
     &              G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,                     &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,ENERGYLIMS)             
            
            bestvalpolyboth(:,Q)=polypboth(:)
            bestxipolyboth=polyxiboth
        
            IF (NOTHING) THEN

                g_VMC_PolyExcitFromWeight1=origpolyboth(1)
                g_VMC_PolyExcitFromWeight2=origpolyboth(2)
                g_VMC_PolyExcitToWeight1=origpolyboth(3)
                g_VMC_PolyExcitToWeight2=origpolyboth(4)
                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,F16.12,A,F6.3,A,F6.3,A,F6.3,A,F6.3)") "Optimum SigmaA POLYEXCITWEIGHT found to be ", polypboth(1), "nA POLYEXCITWEIGHT as ", polypboth(2),           &
     &           ", sigmaB POLYEXCITWEIGHT as ", polypboth(3), ", and nB POLYEXCITWEIGHT as ", polypboth(4),", but using values of ",     &
     &                    origpolyboth(1), ", ", origpolyboth(2), ", ", origpolyboth(3), ", and ", origpolyboth(4)

            ELSE
                g_VMC_PolyExcitFromWeight1=polypboth(1)
                g_VMC_PolyExcitFromWeight2=polypboth(2)
                g_VMC_PolyExcitToWeight1=polypboth(3)
                g_VMC_PolyExcitToWeight2=polypboth(4)
                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12)") "POLYEXCITWEIGHTING parameters optimised to SigmaA=", polypboth(1), ", nA=",polypboth(2),", sigmaB=",polypboth(3), ", and nB=",polypboth(4)
            ENDIF
        ENDIF
            
            !If PolyExcit is specified rather than Excitweighting, we want to search for the minimum in a 4D space, where now we have two different excitto weighting parameters
        IF(TPOLYEXCIT.and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "A, sigma and n POLYEXCITWEIGHTING: (Vertex level, Iteration number, Parameter values, Expected Variance)"
            ENDIF

            !Initial bracketing
            IF((bestvaluespoly(1,(Q-1)).eq.0.D0).and.(bestvaluespoly(2,(Q-1)).eq.0.D0).and.(bestvaluespoly(3,(Q-1)).eq.0.D0)) THEN
                polyp=(/ 0.1,0,1 /)     !Initial A, sigma and n.
                polyxi=RESHAPE( (/ 1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0 /), (/ 3, 3 /) )!Initial directions - unit vectors
            ELSE
                !Choose values which the previous vertex level found as optimum
!                polyp=(/ bestvaluespoly(1,(Q-1)),bestvaluespoly(2,(Q-1)),bestvaluespoly(3,(Q-1)) /)
                polyp=bestvaluespoly(:,(Q-1)) 
                polyxi=bestxipoly
            ENDIF
            n=3                 !Dimensions
            GIDHO=4

            originalpolya=g_VMC_ExcitFromWeight
            originalpolyb1=g_VMC_PolyExcitToWeight1
            originalpolyb2=g_VMC_PolyExcitToWeight2
        
            CALL POWELL(polyp,polyxi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
     &              G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,              &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,ENERGYLIMS)
            
            bestvaluespoly(:,Q)=polyp(:)
            bestxipoly=polyxi
            
            IF (NOTHING) THEN

                g_VMC_ExcitFromWeight=originalpolya
                g_VMC_PolyExcitToWeight1=originalpolyb1
                g_VMC_PolyExcitToWeight2=originalpolyb2
                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,F6.3,A,F6.3,A,F6.3)") "Optimum A POLYEXCITWEIGHT found to be ", polyp(1), ", sigma POLYEXCITWEIGHT as ", polyp(2), ", and n POLYEXCITWEIGHT as ", polyp(3),", but using values of ",     &
     &           originalpolya, ", ", originalpolyb1, ", and ", originalpolyb2

            ELSE
                g_VMC_ExcitFromWeight=polyp(1)
                g_VMC_PolyExcitToWeight1=polyp(2)
                g_VMC_PolyExcitToWeight2=polyp(3)
                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12)") "POLYEXCITWEIGHTING parameters optimised to A=", polyp(1), ", sigma=",polyp(2), ", and n=",polyp(3)

            ENDIF
        ENDIF
            
        !Looking for a & b parameters
        IF ((.not.TPOLYEXCIT).and.(.not.TPOLYEXCITBOTH).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "A and B EXCITWEIGHTING: (Vertex level, Iteration number, Parameter Values, Expected Variance)"
            ENDIF
            
            !Initial bracketing
            IF((bestvalues(1,(Q-1)).eq.0.D0).and.(bestvalues(2,(Q-1)).eq.0.D0)) THEN
                p=(/ 0.1,0.5 /)     !Initial a and b values
                xi= RESHAPE( (/ 1.D0, 0.D0, 0.D0, 1.D0 /), (/ 2, 2 /) ) !Initial directions
            ELSE
                !Choose values which the previous vertex level found as optimum
!                p=(/ bestvalues(1,(Q-1)),bestvalues(2,(Q-1)) /)
                p=bestvalues(:,(Q-1))
                xi=bestxi
            ENDIF
            n=2                 !Dimensions
            GIDHO=1             !To tell brentalgo that we are looking at a & b parameters

            originala=g_VMC_ExcitFromWeight
            originalb=g_VMC_ExcitToWeight
            
            CALL POWELL(p,xi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,          &
     &              G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,           &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,ENERGYLIMS)

            bestvalues(:,Q)=p(:)
            bestxi=xi

            IF (NOTHING) THEN
                
                g_VMC_ExcitFromWeight=originala
                g_VMC_ExcitToWeight=originalb
                WRITE(6,"(A,F16.12,A,F16.12,A,F6.3,A,F6.3)") "Optimum A EXCITWEIGHT found to be ", p(1), ", and B EXCITWEIGHT as ", p(2), ", but using values of ", originala, ", and ", originalb
                
            ELSE
                g_VMC_ExcitFromWeight=p(1)
                g_VMC_ExcitToWeight=p(2)
                WRITE(6,"(A,F16.12,A,F16.12)") "EXCITWEIGHTING parameters optimised to A=", p(1), " and B=",p(2)
            
            ENDIF
        ENDIF
            
        !Looking for importance parameter
        IF (PRE_TAYLOG(3,Q).or.PRE_TAYLOG(4,Q)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "Importance Parameter: (Vertex level, Iteration number, Importance Parameter, Expected Variance)"
            ENDIF
            
            !Initial bracketing guess
            ax=0.7
            bx=0.9
!           cx=0.99 

            GIDHO=2
            originalimport=G_VMC_PI
            
            !Ensuring correct bracketing
            CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,      &
                    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                          &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)
            
            IF (TLOGP) THEN
                WRITE(31,"(A,F16.12,A,F16.12)") "From mnbrak routine, minimum is between ", ax, " and ", bx
            ENDIF

                    
            CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
                       G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,                  &
                       TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,fb,ENERGYLIMS)

            IF (PRE_TAYLOG(3,Q)) THEN
                G_VMC_PI=originalimport
                WRITE(6,"(A,F15.12,A,F5.3)") "Optimum importance parameter found to be ", xmin, ", but using ", originalimport
            ELSEIF (PRE_TAYLOG(4,Q)) THEN
                G_VMC_PI=xmin
                WRITE(6,"(A,F15.12)") "Importance parameter optimised to", xmin

            ENDIF
        ENDIF
        
        !Looking for D Excitweighting parameter
        IF (PRE_TAYLOG(5,Q).or.PRE_TAYLOG(6,Q)) THEN

            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "D Parameter: (Vertex level, Iteration number, D Weighting parameter, Expected Variance)"
            ENDIF
            
            !Initial bracketing guess
            ax=-0.8
            bx=-1.0
            
            GIDHO=5
            originald=g_VMC_ExcitToWeight2

            CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,  &
                    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                      &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)
            IF (TLOGP) THEN
                WRITE(31,"(A,F16.12,A,F16.12)") "From mnbrak routine, minimum is between ", ax, " and ", cx
            ENDIF

            CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,                          &
     &                  RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &                  NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,fb,ENERGYLIMS)
            
            IF (PRE_TAYLOG(5,Q)) THEN
                g_VMC_ExcitToWeight2=originald
                WRITE(6,"(A,F15.12,A,F5.3)") "Optimum D parameter found to be ", xmin, ", but using ", originald
            ELSEIF (PRE_TAYLOG(6,Q)) THEN
                g_VMC_ExcitToWeight2=xmin
                WRITE(6,"(A,F15.12)") "D parameter optimised to", xmin
            ENDIF
        ENDIF
        
        !Looking for C Excitweighting parameter
        IF (PRE_TAYLOG(2,Q).or.(PRE_TAYREAL(1,Q).gt.0).or.NOTHING) THEN
            
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "U Parameter: (Vertex level, Iteration number, C Weighting parameter, Expected Variance)"
            ENDIF
            
         !Initial bracketing
            ax=20
            bx=30
!            cx=90

            GIDHO=3
            originalc=G_VMC_EXCITWEIGHT
            VARSUM=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
     &          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
            ZEROVAR=VARSUM
            
            IF (TLOGP) THEN
                WRITE(31,"(2I3,2G25.16)")  Q, 0, 0.D0, VARSUM
                CALL FLUSH(31)
            ENDIF
            
            !To ensure correct bracketing
            CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,  &
                    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                      &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)
            
            IF (TLOGP) THEN
                WRITE(31,"(A,F16.12,A,F13.9)") "From mnbrak routine, minimum is between ", ax, " and ", cx
            ENDIF
                    
            CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,                          &
     &                  RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &                  NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,fb,ENERGYLIMS)
           
            IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).gt.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN
                G_VMC_EXCITWEIGHT=xmin
                WRITE(6,"(A,F17.12)") "Optimum U weighting outside UEPSILON bounds, so using C=",xmin
            ELSE IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).le.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN
                G_VMC_EXCITWEIGHT=0.D0
                WRITE(6,"(A)") "Optimum U weighting within UEPSILON bounds, so using C=0"
            ELSE
                G_VMC_EXCITWEIGHT=originalc
                WRITE(6,"(A,G20.12,A,F14.8)") "Optimum U weighting found to be", xmin, ", but using", originalc
            END IF
        END IF
    ENDDO
    WRITE(6,*) ""
    
    IF (TLOGP) THEN
        CLOSE(31)
    ENDIF
    
    RETURN
END SUBROUTINE GETVARS

FUNCTION VARIANCEAB(pointab,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
           FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,            &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)

    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    REAL*8 pointab(*),VARIANCEAB,MCPATHSPRE
    EXTERNAL MCPATHSPRE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,L,LT,K,D
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX),GIDHO
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,G
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    
    G=.false.
    
    SELECT CASE (GIDHO)
    CASE(1)
        g_VMC_ExcitFromWeight=pointab(1)
        g_VMC_ExcitToWeight=pointab(2)
    CASE(4)
        g_VMC_ExcitFromWeight=pointab(1)
        g_VMC_PolyExcitToWeight1=pointab(2)
        g_VMC_PolyExcitToWeight2=pointab(3)
        IF ((pointab(2).gt.ENERGYLIMS(2)).or.(pointab(2).lt.ENERGYLIMS(1))) G=.true.
    CASE(7)
        g_VMC_PolyExcitFromWeight1=pointab(1)
        g_VMC_PolyExcitFromWeight2=pointab(2)
        g_VMC_PolyExcitToWeight1=pointab(3)
        g_VMC_PolyExcitToWeight2=pointab(4)
        IF ((pointab(1).lt.ENERGYLIMS(1)).or.(pointab(1).gt.ENERGYLIMS(2))) G=.true.
        IF ((pointab(3).lt.ENERGYLIMS(1)).or.(pointab(3).gt.ENERGYLIMS(2))) G=.true.
    ENDSELECT
    
    IF (G) THEN
        VARIANCEAB=HUGE(VARIANCEAB)
    ELSE
    
    VARIANCEAB=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,    &
                    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,  &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)

    ENDIF
                    
END FUNCTION VARIANCEAB
           

FUNCTION MCPATHSPRE(point,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,       &
     &        FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,        &
     &        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)

    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,L,LT,K
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,GIDHO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR
    REAL*8 NTOTAL,BETA,ECORE,ALAT(3),RHOEPS,DBETA,VARSUM,point,MCPATHSPRE
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) DLWDB,TOTAL,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    ! Saved values only go up to a vertex level of 6 - increase values if we want to go to higher vertex levels
    REAL*8 NTOTSAV(1:7),ENERGYLIMS(2)
    TYPE(HDElement) DLWSAV(1:7),TOTSAV(1:7)
    DATA NTOTSAV/7*1.D0/
    DATA TOTSAV/7*HDElement(1.D0)/
    DATA DLWSAV/7*HDElement(0.D0)/
    SAVE NTOTSAV,TOTSAV,DLWSAV
    
    SELECT CASE (GIDHO)
    !Importance
    CASE(2)
        G_VMC_PI=point
    !U weighting
    CASE(3)
        G_VMC_EXCITWEIGHT=point
    !D Weighting
    CASE(5)
        g_VMC_ExcitToWeight2=point
    !NONE specified
    CASE(6)
        G_VMC_EXCITWEIGHT=point
    ENDSELECT
    
    IF(K.eq.2) THEN
        DLWSAV(1)=HIJS(0)
    ENDIF
        
!    NTOTAL=1.D0
!    TOTAL=1.D0
!    DLWDB=HIJS(0)
!OPEN(PRECALC)
!    do D=2,K
        L=0
        LT=0
        VARSUM=0.D0
        DLWDB2=0.D0
        METH=pre_TAY(1,K)   !Method for vertex level
        CYCLES=pre_TAY(2,K)
        IF(METH.eq.-8) then !Full Rho-diag method
!   This code generates excitations on the fly, and gets
!            F(K)=FMCPR3B(NI,BETA,I_P,IPATH,K,NEL,                                &
!     &          NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,       &
!     &          RHOEPS,0,RHOII,RHOIJ,CNWHTAY,I_CHMAX,LOCTAB,                     &
!     &          ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,            &
!     &          MP2E,NTOTAL)
        ELSEIF(METH.eq.-20) then !Full H-diag method
            EREF=DLWSAV(K-1)/TOTSAV(K-1)

            !EREF=DLWDB/TOTAL

            F(K)=FMCPR3B2(NI,BETA,I_P,IPATH,K,NEL,                               &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,         &
     &         RHOEPS,0,RHOIJ,0,METH,LOCTAB,                                     &
     &         0,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,                    &
     &         MP2E,NTOTSAV(K-1),PREIV_MAX,EREF,VARSUM,PREVAR,TOTSAV(K-1))

!        ELSEIF(METH.eq.-7.or.METH.eq.-19) THEN

            ENDIF
            TOTSAV(K)=TOTSAV(K-1)+F(K)
            DLWSAV(K)=DLWSAV(K-1)+DLWDB2
!            WRITE(6,*) K,TOTSAV(K),DLWSAV(K)
!        ENDDO
        MCPATHSPRE=VARSUM

END FUNCTION MCPATHSPRE

SUBROUTINE BRENTALGO(brent,ax,bx,cx,fun,tol,xmin,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,    &
     &                 G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,   &
     &                 RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,       &
     &                 NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,INITFUNC,ENERGYLIMS)
    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,K,L,LT
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,GIDHO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM,fun,INITFUNC
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER ITMAX
    EXTERNAL fun
    REAL*8 brent,ax,bx,cx,tol,xmin,CGOLD,ZEPS,ENERGYLIMS(2)
    PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
        !  Given a function f, and given a bracketing triplet of abscissas ax,bx,cx (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolated the minimum to a fractional precision of about tol using Brents method.
        !  The abscissa od the minimum is returned as xmin, and the minimum function value is returned as brent, the returned function value.
        !  Parameters: Maximum allowed number of iterations; golden ratio; and a small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero.
    INTEGER iter
    REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(ax,cx)        !a and b must be in ascending order,
                        !though the input abscissas need not be.
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.D0

    write(31,*) "BRENT ALGO STARTING"
    IF (INITFUNC.eq.0.D0) THEN
        WRITE(31,*) "INITFUNC EQUAL 0.D0 - redo initial point"
        VARSUM=fun(x,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                &
              FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
    ELSE
        VARSUM=INITFUNC
    ENDIF
    
    ! Vertex level, iteration number, C weight value, Varsumnu value
    IF(INITFUNC.eq.0.D0) THEN
        IF (TLOGP.and.(GIDHO.ne.1).and.(GIDHO.ne.4).and.(GIDHO.ne.7)) THEN
            WRITE(31,"(2I3,2G25.16)")  K, 0, x, VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.AND.(GIDHO.EQ.1)) THEN
            WRITE(31,"(I3,A,3G25.16)")  K, " brent  0", g_VMC_ExcitFromWeight, g_VMC_ExcitToWeight, VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.and.(GIDHO.eq.4)) THEN
            WRITE(31,"(I3,A,4G25.16)")  K, " brent  0", g_VMC_ExcitFromWeight, g_VMC_PolyExcitToWeight1, g_VMC_PolyExcitToWeight2, VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.and.(GIDHO.eq.7)) THEN
            WRITE(31,"(I3,A,5G25.16)")  K, " brent  0", g_VMC_PolyExcitFromWeight1,g_VMC_PolyExcitFromWeight2,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2,VARSUM
            CALL FLUSH(31)
        ENDIF
    ENDIF
        
    fx=VARSUM
    fv=fx
    fw=fx
    do iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.D0*tol1
        if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.D0*(q-r)
            if(q.gt.0.D0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.             &
              p.ge.q*(b-x)) goto 1  !The above conditions determine the acceptability of the parabolic fit
            d=p/q                    ! Take parabolic step
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2                   ! Skip over golden section step
        end if
     1  if(x.ge.xm) then        ! Golden section step, which we take into
            e=a-x               ! the larger of the two segments
        else
            e=b-x
        endif
        d=CGOLD*e               ! Take the golden section step
     2  if(abs(d).ge.tol1) then ! Arrive here with d computed either from 
            u=x+d               ! parabolic fit, or golden section
        else
            u=x+sign(tol1,d)
        endif
           
        VARSUM=fun(u,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                    &
                  FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
             IF (TLOGP.and.(GIDHO.eq.7)) THEN
                 WRITE(31,"(I3,A,I3,5G25.16)")  K, " brent", iter, g_VMC_PolyExcitFromWeight1,g_VMC_PolyExcitFromWeight2, g_VMC_PolyExcitToWeight1, g_VMC_PolyExcitToWeight2, VARSUM
                 CALL FLUSH(31)
             ENDIF
             IF (TLOGP.and.(GIDHO.eq.4)) THEN
                 WRITE(31,"(I3,A,I3,4G25.16)")  K, " brent", iter, g_VMC_ExcitFromWeight, g_VMC_PolyExcitToWeight1, g_VMC_PolyExcitToWeight2, VARSUM
                 CALL FLUSH(31)
             ENDIF
             IF (TLOGP.and.(GIDHO.eq.1)) THEN
                 WRITE(31,"(I3,A,I3,3G25.16)")  K, " brent", iter, g_VMC_ExcitFromWeight, g_VMC_ExcitToWeight, VARSUM
                 CALL FLUSH(31)
             ENDIF                                  
             IF (TLOGP.and.(GIDHO.ne.1).and.(GIDHO.ne.4).and.(GIDHO.ne.7)) THEN
                 WRITE(31,"(2I3,2G25.16)")  K, iter, u, VARSUM
                 CALL FLUSH(31)
             ENDIF
        
        fu=VARSUM                 ! This is the one function evaluation per iteration
        if(fu.le.fx) then       ! What to do with function
            if(u.ge.x) then
                a=x
            else
                b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if(u.lt.x) then
                a=u
            else
                b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                v=u
                fv=fu
            endif
        endif
        enddo
        write (6,*) "***Brent algorithm exceeds maximum iterations of ", ITMAX ," Correct values may not have been found"
        write (31,*) "***Brent algorithm exceeds maximum iterations of ", ITMAX ," Correct values may not have been found"
     3  xmin=x
        brent=fx
        write(31,*) "BRENT ALGO FINISHED"
    return
END SUBROUTINE BRENTALGO

SUBROUTINE POWELL(p,xi,n,np,ftol,iter,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,                   &
     &        G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,        &
     &        TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,ENERGYLIMS)

    
    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER iter,n,np,NMAX,ITMAX,NEL,I_P,BRR(*),NMSH(*),NTAY,L,LT,NMAXI
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX),Q
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,GIDHO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,TLOGP
    REAL*8 fret,ftol,p(np),xi(np,np),varianceab,BETA,ECORE,NTOTAL,ALAT(3),RHOEPS
    REAL*8 DBETA,VARSUM,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    EXTERNAL varianceab 
    PARAMETER (NMAXI=4,ITMAX=200)
!    USES func,linmin
!       Minimisation of a function 'func' of n variables. (func is not an argument, it is a fixed function name.)  
!Input consists of an initial starting point p(1:n); an initial matrix xi(1:n,1:n) with physical dimensions np by np, and whose columns contain the initial set of directions (usually the n unit vectors). 
!'ftol' is the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness.
!On output, p is set to the best point found, xi is the then current direction set, fret in the returned function value at p, and iter in the number of iterations taken.  The routine linmin is used.
!       Parameters:  Maximum expected value of n, and maximum allowed iterations.
    INTEGER i,ibig,j
    REAL*8 del,fp,fptt,t,pt(NMAXI),ptt(NMAXI),xit(NMAXI)
    fret=varianceab(p,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
           FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE, &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
    IF (TLOGP.and.(GIDHO.eq.7)) THEN
        WRITE(31,"(2I3,5G25.16)") Q, 0, p(1), p(2), p(3), p(4), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.(GIDHO.eq.4)) THEN
        WRITE(31,"(2I3,4G25.16)") Q, 0, p(1), p(2), p(3), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.(GIDHO.eq.1)) THEN
        WRITE(31,"(2I3,3G25.16)") Q, 0, p(1), p(2), fret
        CALL FLUSH(31)
    ENDIF
    do j=1,n                    !Save the initial point
        pt(j)=p(j)
    enddo
    iter=0
1   iter=iter+1
    fp=fret
    ibig=0
    del=0.D0                    !Will be the biggest function decrease
    do i=1,n                    !In each iteration, loop over all directions in the set.
        do j=1,n                !Copy the direction...
            xit(j)=xi(j,i)
        enddo
        fptt=fret
        call linmin(p,xit,n,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,          &
     &    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                    &
     &    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)   !Minimise along it...
        IF (TLOGP.and.(GIDHO.eq.1)) THEN
            WRITE(31,"(2I3,3G25.16)") Q, iter, p(1), p(2), fret
            CALL FLUSH(31)
        ELSEIF (TLOGP.and.(GIDHO.eq.7)) THEN
            WRITE(31,"(2I3,5G25.16)") Q, iter, p(1), p(2), p(3), p(4), fret
            CALL FLUSH(31)
        ELSEIF (TLOGP.and.(GIDHO.eq.4)) THEN
            WRITE(31,"(2I3,4G25.16)") Q, iter, p(1), p(2), p(3), fret
            CALL FLUSH(31)
        ENDIF
        if(abs(fptt-fret).gt.del) then      !And record it if it is the largest decrease so far.
            del=abs(fptt-fret)
            ibig=i
        endif
    enddo
    if(2.D0*abs(fp-fret).le.ftol*(abs(fp)+abs(fret))) return    !Termination criterion
    if(iter.eq.ITMAX) then
        write (6,*) "***Powell algorithm exceeds maximum iterations of ",ITMAX ," Correct values may not have been found"
        write (31,*) "***Powell algorithm exceeds maximum iterations of ",ITMAX ," Correct values may not have been found"
    endif
    do j=1,n                    !construct the extrapolated point and the average direction moved
        ptt(j)=2.D0*p(j)-pt(j)  !Save the old starting point.
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
    enddo
        
    fptt=varianceab(ptt,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
           FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,           &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
    !Function value at extrapolated point
    if (fptt.ge.fp) goto 1      !One reason not to use new direction
    t=2.D0*(fp-2.D0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
    if (t.ge.0.D0) goto 1       !Other reason not to use new direction
    CALL linmin(p,xit,n,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,          &
     &    FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                &
     &    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS) !Move to the minimum of the new direction,
    IF (TLOGP.and.(GIDHO.eq.1)) THEN
        WRITE(31,"(2I3,3G25.16)") Q, iter, p(1), p(2), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.(GIDHO.eq.7)) THEN
        WRITE(31,"(2I3,5G25.16)") Q, iter, p(1), p(2), p(3), p(4), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.(GIDHO.eq.4)) THEN
        WRITE(31,"(2I3,4G25.16)") Q, iter, p(1), p(2), p(3), fret
        CALL FLUSH(31)
    ENDIF
    do j=1,n                    !and save the new direction
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
    enddo
    goto 1                      !Back for another iteration
END SUBROUTINE POWELL

SUBROUTINE linmin(p,xi,n,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,         &
     &       FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,             &
     &       DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)
    
    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),NMSH(*),NMAX,NTAY,Q,L,LT,GIDHO
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER n,NMAXI
    REAL*8 fret,p(n),xi(n),TOL,pl(n),ph(n)
    PARAMETER (NMAXI=4, TOL=1.D-03)      !Maximum anticipated n
    !USES BRENTALGO,F1DIM
    !Given an n-dimensional point p(1:n) and an n-dimensional direction xi(1:n), moves and resets p to where the function func(p) takes on a minimum along the direction xi from p, and replaces xi by the actual vector displacement that p has moved.
    !Also returns as fret the value of func at the returned location p.  This is actually all accomplished by calling the routine BRENTALGO (and mnbrack?)
    INTEGER j,ncom
    REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NMAXI),xicom(NMAXI),brent
    COMMON /f1com/ pcom,xicom,ncom
    EXTERNAL f1dim
    ncom=n                  !Set up the common block
    do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
    enddo
    ax=0.D0
    xx=5.D-01
    
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
     &          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,              &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)

    call BRENTALGO(fret,ax,xx,bx,f1dim,TOL,xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &       G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,              &
     &       RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &       NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP,fx,ENERGYLIMS)
     do j=1,n                !Construct the vector results to return
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
    enddo
    return
    END SUBROUTINE linmin

FUNCTION f1dim(x,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
     &          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)      

    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,L,LT,K,D,Q
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,GIDHO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,varianceab
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER NMAXI
    REAL*8 f1dim,func,x,ENERGYLIMS(2)
    PARAMETER (NMAXI=4)
    EXTERNAL varianceab
    !USES func
!Used by linmin as the function passed to mnbrak and brentalgo
    INTEGER j,ncom
    REAL*8 pcom(NMAXI),xicom(NMAXI),xt(NMAXI)
    COMMON /f1com/ pcom,xicom,ncom
    do j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
    enddo
    f1dim=varianceab(xt,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
              FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,        &
              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
    return
END FUNCTION f1dim
   
SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
               FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,TLOGP,ENERGYLIMS)
    
    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,MINI
    PARAMETER (GOLD=1.618034, GLIMIT=100.D0,MINI=1.D-20)
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,L,LT,Q
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX),GIDHO
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,t
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,savedax,savedbx
    REAL*8 ENERGYLIMS(2)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    EXTERNAL func
    !Given a function 'func', and given distinct initial points ax and bx, this routine searches in the downhill direction (defined by the function as evaluated at the initial points) and returns new points ax,bx,cx that bracket a minimum of the function.
    !Also returned are the function values at the three points, fa, fb, fc.
    !Parameters: GOLD in the default ratio by which successive intervals are magnified; GLIMIT is the maximum magnification allowed for a parabolic-fit step.
    REAL*8 dum,fu,qu,r,u,ulim
!    WRITE(31,*) "MNBRAK ALGO STARTING"
2   savedax=ax
    savedbx=bx
    fa= func(ax,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
               FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
!    IF (TLOGP.and.(GIDHO.eq.4)) THEN
!        WRITE(31,"(2I3,4G25.16)") Q, 0, g_VMC_ExcitFromWeight,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2, fa
!        CALL FLUSH(31)
!    ENDIF
    fb= func(bx,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
               FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
!    IF (TLOGP.and.(GIDHO.eq.4)) THEN
!        WRITE(31,"(2I3,4G25.16)") Q, 0, g_VMC_ExcitFromWeight,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2, fb
!        CALL FLUSH(31)
!    ENDIF
    IF(fb.gt.fa) THEN      !Switch roles of a & b so that we can go downhill in direction from a to b
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
    ENDIF
    cx=bx+GOLD*(bx-ax)     !First guess for c
    fc=func(cx,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,            &
          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
    IF (TLOGP.and.(GIDHO.eq.4)) THEN
!        WRITE(31,"(2I3,4G25.16)") Q, 0, g_VMC_ExcitFromWeight,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2, fc
!        CALL FLUSH(31)
!    ENDIF
1      IF (fb.ge.fc) THEN      !"do while": keep returning here until we bracket.
            IF(fb.eq.fc) THEN
                t=t+1
                IF(t.gt.3) THEN
                    WRITE(31,*) "***mnbrak routine stuck - restarting with different values***"
                    ax=savedax-1.D0
                    bx=savedbx-1.D0
                    GOTO 2
                ENDIF
            ENDIF
            IF(fb.ne.fc) t=0
            r=(bx-ax)*(fb-fc)  !Compute u by parabolic extrapolation extrapolation from a,b,c
            qu=(bx-cx)*(fb-fa)  !TINY is used to prevent possible division by zero
            u=bx-((bx-cx)*qu-(bx-ax)*r)/(2.D0*sign(max(abs(qu-r),MINI),qu-r))
            ulim=bx+GLIMIT*(cx-bx)   !We won't go father than this.  Test various possibilities:
            IF((bx-u)*(u-cx).gt.0.D0) THEN  !Parabolic u is between b and c: try it
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
     
                IF (fu.lt.fc) THEN      !Got a minimum between b and c
                    ax=bx
                    fa=fb
                    bx=u
                    fb=fu
                    return
                ELSEIF (fu.gt.fb) THEN  !Got a minimum between a and u
                    cx=u
                    fc=fu
                    return
                ENDIF
                u=cx+GOLD*(cx-bx)       !Parabolic fit was no use. Use default magnification
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
!                IF (TLOGP.and.(GIDHO.eq.4)) THEN
!                    WRITE(31,"(2I3,4G25.16)") Q, 0, g_VMC_ExcitFromWeight,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2, fu
!                    CALL FLUSH(31)
!                ENDIF
            ELSEIF((cx-u)*(u-ulim).gt.0.D0) THEN    !Parabolic fit is between c and its allowed limit
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
                IF(fu.lt.fc) THEN
                    bx=cx
                    cx=u
                    u=cx+GOLD*(cx-bx)
                    fb=fc
                    fc=fu
                    fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
                ENDIF
            ELSEIF((u-ulim)*(ulim-cx).ge.0.D0) THEN     !Limit parabolic u to maximum allowed value
                u=ulim
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
            ELSE        !Reject parabolic u, use default magnification
                u=cx+GOLD*(cx-bx)
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,PREVAR,GIDHO,ENERGYLIMS)
!                IF (TLOGP.and.(GIDHO.eq.4)) THEN
!                    WRITE(31,"(2I3,4G25.16)") Q, 0, g_VMC_ExcitFromWeight,g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2, fu
!                    CALL FLUSH(31)
!                ENDIF
            ENDIF
            ax=bx           !Eliminate oldest point and continue
            bx=cx
            cx=u
            fa=fb
            fb=fc
            fc=fu
            goto 1
        ENDIF
    RETURN
END SUBROUTINE mnbrak
