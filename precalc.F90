module PreCalc

    USE Calc , only : G_VMC_PI,G_VMC_EXCITWEIGHT,G_VMC_EXCITWEIGHTS,G_VMC_SEED,     &
     &      CUR_VERT,EXCITFUNCS
    USE PRECALCREAD , only : PREIV_MAX,TOTALERROR,PRE_TAYREAL,MEMSAV,PRE_TAYLOG,PRE_TAY, &
     &      USEVAR,TRUECYCLES,TGRIDVAR,GRIDVARPAR,TLINEVAR,LINEVARPAR
    
    Use Integrals, only : ChemPot
    REAL*8, POINTER, DIMENSION(:) :: PGENLIST
    INTEGER, POINTER, DIMENSION(:) :: NMEM
    REAL*8, POINTER, DIMENSION(:,:) :: GRAPHPARAMS
    INTEGER, POINTER, DIMENSION(:,:,:) :: GRAPHS
#if defined(POINTER8)
    INTEGER*8, POINTER, DIMENSION(:,:) :: PVERTMEMS
#else
    INTEGER, POINTER, DIMENSION(:,:) :: PVERTMEMS
#endif
    SAVE PGENLIST,NMEM,GRAPHPARAMS,GRAPHS,PVERTMEMS
   
    contains

SUBROUTINE GETVARS(NI,BETA,I_P,IPATH,I,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,              &
     &         FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
     &         DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,TLOGP,KSYM,NWHTAY,I_VMAX)

     USE HElem
     use System, only: BasisFN
     IMPLICIT NONE
     TYPE(BasisFN) G1(*)
     INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),I,L,LT,Q
     INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX),NBASIS
     INTEGER IPATH(NEL,0:PREIV_MAX),GIDHO,n,gg,zz
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
     INTEGER iters,r,s,rr,KSYM(5),UNITNO,vv,kk,cc,I_VMAX,NWHTAY(3,I_VMAX)
     COMPLEX*16 FCK(*)
     REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM
     REAL*8 bestvals(6,PREIV_MAX)!,originalvals(6)
     REAL*8 ax,bx,cx,minvar,xmin,originalc,originalimport
     REAL*8 ZEROVAR,p(2),pl(2),ph(2),xi(2,2)
     REAL*8 fret,fa,fb,fc!,distance
     REAL*8 polyp(3),polyxi(3,3),SUMSD
     REAL*8 bestxipolyboth(4,4),polypboth(4)
     REAL*8 bestxipoly(3,3),bestxi(2,2),polyxiboth(4,4)
     REAL*8 ORBENERGY,ENERGYLIMS(2),xxx,g_VMC_FINAL(6,2:10)
     REAL*8 G_VMC_EXCITFINAL(2:10),VARIANCES(2:PREIV_MAX)
     LOGICAL TSYM,NOTHING,TLOGP,DEALLOC,LOWNEL,check
     TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
     TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
     TYPE(HElement) HIJS(0:PREIV_MAX)
     TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
     CHARACTER(len=12) :: abstr
!     REAL*8 MCPATHSPRE
!     EXTERNAL MCPATHSPRE
     
     DEALLOC=.FALSE.
     LOWNEL=.FALSE.
     
     IF (EXCITFUNCS(2).or.EXCITFUNCS(3)) THEN
        ENERGYLIMS(1)=ORBENERGY(NMAX,1)
        ENERGYLIMS(2)=ORBENERGY(NMAX,NBASIS)
     ENDIF
     
     IF(NEL.LE.2) LOWNEL=.TRUE.
     
!     do b=1,nbasis
!        energy=ORBENERGY(NMAX,b)
!        WRITE(6,*) "ENERGIES ARE", b, energy
!     enddo
!     CALL FLUSH(6)
     bestvals=0.D0
     
     !Only open PRECALC file if logging option on
     IF (TLOGP) THEN
        OPEN(31,FILE="PRECALC",STATUS="UNKNOWN")
        IF (EXCITFUNCS(2).or.EXCITFUNCS(3)) THEN
            WRITE(31,"(A,G25.16,A,G25.16)") "Energy limits on sigma given by:",ENERGYLIMS(1),", and ",ENERGYLIMS(2)
        ENDIF
        IF(LOWNEL) THEN
            WRITE(31,"(A)") "Low number of active electrons: All parameters may not be able to be optimised"
        ENDIF
     ENDIF
      
     !Loop over vertex levels to look at
     DO Q=I,PREIV_MAX
         
     !if usevar has vertex levels specified which are larger than I_VMAX, then set them to zero (caused if USE has no arguments specified)
        DO rr=1,8
            IF(USEVAR(Q,rr).gt.I_VMAX) USEVAR(Q,rr)=0
        ENDDO
     
     
         CUR_VERT=Q
         NOTHING=.FALSE.
     
         ! IF NONE IS SPECIFIED
         IF (PRE_TAY(3,Q).eq.0) THEN
              
             GIDHO=6
             VARSUM=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
     &                  FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,  &
     &                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,&
     &                  KSYM)
             
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
        IF(EXCITFUNCS(3).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "SigmaA, nA, sigmaB and nB POLYEXCITWEIGHTING: (Vertex level, Iteration number, Parameter values, Expected Variance)"
            ENDIF

            !Initial bracketing
            IF((bestvals(1,(Q-1)).eq.0.D0).and.(bestvals(2,(Q-1)).eq.0.D0).and.(bestvals(3,(Q-1)).eq.0.D0).and.(bestvals(4,(Q-1)).eq.0.D0)) THEN

                polypboth=(/ 0,1,0,1 /)     !Initial sigmaA, nA, sigmaB and nB.
                polyxiboth=RESHAPE( (/ 1.D0,((0.D0,r=1,4),1.D0,s=1,3) /), (/ 4, 4 /) )!Initial directions - unit vectors
            ELSE
                !Choose values which the previous vertex level found as optimum
                polypboth=bestvals(1:4,(Q-1))  
                polyxiboth=bestxipolyboth
            ENDIF
            n=4                 !Dimensions
            GIDHO=7
            
            CALL POWELL(polypboth,polyxiboth,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,    & 
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,                     &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM,LOWNEL)             
            
            bestvals(1:4,Q)=polypboth(:)
            bestxipolyboth=polyxiboth
            
            VARIANCES(Q)=fret
            
            IF (NOTHING) THEN

                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,F16.12,A,I3,A)") "Optimum SigmaA POLYEXCITWEIGHT found to be ", polypboth(1), "nA POLYEXCITWEIGHT as ", polypboth(2),           &
     &           ", sigmaB POLYEXCITWEIGHT as ", polypboth(3), ", and nB POLYEXCITWEIGHT as ", polypboth(4),", at vertex level ",     &
     &                    Q, ", but not using these values"

            ELSE
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        g_VMC_FINAL(1:4,USEVAR(Q,r))=polypboth(:)
                        WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,I3)") "POLYEXCITWEIGHTING parameters optimised to SigmaA=", polypboth(1), ", nA=",polypboth(2),", sigmaB=",polypboth(3), ", and nB=",polypboth(4)," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
            
            !If PolyExcit is specified rather than Excitweighting, we want to search for the minimum in a 4D space, where now we have two different excitto weighting parameters
        IF(EXCITFUNCS(2).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "A, sigma and n POLYEXCITWEIGHTING: (Vertex level, Iteration number, Parameter values, Expected Variance)"
            ENDIF

            !Initial bracketing
            IF((bestvals(1,(Q-1)).eq.0.D0).and.(bestvals(2,(Q-1)).eq.0.D0).and.(bestvals(3,(Q-1)).eq.0.D0)) THEN
                polyp=(/ 0.1,0.0,1.0 /)     !Initial A, sigma and n.
                polyxi=RESHAPE( (/ 1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0 /), (/ 3, 3 /) )!Initial directions - unit vectors
!                polyxi=RESHAPE( (/ 1.D0,0.D0,0.D0,0.D0,0.D0,1.D0,0.D0,1.D0,0.D0 /), (/ 3, 3 /) )!Initial directions -find n before sigma 
            ELSE
                !Choose values which the previous vertex level found as optimum
                polyp=bestvals(1:3,(Q-1))
                polyxi=bestxipoly
            ENDIF
            n=3                 !Dimensions
            GIDHO=4

            CALL POWELL(polyp,polyxi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,              &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,  &
     &              KSYM,LOWNEL)
            
            bestvals(1:3,Q)=polyp(:)
            bestxipoly=polyxi
        
            VARIANCES(Q)=fret
            
            IF (NOTHING) THEN

                WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,I3,A)") "Optimum A POLYEXCITWEIGHT found to be ", polyp(1), ", sigma POLYEXCITWEIGHT as ", polyp(2), ", and n POLYEXCITWEIGHT as ", polyp(3),", at vertex level ", Q,     &
     &           ", but not using these values"

            ELSE
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        g_VMC_FINAL(1:3,USEVAR(Q,r))=polyp(:)
                        WRITE(6,"(A,F16.12,A,F16.12,A,F16.12,A,I3)") "POLYEXCITWEIGHTING parameters optimised to A=", polyp(1), ", sigma=",polyp(2), ", and n=",polyp(3)," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
       
!*******************************************************************! - Looking at the Chemical potential with three adjustable parameters - two polynomials to make up the 'from' orbital choice

        IF (EXCITFUNCS(5).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "An, An' and Bn CHEMPOT-TWOFROM: (Vertex level, Iteration number, Parameter Values, Expected Variance)"
            ENDIF

            IF((bestvals(1,(Q-1)).eq.0.D0).and.(bestvals(2,(Q-1)).eq.0.D0).and.(bestvals(3,(Q-1)).eq.0.D0)) THEN
                polyp=(/ 1,1,1 /)    !Initial values
                polyxi=RESHAPE( (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/ 3, 3 /) ) !Initial Direction
            ELSE
                polyp=bestvals(1:3,(Q-1))
                polyxi=bestxipoly
            ENDIF
            n=3
            GIDHO=9

            CALL POWELL(polyp,polyxi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,                    &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,        &
     &              KSYM,LOWNEL)
        
            bestvals(1:3,Q)=polyp(:)
            bestxipoly=polyxi
            VARIANCES(Q)=fret

            IF(NOTHING) THEN
                
                WRITE(6,"(A,F17.12,A,F17.12,A,F17.12,A,I3,A)") "Optimum An CHEMPOT-TWOFROM parameters found to be ", polyp(1), " , An' as ", polyp(2), " , and Bn as ", polyp(3), " , at vertex level ", Q, " , but not using these values"
            ELSE
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        g_VMC_FINAL(1:3,USEVAR(Q,r))=polyp(:)
                        WRITE(6,"(A,F17.12,A,F17.12,A,F17.12,A,I3)") "CHEMPOT-TWOFROM parameters optimised to An= ", polyp(1), " , An'= ",polyp(2)," and Bn= ",polyp(3)," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            ENDIF
        ENDIF

!***********************************************************************************! - CHEMPOTWEIGHTING FUNCTION
        
        !Looking through the parameters for a polynomial function for occupied and virtual orbitals, with a cut-off at the chemical potential
        IF (EXCITFUNCS(4).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF(.not.TGRIDVAR(Q)) THEN    
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "An and Bn CHEMPOTWEIGHTING: (Vertex level, Iteration number, Parameter Values, Expected Variance)"
            ENDIF
                
            !Initial bracketing
            IF((bestvals(1,(Q-1)).eq.0.D0).and.(bestvals(2,(Q-1)).eq.0.D0)) THEN
                p=(/ 0.5,1.5 /)     !Initial values
                xi=RESHAPE( (/ 1.D0, 0.D0, 0.D0, 1.D0 /), (/ 2, 2 /) ) !Initial directions
            
            ELSE        !Choose values which previous vertex levels found were optimum
               p=bestvals(1:2,(Q-1))
               xi=bestxi
           ENDIF
           n=2
           GIDHO=8
           
            CALL POWELL(p,xi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,          &
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,           &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,&
     &              KSYM,LOWNEL)
        
            bestvals(1:2,Q)=p(:)
            bestxi=xi

            VARIANCES(Q)=fret

            IF (NOTHING) THEN
             
                WRITE(6,"(A,F16.12,A,F16.12,A,I3,A)") "Optimum An CHEMPOTWEIGHTING found to be ", p(1), ", and Bn CHEMPOTWEIGHTING as ", p(2), ", at vertex level ", Q, ",but not using these values"
            ELSE
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        
                        g_VMC_FINAL(1:2,USEVAR(Q,r))=p(:)
                
                        WRITE(6,"(A,F16.12,A,F16.12,A,I3)") "CHEMPOTWEIGHTING parameters optimised to An= ", p(1), " and Bn= ",p(2)," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            ENDIF
        !to print out landscape    
        ELSE

            abstr=''
            write (abstr,'(I1)') Q
            abstr='GRIDVAR-'//abstr
            n=2
            GIDHO=8
            UNITNO=100+Q
            OPEN(UNITNO,FILE=abstr,STATUS="UNKNOWN")
            CALL MAKEGRID(NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                     &
                         FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,         &
                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM,UNITNO)
            CLOSE(UNITNO)
        ENDIF
        ENDIF
        
        !Looking for a & b parameters
        IF ((EXCITFUNCS(1)).and.(PRE_TAYLOG(1,Q).or.NOTHING)) THEN
            IF(.not.TGRIDVAR(Q)) THEN 
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "A and B EXCITWEIGHTING: (Vertex level, Iteration number, Parameter Values, Expected Variance)"
            ENDIF
           
            !Initial bracketing
            IF((bestvals(1,(Q-1)).eq.0.D0).and.(bestvals(2,(Q-1)).eq.0.D0)) THEN
                p=(/ 0.5,0.5 /)     !Initial a and b values
                xi= RESHAPE( (/ 1.D0, 0.D0, 0.D0, 1.D0 /), (/ 2, 2 /) ) !Initial directions
            ELSE
                !Choose values which the previous vertex level found as optimum
                p=bestvals(1:2,(Q-1))
                xi=bestxi
            ENDIF
            n=2                 !Dimensions
            GIDHO=1             !To tell brentalgo that we are looking at a & b parameters

            CALL POWELL(p,xi,n,n,pre_TAYREAL(2,Q),iters,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,          &
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,           &
     &              TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,&
     &              KSYM,LOWNEL)

            bestvals(1:2,Q)=p(:)
            bestxi=xi

            VARIANCES(Q)=fret

            IF (NOTHING) THEN
                
                WRITE(6,"(A,F16.12,A,F16.12,A,I3,A)") "Optimum A EXCITWEIGHT found to be ", p(1), ", and B EXCITWEIGHT as ", p(2), ", at vertex level ", Q, ", but not using these values"
                
            ELSE
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        
                        g_VMC_FINAL(1:2,USEVAR(Q,r))=p(:)
                        WRITE(6,"(A,F16.12,A,F16.12,A,I3)") "EXCITWEIGHTING parameters optimised to A=", p(1), " and B=",p(2)," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            
            ENDIF
            
        !Print out landscape
        ELSE
            abstr=''
            write (abstr,'(I1)') Q
            abstr='GRIDVAR-'//abstr
            n=2
            GIDHO=1
            UNITNO=100+Q
            OPEN(UNITNO,FILE=abstr,STATUS="UNKNOWN")
            CALL MAKEGRID(NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                     &
                         FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,         &
                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM,UNITNO)
            CLOSE(UNITNO)
        ENDIF
!            distance=SQRT((ABS(xxxx-g_VMC_ExcitWeights(1))**2)+(ABS(xxxx-g_VMC_ExcitWeights(2))**2))
!            WRITE(6,*) "DISTANCE AWAY FROM EXACT VALUE IS: ", distance
            
        ENDIF
            
        !Looking for importance parameter
        IF (PRE_TAYLOG(3,Q).or.PRE_TAYLOG(4,Q)) THEN
            IF(.NOT.TLINEVAR(Q)) THEN
            
                IF (TLOGP) THEN
                    WRITE(31,*) ""
                    WRITE(31,"(A)") "Importance Parameter: (Vertex level, Iteration number, Importance Parameter, Expected Variance)"
                ENDIF
       
            !Reset 'FIRST(K)' in MCPATHSPRE     
                xxx=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
     &                  FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,       &
     &                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,155,ENERGYLIMS,KSYM)
     
                !Initial bracketing guess
                ax=0.7
                bx=0.9
!               cx=0.99 

                GIDHO=2
                originalimport=G_VMC_PI
            
                !Ensuring correct bracketing
                CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,      &
                        FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                          &
                        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)
            
                IF (TLOGP) THEN
                    WRITE(31,"(A,F16.12,A,F16.12)") "From mnbrak routine, minimum is between ", ax, " and ", bx
                ENDIF

                CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
                           G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,                  &
                           TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,fb,ENERGYLIMS,   &
                           KSYM)

                VARIANCES(Q)=minvar

                IF (PRE_TAYLOG(3,Q)) THEN
                    G_VMC_PI=originalimport
                    WRITE(6,"(A,F15.12,A,F5.3)") "Optimum importance parameter found to be ", xmin, ", but using ", originalimport
                ELSEIF (PRE_TAYLOG(4,Q)) THEN
                    G_VMC_PI=xmin
                    WRITE(6,"(A,F15.12)") "Importance parameter optimised to", xmin

                ENDIF
            ELSE

                abstr=''
                write (abstr,'(I1)') Q
                abstr='LINEVAR-'//abstr
                n=1
                GIDHO=2
                UNITNO=150+Q
                OPEN(UNITNO,FILE=abstr,STATUS="UNKNOWN")
                CALL MAKEGRID(NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                     &
                             FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,         &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM,UNITNO)
                CLOSE(UNITNO)
            ENDIF
        ENDIF
        
        !Looking for D Excitweighting parameter
        IF (PRE_TAYLOG(5,Q).or.PRE_TAYLOG(6,Q)) THEN

            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "D Parameter: (Vertex level, Iteration number, D Weighting parameter, Expected Variance)"
            ENDIF
            
            !Reset 'FIRST(K)' in MCPATHSPRE     
        xxx=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
     &          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,       &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,155,ENERGYLIMS,KSYM)
            
            !Initial bracketing guess
            ax=-0.8
            bx=-1.0
            
            GIDHO=5

            CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,  &
                    FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                      &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)
            IF (TLOGP) THEN
                WRITE(31,"(A,F16.12,A,F16.12)") "From mnbrak routine, minimum is between ", ax, " and ", cx
            ENDIF

            CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &                  G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                          &
     &                  RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &                  NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,fb,ENERGYLIMS,KSYM)
            
            VARIANCES(Q)=minvar
            
            IF (PRE_TAYLOG(5,Q)) THEN
                WRITE(6,"(A,F15.12,A,F5.3)") "Optimum D parameter found to be ", xmin, ", but not using this value"
            ELSEIF (PRE_TAYLOG(6,Q)) THEN
                DO r=1,I_VMAX
                    IF(USEVAR(Q,r).ne.0) THEN
                        
                        g_VMC_FINAL(3,USEVAR(Q,r))=xmin
                        WRITE(6,"(A,F15.12,A,I3)") "D parameter optimised to ", xmin," for vertex level ",USEVAR(Q,r)
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
        
        !Looking for C Excitweighting parameter
        IF (PRE_TAYLOG(2,Q).or.(PRE_TAYREAL(1,Q).gt.0).or.NOTHING) THEN
            IF(.not.TLINEVAR(Q)) THEN    
            
                IF (TLOGP) THEN
                    WRITE(31,*) ""
                    WRITE(31,"(A)") "U Parameter: (Vertex level, Iteration number, C Weighting parameter, Expected Variance)"
                ENDIF
           
                !Reset 'FIRST(K)' in MCPATHSPRE     
            xxx=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
     &              FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,       &
     &              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,155,ENERGYLIMS,KSYM)
        
                !Initial bracketing
                ax=20
                bx=30
!               cx=90

                GIDHO=3
                VARSUM=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
     &              FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                 &
     &              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
                ZEROVAR=VARSUM
            
                IF (TLOGP) THEN
                    WRITE(31,"(2I3,2G25.16)")  Q, 0, 0.D0, VARSUM
                    CALL FLUSH(31)
                ENDIF
            
                !To ensure correct bracketing
                CALL mnbrak(ax,bx,cx,fa,fb,fc,MCPATHSPRE,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,  &
                        FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                      &
                        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)
            
                IF (TLOGP) THEN
                    WRITE(31,"(A,F16.12,A,F13.9)") "From mnbrak routine, minimum is between ", ax, " and ", cx
                ENDIF
                    
                CALL BRENTALGO(minvar,ax,bx,cx,MCPATHSPRE,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &                                 G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                          &
     &                                 RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &                                  NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,fb,ENERGYLIMS,KSYM)
           
                VARIANCES(Q)=minvar
     
                IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).gt.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN

                    DO r=1,I_VMAX
                        IF(USEVAR(Q,r).ne.0) THEN
                            G_VMC_EXCITFINAL(USEVAR(Q,r))=xmin
                            WRITE(6,"(A,F17.12,A,I3)") "Optimum U weighting outside UEPSILON bounds, so using C=",xmin," for vertex level ",USEVAR(Q,r)
                        ENDIF
                    ENDDO
                
                ELSE IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).le.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN
                
                    DO r=1,I_VMAX
                        IF(USEVAR(Q,r).ne.0) THEN
                    
                            G_VMC_EXCITFINAL(USEVAR(Q,r))=0.D0
                            WRITE(6,"(A,I3)") "Optimum U weighting within UEPSILON bounds, so using C=0, for vertex level ",USEVAR(Q,r)
                        ENDIF
                    ENDDO
                ELSE
                
                    WRITE(6,"(A,G20.12,A,I3,A)") "Optimum U weighting found to be", xmin, " for vertex level ",Q, " ,but not using this value"
                END IF
        
            ELSE
                abstr=''
                write (abstr,'(I1)') Q
                abstr='LINEVAR-'//abstr
                n=1
                GIDHO=3
                UNITNO=151+Q
                OPEN(UNITNO,FILE=abstr,STATUS="UNKNOWN")
                CALL MAKEGRID(NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                     &
                             FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,         &
                             DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM,UNITNO)
                CLOSE(UNITNO)
            ENDIF
            
        END IF
    
    CALL FLUSH(6)


    !Save & print variances for each precalc vertex level!
    
    !End of vertex level
    ENDDO
    WRITE(6,*) ""
    
    IF (TLOGP) THEN
        WRITE(31,*) ""
        WRITE(31,"(A)") "Calculated expected variances for the vertex levels are:"
        DO zz=2,PREIV_MAX
            WRITE(31,"(A,I3,G25.16)") "Expected final variance for vertex ", zz, VARIANCES(zz)
        ENDDO
        CLOSE(31)
    ENDIF
    
    g_VMC_ExcitWeights(:,2:10)=G_VMC_FINAL
    g_VMC_EXCITWEIGHT(2:10)=G_VMC_EXCITFINAL
   
    do vv=2,I_VMAX
        check=.false.
        do kk=2,preIV_MAX
            do cc=1,I_VMAX
                IF(USEVAR(kk,cc).eq.vv) check=.true.
            enddo
        enddo
    !If use isn't specified for a vertex level, use the values given in the input file
        IF(.not.check) THEN
            g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
            G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
            IF(((TRUECYCLES.ne.0).or.(TOTALERROR.ne.0.D0)).and.(NWHTAY(1,vv).eq.-19)) THEN
                WRITE(6,*) "***VERTEX WEIGHTING PRECALCULATION not possible, since USE statement not attributed to all true MC levels***"
                TRUECYCLES=0
            ENDIF
        ENDIF
    enddo
    
    !calculate vertex level splitting!
    IF((TRUECYCLES.ne.0).or.(TOTALERROR.ne.0.D0)) THEN
        do vv=2,preIV_MAX
            do kk=1,I_VMAX
                IF(USEVAR(vv,kk).ne.0) THEN
                    gg=USEVAR(vv,kk)
                    IF(NWHTAY(1,gg).eq.-19) THEN
                        SUMSD=SUMSD+sqrt(VARIANCES(vv))
                        WRITE(6,"(A,I3,A,G20.12,A,I3)") "Expected MC variance for ",gg," vertices, is ", VARIANCES(vv)," taken from precalc level ",vv
                    ENDIF
                ENDIF
            enddo
        enddo
    ENDIF

    IF(TOTALERROR.ne.0.D0) THEN
        TRUECYCLES=(SUMSD/TOTALERROR)**2
    ENDIF
    
    IF(TRUECYCLES.ne.0) THEN
        do vv=2,preIV_MAX
            do kk=1,I_VMAX
                IF(USEVAR(vv,kk).ne.0) THEN
                    gg=USEVAR(vv,kk)
                    IF(NWHTAY(1,gg).eq.-19) THEN
                        NWHTAY(2,gg)=NINT(TRUECYCLES*(sqrt(VARIANCES(vv)))/SUMSD)
                        WRITE(6,"(A,I2,A,I9,A)") "Vertex weighting optimised, so that vertex level ", gg," has ",NWHTAY(2,gg)," cycles"
                    ENDIF
                ENDIF
            enddo
        enddo
    ENDIF

    !Deallocate arrays
    DO rr=1,7
        IF((pre_TAY(1,rr).eq.-7).or.(pre_TAY(1,rr).eq.-19)) DEALLOC=.true.
    ENDDO
    IF(DEALLOC) THEN
        xxx=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
     &          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,       &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,154,ENERGYLIMS,KSYM)
    ENDIF
    
!    write(6,*) g_VMC_Excitweights
!    call flush(6)
    
    RETURN
END SUBROUTINE GETVARS

FUNCTION VARIANCEAB(pointab,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,        &
           FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,            &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)

    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    REAL*8 pointab(*),VARIANCEAB
!    EXTERNAL MCPATHSPRE
!    REAL*8 MCPATHSPRE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),L,LT,K,D
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX),GIDHO
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,G
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    
    G=.false.
    
    SELECT CASE (GIDHO)
    CASE(1)
        g_VMC_ExcitWeights(1:2,K)=pointab(1:2)
    CASE(8)
        g_VMC_ExcitWeights(1:2,K)=pointab(1:2)
    CASE(9)
        g_VMC_ExcitWeights(1:3,K)=pointab(1:3)
    CASE(4)
        g_VMC_ExcitWeights(1:3,K)=pointab(1:3)
        IF ((pointab(2).gt.ENERGYLIMS(2)).or.(pointab(2).lt.ENERGYLIMS(1))) G=.true.
    CASE(7)
        g_VMC_ExcitWeights(1:4,K)=pointab(1:4)
        IF ((pointab(1).lt.ENERGYLIMS(1)).or.(pointab(1).gt.ENERGYLIMS(2))) G=.true.
        IF ((pointab(3).lt.ENERGYLIMS(1)).or.(pointab(3).gt.ENERGYLIMS(2))) G=.true.
    ENDSELECT
    
    IF (G) THEN
        VARIANCEAB=HUGE(VARIANCEAB)
    ELSE
    
    VARIANCEAB=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,     &
                    FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,   &
                    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,&
                    KSYM)

    ENDIF
    RETURN
                    
END FUNCTION VARIANCEAB
           

FUNCTION MCPATHSPRE(point,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,         &
              FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,        &
              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)

    USE HElem
    Use Determinants, only: GetHElement2
    use System, only: BasisFN
    Use Logging, only: PrevarLogging
    IMPLICIT NONE
    include 'irat.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),L,LT,K
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX),I,CNWHTAY
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,GIDHO
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5),sss,ierr,ierr2,ierr3,STORE(6),NMEMLEN,INODE2(NEL)
    INTEGER IEXCITS,J,EXCITGEN(0:PREIV_MAX),ierr4,b,dd,bb,aa,DEALLOCYC(2)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,FIRST(2:8)
    REAL*8 NTOTAL,BETA,ECORE,ALAT(3),RHOEPS,DBETA,VARSUM,point,MCPATHSPRE
    REAL*8 PROB,SumX,SumY,SumXsq,SumYsq,SumXY,ORIGIMPORT
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) DLWDB,TOTAL,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX),TOTAL2
    ! Saved values only go up to a vertex level of 6 - increase values if we want to go to higher vertex levels
    REAL*8 NTOTSAV(1:7),ENERGYLIMS(2),PGR,SUMPGEN,NTOTAL2
    TYPE(HDElement) DLWSAV(1:7),TOTSAV(1:7),RH,DLWDBCORE,WCORE,FF,OWEIGHT,OETILDE,FMCPR4D2
    REAL*8 INWI,PFAC,OPROB,VarX,VarY,Covar,NORMALISE
    DATA NTOTSAV/7*1.D0/
    DATA TOTSAV/7*HDElement(1.D0)/
    DATA DLWSAV/7*HDElement(0.D0)/
    DATA FIRST/7*.TRUE. /
    SAVE NTOTSAV,TOTSAV,DLWSAV
    SAVE FIRST,NORMALISE,DEALLOCYC
    LOGICAL LISNAN
    INTEGER ISEED,aaa,t,tt
    INTEGER I_OCLS,ITREE,ILOGGING,I_OVCUR,IACC
    REAL*8 ORIGEXCITWEIGHTS(6),ORIGEXCITWEIGHT
    REAL*8 XIJ(0:PREIV_MAX-1,0:PREIV_MAX-1)
   
    SELECT CASE (GIDHO)
    !Importance
    CASE(2)
        IF((point.GE.1.D0).or.(point.LE.0.D0)) THEN
            MCPATHSPRE=HUGE(MCPATHSPRE)
            RETURN
        ELSE
            G_VMC_PI=point
        ENDIF
    !U weighting
    CASE(3)
        G_VMC_EXCITWEIGHT(K)=point
    !D Weighting
    CASE(5)
        g_VMC_ExcitWeights(3,K)=point
    !NONE specified
    CASE(6)
        G_VMC_EXCITWEIGHT(K)=point
    CASE(154)
        IF(ASSOCIATED(GRAPHPARAMS)) THEN
            CALL MemDealloc(GRAPHPARAMS)
            DEALLOCATE(GRAPHPARAMS)
        ENDIF
        IF(ASSOCIATED(GRAPHS)) THEN
            CALL MemDealloc(GRAPHS)
            DEALLOCATE(GRAPHS)
        ENDIF
        IF(ASSOCIATED(NMEM)) THEN
            CALL MemDealloc(NMEM)
            DEALLOCATE(NMEM)
        ENDIF
        IF(ASSOCIATED(PVERTMEMS)) THEN
!            DO t=1,DEALLOCYC(2)
!                CALL FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,K,XIJ,          &
!     &              NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,   &
!     &              NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,0,154,         &
!     &              PVERTMEMS(:,t))
!            ENDDO
            CALL MemDealloc(PVERTMEMS)
            DO t=1,DEALLOCYC(2)
                DO tt=1,(DEALLOCYC(1)-1)
                    CALL FREEM(PVERTMEMS(tt,t))
                ENDDO
            ENDDO
            DEALLOCATE(PVERTMEMS)
        ENDIF
        IF(ASSOCIATED(PGENLIST)) THEN
            CALL MemDealloc(PGENLIST)
            DEALLOCATE(PGENLIST)
        ENDIF
        RETURN
    CASE(155)
        FIRST(K)=.TRUE.
        RETURN
    ENDSELECT
    
    IF((K.eq.2).and.(FIRST(2))) THEN
        DLWSAV(1)=HIJS(0)
    ENDIF
        
!    do D=2,K
        NTOTAL2=NTOTSAV(K-1)
        TOTAL2=TOTSAV(K-1)%v
        L=0
        LT=0
        VARSUM=0.D0
        DLWDB2=0.D0
        METH=pre_TAY(1,K)   !Method for vertex level
        CYCLES=pre_TAY(2,K)
        IF(METH.eq.-8) then !Full Rho-diag method
            CNWHTAY=0   !options for disallowing certain connections/graphs
            EREF=DLWSAV(K-1)/TOTSAV(K-1)
            
            F(K)=FMCPR3B(NI,BETA,I_P,IPATH,K,NEL,                                &
     &          NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,       &
     &          RHOEPS,0,RHOII,RHOIJ,CNWHTAY,METH,LOCTAB,                        &
     &          PreVarLogging,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,                   &
     &          MP2E,NTOTAL,EREF,VARSUM,TOTAL2)
        
        
        ELSEIF(METH.eq.-7) THEN !RHO-diag MC precalc
            STOP 'Rho-Diag not yet working for MC precalc'
        ELSEIF(METH.eq.-20) then !Full H-diag method
            EREF=DLWSAV(K-1)/TOTSAV(K-1)

            F(K)=FMCPR3B2(NI,BETA,I_P,IPATH,K,NEL,                               &
               NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
               RHOEPS,0,RHOIJ,0,METH,LOCTAB,                                     &
               PreVarLogging,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,                    &
               MP2E,NTOTAL2,PREIV_MAX,EREF,VARSUM,TOTAL2)

               
        ELSEIF(METH.eq.-19) THEN
           
            SumX=0.D0
            SumY=0.D0
            SumXsq=0.D0
            SumYsq=0.D0
            SumXY=0.D0
            VarX=0.D0
            VarY=0.D0
            Covar=0.D0
        
            !Will only need to run everytime a new vertex MC level is used
            !Deallocate the arrays if already been used by previous vertex level
            IF(FIRST(K)) THEN
                IF (ASSOCIATED(PVERTMEMS)) THEN
                    DO t=1,DEALLOCYC(2)
                        DO tt=1,(DEALLOCYC(1)-1)
                            CALL FREEM(PVERTMEMS(tt,t))
                        ENDDO
                    ENDDO
!                    DO t=1,DEALLOCYC(2)
                
!                        CALL FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,K,XIJ,          &
!                           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,    &
!                           NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,0,154,          &
!                           PVERTMEMS(:,t))
!                    ENDDO

                ENDIF
                
                IF (ASSOCIATED(PVERTMEMS)) THEN
                    CALL MemDealloc(PVERTMEMS)
                    DEALLOCATE(PVERTMEMS)
                ENDIF
                IF (ASSOCIATED(GRAPHS)) THEN
                    CALL MemDealloc(GRAPHS)
                    DEALLOCATE(GRAPHS)
                ENDIF
                IF (ASSOCIATED(PGENLIST)) THEN
                    CALL MemDealloc(PGENLIST)
                    DEALLOCATE(PGENLIST)
                ENDIF
                IF (ASSOCIATED(GRAPHPARAMS)) THEN
                    CALL MemDealloc(GRAPHPARAMS)
                    DEALLOCATE(GRAPHPARAMS)
                ENDIF
                IF (ASSOCIATED(NMEM)) THEN
                    CALL MemDealloc(NMEM)
                    DEALLOCATE(NMEM)
                ENDIF
                IF(MEMSAV(K)) THEN
                    ALLOCATE(PVERTMEMS(0:K,CYCLES),STAT=ierr4)
                    CALL MemAlloc(ierr4,PVERTMEMS,(K+1)*CYCLES/IRAT,'PRECALC_PVERTMEMS')
                    ALLOCATE(GRAPHS(NEL,0:K,CYCLES),STAT=ierr)
                    CALL MemAlloc(ierr,GRAPHS,NEL*(K+1)*CYCLES/IRAT,'PRECALC_GRAPHS')
                    ALLOCATE(GRAPHPARAMS(3,CYCLES),STAT=ierr2)
                    CALL MemAlloc(ierr2,GRAPHPARAMS,3*CYCLES,'PRECALC_GRAPHPARAMS')
                    ALLOCATE(PGENLIST(CYCLES),STAT=ierr2)
                    CALL MemAlloc(ierr2,PGENLIST,CYCLES,'PRECALC_GRAPHPARAMS')
                    DEALLOCYC=(/ K,CYCLES /)
                ENDIF

                CALL ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
                CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
                RHOII(0)=RH
                RHOIJ(0,0)=RHOII(0)
                HIJS(0)=GETHELEMENT2(NI,NI,NEL,NBASISMAX,                   &
     &          G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,0,ECORE)
! These variables have been passed in and represent the precalculated values for 1..I_VMIN-1 vertices
!Do these variable want to stay the same no matter what vertex level MC we are on?
                CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
!Setup the spin excit generator
                STORE(1)=0
                CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,.TRUE.,NMEMLEN,INODE2,I,0,STORE,3)
                ALLOCATE(NMEM(NMEMLEN),STAT=ierr3)
                CALL MemAlloc(ierr3,NMEM,NMEMLEN,'PRECALC_NMEM')
                NMEM(1)=0
                CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,.TRUE.,NMEM,INODE2,I,0,STORE,3)
!    Count the excitations (and generate a random one which we throw)
                ISEED=G_VMC_SEED
                CALL GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,NMEM,INODE2,ISEED,IEXCITS,0,UMAT,NMAX,PGR)
            
            ENDIF
                
                IF((FIRST(K).and.MEMSAV(K)).or.(.not.MEMSAV(K))) THEN 
                
                    !Now generate the graphs needed at that vertex
                    TOTAL=1.D0
                    DLWDBCORE=DLWSAV(K-1)
                    WCORE=TOTSAV(K-1)
                    !EREF IS 1v DLWDB2
                    EREF=DLWDBCORE/WCORE
                    OETILDE=EREF
                    DLWDB2=0.D0
                    TSYM=.false.
                    ILOGGING=0
                    I_OVCUR=1
                    I_OCLS=0
                    ITREE=1
                    OPROB=1.D0 !come from 1v graph
                    ISEED=g_VMC_SEED
                    PFAC=1.D0
                !*** Would want to use gidho here to determine which parameters we want to unbias for ***
!Unbias graph generation (still relative probabilities of generating different classes)
                    IF(GIDHO.eq.3) THEN
                        !biasing against U parameter
                        ORIGEXCITWEIGHT=G_VMC_EXCITWEIGHT(K)
                        G_VMC_EXCITWEIGHT(K)=0.D0
                    ELSEIF(GIDHO.eq.2) THEN 
                        !Biasing against importance parameter
                        ORIGIMPORT=G_VMC_PI
                        G_VMC_PI=-1.D0             !If Importance parameter set to -1, this means that we want an equal weighting
                    ELSE
                        !biasing against other excitation parameters
                        ORIGEXCITWEIGHTS(:)=g_VMC_ExcitWeights(:,K)
                        g_VMC_ExcitWeights(:,K)=0.D0
                    ENDIF
                    NORMALISE=0.D0

                    !If MEMSAV not on, then don't save excitation generators

!                OPEN(43,FILE="MCGRAPHS",STATUS="UNKNOWN")
                    DO b=1,CYCLES
                        OWEIGHT=0.D0
                        EXCITGEN=0
                        INWI=0.D0
                        !DLWDB2 is just the energy, OETILDE is energy*weight
                        FF=FMCPR4D2(NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,                            &
                                              NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,CYCLES,METH,                   &
                                             PreVarLOGGING,TSYM,ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,OETILDE,OPROB,I_OVCUR,&
                                              I_OCLS,ITREE,OWEIGHT,PFAC,IACC,INWI,K,EXCITGEN(0:K))
!                        WRITE(60,*) IPATH(:,0),IPATH(:,1),IPATH(:,2),OPROB
!                        WRITE(60,*) OPROB
!                   WRITE(43,("I3,3G25.16")) K, DLWDB2, OETILDE, OWEIGHT
!                   WRITE(43,*) DLWDB2, OETILDE/OWEIGHT
!                   WRITE(43,*) DLWDB2
!                   CALL WRITEPATH(43,IPATH(:,0:K),K,NEL,.TRUE.)
                        NORMALISE=NORMALISE+OPROB

                   IF(MEMSAV(K)) THEN
                        GRAPHS(:,:,b)=IPATH(:,0:K)
                        GRAPHPARAMS(1,b)=OWEIGHT%v
                        GRAPHPARAMS(2,b)=(DLWDB2%v-EREF%v)
                        GRAPHPARAMS(3,b)=OPROB
                        !PVERTMEMS(0,b)=EXCITGEN(0)
                        PVERTMEMS(0:K,b)=EXCITGEN(0:K)
!                   WRITE(6,("5G25.16")) FF%v,OWEIGHT%v,DLWDB2%v,OPROB,OETILDE%v
!                   CALL FLUSH(6)
                   
                    ELSE

                        DLWDB2=DLWDB2-EREF
                        J=0
                        PROB=0.D0
                        IF(GIDHO.eq.3) THEN
                            G_VMC_EXCITWEIGHT(K)=ORIGEXCITWEIGHT
                        ELSEIF(GIDHO.eq.2) THEN
                            G_VMC_PI=ORIGIMPORT
                        ELSE
                            g_VMC_ExcitWeights(:,K)=ORIGEXCITWEIGHTS(:)
                        ENDIF
                        CALL CalcWriteGraphPGen(J,IPATH,K,NEl,LOCTAB,G1,               &
                                  NBASISMAX,UMat,NMAX,NBASIS,PROB,EXCITGEN(0:K))
                   
                        SumX=SumX+((OWEIGHT%v*DLWDB2%v)/OPROB)
                        SumY=SumY+((OWEIGHT%v+(WCORE%v*PROB))/OPROB)
                        SumXsq=SumXsq+(((OWEIGHT%v*DLWDB2%v)**2)/(OPROB*PROB))
                        SumYsq=SumYsq+((((OWEIGHT%v/PROB)+WCORE%v)**2)*(PROB/OPROB))
                        SumXY=SumXY+(((OWEIGHT%v/PROB)+WCORE%v)*(OWEIGHT%v*DLWDB2%v/PROB)*    &
                                (PROB/OPROB))

                        IF(GIDHO.eq.3) THEN
                            G_VMC_EXCITWEIGHT(K)=0.D0
                        ELSEIF(GIDHO.eq.2) THEN
                            G_VMC_PI=-1.D0
                        ELSE
                            g_VMC_ExcitWeights(:,K)=0.D0
                        ENDIF
                    CALL FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,K,XIJ,          &
                       NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,    &
                       NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,0,154,          &
                       EXCITGEN)
!                        DO tt=1,K-1
!                            CALL FREEM(EXCITGEN(tt))
!                        ENDDO
                        
                    ENDIF
                    
                ENDDO
                IF(GIDHO.eq.3) THEN
                    G_VMC_EXCITWEIGHT(K)=ORIGEXCITWEIGHT
                ELSEIF(GIDHO.eq.2) THEN
                    G_VMC_PI=ORIGIMPORT
                ELSE
                    g_VMC_ExcitWeights(:,K)=ORIGEXCITWEIGHTS(:)
                ENDIF
            !Choosing whether to look at graphs or not
            ENDIF                
!                CLOSE(43)
            
            !So what is av. Energy for this vertex level??
                          
            IF(MEMSAV(K)) THEN
                
                J=0
            !    
                DO bb=1,CYCLES
                   
            !        
                    PROB=0.D0
                    IPATH(:,0:K)=GRAPHS(:,0:K,bb)
                    EXCITGEN(0:K)=PVERTMEMS(:,bb)
                    CALL  CalcWriteGraphPGen(J,IPATH,K,NEl,LOCTAB,G1,               &
                               NBASISMAX,UMat,NMAX,NBASIS,PROB,EXCITGEN(0:K))
                    
                    PGENLIST(bb)=PROB
            !
                ENDDO
            !
                !CALCULATE VARIANCE
                !Does EREF/WREF want to change as we go to higher MC vertex levels?
                DO aa=1,CYCLES
            !
                    SumX=SumX+((GRAPHPARAMS(1,aa)*GRAPHPARAMS(2,aa))/GRAPHPARAMS(3,aa))
                    SumY=SumY+((GRAPHPARAMS(1,aa)+(WCORE%v*PGENLIST(aa)))/GRAPHPARAMS(3,aa))
                    SumXsq=SumXsq+(((GRAPHPARAMS(1,aa)*GRAPHPARAMS(2,aa))**2)/(GRAPHPARAMS(3,aa)*PGENLIST(aa)))
                    SumYsq=SumYsq+((((GRAPHPARAMS(1,aa)/PGENLIST(aa))+WCORE%v)**2)*(PGENLIST(aa)/GRAPHPARAMS(3,aa)))
                    SumXY=SumXY+(((GRAPHPARAMS(1,aa)/PGENLIST(aa))+WCORE%v)*(GRAPHPARAMS(1,aa)*GRAPHPARAMS(2,aa)/PGENLIST(aa))*    &
                        (PGENLIST(aa)/GRAPHPARAMS(3,aa)))
            !
            
!                 SUMPGENS=SUMPGENS+PGENLIST(aa)
                ENDDO    
            !
            ENDIF
            
            SumX=SumX/(CYCLES+0.D0)
            SumY=SumY/(CYCLES+0.D0)
            SumXsq=SumXsq/(CYCLES+0.D0)
            SumYsq=SumYsq/(CYCLES+0.D0)
            SumXY=SumXY/(CYCLES+0.D0)
            VarX=SumXsq-(SumX**2)
            VarY=SumYsq-(SumY**2)
            Covar=SumXY-(SumX*SumY)
            
!            WRITE(18,*) g_VMC_ExcitWeights(1,K),g_VMC_ExcitWeights(2,K)
!            WRITE(18,*) SumX,SumY,VarX,VarY,Covar

            VARSUM=((SumX/SumY)**2)*((VarX/SumX**2)+(VarY/SumY**2)-(2*Covar/(SumX*SumY)))
            
        ENDIF
        
        IF(FIRST(K)) THEN
            IF((METH.EQ.-20).or.(METH.EQ.-8)) THEN
                
                NTOTSAV(K)=NTOTAL2
                TOTSAV(K)=TOTSAV(K-1)+F(K)
                DLWSAV(K)=DLWSAV(K-1)+DLWDB2
!                WRITE(6,*) g_VMC_Excitweights(:,K)
!                WRITE(6,*) "EXACT EXPECTED VARIANCE: ", VARSUM
            ELSEIF(METH.EQ.-19) THEN
            
!                WRITE(6,*) g_VMC_Excitweights(:,K)
!                WRITE(6,*) "CONVERGENCE: ", VARSUM
                WRITE(6,*) "SUM OF UNBIASED PGENS IN MCPRECALC: ", NORMALISE
                !???
                TOTSAV(K)=TOTSAV(K-1)
                DLWSAV(K)=DLWSAV(K-1)
            ENDIF
        ENDIF
        
        MCPATHSPRE=VARSUM
        FIRST(K)=.FALSE.
        
        RETURN
END FUNCTION MCPATHSPRE

!Now not called as didn't seem to like having an allocatable array passed to it - wanted to have pointer passed to it
SUBROUTINE GETGRAPHS(METH,CYCLES,GRAPHS,GRAPHPARAMS,PVERTMEMS,NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,G1,    &
                         NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,TSYM,           &
                         ECORE,KSYM,DBETA,DLWDB2,HIJS,NMEM,ISEED)

    USE HElem
    use System, only: BasisFN
    use Logging, only: PreVarLogging
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),NI(NEL)
    INTEGER nBasisMax(5,5),NBASIS,IPATH(NEL,0:PREIV_MAX),I_V
    INTEGER IACC,KSYM(5),b
    COMPLEX*16 FCK(*)
    REAL*8 BETA,ECORE,ALAT(3),RHOEPS,DBETA,OPROB
    LOGICAL TSYM
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) DLWDB2,RHOII(0:PREIV_MAX),OETILDE,FMCPR4D2
    TYPE(HDElement) FF,OWEIGHT
    INTEGER GRAPHS(NEL,0:I_V,CYCLES),ISEED,ILOGGING,I_OVCUR
#if defined(POINTER8)
    INTEGER*8 PVERTMEMS(0:I_V,CYCLES)
#else
    INTEGER PVERTMEMS(0:I_V,CYCLES)
#endif
    INTEGER EXCITGEN(0:I_V)
    INTEGER I_OCLS,ITREE
    REAL*8 GRAPHPARAMS(3,CYCLES)
    REAL*8 INWI,PFAC
    REAL*8 NMEM(*),ORIGEXCITWEIGHTS(6)
   
    OETILDE=DLWDB2 
    DLWDB2=0.D0
    TSYM=.false.
    ILOGGING=0
    I_OVCUR=1
    I_OCLS=0
    ITREE=1
    OPROB=1.D0 !come from 1v graph
    ISEED=g_VMC_SEED
    PFAC=1.D0

    !Unbias graph generation (still relative probabilities of generating different classes)
    ORIGEXCITWEIGHTS(:)=g_VMC_ExcitWeights(:,I_V)
    g_VMC_ExcitWeights(:,I_V)=0.D0
    
    WRITE(6,*) "WEIGHT       OWEIGHT      DLWDB2        PGEN         OETILDE"
    DO b=1,CYCLES
        OWEIGHT=0.D0
        EXCITGEN=0
        INWI=0.D0
        FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,                          &
                              NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,CYCLES,METH,                   &
                              PrevarLOGGING,TSYM,ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,OETILDE,OPROB,I_OVCUR,&
                              I_OCLS,ITREE,OWEIGHT,PFAC,IACC,INWI,I_V,EXCITGEN)

        GRAPHS(:,:,b)=IPATH(:,0:I_V)
        GRAPHPARAMS(1,b)=OWEIGHT%v
        GRAPHPARAMS(2,b)=DLWDB2%v
        GRAPHPARAMS(3,b)=OPROB
        PVERTMEMS(:,b)=EXCITGEN
        
        !E/w
        WRITE(6,"(5G25.16)") FF%v,OWEIGHT%v,DLWDB2%v,OPROB,OETILDE%v
        CALL FLUSH(6)
    ENDDO

    g_VMC_ExcitWeights(:,I_V)=ORIGEXCITWEIGHTS(:)

    RETURN
    END SUBROUTINE GETGRAPHS 

SUBROUTINE BRENTALGO(brent,ax,bx,cx,fun,tol,xmin,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,    &
     &                 G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,   &
     &                 RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,       &
     &                 NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,INITFUNC,ENERGYLIMS,KSYM)
    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),K,L,LT
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,GIDHO
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM,fun,INITFUNC
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER ITMAX
!    EXTERNAL fun
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

    IF (TLOGP) write(31,*) "BRENT ALGO STARTING"
    IF (INITFUNC.eq.0.D0) THEN
        IF (TLOGP) WRITE(31,*) "INITFUNC EQUAL 0.D0 - redo initial point"
        VARSUM=fun(x,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                &
              FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
    ELSE
        VARSUM=INITFUNC
    ENDIF
    
    ! Vertex level, iteration number, C weight value, Varsumnu value
    IF(INITFUNC.eq.0.D0) THEN
        IF (TLOGP.and.(GIDHO.ne.1).and.(GIDHO.ne.4).and.(GIDHO.ne.7).and.(GIDHO.ne.8).and.(GIDHO.ne.9)) THEN
            WRITE(31,"(2I3,2G25.16)")  K, 0, x, VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.AND.(GIDHO.EQ.8)) THEN
            WRITE(31,"(I3,A,3G25.16)")  K, " brent  0", g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.AND.(GIDHO.EQ.1)) THEN
            WRITE(31,"(I3,A,3G25.16)")  K, " brent  0", g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.and.((GIDHO.eq.4).or.(GIDHO.eq.9))) THEN
            WRITE(31,"(I3,A,4G25.16)")  K, " brent  0", g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), g_VMC_ExcitWeights(3,K), VARSUM
            CALL FLUSH(31)
        ENDIF
        IF (TLOGP.and.(GIDHO.eq.7)) THEN
            WRITE(31,"(I3,A,5G25.16)")  K, " brent  0", g_VMC_ExcitWeights(1,K),g_VMC_ExcitWeights(2,K),g_VMC_ExcitWeights(3,K),g_VMC_ExcitWeights(4,K),VARSUM
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
                  FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
                  DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
             IF (TLOGP.and.(GIDHO.eq.7)) THEN
                 WRITE(31,"(I3,A,I3,5G25.16)")  K, " brent", iter, g_VMC_ExcitWeights(1,K),g_VMC_ExcitWeights(2,K), g_VMC_ExcitWeights(3,K), g_VMC_ExcitWeights(4,K), VARSUM
                 CALL FLUSH(31)
             ENDIF
             IF (TLOGP.and.((GIDHO.eq.4).or.(GIDHO.eq.9))) THEN
                 WRITE(31,"(I3,A,I3,4G25.16)")  K, " brent", iter, g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), g_VMC_ExcitWeights(3,K), VARSUM
                 CALL FLUSH(31)
             ENDIF
             IF (TLOGP.and.(GIDHO.eq.8)) THEN
                 WRITE(31,"(I3,A,I3,3G25.16)")  K, " brent", iter, g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), VARSUM
                 CALL FLUSH(31)
             ENDIF                                  
             IF (TLOGP.and.(GIDHO.eq.1)) THEN
                 WRITE(31,"(I3,A,I3,3G25.16)")  K, " brent", iter, g_VMC_ExcitWeights(1,K), g_VMC_ExcitWeights(2,K), VARSUM
                 CALL FLUSH(31)
             ENDIF                                  
             IF (TLOGP.and.(GIDHO.ne.1).and.(GIDHO.ne.4).and.(GIDHO.ne.7).and.(GIDHO.ne.8).and.(GIDHO.ne.9)) THEN
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
        IF (TLOGP) write(31,*) "BRENT ALGO FINISHED"
    return
END SUBROUTINE BRENTALGO

SUBROUTINE POWELL(p,xi,n,np,ftol,iter,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,                   &
     &        G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,        &
     &        TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,&
     &        KSYM,LOWNEL)

    
    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER iter,n,np,NMAX,ITMAX,NEL,I_P,BRR(*),NMSH,NTAY(2),L,LT,NMAXI
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX),Q
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,GIDHO
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,TLOGP,LOWNEL
    REAL*8 fret,ftol,p(np),xi(np,np),BETA,ECORE,NTOTAL,ALAT(3),RHOEPS
    REAL*8 DBETA,VARSUM,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
!    REAL*8 varianceab
!    EXTERNAL varianceab 
    PARAMETER (NMAXI=4,ITMAX=200)
!    USES func,linmin
!       Minimisation of a function 'func' of n variables. (func is not an argument, it is a fixed function name.)  
!Input consists of an initial starting point p(1:n); an initial matrix xi(1:n,1:n) with physical dimensions np by np, and whose columns contain the initial set of directions (usually the n unit vectors). 
!'ftol' is the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness.
!On output, p is set to the best point found, xi is the then current direction set, fret in the returned function value at p, and iter in the number of iterations taken.  The routine linmin is used.
!       Parameters:  Maximum expected value of n, and maximum allowed iterations.
    INTEGER i,ibig,j
    REAL*8 del,fp,fptt,t,pt(NMAXI),ptt(NMAXI),xit(NMAXI)
    fret=varianceab(p,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,      &
           FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
           KSYM)
    IF (TLOGP.and.(GIDHO.eq.7)) THEN
        WRITE(31,"(2I3,5G25.16)") Q, 0, p(1), p(2), p(3), p(4), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.((GIDHO.eq.4).or.(GIDHO.eq.9))) THEN
        WRITE(31,"(2I3,4G25.16)") Q, 0, p(1), p(2), p(3), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.((GIDHO.eq.1).or.(GIDHO.eq.8))) THEN
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
        IF(LOWNEL.AND.(BTEST(i,0)).AND.((GIDHO.EQ.1).OR.(GIDHO.EQ.7).OR.(GIDHO.EQ.8))) CYCLE
        IF((GIDHO.eq.9).AND.(i.eq.2).AND.(Q.eq.2)) CYCLE  !At the two-vertex level, changing the excit-from above chempot will not affect variance
        do j=1,n                !Copy the direction...
            xit(j)=xi(j,i)
        enddo
        fptt=fret
        
        call linmin(p,xit,n,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,          &
     &    FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                    &
     &    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)   !Minimise along it...
        IF (TLOGP.and.((GIDHO.eq.1).or.(GIDHO.eq.8))) THEN
            WRITE(31,"(2I3,3G25.16)") Q, iter, p(1), p(2), fret
            CALL FLUSH(31)
        ELSEIF (TLOGP.and.(GIDHO.eq.7)) THEN
            WRITE(31,"(2I3,5G25.16)") Q, iter, p(1), p(2), p(3), p(4), fret
            CALL FLUSH(31)
        ELSEIF (TLOGP.and.((GIDHO.eq.4).or.(GIDHO.eq.9))) THEN
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
           FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,           &
           DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
    !Function value at extrapolated point
    if (fptt.ge.fp) goto 1      !One reason not to use new direction
    t=2.D0*(fp-2.D0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
    if (t.ge.0.D0) goto 1       !Other reason not to use new direction
    CALL linmin(p,xit,n,fret,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,          &
     &    FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                &
     &    DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM) !Move to the minimum of the new direction,
    IF (TLOGP.and.((GIDHO.eq.1).or.(GIDHO.eq.8))) THEN
        WRITE(31,"(2I3,3G25.16)") Q, iter, p(1), p(2), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.(GIDHO.eq.7)) THEN
        WRITE(31,"(2I3,5G25.16)") Q, iter, p(1), p(2), p(3), p(4), fret
        CALL FLUSH(31)
    ELSEIF (TLOGP.and.((GIDHO.eq.4).or.(GIDHO.eq.9))) THEN
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
     &       FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,             &
     &       DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)
    
    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),NMSH,NMAX,NTAY(2),Q,L,LT,GIDHO
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM,ENERGYLIMS(2)
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
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
!    EXTERNAL f1dim
!    REAL*8 f1dim    
    ncom=n                  !Set up the common block
    do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
    enddo
    SELECT CASE (GIDHO)
    CASE(8)
        ax=0.D0
        xx=2.D0
    CASE DEFAULT
        ax=0.D0
        xx=5.D-01
    END SELECT
    
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
     &          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,              &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)

    call BRENTALGO(fret,ax,xx,bx,f1dim,TOL,xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &       G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,              &
     &       RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &       NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,fx,ENERGYLIMS,KSYM)
     do j=1,n                !Construct the vector results to return
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
    enddo
    return
    END SUBROUTINE linmin

FUNCTION f1dim(x,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
     &          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM) 

    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),L,LT,K,D,Q
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,GIDHO
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,LISNAN
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER NMAXI
    REAL*8 f1dim,func,x,ENERGYLIMS(2)
    PARAMETER (NMAXI=4)
!    EXTERNAL varianceab
!    REAL*8 varianceab
    !USES func
!Used by linmin as the function passed to mnbrak and brentalgo
    INTEGER j,ncom
    REAL*8 pcom(NMAXI),xicom(NMAXI),xt(NMAXI)
    COMMON /f1com/ pcom,xicom,ncom
    do j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
!        IF(LISNAN(xt(j))) THEN
!            WRITE(6,*) xt(j)
!            WRITE(6,*) pcom(j),x,xicom(j)
!            CALL FLUSH(6)
!        ENDIF
    enddo
    f1dim=varianceab(xt,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
              FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,        &
              DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
    return
END FUNCTION f1dim
   
SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,   &
               FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,TLOGP,ENERGYLIMS,KSYM)
    
    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,MINI
    PARAMETER (GOLD=1.618034, GLIMIT=100.D0,MINI=1.D-20)
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH,NMAX,NTAY(2),L,LT,Q
    INTEGER NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX),GIDHO
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,t
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,savedax,savedbx
    REAL*8 ENERGYLIMS(2)
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
!    EXTERNAL func
    !Given a function 'func', and given distinct initial points ax and bx, this routine searches in the downhill direction (defined by the function as evaluated at the initial points) and returns new points ax,bx,cx that bracket a minimum of the function.
    !Also returned are the function values at the three points, fa, fb, fc.
    !Parameters: GOLD in the default ratio by which successive intervals are magnified; GLIMIT is the maximum magnification allowed for a parabolic-fit step.
    REAL*8 dum,fu,qu,r,u,ulim
!    WRITE(31,*) "MNBRAK ALGO STARTING"
2   savedax=ax
    savedbx=bx
    fa= func(ax,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
               FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,   &
               KSYM)
    fb= func(bx,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,                  &
               FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
               DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,   &
               KSYM)
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
          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
          KSYM)
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
                      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                      KSYM)
     
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
                      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                      KSYM)
            ELSEIF((cx-u)*(u-ulim).gt.0.D0) THEN    !Parabolic fit is between c and its allowed limit
                
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                      KSYM)
                IF(fu.lt.fc) THEN
                    bx=cx
                    cx=u
                    u=cx+GOLD*(cx-bx)
                    fb=fc
                    fc=fu
                    fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                          KSYM)
                ENDIF
            ELSEIF((u-ulim)*(ulim-cx).ge.0.D0) THEN     !Limit parabolic u to maximum allowed value
                u=ulim
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                      KSYM)
            ELSE        !Reject parabolic u, use default magnification
                u=cx+GOLD*(cx-bx)
                fu=func(u,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,             &
                      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,    &
                      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS, &
                      KSYM)
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

SUBROUTINE MAKEGRID(NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,         &
        FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,        &
        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM,UNITNO)

    USE HElem
    use System, only: BasisFN
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),NMSH,NMAX,NTAY(2),K
    INTEGER L,LT,NI(NEL),nBasisMax(5,5),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),NBASIS,GIDHO
#if defined(POINTER8)
     INTEGER*8 LOCTAB(3,PREIV_MAX)
#else
     INTEGER LOCTAB(3,PREIV_MAX)
#endif
    INTEGER KSYM(5)
    INTEGER IEXCITS,UNITNO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM
    REAL*8 NTOTAL,BETA,ECORE,ALAT(3),RHOEPS,DBETA,VARSUM,A,B
    TYPE(HElement) UMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) DLWDB,TOTAL,DLWDB2,EREF
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    ! Saved values only go up to a vertex level of 6 - increase values if we want to go to higher vertex levels
    REAL*8 ENERGYLIMS(2),origvals(2),p(2),origval
!    EXTERNAL MCPATHSPRE
!    REAL*8 MCPATHSPRE
IF(TLINEVAR(K)) THEN
    IF(GIDHO.eq.2) THEN
        origval=g_VMC_PI
    ELSEIF(GIDHO.eq.3) THEN
        origval=G_VMC_EXCITWEIGHT(K)
    ENDIF

    DO A=LINEVARPAR(K,1),LINEVARPAR(K,2),LINEVARPAR(K,3)
        
        VARSUM=MCPATHSPRE(A,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,              &
     &      FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                 &
     &      DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
        WRITE(UNITNO,"(F9.3,G25.16)") A,VARSUM
        CALL FLUSH(UNITNO)
    ENDDO
    
    IF(GIDHO.eq.2) THEN
        g_VMC_PI=origval
    ELSEIF(GIDHO.eq.3) THEN
        G_VMC_EXCITWEIGHT(K)=origval
    ENDIF
    
ELSE
    origvals(:)=g_VMC_ExcitWeights(1:2,K)
    
    DO A=GRIDVARPAR(K,1),GRIDVARPAR(K,2),GRIDVARPAR(K,3)
        DO B=GRIDVARPAR(K,4),GRIDVARPAR(K,5),GRIDVARPAR(K,6)
            p=(/ A,B /)     !Initial a and b values
            g_VMC_ExcitWeights(1:2,K)=p(:)
        
            VARSUM=MCPATHSPRE(0.D0,NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
     &          FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,                 &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,ENERGYLIMS,KSYM)
            WRITE(UNITNO,"(2F9.3,G25.16)") A,B,VARSUM
            CALL FLUSH(UNITNO)
        ENDDO
        WRITE(UNITNO,*) ""
    ENDDO
    g_VMC_ExcitWeights(1:2,K)=origvals(:)
ENDIF

ENDSUBROUTINE MAKEGRID

end module PreCalc
