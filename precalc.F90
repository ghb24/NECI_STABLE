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
     INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),GIDHO
     COMPLEX*16 FCK(*)
     REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM
     REAL*8 ax,bx,cx,minvar,xmin,originalc,originalimport
     REAL*8 ZEROVAR
     LOGICAL TSYM,PREVAR,NOTHING,TLOGP
     TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
     TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
     TYPE(HElement) HIJS(0:PREIV_MAX)
     TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
     
     !Only open PRECALC file if logging option on
     IF (TLOGP) THEN
        OPEN(31,FILE="PRECALC",STATUS="UNKNOWN")
     ENDIF
      
     !Loop over vertex levels to look at
     DO Q=I,PREIV_MAX
         
         NOTHING=.FALSE.
     
         ! IF NONE IS SPECIFIED
         IF (PRE_TAY(3,Q).eq.0) CYCLE
     
         WRITE(6,*) ""
         WRITE(6,"(A,I2,A)") "For a vertex level of", Q, ", PRECALC finds:"
     !IF ... do for a&b
!     IF(PRE_TAYLOG(3,K)) THEN   !wants to use a&b calculated earlier (k should be more than 3)
         !Find import
!     ENDIF
        
        !IF NOTHING SPECIFIED
        IF (PRE_TAYLOG(1,Q).and.PRE_TAYLOG(2,Q).and.PRE_TAYLOG(3,Q).and.PRE_TAYLOG(4,Q).and.            &
     &      (PRE_TAYREAL(1,Q).eq.0.D0).and.(PRE_TAYREAL(2,Q).eq.1.D-01).and.(PRE_TAY(3,Q).ne.0)) THEN
            NOTHING=.TRUE.
        ENDIF

        !Looking for importance parameter
        IF (PRE_TAYLOG(3,Q).or.PRE_TAYLOG(4,Q)) THEN
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "Importance Parameter: (Vertex level, Iteration number, Importance Parameter, Expected Variance)"
            ENDIF
            
            !Initial bracketing
            ax=0.6
            bx=0.9
            cx=0.99

            GIDHO=2
            originalimport=G_VMC_PI
            CALL BRENTALGO(minvar,ax,bx,cx,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,     &
                       G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,       &
                       TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP)

            IF (PRE_TAYLOG(3,Q)) THEN
                G_VMC_PI=originalimport
                WRITE(6,"(A,F15.12,A,F5.3)") "Optimum importance parameter found to be ", xmin, ", but using ", originalimport
            ELSEIF (PRE_TAYLOG(4,Q)) THEN
                G_VMC_PI=xmin
                WRITE(6,"(A,F15.12)") "Importance parameter optimised to", xmin

            ENDIF
        ENDIF

        !Looking for C Excitweighting parameter
        IF (PRE_TAYLOG(2,Q).or.(PRE_TAYREAL(1,Q).gt.0).or.NOTHING) THEN
            
            IF (TLOGP) THEN
                WRITE(31,*) ""
                WRITE(31,"(A)") "U Parameter: (Vertex level, Iteration number, C Weighting parameter, Expected Variance)"
            ENDIF
            
         !Initial bracketing
            ax=0
            bx=30
            cx=90

            GIDHO=3
            originalc=G_VMC_EXCITWEIGHT
            G_VMC_EXCITWEIGHT=0.D0
            CALL MCPATHSPRE(NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
     &          FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,          &
     &          DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,VARSUM,PREVAR)
            ZEROVAR=VARSUM
            
            IF (TLOGP) THEN
                WRITE(31,"(2I3,2G25.16)")  Q, 0, 0.D0, VARSUM
                CALL FLUSH(31)
            ENDIF
            
            CALL BRENTALGO(minvar,ax,bx,cx,pre_TAYREAL(2,Q),xmin,NI,BETA,I_P,IPATH,Q,NEL,NBASISMAX,      &
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,                          &
     &                  RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,                  &
     &                  NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP)
           
            IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).gt.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN
                G_VMC_EXCITWEIGHT=xmin
                WRITE(6,"(A,G25.12)") "Optimum U weighting outside UEPSILON bounds, so using C=",xmin
            ELSE IF ((PRE_TAYLOG(1,Q).or.(PRE_TAYREAL(1,Q).gt.0)).and.((zerovar-minvar).le.PRE_TAYREAL(1,Q)).and..not.NOTHING) THEN
                G_VMC_EXCITWEIGHT=0
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

!NEED ROUTINE TO GET INITIAL THREE POINTS

SUBROUTINE MCPATHSPRE(NI,BETA,I_P,IPATH,K,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,           &
     &        FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
     &        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,VARSUM,PREVAR)

    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,L,LT,K,D
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM
!    REAL*8 VARIANCE(2,2:PREIV_MAX)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    
    NTOTAL=1.D0
    TOTAL=1.D0
    DLWDB=HIJS(0)
!OPEN(PRECALC)
    do D=2,K
        L=0
        LT=0
        VARSUM=0.D0
        DLWDB2=0.D0
        METH=pre_TAY(1,D)   !Method for vertex level
        CYCLES=pre_TAY(2,D)
        IF(METH.eq.-8) then !Full Rho-diag method
!   This code generates excitations on the fly, and gets
!            F(D)=FMCPR3B(NI,BETA,I_P,IPATH,D,NEL,                                &
!     &          NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,       &
!     &          RHOEPS,0,RHOII,RHOIJ,CNWHTAY,I_CHMAX,LOCTAB,                     &
!     &          ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,            &
!     &          MP2E,NTOTAL)
        ELSEIF(METH.eq.-20) then !Full H-diag method
            EREF=DLWDB/TOTAL

            F(K)=FMCPR3B2(NI,BETA,I_P,IPATH,D,NEL,                               &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,         &
     &         RHOEPS,0,RHOIJ,0,METH,LOCTAB,                                     &
     &         0,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,                    &
     &         MP2E,NTOTAL,PREIV_MAX,EREF,VARSUM,PREVAR,TOTAL)

!        ELSEIF(METH.eq.-7.or.METH.eq.-19) THEN

        ENDIF
        TOTAL=TOTAL+F(D)
        DLWDB=DLWDB+DLWDB2
        ENDDO

            RETURN
        END SUBROUTINE MCPATHSPRE

SUBROUTINE BRENTALGO(brent,ax,bx,cx,tol,xmin,NI,BETA,I_P,IPATH,I,NEL,NBASISMAX,        &
     &                 G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,   &
     &                 RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,       &
     &                 NTOTAL,DLWDB,TOTAL,GIDHO,PREVAR,TLOGP)
    USE HElement
    IMPLICIT NONE
    include 'vmc.inc'
    include 'basis.inc'
    TYPE(BasisFN) G1(*)
    INTEGER NEL,I_P,BRR(*),METH,CYCLES,NMSH(*),NMAX,NTAY,I,L,LT
    INTEGER NI(NEL),NBASISMAX(5,2),IFRZ(0:NBASIS,PREIV_MAX)
    INTEGER IPATH(NEL,0:PREIV_MAX),LOCTAB(3,PREIV_MAX),NBASIS,GIDHO
    COMPLEX*16 FCK(*)
    LOGICAL TSYM,PREVAR,TLOGP
    REAL*8 BETA,ECORE,NTOTAL,ALAT(3),RHOEPS,DBETA,VARSUM
!    REAL*8 VARIANCE(2,2:PREIV_MAX)
    TYPE(HElement) UMat(*),TMat(*),RHOIJ(0:PREIV_MAX,0:PREIV_MAX)
    TYPE(HDElement) TOTAL,DLWDB,DLWDB2,EREF,FMCPR3B,FMCPR3B2,F(2:PREIV_MAX)
    TYPE(HElement) HIJS(0:PREIV_MAX)
    TYPE(HDElement) MP2E(2:PREIV_MAX),RHOII(0:PREIV_MAX)
    INTEGER ITMAX
    REAL*8 brent,ax,bx,cx,tol,xmin,CGOLD,ZEPS
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
    SELECT CASE (GIDHO)

    ! Looking through importance variable
    CASE (2)
        G_VMC_PI=x
        
    ! Looking through C excitweighting variable
    CASE (3)
        G_VMC_EXCITWEIGHT=x
    END SELECT
    call MCPATHSPRE(NI,BETA,I_P,IPATH,I,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,            &
    &        FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
    &        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,VARSUM,PREVAR)
    ! Vertex level, iteration number, C weight value, Varsumnu value
    IF (TLOGP) THEN
        WRITE(31,"(2I3,2G25.16)")  I, 0, x, VARSUM
        CALL FLUSH(31)
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
     &         p.ge.q*(b-x)) goto 1  !The above conditions determine the acceptability of the parabolic fit
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
        SELECT CASE (GIDHO)
        ! Looking through Importance variable
        CASE (2)
            G_VMC_PI=u
        
        ! Looking through C excitweighting variable
        CASE (3)
            G_VMC_EXCITWEIGHT=u
        END SELECT
        VARSUM=0.D0
        CALL MCPATHSPRE(NI,BETA,I_P,IPATH,I,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,            &
        &        FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,      &
        &        DBETA,DLWDB2,HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,VARSUM,PREVAR)
        
        ! Vertex level, iteration number, value, Varsumnu value
        IF (TLOGP) THEN
            WRITE(31,"(2I3,2G25.16)")  I, iter, u, VARSUM
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
        write (6,*) "***Brent algorithm exceeds maximum iterations of", ITMAX ,"Correct values may not have been found"
        write (31,*) "***Brent algorithm exceeds maximum iterations of", ITMAX ,"Correct values may not have been found"
     3  xmin=x
        brent=fx
        return
        END SUBROUTINE BRENTALGO
