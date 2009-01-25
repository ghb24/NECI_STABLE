MODULE Logging

    IMPLICIT NONE
    Save

    INTEGER ILOGGING,ILOGGINGDef,iGlobalTimerLevel,nPrintTimer,G_VMC_LOGCOUNT
    INTEGER HFLOGLEVEL,iWritePopsEvery
    INTEGER PreVarLogging,WavevectorPrint,NoHistBins
    REAL*8 MaxHistE
    LOGICAL TDistrib,TPopsFile,TCalcWavevector,TDetPops
    LOGICAL TZeroProjE,TWriteDetE,TAutoCorr
    INTEGER NoACDets(2:4),iPopsPartEvery

    contains

    subroutine SetLogDefaults()
      != Set defaults for Logging data items.

      use default_sets
      implicit none

      NoACDets(:)=0
      TAutoCorr=.false.
      MaxHistE=50.D0
      NoHistBins=200
      iWritePopsEvery=100000
      TCalcWavevector=.false.
      WavevectorPrint=100
      TPopsFile=.true.
      TDistrib=.false.
      ILOGGINGDef=0
      iGlobalTimerLevel=40
      nPrintTimer=10
      HFLOGLEVEL=0
      PreVarLogging=0
      TDetPops=.false.
      TZeroProjE=.false.
      TWriteDetE=.false.
      iPopsPartEvery=1

! Feb08 defaults
      IF(Feb08) THEN
          !Mcpaths set
          ILOGGINGDef=2
      ENDIF

    end subroutine SetLogDefaults



    SUBROUTINE LogReadInput()
      USE input
      IMPLICIT NONE
      LOGICAL eof
      INTEGER :: i
      CHARACTER (LEN=100) w

      ILogging=iLoggingDef

      logging: do
        call read_line(eof)
        if (eof) then
            exit
        end if
        call readu(w)
        select case(w)
        case("AUTOCORR")
!This is a Parallel FCIMC option - it will calculate the largest weight MP1 determinants and histogramm them
!HF Determinant is always histogrammed. NoACDets(2) is number of doubles. NoACDets(3) is number of triples and NoACDets(4) is 
!number of quads to histogram.
            TAutoCorr=.true.
            do i=2,4
                IF(item.lt.nitems) call readi(NoACDets(i))
            enddo
        case("DETPOPS")
!This option no longer works...
            TDetPops=.true.
        case("DISTRIBS")
            TDistrib=.true.
        case("POPSFILE")
! This is so that the determinants at the end of the MC run are written
! out, to enable them to be read back in using READPOPS in the Calc section,
! if you want to restart the simulation at a later date.  !iWritePopsEvery
! will write the configuration of particles out each time the iteration
! passes that many.
            TPopsFile=.true.
            IF(item.lt.nitems) THEN
                call readi(iWritePopsEvery)
                IF(iWritePopsEvery.lt.0) THEN
!If a negative argument is supplied to iWritePopsEvery, then the POPSFILE will never be written out, even at the end of a simulation.
                    TPopsFile=.false.
                ENDIF
            ENDIF
        case("REDUCEDPOPSFILE")
!A reduced popsfile works in exactly the same way as a normal popsfile, but only every iPopsPartEvery particle is printed out.
            TPopsFile=.true.
            call readi(iWritePopsEvery)
            call readi(iPopsPartEvery)
        case("WRITEDETE")
!This logging option will write out the energies of all determinants which have been spawned at in the simulation
! The two input options are the number of bins, and the maximum determinant energy to be histogrammed.
            TWriteDetE=.true.
            IF(item.lt.nitems) call readi(NoHistBins)
            IF(item.lt.nitems) call readf(MaxHistE)
        case("ZEROPROJE")
! This is for FCIMC when reading in from a POPSFILE. If this is on, then the energy 
! estimator will be restarted.
            TZeroProjE=.true.
        case("WAVEVECTORPRINT")
! This is for FCIMC - if on, it will calculate the exact eigenvector and
! values initially, and then print out the running wavevector every
! WavevectorPrint MC steps. However, this is slower.
            TCalcWavevector=.true.
            call readi(WavevectorPrint)
        case("MCPATHS")
            ILOGGING = IOR(ILOGGING,2**1)
        case("BLOCKING")
            ILOGGING = IOR(ILOGGING,2**13)
        case("PREVAR")
            ILOGGING = IOR(ILOGGING,2**14)
        case("FMCPR")
!  We log the value
            ILOGGING = IOR(ILOGGING,2**0)
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
               case("LABEL")
                   ILOGGING = IOR(ILOGGING,2**2)
               case("RHO")
                   ILOGGING = IOR(ILOGGING,2**3)
               case("1000")
                   ILOGGING = IOR(ILOGGING,2**9)
               case("EXCITATION")
                   ILOGGING = IOR(ILOGGING,2**12)
               case("XIJ")
                   ILOGGING = IOR(ILOGGING,2**6)
               case("")
                   ILOGGING = IOR(ILOGGING,2**2)
               case default
                  CALL report("Logging keyword FMCPR "//trim(w)       &
     &               //" not recognised",.true.)
               end select
            enddo
        case("CALCPATH")
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
               case("LABEL")
                   ILOGGING = IOR(ILOGGING,2**4)
               case("RHO")
                   ILOGGING = IOR(ILOGGING,2**5)
               case("")
                   ILOGGING = IOR(ILOGGING,2**4)
               case default
                  CALL report("Logging keyword CALCPATH "//trim(w)    &
     &               //" not recognised",.true.)
               end select
            enddo
        case("XIJ")
            ILOGGING = IOR(ILOGGING,2**6)
        case("HAMILTONIAN")
            ILOGGING = IOR(ILOGGING,2**7)
        case("PSI")
            ILOGGING = IOR(ILOGGING,2**8)
        case("TIMING")
            do while(item.lt.nitems)
                call readu(w)
                select case(w)
                case("LEVEL")
                    call readi(iGlobalTimerLevel)
                case("PRINT")
                    call readi(nPrintTimer)
                case default
                    call reread(-1)
                    call readi(iGlobalTimerLevel)
                end select
            end do
        case("VERTEX")
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
              ! case("1000")
              !     ILOGGING = IOR(ILOGGING,2**9)
               case("EVERY")
                   ILOGGING = IOR(ILOGGING,2**10)
               case default
                  call reread(-1)
                  call geti(G_VMC_LOGCOUNT)
                  ILOGGING = IOR(ILOGGING,2**9)
!                  CALL report("Logging keyword VERTEX "//trim(w)    &
!     &               //" not recognised",.true.)
               end select
            end do
        case("HFBASIS")
            ILOGGING = IOR(ILOGGING,2**11)
        case("HFLOGLEVEL")
            call geti(HFLOGLEVEL)
        case("SAVEPREVARLOGGING")
             PreVarLogging=iLogging
             iLogging=iLoggingDef
        case("ENDLOG")
            exit logging
        case default
           CALL report("Logging keyword "//trim(w)//" not recognised",.true.)
        end select
      end do logging
    END SUBROUTINE LogReadInput

END MODULE Logging
