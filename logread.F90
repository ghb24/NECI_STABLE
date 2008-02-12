      MODULE LOGREAD
        USE input
        USE sysread , only : defaults,Feb08
        IMPLICIT NONE

        INTEGER ILOGGING,iGlobalTimerLevel,G_VMC_LOGCOUNT

        contains

        SUBROUTINE readinputlog()
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
      
      !Logging defaults
      ILOGGING=0
      iGlobalTimerLevel=40

! Feb08 defaults
      IF(Feb08) THEN
          !Mcpaths set
          ILOGGING=2
      ENDIF
      
        logging: do
          call read_line(eof)
          if (eof) then
              exit
          end if
          call readu(w)
          select case(w)
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
     &                 //" not recognised",.true.)
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
     &                 //" not recognised",.true.)
                 end select
              enddo
          case("XIJ")
              ILOGGING = IOR(ILOGGING,2**6)
          case("HAMILTONIAN")
              ILOGGING = IOR(ILOGGING,2**7)
          case("PSI")
              ILOGGING = IOR(ILOGGING,2**8)
          case("TIMING")
              call readi(iGlobalTimerLevel)
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
!                    CALL report("Logging keyword VERTEX "//trim(w)    &
!     &                 //" not recognised",.true.)
                 end select
              end do
          case("HFBASIS")
              ILOGGING = IOR(ILOGGING,2**11)
          case("ENDLOG")
              exit logging
          case default
             CALL report("Logging keyword "//trim(w)                    &
     &                 //" not recognised",.true.)
          end select
        end do logging
        END SUBROUTINE 

      END MODULE LOGREAD
