
      MODULE PRECALCREAD
        USE input
        USE sysread , only : defaults,Feb08
        IMPLICIT NONE

        LOGICAL TPREVAR,PRE_TAYLOG(6,10),TLINEVAR(2:10),TGRIDVAR(2:10)
        LOGICAL MEMSAV(2:10)
        INTEGER PREIV_MAX,PRE_TAY(3,10),TRUECYCLES,USEVAR(2:10,8)
        REAL*8 PRE_TAYREAL(2,10),TOTALERROR,PREWEIGHTEPS
        REAL*8 GRIDVARPAR(2:10,6),LINEVARPAR(2:10,6)

        contains

        SUBROUTINE readinputprecalc()
        IMPLICIT NONE
        LOGICAL eof
        INTEGER zz,z,k
        CHARACTER (LEN=100) w
        
! Precalc defaults      
          TPREVAR=.false.
          PREIV_MAX=0      
          PRE_TAY(:,:)=0
          PRE_TAYLOG(:,:)=.false.
          PRE_TAYREAL(:,:)=0.D0
          PRE_TAYREAL(2,:)=0.1      !TOLERANCE   
          TLINEVAR(:)=.false.
          TGRIDVAR(:)=.false.
          TRUECYCLES=0
          TOTALERROR=0.D0
          PREWEIGHTEPS=1.D-8
          USEVAR(:,:)=0
          MEMSAV(:)=.false.
          GRIDVARPAR(:,:)=0.D0
          LINEVARPAR(:,:)=0.D0

!Feb08 defaults...
          IF(Feb08) THEN
              !...add defaults...
          ENDIF

        preIV_MAX=1
!        preIH=-10
        precalc: do
           call read_line(eof)
           if (eof) then
               call report("No arguments for PRECALC block",.true.)
           end if
           call readu(w)
           select case(trim(w))
           case("VERTEX")
               preIV_MAX=preIV_MAX+1
               pre_TAY(3,preIV_MAX)=preIV_MAX
               call inpgetprecalc(pre_TAY(1,preIV_MAX))
            
           case("TOLERANCE")
               call readf(pre_TAYREAL(2,preIV_MAX))
               if ( pre_TAYREAL(2,preIV_MAX) .le. 0.D0 ) then
                   call report(trim(w)//" must take a value "           &
     &                  //"more than 0.D0",.true.)
               endif
           case("LINEVAR")
               TLINEVAR(preIV_MAX)=.true.
               IF(pre_TAYLOG(4,preIV_MAX)) THEN
                   call report(trim(w)//" can not be applied"           &
     &                 //" in conjunction with USEIMPORT",.true.)
               ENDIF
               do zz=1,3
                   call readf(LINEVARPAR(preIV_MAX,zz))
               enddo
           case("GRIDVAR")
               TGRIDVAR(preIV_MAX)=.true.
               IF(pre_TAYLOG(1,preIV_MAX)) THEN
                   call report(trim(w)//" can not be applied"           &
     &                     //" in conjunction with USE",.true.)
               ENDIF
               do zz=1,6
                   call readf(GRIDVARPAR(preIV_MAX,zz))
               enddo
           case("TOTALERROR")
               call readf(TOTALERROR)
           case("TRUECYCLES")
               call readi(TRUECYCLES)
           case("PREGRAPHEPSILON")
               call readf(PREWEIGHTEPS)
           case("USE")
                pre_TAYLOG(1,preIV_MAX)=.true.    
                IF(TGRIDVAR(preIV_MAX)) THEN
                    call report(trim(w)//" can not be applied"          &
     &                      //" in conjunction with GRIDVAR",.true.)
                ENDIF
                !USE is specified on its own - all vertex levels take this value
                IF(nitems.eq.1) THEN
                    do k=2,10
                        IF(k.ne.preIV_MAX) USEVAR(k,:)=0
                        IF(k.eq.preIV_MAX) THEN
                            do z=1,8
                                USEVAR(preIV_MAX,z)=z+1
                            ENDDO
                        ENDIF
                    enddo
                ELSE
                call readi(USEVAR(preIV_MAX,1))
                do k=2,8
                    if(item.lt.nitems) THEN
                        call readi(USEVAR(preIV_MAX,k))
                
                    do z=1,(preIV_MAX-1)
                        do zz=1,8
                       if (USEVAR(z,zz).eq.USEVAR(preIV_MAX,k)) then
                       call report(trim(w)//" can only be applied"      &
     &                      //" to a vertex level once",.true.)
                        endif
                        enddo
                    enddo
                
                    endif
                enddo
                ENDIF
           case("MEMORISE")
               MEMSAV(preIV_MAX)=.true.
               if ((pre_TAY(1,preIV_MAX).ne.-7).and.                    &
     &              (pre_TAY(1,preIV_MAX).ne.-19)) then
                    call report(trim(w)//" not yet valid for full"      &
     &                  //" sum precalc",.true.)
               endif
           case("NONE")
                pre_TAY(3,preIV_MAX)=0
           case("CYCLES")
               call readi(pre_TAY(2,preIV_MAX))
               if ((pre_TAY(1,preIV_MAX).ne. -7).and.                   &
     &                  (pre_TAY(1,preIV_MAX).ne. -19)) then
                  call report(trim(w)//" only valid for MC"             &
     &                 //" method",.true.)
               end if
           case("UEPSILON")
               call readf(pre_TAYREAL(1,preIV_MAX))
           case("FINDC")
               pre_TAYLOG(2,preIV_MAX)=.true.
           case("FINDD")
               pre_TAYLOG(5,preIV_MAX)=.true.
           case("USED")
               pre_TAYLOG(6,preIV_MAX)=.true.
           case("FINDIMPORT")
               if ( preIV_MAX .lt. 3 ) then
                   call report(trim(w)//" only valid for vertex "       &
     &                  //"levels of 3 or higher",.true.)
               endif
               pre_TAYLOG(3,preIV_MAX)=.true.
           case("USEIMPORT")
               if ( preIV_MAX .lt. 3 ) then
                   call report(trim(w)//" only valid for vertex "       &
     &                  //"levels of 3 or higher",.true.)
               endif
               pre_TAYLOG(4,preIV_MAX)=.true.
                   do zz=1,(preIV_MAX-1)
                      if ( pre_TAYLOG(4,zz) ) then
                         call report(trim(w)//" can only be applied"    &
     &                       //" to a single vertex level",.true.)
                     endif
                   enddo
           case("ENDPRECALC")
               exit precalc
           case default
               call report ("Keyword "//trim(w)//                       &
     &              " not recognised",.true.)
           end select
           end do precalc

       END SUBROUTINE
       
       END MODULE precalcread
