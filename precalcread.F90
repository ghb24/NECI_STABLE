
      MODULE PRECALCREAD
        IMPLICIT NONE

!TPREVAR is set to true when the precalc routines have yet to run.
        LOGICAL TPREVAR

!PRE_TAYLOG(1,precalc_vertex_level)='use' those sets of optimised parameters
!PRE_TAYLOG(2,...)=Does the 'c' (matrix element) parameter want to be optimised?
!PRE_TAYLOG(3,...)=Does the 'importance' parameter want to be optimised?
!PRE_TAYLOG(4,...)=Does the 'importance' parameter want to be used?
!PRE_TAYLOG(5,...)=Does the 'd' (inverse energy change) parameter want to be optimised?
!PRE_TAYLOG(6,...)=Does the 'd' parameter want to be used?
        LOGICAL PRE_TAYLOG(6,10)

!TLINEVAR/GRIDVAR will create a surface of expected variances for the vertex levels indicated for 1D/2D probability distributions respectivly
        LOGICAL TLINEVAR(2:10),TGRIDVAR(2:10)

!MEMSAV will hold the excitations for MC-Precalc in memory for that vertex level - memory/speed trade-off
        LOGICAL MEMSAV(2:10)

!PRE_TAY(1,...)=The method to use for the precalc at that vertex level
!PRE_TAY(2,...)=The number of cycles to use for MC-precalc optimisation at that level
!PRE_TAY(3,...)=The precalc vertex level
!PREIV_MAX = The largest vertex level to be optimised in precalc
        INTEGER PREIV_MAX,PRE_TAY(3,10)
        
!TRUECYCLES are the cycles needed in the main program to achieve an error specified by TOTALERROR
        INTEGER TRUECYCLES
        
!USEVAR indicates which precalc optimised values at each vertex level are to be used for each vertex level in the main program.
        INTEGER USEVAR(2:10,8)

!PRE_TAYREAL(1,...)= 'c epsilon' - if the difference between the optimised 'c' parameter and 0.D0 is less than this, then just set it to zero.
!PRE_TAYREAL(2,...) is the tolerance for convergence onto an optimised parameter in terms of the change in the brent algorithm for the parameter for each iteration.
!TOTALERROR is a specified desired error in the main program, from which the number of cycles needed will be approximated.
!PREWEIGHTEPS means a graph will be discarded in the MC-precalc if its weight is less than this.
        REAL*8 PRE_TAYREAL(2,10),TOTALERROR,PREWEIGHTEPS

!GRID/LINEVARPAR gives the parameters for the search when TGRID/LINEVAR is set.
        REAL*8 GRIDVARPAR(2:10,6),LINEVARPAR(2:10,6)

        contains

        SUBROUTINE PrecalcReadInput()
        USE input
        use default_sets
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

       END SUBROUTINE PrecalcReadInput
      
      subroutine inpgetprecalc(preIH)
         use input
         use SystemData , only : TUSEBRILLOUIN
         implicit none
         integer preIH
         CHARACTER(LEN=16) w
                do while ( item .lt. nitems )
                  call readu(w)
                  select case(w)
                  case("HDIAG")
                      call readu(w)
                      select case(w)
                      case("FULL")
                          preIH=-20
                      case("MC")
                          preIH=-19
                          IF(.NOT.TUSEBRILLOUIN) THEN
                              write(6,*) "Warning  USEBRILLOUINTHEOREM" &
     &                       //" might need to be specified in system " &
     &                       //"block to use MC-PRECALC"
                          ENDIF
                      case default
                         call report("Error - must specify FULL"        &
     &                   //" or MC after HDIAG in PRECALC block",.true.)
                      end select
                   case("RHODIAG")
                         call readu(w)
                         select case(w)
                         case("FULL")
                             preIH=-8
                         case("MC")
                             preIH=-7
                         case default
                           call report("Error - must specify FULL or "  &
     &                     //"MC after RHODIAG in PRECALC block",.true.)
                         end select
                     case default
                           call report("Keyword error with "//trim(w),  &
     &                     .true.)
                     end select
                 end do
        end subroutine

       
       END MODULE precalcread
