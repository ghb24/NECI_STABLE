! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 

      SUBROUTINE LOGNAN(NI,NEL,BETA,ITIME)
      use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,ITIME,NI(NEL),K,I,neci_system
         real(dp) BETA
         CHARACTER(40) STR
         WRITE(STR,*) 'mkdir "nan',ITIME,'"'
         I=neci_SYSTEM(STR)
         WRITE(STR,*) 'cp * "nan',ITIME,'"'
         I=neci_SYSTEM(STR)
         WRITE(STR,*) 'nan',ITIME,"/LOG"
         OPEN(12,FILE=STR,STATUS="UNKNOWN")
         WRITE(12,*) "WARNING: nan found at time",ITIME
         WRITE(12,*) "Beta=",BETA
         WRITE(12,"(A)",advance='no') "  nan det=("
         DO K=1,NEL
            WRITE(12,"(I3,A)",advance='no') NI(K),","
         ENDDO
         WRITE(12,"(A)",advance='no') "),"
         CLOSE(12)
         RETURN

      END


      SUBROUTINE WRITE_XMATRIX(IUNIT,RHOIJ,I_V)
      use constants, only: dp
         IMPLICIT NONE
         INTEGER I_V,J,K,IUNIT
         real(dp) RHOIJ(0:I_V-1,0:I_V-1)      
         WRITE(IUNIT,"(A)",advance='no') "{"
         DO J=0,I_V-1
            DO K=0,I_V-1
               WRITE(IUNIT,"(E25.16,A)",advance='no') RHOIJ(J,K),","
            ENDDO
            WRITE(IUNIT,"(A)",advance='no') "|"
         ENDDO
         WRITE(IUNIT,*) "}"
         RETURN
      END

