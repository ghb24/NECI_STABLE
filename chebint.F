! JSS 02 Sept 08:
! This file has not been touch since the initial import into svn in 2006.  It is not currently used anywhere.
! I don't know where it came from.  Candidate for removal?

C.. Chebyshev Integration using code from Chris Nex
C.. 

c Received: from [192.88.108.4] by finnis (5.52/901028.PI)
c  (for finnis) id AA01385; Thu, 19 Sep 91 15:10:38 MDT
c Received: by ixgate.gmd.de id AA06336; Thu, 19 Sep 91 15:10:26 +0200
c Date: 19 Sep 91  9:13 
c From: Chris Nex <CMMN@phoenix.cambridge.ac.uk>
c To: finnis
c Message-Id: <A4850D7BAA9FA630@UK.AC.CAMBRIDGE.PHOENIX>
c Subject: qfint

      SUBROUTINE QFINT(A,B,F,FINT,ERR,E,MMAX,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A,B,F,FINT,ERR,E
      REAL*8 ETA,SIGMAM
      INTEGER MM,IC,MMAX,IP
      EXTERNAL F
      COMMON /BLKI/XV(511),WV(511)
      IF(IP.GT.0) WRITE(IP,1)
1     FORMAT(27H NO. OF FN. CALLS  INTEGRAL,15X,5HERROR) 
      RM=4.0
      IC=1
      M=4
      J=1
      ESM=0.0
2     MM=M-1
      SIGMAM=0.0
      DO 3 I=1,MM,IC
      ETA=(B-A)*XV(J)+A
      SIGMAM=SIGMAM+WV(J)*F(ETA)
3     J=J+1
      SIGMAM=SIGMAM*2.0*(B-A)/RM
      E=ABS(ESM-SIGMAM)*0.5
      ESM=(ESM+SIGMAM)*0.5
      IF(IP.GT.0)WRITE(IP,5)MM,ESM,E
c      WRITE(6,*) ' ESM,MM,E ' , ESM,MM,E
5     FORMAT(I12,4X,3E20.12)
      M=2*M
      RM=2.0*RM
      GOTO (4,6),IC

4     IC=2
      GOTO 2
6     IF(E-ERR)8,8,7
7     IF(M-MMAX-1)2,2,8
8     FINT=ESM
      RETURN
      END

c Received: from [192.88.108.4] by finnis (5.52/901028.PI)
c  (for finnis) id AA01386; Thu, 19 Sep 91 15:10:39 MDT
c Received: by ixgate.gmd.de id AA06340; Thu, 19 Sep 91 15:10:30 +0200
c Date: 19 Sep 91  9:13 
c From: Chris Nex <CMMN@phoenix.cambridge.ac.uk>
c To: finnis
c Message-Id: <A4850D9560007070@UK.AC.CAMBRIDGE.PHOENIX>
c Subject: inital
c
      SUBROUTINE INITAL(MMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RM,M,THETA
      INTEGER J,IC
      COMMON /BLKI/XV(511),WV(511)
      PI2=6.28318530718D0
      J=1
      RM=4.0
      IC=1
      M=4
2     MM=M-1
      DO 3 I=1,MM,IC
      RI=FLOAT(I)
      THETA=RI*PI2/RM
      XV(J)=RI/RM-SIN(THETA)/PI2
      WV(J)=1.0-COS(THETA)
3     J=J+1
      M=2*M
      RM=2.0*RM
      GOTO(4,6),IC
4     IC=2
      GOTO 2
6     IF(M-MMAX-1)2,2,8
8     CONTINUE
      RETURN
      END


