C     ==================================================================
C     == Data passed from CPMD.                                       ==
C     ==================================================================
C     == NSTATES:   number of states.                                 ==
C     == NROTOP:    number of rotational symmetry operations.         ==
C     == NTRANSLAT: number of translational symmetry operations.      ==
C     == ROT_CHAR:  characters of (the degenerate sets of) the        ==
C     ==            wavefunctions.                                    ==
C     == TRANS_CHAR:characters of (the degenerate sets of) the        ==
C     ==            wavefunctions.  If the translational irreducible  ==
C     ==            representations are real, then these are the      ==
C     ==            representations of the (1D) symmetry states.      ==
C     == GROUP_CHAR:characters due to the operations of the full      ==
C     ==            space group (or at least a large subgroup thereof.==
C     == EIGENVALUES:eigenvalues of the states, in rydberg.           ==
C     == XI:        finite size correction.                           ==
C     == BIGCELL:   true if the cell is non-primitive.                ==
C     == ROT/TRANS_LABEL : label characters of states by unique       ==
C     ==            representations.                                  ==
C     == IMPROPER_ROT: true if operation is improper.                 ==
C     == UMAT2D : Stores the integrals of the form <ij|ij> and <ij|ji>==
C     == <ij|ij> is stored in the upper diagaonal, <ij|ji> in the     ==
C     == off-diagonal elements of the lower triangle.                 ==
C     == K_VECTORS: Stores the momentum of each translational label   ==
C     == PROJECT:  True if states are projected according to          ==
C     == translational symmetry.                                      ==
C     == TSPACEGP:  True if using all space group operations, rather  ==
C     ==            than considering rotational and translaitonal     ==
C     ==            symmetries separately.                            ==
C     ==================================================================
      INTEGER NSTATES,NROTOP,NTRANSLAT,ROT_LABEL(NSTATES),
     &        TRANS_LABEL(NSTATES),NCOMPLETE
      REAL*8  ROT_CHAR(NSTATES,NROTOP),EIGENVALUES(NSTATES),XI
      COMPLEX*16  TRANS_CHAR(NSTATES,NTRANSLAT),
     &            GROUP_CHAR(NTRANSLAT*NROTOP,NSTATES),
     &            TRANS_CHAR_TABLE(NTRANSLAT,NTRANSLAT)
      INTEGER K_VECTORS(3,NTRANSLAT),CELLSIZES(3)
      LOGICAL BIGCELL,IMPROPER_OP(NROTOP),PROJECT,TSPACEGP
      REAL*8 EIONION,DEGENTOL
      POINTER (IP_EIGENVALUES,EIGENVALUES),(IP_ROT_CHAR,ROT_CHAR),
     &        (IP_TRANS_CHAR,TRANS_CHAR),(IP_ROT_LABEL,ROT_LABEL),
     &        (IP_TRANS_LABEL,TRANS_LABEL),(IP_IMPROPER_OP,IMPROPER_OP),
     &        (IP_K_VECTORS,K_VECTORS),(IP_GROUP_CHAR,GROUP_CHAR),
     &        (IP_TRANS_CHAR_TABLE,TRANS_CHAR_TABLE)
      COMMON /CPMDDATA/ NSTATES,NROTOP,NTRANSLAT,IP_ROT_CHAR,
     &                  IP_EIGENVALUES,XI,IP_TRANS_CHAR,BIGCELL,
     &                  IP_ROT_LABEL,IP_TRANS_LABEL,IP_IMPROPER_OP,
     &                  EIONION,NCOMPLETE,DEGENTOL,TSPACEGP,
     &                  IP_K_VECTORS,IP_GROUP_CHAR,IP_TRANS_CHAR_TABLE
C     ==================================================================
