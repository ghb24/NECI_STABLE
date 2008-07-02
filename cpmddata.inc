!     ==================================================================
!     == Data passed from CPMD.                                       ==
!     ==================================================================
!     == NSTATES:   number of states.                                 ==
!     == NROTOP:    number of rotational symmetry operations.         ==
!     == NTRANSLAT: number of translational symmetry operations.      ==
!     == NSpGpOp:   number of space group operations of the cell.     ==
!     == ROT_CHAR:  characters of (the degenerate sets of) the        ==
!     ==            wavefunctions.                                    ==
!     == TRANS_CHAR:characters of (the degenerate sets of) the        ==
!     ==            wavefunctions.  If the translational irreducible  ==
!     ==            representations are real, then these are the      ==
!     ==            representations of the (1D) symmetry states.      ==
!     == GROUP_CHAR:characters due to the operations of the full      ==
!     ==            space group (or at least a large subgroup thereof.==
!     == EIGENVALUES:eigenvalues of the states, in rydberg.           ==
!     == XI:        finite size correction.                           ==
!     == BIGCELL:   true if the cell is non-primitive.                ==
!     == ROT/TRANS_LABEL : label characters of states by unique       ==
!     ==            representations.                                  ==
!     == IMPROPER_ROT: true if operation is improper.                 ==
!     == UMAT2D : Stores the integrals of the form <ij|ij> and <ij|ji>==
!     == <ij|ij> is stored in the upper diagaonal, <ij|ji> in the     ==
!     == off-diagonal elements of the lower triangle.                 ==
!     == K_VECTORS: Stores the momentum of each translational label   ==
!     == PROJECT:  True if states are projected according to          ==
!     == translational symmetry.                                      ==
!     == TSPACEGP:  True if using all space group operations, rather  ==
!     ==            than considering rotational and translaitonal     ==
!     ==            symmetries separately.                            ==
!     == NKPOINTS:  Number of k-points.                               ==
!     == KPNT_CHAR: Representation of a state with a given k-point    ==
!     ==            the translational symmetry operations of the      ==
!     ==            equivalent supercell.                             ==
!     == NKPTRANS:  Number of translational operations of supercell   ==
!     ==            equivalent to the k-point mesh.                   ==
!     == NKPS: Number of k-points (excluding symmetry equivalent ones.)=
!     ==================================================================
      INTEGER NSTATES,NROTOP,NTRANSLAT,ROT_LABEL(NSTATES),              &
     &        TRANS_LABEL(NSTATES),NCOMPLETE,NKPS,NKPTRANS,NSpGpOp      
      REAL*8  ROT_CHAR(NSTATES,NROTOP),EIGENVALUES(NSTATES),XI,         &
     &        PIInt(NStates)
      COMPLEX*16  TRANS_CHAR(NSTATES,NTRANSLAT),                        &
     &            GROUP_CHAR(nSpGpOp,NSTATES),                          &
     &            TRANS_CHAR_TABLE(NTRANSLAT,NTRANSLAT),                &
     &            KPNT_CHAR(NKPTRANS,NKPS)
      INTEGER K_VECTORS(3,NTRANSLAT),CELLSIZES(3)
      LOGICAL BIGCELL,IMPROPER_OP(NROTOP*NTRANSLAT),PROJECT,TSPACEGP
      REAL*8 EIONION,DEGENTOL
      POINTER (IP_EIGENVALUES,EIGENVALUES),(IP_ROT_CHAR,ROT_CHAR),      &
     &        (IP_TRANS_CHAR,TRANS_CHAR),(IP_ROT_LABEL,ROT_LABEL),      &
     &        (IP_TRANS_LABEL,TRANS_LABEL),(IP_IMPROPER_OP,IMPROPER_OP),&
     &        (IP_K_VECTORS,K_VECTORS),(IP_GROUP_CHAR,GROUP_CHAR),      &
     &        (IP_TRANS_CHAR_TABLE,TRANS_CHAR_TABLE),                   &
     &        (IP_KPNT_CHAR,KPNT_CHAR),(IP_PIInt,PIInt)
      LOGICAL TKP
      POINTER (IP_KPNTINDEXNECI,KPNTIND)
      INTEGER KPNTIND(2*nkps*nStates)
      COMPLEX*16 PICorr(NSTATES)
      POINTER    (IP_PICorr,PICorr)
      COMMON /CPMDDATA/ NSTATES,NROTOP,NTRANSLAT,NKPTRANS,IP_ROT_CHAR,  &
     &                  IP_EIGENVALUES,XI,IP_TRANS_CHAR,BIGCELL,NSpGpOp,&
     &                  IP_ROT_LABEL,IP_TRANS_LABEL,IP_IMPROPER_OP,     &
     &                 EIONION,NCOMPLETE,DEGENTOL,TSPACEGP,IP_KPNT_CHAR,&
     &                 IP_K_VECTORS,IP_GROUP_CHAR,IP_TRANS_CHAR_TABLE,  &
     &                 TKP,NKPS,IP_KPNTINDEXNECI,IP_PICorr,IP_PIInt
!     ==================================================================
