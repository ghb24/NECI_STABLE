      REAL*8 G_VMC_PI,G_VMC_FAC,G_VMC_EXCITWEIGHT
      REAL*8 g_VMC_ExcitFromWeight,g_VMC_ExcitToWeight 
      REAL*8 g_VMC_ExcitToWeight2
      REAL*8 g_VMC_PolyExcitToWeight1,g_VMC_PolyExcitToWeight2
      INTEGER G_VMC_SEED, G_VMC_LOGCOUNT,PREIV_MAX
      REAL*8 g_MultiWeight(0:10)
      REAL*8 PRE_TAYREAL(2,10)
      LOGICAL PRE_TAYLOG(6,10)
      INTEGER PRE_TAY(3,10)
      LOGICAL TVARCALC(0:10)
      LOGICAL TMPTHEORY
      LOGICAL TPOLYEXCIT
      LOGICAL TVVDISALLOW
      LOGICAL TMCDIRECTSUM
      COMMON /VMC/ G_VMC_PI,G_VMC_SEED,G_VMC_FAC,                       &
     &   G_VMC_EXCITWEIGHT,g_VMC_ExcitFromWeight,                       &
     &   g_VMC_ExcitToWeight,g_MultiWeight,TMPTHEORY,                   &
     &   TVVDISALLOW,TMCDIRECTSUM,G_VMC_LOGCOUNT,                       &
     &   TVARCALC, TPOLYEXCIT,PRE_TAYREAL, PRE_TAYLOG, PRE_TAY,         &
     &   PREIV_MAX,g_VMC_ExcitToWeight2,g_VMC_PolyExcitToWeight1,       &
     &   g_VMC_PolyExcitToWeight2
