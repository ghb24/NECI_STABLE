      REAL*8 G_VMC_PI,G_VMC_FAC,G_VMC_EXCITWEIGHT
      REAL*8 g_VMC_ExcitFromWeight,g_VMC_ExcitToWeight
      INTEGER G_VMC_SEED
      REAL*8 g_MultiWeight(0:10)
      COMMON /VMC/ G_VMC_PI,G_VMC_SEED,G_VMC_FAC,                       &
     &   G_VMC_EXCITWEIGHT,g_VMC_ExcitFromWeight,                       &
     &   g_VMC_ExcitToWeight,g_MultiWeight
