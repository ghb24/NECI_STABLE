      INTEGER NUHFDET(1000)
      LOGICAL  TUSEBRILLOUIN
      INTEGER nActiveSpace(2)
      LOGICAL TSTAR
      LOGICAL tStoreAsExcitations
      REAL*8  E0HFDET
      LOGICAL TQUADRHO
      LOGICAL TEXPRHO
!  True if Epstein-Nesbet PT rather than Moller-Plesset
      LOGICAL tENPT  
      COMMON /UHFDET/ NUHFDET,TUSEBRILLOUIN,nActiveSpace,TSTAR,         &
     &   tStoreAsExcitations,E0HFDET,tENPT,TQUADRHO,TEXPRHO
      
