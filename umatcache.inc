      TYPE(HElement) UMATCACHE(0:NTYPES,NSLOTS,NPAIRS)
      INTEGER UMATLABELS(NSLOTS,NPAIRS)
      INTEGER NSLOTS,NPAIRS,NSTATES,NTYPES
      INTEGER NSLOTSINIT
      INTEGER NHITS,NMISSES
      POINTER (IP_UMATCACHE,UMATCACHE)
      POINTER (IP_UMATLABELS,UMATLABELS)
      TYPE(HElement) UMAT2D(NSTATES,NSTATES)
      POINTER (IP_UMAT2D,UMAT2D)
      LOGICAL TUMAT2D
      POINTER (IP_TRANSTABLE,TRANSTABLE)
      INTEGER TRANSTABLE(NSTATES)
!.. UMATCACHEFLAG=0 is normal operation
!.. UMATCACHEFLAG=1 means cache from bottom up
!..      This is useful when lots of sequential pieces of data
!..    are being sought.
      INTEGER UMATCACHEFLAG
         integer nAuxBasis,nBasisPairs
         
         Real*8 DFCoeffs(nAuxBasis,nBasisPairs)
         Real*8 DFInts(nAuxBasis,nBasisPairs)
         Pointer (IP_DFCoeffs,DFCoeffs)
         Pointer (IP_DFInts,DFInts)
      COMMON /UMATCACHE/ NSLOTS,NPAIRS,IP_UMATCACHE,IP_UMATLABELS,      &
     &   NSTATES,NSLOTSINIT,UMATCACHEFLAG,NHITS,NMISSES,IP_UMAT2D,      &
     &   TUMAT2D,IP_TRANSTABLE,nAuxBasis, nBasisPairs,IP_DFCoeffs,      &
     &   IP_DFInts,NTYPES
