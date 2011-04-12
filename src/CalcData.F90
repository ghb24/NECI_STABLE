module CalcData

    use constants, only: dp
    use MemoryManager, only: TagIntType
    implicit none

    save

LOGICAL :: TSTAR,TTROT,TMCExcits,TGrowInitGraph
LOGICAL :: TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
LOGICAL :: TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER,tMP2Standalone
LOGICAL :: EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
LOGICAL :: TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
LOGICAL :: TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
LOGICAL :: TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL ::THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets,tDirectAnnihil,tCCMC
LOGICAL :: TFullUnbias,TNoAnnihil,tStartMP1
LOGICAL :: TRhoElems,TReturnPathMC,TSignShift
LOGICAL :: THFRetBias,TProjEMP2,TFixParticleSign
LOGICAL :: TStartSinglePart,TRegenExcitgens
LOGICAL :: TUnbiasPGeninProjE, tCheckHighestPopOnce
LOGICAL :: tCheckHighestPop,tRestartHighPop,tChangeProjEDet
LOGICAL :: tRotoAnnihil,tRegenDiagHEls,tSpawnAsDet,tFindGroundDet
LOGICAL :: tTruncCAS,tTruncInitiator,tDelayTruncInit,tAddtoInitiator    !Truncation the FCIMC excitation space by CAS
LOGICAL :: tInitIncDoubs,tWalkContGrow,tAnnihilatebyRange,tRetestAddtoInit
logical :: tReadPopsRestart, tReadPopsChangeRef

logical :: tStartCAS    !Start FCIMC dynamic with walkers distributed according to CAS diag.

! Base hash values only on spatial orbitals
! --> All dets with same spatial structure on the same processor.
logical :: tSpatialOnlyHash

! Do we allow walkers to survive (in the initiator approx.) if a determinant
! with the same spatial configuration is an initiator?
logical :: tSpawnSpatialInit

!These options mean that only initiators can spawn walkers.
!tSpawn_Only_Init_Grow means that this option is removed once variable shift is entered.
logical :: tSpawn_Only_Init,tSpawn_Only_Init_Grow

! Do we truncate spawning based on the number of unpaired electrons
logical :: tTruncNOpen
integer :: trunc_nopen_max

logical :: tMaxBloom    !If this is on, then we only print out a bloom warning if it is the biggest to date.

INTEGER :: NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits,IterTruncInit,InitiatorWalkNo,NShiftEquilSteps
INTEGER :: NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED,HApp,iFullSpaceIter
INTEGER :: IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets,iDetGroup
INTEGER :: CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
INTEGER :: InitWalkers,NMCyc,StepsSft,CLMax
INTEGER :: NEquilSteps,InitialPart
INTEGER :: OccCASorbs,VirtCASorbs,iAnnInterval
integer :: iPopsFileNoRead, iPopsFileNoWrite,iWeightPopRead,iRestartWalkNum
INTEGER *8 :: MaxNoatHF,HFPopThresh

REAL*8 :: g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
REAL*8 :: G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10)
REAL*8 :: BETAP,RHOEPSILON,DBETA,STARCONV,GraphBias
REAL*8 :: GrowGraphsExpo,Tau,SftDamp,ScaleWalkers
REAL*8 :: GrowMaxFactor,CullFactor,PRet,FracLargerDet
REAL*8 :: MemoryFacPart,MemoryFacAnnihil
REAL*8 :: MemoryFacSpawn,SinglesBias,TauFactor,StepsSftImag

real(dp) :: MemoryFacInit

real(dp), target :: DiagSft

REAL*8 :: GraphEpsilon
REAL*8 :: PGenEpsilon


!// additional from NECI.F
INTEGER, Allocatable :: MCDet(:)
INTEGER(TagIntType) :: tagMCDet=0
REAL*8 :: RHOEPS ! calculated from RHOEPSILON

!// set if we include no triple-excitations as the 3rd vertex in 3+ vertex graphs.
LOGICAL :: lNoTriples

LOGICAL tUseProcsAsNodes  !Set if we treat each processor as its own node.
INTEGER iLogicalNodeSize  !An alternative to the above, create logical nodes of at most this size.
                          ! 0 means use physical nodes.

end module CalcData
