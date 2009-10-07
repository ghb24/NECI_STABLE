module CalcData

implicit none

save

LOGICAL :: TSTAR,TTROT,TMCExcitSpace,TGrowInitGraph
LOGICAL :: TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
LOGICAL :: TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER,tMP2Standalone
LOGICAL :: EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
LOGICAL :: TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
LOGICAL :: TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
LOGICAL :: TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips,tConstructNOs
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL ::THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets,tDirectAnnihil,tCCMC
LOGICAL :: TStartMP1,TNoBirth,TDiffuse,TFlipTau,TExtraPartDiff,tMultipleDetsSpawn
LOGICAL :: TFullUnbias,TNodalCutoff,TNoAnnihil,TMCDiffusion,tPrintDominant
LOGICAL :: TRhoElems,TReturnPathMC,TResumFCIMC,TSignShift,tFindGuide,tUseGuide
LOGICAL :: THFRetBias,TExcludeRandGuide,TProjEMP2,TFixParticleSign,tSpawnDominant
LOGICAL :: TStartSinglePart,TRegenExcitgens,TFixShiftDoubs,tMinorDetsStar
LOGICAL :: TUnbiasPGeninProjE,TAnnihilonproc,TFixShiftShell,tNoReturnStarDets,tNoDomSpinCoup
LOGICAL :: TFixShiftKii,tMagnetize,tSymmetricField,tFixCASShift,tAllSpawnStarDets
LOGICAL :: TDistAnnihil,TLocalAnnihilation,tGlobalSftCng,tAnnihilatebyrange
LOGICAL :: tRotoAnnihil,tRegenDiagHEls,tSpawnAsDet,tFindGroundDet,tStarOrbs,tHighExcitsSing
LOGICAL :: tTruncCAS    !Truncation the FCIMC excitation space by CAS

INTEGER :: NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits
INTEGER :: NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED,HApp,iFullSpaceIter
INTEGER :: IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets,iDetGroup
INTEGER :: CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
INTEGER :: InitWalkers,NMCyc,StepsSft,FlipTauCyc,CLMax,MaxExcDom,MinExcDom
INTEGER :: RhoApp,NEquilSteps,ShellFix,NoMagDets,iGuideDets,iNoDominantDets
INTEGER :: OccCASorbs,VirtCASorbs,iStarOrbs,iHighExcitsSing,iInitGuideParts,iAnnInterval


REAL*8 :: g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
REAL*8 :: G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10)
REAL*8 :: BETAP,RHOEPSILON,DBETA(3),STARCONV,GraphBias
REAL*8 :: GrowGraphsExpo,DiagSft,Tau,SftDamp,ScaleWalkers
REAL*8 :: GrowMaxFactor,CullFactor,Lambda,NodalCutoff,PRet
REAL*8 :: FixShift,MemoryFacPart,MemoryFacAnnihil,FixedKiiCutoff
REAL*8 :: BField,MemoryFacSpawn,SinglesBias



!// additional from NECI.F
INTEGER, Allocatable :: MCDet(:)
INTEGER :: tagMCDet=0
REAL*8 :: RHOEPS ! calculated from RHOEPSILON

!// set if we include no triple-excitations as the 3rd vertex in 3+ vertex graphs.
LOGICAL :: lNoTriples

LOGICAL tFCIMCSerial

end module CalcData
