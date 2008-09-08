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
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL :: THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets
LOGICAL :: TStartMP1,TNoBirth,TDiffuse,TFlipTau,TExtraPartDiff
LOGICAL :: TFullUnbias,TNodalCutoff,TNoAnnihil,TMCDiffusion
LOGICAL :: TRhoElems,TReturnPathMC,TResumFCIMC,TSignShift
LOGICAL :: THFRetBias,TExcludeRandGuide,TProjEMP2,TFixParticleSign
LOGICAL :: TStartSinglePart,TRegenExcitgens

INTEGER :: NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits
INTEGER :: NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED,HApp
INTEGER :: IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets
INTEGER :: CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
INTEGER :: InitWalkers,NMCyc,StepsSft,FlipTauCyc,CLMax
INTEGER :: RhoApp,NEquilSteps,MemoryFac


REAL*8 :: g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
REAL*8 :: G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10)
REAL*8 :: BETAP,RHOEPSILON,DBETA(3),STARCONV,GraphBias
REAL*8 :: GrowGraphsExpo,DiagSft,Tau,SftDamp,ScaleWalkers
REAL*8 :: GrowMaxFactor,CullFactor,Lambda,NodalCutoff,PRet



!// additional from NECI.F
INTEGER, Allocatable :: MCDet(:)
INTEGER :: tagMCDet=0
REAL*8 :: RHOEPS ! calculated from RHOEPSILON

!// set if we include no triple-excitations as the 3rd vertex in 3+ vertex graphs.
LOGICAL :: lNoTriples

end module CalcData
