module CCMCData
   save
   real*8   dT1SqCuml
   logical  tExactCluster  ! Go through all combinations of det
   logical  tCCMCFCI       ! Run CCMC code without excitation clusters, recovering the FCIMC result
   logical  tAmplitudes    ! Use real numbers to indicate the amplitudes rather than stochastically sampling
end module CCMCData
